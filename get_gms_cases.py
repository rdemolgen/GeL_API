import os
import requests
import argparse
import csv
import jellypy.pyCIPAPI.auth as auth
import jellypy.pyCIPAPI.config as config
import jellypy.pyCIPAPI.interpretation_requests as irs
import jellypy.pyCIPAPI.opencga as opencga
from datetime import datetime, date
from protocols.reports_6_0_1 import InterpretedGenome, InterpretationRequestRD, CancerInterpretationRequest, ClinicalReport
from protocols.participant_1_2_0 import Referral, ReferralTest


## Define settings ##
testing_on=True
category='gms'
interpreter_organisation_name='South West Genomic Laboratory Hub'
sites=["69A50","RA7"]
assembly='GRCh38'
outdir = '/mnt/data3/genome_sequencing/GEL_GMS/'

def get_interpretation_request(case):
    """
    Return the interpretation request for a case
    """
    ir_id, ir_version = case['interpretation_request_id'].split('-')
    ir = irs.get_interpretation_request_json(ir_id=ir_id, ir_version=ir_version, testing_on=testing_on)
    return ir

def get_ir_obj(ir):
    """
    Return the json data for the interpretation request
    """
    json_request_data = ir["interpretation_request_data"]["json_request"]
    # return rare disease cases only
    if ir["sample_type"] == 'raredisease':
        ir_obj = InterpretationRequestRD.fromJsonDict(json_request_data)
    return ir_obj

def get_referral_data(ir):
    """
    Use the referral protocol to return the referral data for a case
    """
    try:
        referral_obj = Referral.fromJsonDict(ir['referral']['referral_data'])
        # get the test order date
        test_order_date_dict = referral_obj.referralTests[0].referralTestOrderingDate.toJsonDict()
        test_order_date = date(test_order_date_dict['year'], test_order_date_dict["month"], test_order_date_dict["day"])

        for member in referral_obj.pedigree.members:
            if member.isProband:
                participant_id = member.participantId
                participant_uid = member.participantUid       
                sex = member.sex.lower()
                hpo_list = []
                for hpo_term in member.hpoTermList:
                    hpo_list.append(hpo_term.term)
                hpo_terms = hpo_list

        referral_data = {
            'referral_id': referral_obj.referralId,
            'participant_id': participant_id,
            'hpo_terms': hpo_terms,
            'requester': referral_obj.requester.organisationName,
            'interpreterName': referral_obj.referralTests[0].interpreter.organisationName,
            'interpreterCode': referral_obj.referralTests[0].interpreter.organisationCode,
            'clinical_indication': referral_obj.clinicalIndication.clinicalIndicationFullName,
            'clinical_indication_code': referral_obj.clinicalIndication.clinicalIndicationCode,
            'test_order_date': test_order_date,
        }

    
    except:
        #print("No referral data")
        referral_data = {
            'referral_id': "No referral data",
            'participant_id': "unknown",
            'requester': "unknown",
            'interpreterName': "unknown",
            'clinical_indication': "unknown",
            'test_order_date': "unknown",
        }
    
    return referral_data

def get_family_details(ir_obj):
    """
    Get the sequenced and affected family memebrs for a case.
    """

    family = []

    for member in ir_obj.pedigree.members:
        if member.isProband:
            family.append({
                'participant_id':member.participantId, 
                'sampleId':member.samples[0].sampleId, 
                'relation':'Proband',
                'sex': member.sex,
                'affected':member.affectionStatus})
        if (not member.isProband) and (member.samples):
            if member.additionalInformation:
                family.append({
                    'participant_id':member.participantId, 
                    'sampleId':member.samples[0].sampleId, 
                    'relation':member.additionalInformation['relation_to_proband'], 
                    'sex': member.sex,
                    'affected':member.affectionStatus})

    return family

def output_sample_data(case_data):
    """
    Append sample member data to text file for running pipeline
    """

    out_file = outdir + 'sample_details.tsv'

    sample_details = []

    # reformat case_data to sample_details format
    for member in case_data['members']:
        if member['relation'] == 'Proband':
            sample_details.append({
                'sampleId': member['sampleId'],
                'referral_id': case_data['referral_id'],
                'participant_id': member['participant_id'],
                'sex': member['sex'],
                'relation': member['relation'],
                'clinical_indication': case_data['clinical_indication']
            })

        else:
            sample_details.append({
                'sampleId': member['sampleId'],
                'referral_id': case_data['referral_id'],
                'participant_id': member['participant_id'],
                'sex': member['sex'],
                'relation': member['affected'].capitalize() + ' ' + member['relation'],
                'clinical_indication': ''
            })

    with open(out_file, 'a+') as csvfile:
        fieldnames = [
            'sampleId',
            'referral_id',
            'participant_id',
            'sex',
            'relation',
            'clinical_indication'
        ]

        # remove sample from list of details to be added if already exists in file
        csvfile.seek(0)
        for line in csvfile.readlines():
            for sample in sample_details:
                if sample['sampleId'] in line:
                    sample_details = [i for i in sample_details if not (i['sampleId'] == sample['sampleId'])]
                    #print(sample_details)
        
        # write sample details list to file
        writer = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)        
        writer.writerows(sample_details)


def output_case_data(cip_api_data):
    """
    Write a summary of case details to a file
    """
    todays_date = str(datetime.today().date())
    out_file = outdir + 'pilot_cases_' + todays_date + '.csv'

    with open(out_file, 'w') as csvfile:
        fieldnames = [
            'referral_id',
            'status',
            'last_modified',
            'case_id',
            'proband_id',
            'proband',
            'mother',
            'father',
            'family_structure',
            'affected',
            'clinical_indication',
            'hpo_terms',     
            'requester',
            'interpreter_name',
            'interpreter_code',
            'test_order_date'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(cip_api_data)


def get_ir_list():
    """
    Get GEL family information
    """
    
    cip_api_data = []
    
    irlist_all = irs.get_interpretation_request_list(sample_type='raredisease', testing_on=testing_on, interpreter_organisation_name=interpreter_organisation_name)

    for case in irlist_all:
        # get interpretation request for each case
        ir = get_interpretation_request(case)
        # get json request data for each case    
        ir_obj = get_ir_obj(ir)
        # get referral data for each case
        referral_data = get_referral_data(ir)
        
        # define family details
        family = get_family_details(ir_obj)

        if len(family) == 1:
            family_structure = "Singleton"
        elif (len(family) == 3) and ((set([member['relation'] for member in family]) == set(['Proband', 'Mother', 'Father'])) == True):
            family_structure = "Trio"
        else:
            family_structure = ', '.join([member['relation'] for member in family])

        affected = ', '.join([member['relation'] for member in family if member['affected'] == 'AFFECTED'])

        # get the sample Ids for the trio members
        if family_structure == 'Trio':
            mother_sampleId = ''.join([member['sampleId'] for member in family if member['relation'] == 'Mother' ])
            father_sampleId = ''.join([member['sampleId'] for member in family if member['relation'] == 'Father' ])
        else:
            mother_sampleId = ''
            father_sampleId = ''

        proband_sampleId = ''.join([member['sampleId'] for member in family if member['relation'] == 'Proband' ])

        case_data = {
            'referral_id': case['family_id'],
            'status': case['last_status'],
            'last_modified': datetime.strptime(case['last_modified'], '%Y-%m-%dT%H:%M:%S.%fZ').date(),
            'case_id': case['case_id'],
            'proband_id': case['proband'],
            'proband': proband_sampleId,
            'mother': mother_sampleId,
            'father': father_sampleId,
            'members': family,
            'family_structure': family_structure,
            'affected': affected,
            'clinical_indication': referral_data['clinical_indication'],
            'hpo_terms': referral_data['hpo_terms'],     
            'requester': referral_data['requester'],
            'interpreter_name': referral_data['interpreterName'],
            'interpreter_code': referral_data['interpreterCode'],
            'test_order_date': referral_data['test_order_date']
        }

        # write sample data to sample_details file for pipeline input
        output_sample_data(case_data)
        
        cip_api_data.append(case_data)

    # write all cases to file
    output_case_data(cip_api_data)


if __name__ == '__main__':
    get_ir_list()

    


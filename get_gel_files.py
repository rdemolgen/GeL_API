import os
import requests
import argparse
import csv
import jellypy.pyCIPAPI.auth as auth
import jellypy.pyCIPAPI.config as config
import jellypy.pyCIPAPI.interpretation_requests as irs
import jellypy.pyCIPAPI.opencga as opencga
from datetime import datetime
from protocols.reports_6_0_1 import InterpretedGenome, InterpretationRequestRD, CancerInterpretationRequest, ClinicalReport
from protocols.participant_1_2_0 import Referral, ReferralTest


testing_on=True
category='gms'
interpreter_organisation_name='South West Genomic Laboratory Hub'
assembly='GRCh38'
base_path = '/mnt/data3/genome_sequencing/GEL_GMS/'


def get_interpretation_request(case):
    """
    Return the interpretation request for a case
    """
    ir_id, ir_version = case['interpretation_request_id'].split('-')
    ir = irs.get_interpretation_request_json(ir_id=ir_id, ir_version=ir_version, testing_on=testing_on)
    return ir


def setup_folder(family_id):
    """
    Make a new directory for the GEL family
    """
    family_folder = os.path.join(base_path, family_id)
    if not os.path.exists(family_folder):
        os.makedirs(family_folder)
    return family_folder


def get_gel_crams(ir, family_folder):
    """
    Get the file name of the bams from the ir json data
    """
    json_request_data = ir["interpretation_request_data"]["json_request"]
    # find the file names of the bams
    for bam in json_request_data["bams"]:
        file_name = bam["uriFile"].split('/')[-1]
        download_files(file_name, family_folder)


def download_files(file_name, family_folder):
    """
    Download the bam files from OpenCGA
    """
    # get the study ID
    #study_id = opencga.get_study_id(study_type="raredisease", assembly=assembly)
    
    study_id = "1053593329"
    # get the file id for the crams
    file_id = opencga.find_file_id(study_id=study_id, file_format="CRAM", file_name=file_name)
    # download the cram files and output to the family folder
    print("Downloading files...")
    opencga.download_file(file_id=file_id, study_id=study_id, file_name=file_name, download_folder=family_folder)


def get_files(case):
    """
    Find the interpretation request based on the family_id and use that information to download the cram files from GEL
    """

    irlist = irs.get_interpretation_request_list(sample_type='raredisease', testing_on=testing_on, interpreter_organisation_name=interpreter_organisation_name, family_id=case)

    for case in irlist:
        case_id = case['case_id']
        participant_id = case['proband']
        status = case['last_status']
        last_modified = datetime.strptime(case['last_modified'], '%Y-%m-%dT%H:%M:%S.%fZ').date()
        sample_type = case['sample_type']
        family_id = case['family_id']
        sites = ', '.join(case['sites'])

        # return the interpretation request for the case
        ir = get_interpretation_request(case)
        # create a new folder for each family
        family_folder = setup_folder(family_id)
        # download the files from OpenCGA
        get_gel_crams(ir, family_folder)



if __name__ == '__main__':
    with open('cases_to_process', 'r') as inFile:
        for case in inFile:
            get_files(case)

    


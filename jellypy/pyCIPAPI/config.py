#!/usr/bin/env python

# Configuration file for setting common variables to avoid hard-coding them in code:

# Set to true to use Active Directory authentication, or false to use legacy LDAP authentication
use_active_directory = True

# CIP-API AD authentication URLs
live_100K_auth_url = 'https://login.microsoftonline.com/0a99a061-37d0-475e-aa91-f497b83269b2/oauth2/token'
beta_testing_auth_url = 'https://login.microsoftonline.com/afee026d-8f37-400e-8869-72d9124873e4/oauth2/token'

# CIP-API base URLs for live data and beta testing:
live_100k_data_base_url = 'https://cipapi.genomicsengland.nhs.uk/api/2/'
beta_testing_base_url = 'https://cipapi-gms-beta.genomicsengland.nhs.uk/api/2/' #uat
# sit beta_testing_base_url = 'https://cipapi-beta.genomicsengland.co.uk/api/2/'

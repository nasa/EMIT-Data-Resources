# -*- coding: utf-8 -*-
"""
---------------------------------------------------------------------------------------------------
 How to Set Up Direct Access the LP DAAC Data Pool with Python
 The following Python code will configure a netrc profile that will allow users to download data
 from an Earthdata Login enabled server. See README for additional information
---------------------------------------------------------------------------------------------------
 Author: Cole Krehbiel
 Last Updated: 11/20/2018
-------------------------------------------------------------------------------
"""
# Load necessary packages into Python
from netrc import netrc
from subprocess import Popen
from getpass import getpass
import os

# -----------------------------------AUTHENTICATION CONFIGURATION-------------------------------- #
urs = 'urs.earthdata.nasa.gov'    # Earthdata URL to call for authentication
prompts = ['Enter NASA Earthdata Login Username \n(or create an account at urs.earthdata.nasa.gov): ',
           'Enter NASA Earthdata Login Password: ']

# Determine if netrc file exists, and if so, if it includes NASA Earthdata Login Credentials
try:
    netrcDir = os.path.expanduser("~/.netrc")
    netrc(netrcDir).authenticators(urs)[0]

# Below, create a netrc file and prompt user for NASA Earthdata Login Username and Password
except FileNotFoundError:
    homeDir = os.path.expanduser("~")
    Popen('touch {0}.netrc | chmod og-rw {0}.netrc | echo machine {1} >> {0}.netrc'.format(homeDir + os.sep, urs), shell=True)
    Popen('echo login {} >> {}.netrc'.format(getpass(prompt=prompts[0]), homeDir + os.sep), shell=True)
    Popen('echo password {} >> {}.netrc'.format(getpass(prompt=prompts[1]), homeDir + os.sep), shell=True)

# Determine OS and edit netrc file if it exists but is not set up for NASA Earthdata Login
except TypeError:
    homeDir = os.path.expanduser("~")
    Popen('echo machine {1} >> {0}.netrc'.format(homeDir + os.sep, urs), shell=True)
    Popen('echo login {} >> {}.netrc'.format(getpass(prompt=prompts[0]), homeDir + os.sep), shell=True)
    Popen('echo password {} >> {}.netrc'.format(getpass(prompt=prompts[1]), homeDir + os.sep), shell=True)

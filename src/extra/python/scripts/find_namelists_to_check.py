import subprocess
import os
import glob
from pathlib import Path
import pdb

#A script to find the fortran files within Isca's src directory
#that include namelists, and to check if namelist checking is done.

#find the location of the source code
GFDL_BASE = os.environ['GFDL_BASE']

#setup some output dictionaries and lists
fortran_file_dict = {}
includes_namelist_dict = {}
includes_check_namelist_dict = {}
n_check_namelist_dict = {}
if_def_internal_nml_dict={}

files_with_namelists = []
namelists_to_flag = []
namelists_to_flag_possible = []


#find ALL of the fortran files within GFDL_BASE/src directory
for path in Path(f'{GFDL_BASE}/src/').rglob('*.*90'):
    #exclude files with ._ at the start
    if path.name[0:2]!='._':
        #add all the remaining files to a dictionary
        fortran_file_dict[path.name] = path
    
#go through each file and check if it contains a namelist, and if it does namelist checking    
for file_name in fortran_file_dict.keys():
    file_path = fortran_file_dict[file_name]

    #initialise some of the checking variables
    namelist_in_file=False
    check_namelist_in_file=False
    number_of_checks=0
    if_def_internal_nml_in_file=False
    #open each of the fortran files
    with open(file_path, 'r') as read_obj:
        for line in read_obj:
            #check if it contains a namelist
            if 'namelist /' in line and not namelist_in_file:
                namelist_in_file=True
            # does it contain the check_nml_error command? 
            if 'check_nml_error' in line and not check_namelist_in_file:
                check_namelist_in_file=True
            # count how many times this string is mentioned
            if 'check_nml_error' in line:
                number_of_checks=number_of_checks+1
            #check if there's more than one type of namelist reading available
            if '#ifdef INTERNAL_FILE_NML' in line and not if_def_internal_nml_in_file:
                if_def_internal_nml_in_file=True                

    #make a list of those files that do have a namelist     
    if namelist_in_file:
        files_with_namelists.append(file_name)

    #make a list of those files that do have a namelist but don't do checking       
    if namelist_in_file and not check_namelist_in_file:
        namelists_to_flag.append(file_name)

    #making a list of files that have namelists, that read them in more than one way, and have fewer than 3 mentions of check_nml_error. This is to catch cases where there is some namelist checking taking place, but it's not on all the methods of namelist reading.
    if namelist_in_file and if_def_internal_nml_in_file and number_of_checks<3:
        namelists_to_flag_possible.append(file_name)

    #keep a record of the files that include a namelist
    includes_namelist_dict[file_name]=namelist_in_file

    #keep a record of the files that do and don't do namelist checking
    includes_check_namelist_dict[file_name]=check_namelist_in_file

    #keep a record of the number of checks taking place
    n_check_namelist_dict[file_name] = number_of_checks

#create a list of files that appear in namelists_to_flag_possible
list_of_filepaths_to_check = [str(fortran_file_dict[path]) for path in namelists_to_flag_possible]

#print the number of checks
print([n_check_namelist_dict[path] for path in namelists_to_flag_possible])

#print the list of files
print(namelists_to_flag_possible)

#print their directories
print(list_of_filepaths_to_check)
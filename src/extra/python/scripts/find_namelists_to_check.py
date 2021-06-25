import subprocess
import os
import glob
from pathlib import Path
import pdb

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
            if 'check_nml_error' in line:
                number_of_checks=number_of_checks+1
            if '#ifdef INTERNAL_FILE_NML' in line and not if_def_internal_nml_in_file:
                if_def_internal_nml_in_file=True                

    #make a list of those files that do have a namelist but don't do checking            
    if namelist_in_file:
        files_with_namelists.append(file_name)

    if namelist_in_file and not check_namelist_in_file:
        namelists_to_flag.append(file_name)

    if namelist_in_file and if_def_internal_nml_in_file and number_of_checks<3:
        namelists_to_flag_possible.append(file_name)

    #keep a record of the files that include a namelist
    includes_namelist_dict[file_name]=namelist_in_file
    #keep a record of the files that do and don't do namelist checking
    includes_check_namelist_dict[file_name]=check_namelist_in_file

    n_check_namelist_dict[file_name] = number_of_checks

list_of_filepaths_to_check = [str(fortran_file_dict[path]) for path in namelists_to_flag_possible]

print([n_check_namelist_dict[path] for path in namelists_to_flag_possible])

print(namelists_to_flag_possible)
print(list_of_filepaths_to_check)
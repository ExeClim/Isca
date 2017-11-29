"""
Script for comparing the results of test cases with two different commit IDs.
Purpose is to make sure that any new commits keep the results of the test cases the same,
or only change the test cases it expects to (e.g. a bug fix will change the result).

When you submit a new pull request, please run this test and report the results in the pull request.
"""
import numpy as np
from isca import Experiment, IscaCodeBase, FailedRunError, GFDL_BASE, DiagTable
from isca.util import exp_progress
import xarray as xar
import pdb
import numpy as np
import sys

def get_nml_diag(test_case_name):
   """Gets the appropriate namelist and input files from each of the test case scripts in the test_cases folder"""
   
    if 'axisymmetric' in test_case_name:
        sys.path.insert(0, GFDL_BASE+'exp/test_cases/axisymmetric/')     
        from axisymmetric_test_case import namelist as nml_out
        from axisymmetric_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles        

    if 'bucket_model' in test_case_name:
        sys.path.insert(0, GFDL_BASE+'exp/test_cases/bucket_hydrology/')     
        from bucket_model_test_case import namelist as nml_out   
        from bucket_model_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles
             
    if 'frierson' in test_case_name:
        sys.path.insert(0, GFDL_BASE+'exp/test_cases/frierson/')     
        from frierson_test_case import namelist as nml_out
        from frierson_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles                 
        
    if 'giant_planet' in test_case_name:
        sys.path.insert(0, GFDL_BASE+'exp/test_cases/giant_planet/')     
        from giant_planet_test_case import namelist as nml_out
        from giant_planet_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles           

    if 'held_suarez' in test_case_name:
        sys.path.insert(0, GFDL_BASE+'exp/test_cases/held_suarez/')
        from held_suarez_test_case import namelist as nml_out
        from held_suarez_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles           

    if 'MiMA' in test_case_name:
        sys.path.insert(0, GFDL_BASE+'exp/test_cases/MiMA/')
        from MiMA_test_case import namelist as nml_out
        from MiMA_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles                      
        
    if 'realistic_continents_fixed_sst' in test_case_name:
        sys.path.insert(0, GFDL_BASE+'exp/test_cases/realistic_continents/')
        from realistic_continents_fixed_sst_test_case import namelist as nml_out
        from realistic_continents_fixed_sst_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles           

    if 'realistic_continents_variable_qflux' in test_case_name:
        sys.path.insert(0, GFDL_BASE+'exp/test_cases/realistic_continents/')
        from realistic_continents_variable_qflux_test_case import namelist as nml_out
        from realistic_continents_variable_qflux_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles           

    if 'top_down_test' in test_case_name:
        sys.path.insert(0, GFDL_BASE+'exp/test_cases/top_down_test/')
        from top_down_test import namelist as nml_out
        from top_down_test import exp as exp_temp
        input_files = exp_temp.inputfiles           

    if 'variable_co2_grey' in test_case_name:
        sys.path.insert(0, GFDL_BASE+'exp/test_cases/variable_co2_concentration/')
        from variable_co2_grey import namelist as nml_out
        from variable_co2_grey import exp as exp_temp
        input_files = exp_temp.inputfiles           

    if 'variable_co2_rrtm' in test_case_name:
        sys.path.insert(0, GFDL_BASE+'exp/test_cases/variable_co2_concentration/')
        from variable_co2_rrtm import namelist as nml_out
        from variable_co2_rrtm import exp as exp_temp
        input_files = exp_temp.inputfiles           
                 
    return nml_out, input_files  

def define_simple_diag_table():
    """Defines a simple diag table for the test cases."""
    
    diag = DiagTable()
    diag.add_file('atmos_daily', 1, 'days', time_units='days')

    #Tell model which diagnostics to write
    diag.add_field('dynamics', 'ps', time_avg=True)
    diag.add_field('dynamics', 'bk')
    diag.add_field('dynamics', 'pk')
    diag.add_field('dynamics', 'ucomp', time_avg=True)
    diag.add_field('dynamics', 'vcomp', time_avg=True)
    diag.add_field('dynamics', 'temp', time_avg=True)
    diag.add_field('dynamics', 'vor', time_avg=True)
    diag.add_field('dynamics', 'div', time_avg=True)

    return diag

def conduct_comparison_on_test_case(base_commit, later_commit, test_case_name, repo_to_use='git@github.com:execlim/Isca', num_cores_to_use=4):
    """Process here is to checkout each commit in turn, compiles it if necessary, uses the appropriate nml for the test
    case under consideration, and runs the code with the two commits in turn. The output is then compared for all variables
    in the diag file. If there are any differences in the output variables then the test classed as a failure."""

    data_dir_dict = {}
    nml_use, input_files_use  = get_nml_diag(test_case_name)
    diag_use = define_simple_diag_table()
    test_pass = True    

    #Do the run for each of the commits in turn
    for s in [base_commit, later_commit]:
        exp_name = test_case_name+'_trip_test_'+s
        cb = IscaCodeBase(repo=repo_to_use, commit=s)
        cb.compile()
        exp = Experiment(exp_name, codebase=cb)        
        exp.namelist = nml_use.copy()
        exp.diag_table = diag_use
        exp.inputfiles = input_files_use

        #Only run for 3 days to keep things short.
        exp.update_namelist({
        'main_nml': {
        'days': 3,
        }})

        try:
            # run with a progress bar
            with exp_progress(exp, description=s) as pbar:
                exp.run(1, use_restart=False, num_cores=num_cores_to_use)

        except FailedRunError as e:
            #If run fails then test automatically fails
            test_pass = False
            continue
     
        data_dir_dict[s] = exp.datadir

    #For each of the diag files defined, compare the output
    for diag_file_entry in diag_use.files.keys():
        base_commit_dataset  = xar.open_dataset(data_dir_dict[base_commit] +'/run0001/'+diag_file_entry+'.nc', decode_times=False)
        later_commit_dataset = xar.open_dataset(data_dir_dict[later_commit]+'/run0001/'+diag_file_entry+'.nc', decode_times=False)
    
        diff = later_commit_dataset - base_commit_dataset
    
        #Check each of the output variables for differences
        for var in diff.data_vars.keys():
            maxval = np.abs(diff[var]).max()
            if maxval !=0.:
                print('Test failed for '+var+' '+str(maxval))
                test_pass = False
    
    if test_pass:    
        print('Test passed for '+test_case_name+'. Commit '+later_commit+' gives the same answer as commit '+base_commit)
        return_test_result = 'pass'
    else:
        print('Test failed for '+test_case_name+'. Commit '+later_commit+' gives a different answer to commit '+base_commit)   
        return_test_result = 'fail'
    
    return return_test_result
    
if __name__=="__main__":
    #Base commit is the earlier commit you want to compare against
    base_commit  = '155661f8c7945049cbac0dcf2019bb17fe7a6a8d'
    #later commit is the newer commit you're wanting to test
    later_commit = 'HEAD'
    
    #List of test cases to check
    exps_to_check = ['axisymmetric', 'bucket_model', 'frierson', 'giant_planet', 'held_suarez', 'MiMA', 'realistic_contients_fixed_sst', 'realistic_continents_variable_qflux', 'top_down_test', 'variable_co2_grey', 'variable_co2_rrtm']

    exp_outcome_dict = {}

    #Run the test on each test case in turn
    for exp_name in exps_to_check:
        exp_outcome_dict[exp_name] = conduct_comparison_on_test_case(base_commit, later_commit, exp_name, num_cores_to_use=4)
    
    #Decide if all tests passed or not
    overall_result = all([ k=='pass' for k in exp_outcome_dict.values() ])
    
    #Print results of each test case in turn, then overall results
    print('Results for all of the test cases ran comparing '+base_commit+' and '+later_commit+' are as follows...')
    for exp_key in exp_outcome_dict.keys():
        print(exp_key, ':', exp_outcome_dict[exp_key])
        
    if overall_result:
        print('Congratulations, all tests have passed')
    else:
        print('Nightmare, some tests have failed')
        
        
"""
Script for comparing the results of test cases with two different commit IDs.
Purpose is to make sure that any new commits keep the results of the test cases the same,
or only change the test cases it expects to (e.g. a bug fix will change the result).

When you submit a new pull request, please run this test and report the results in the pull request.
"""
import numpy as np
from isca import Experiment, IscaCodeBase, SocratesCodeBase, FailedRunError, GFDL_BASE, DiagTable
from isca.util import exp_progress
import xarray as xar
import pdb
import numpy as np
import os
import sys
import f90nml

# Mars test cases (grey_mars, radiative_eq_mars, socrates_mars, socrates_mars_dust) use dt_atmos=110s,
# chosen to divide evenly into a Mars "day" of MARS_DAY_LENGTH_SECONDS=88440s (804 steps) - a calendar
# convention (see e.g. grey_mars_test_case.py) rather than the true Martian sol (rotation_period=88308s
# in those same namelists), chosen so an integer number of days fits into a Mars year. main_nml's
# days/seconds units and DiagTable's 'days' time_units always mean Earth's fixed 86400s/day regardless
# of planet, and 86400 is not a multiple of 110 - so the generic short-test-run overrides below (which
# assume a run length in whole Earth days divides evenly by dt_atmos) need a Mars-specific equivalent.
MARS_DAY_LENGTH_SECONDS = 88440

def is_mars_test_case(test_case_name):
    return any(name in test_case_name for name in ('grey_mars', 'radiative_eq_mars', 'socrates_mars'))

def get_nml_diag(test_case_name):
    """Gets the appropriate namelist and input files from each of the test case scripts in the test_cases folder
    """

    if 'axisymmetric' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/axisymmetric/'))
        from axisymmetric_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles   
        nml_out = exp_temp.namelist    
        codebase_to_use = IscaCodeBase

    if 'bucket_model' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/bucket_hydrology/'))
        from bucket_model_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist     
        codebase_to_use = IscaCodeBase
            
    if 'frierson' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/frierson/'))
        from frierson_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist
        codebase_to_use = IscaCodeBase

    if 'frierson_dry_heating' in test_case_name:
        # Note: 'frierson_dry_heating' also matches the 'frierson' check above - this block
        # runs afterwards and overwrites nml_out/input_files/codebase_to_use, so it wins.
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/frierson_dry_heating/'))
        from frierson_dry_heating_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist
        codebase_to_use = IscaCodeBase

    if 'grey_mars' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/grey_mars/'))
        from grey_mars_test_case import exp as exp_temp
        from isca import GreyCodeBase
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist
        codebase_to_use = GreyCodeBase

    if 'radiative_eq_mars' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/radiative_eq_mars/'))
        from radiative_eq_mars_test_case import exp as exp_temp
        from isca import GreyCodeBase
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist
        codebase_to_use = GreyCodeBase

    if 'socrates_mars' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/socrates_mars/'))
        from socrates_mars_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist
        codebase_to_use = SocratesCodeBase

    if 'socrates_mars_dust' in test_case_name:
        # Note: 'socrates_mars_dust' also matches the 'socrates_mars' check above - this
        # block runs afterwards and overwrites nml_out/input_files/codebase_to_use, so it wins.
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/socrates_mars_dust/'))
        from socrates_mars_dust_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist
        codebase_to_use = SocratesCodeBase

    if 'column_test' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/column_test_case/'))
        from column_test_case import exp as exp_temp
        from isca import ColumnCodeBase
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist
        codebase_to_use = ColumnCodeBase

    if 'giant_planet' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/giant_planet/'))
        from giant_planet_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles   
        nml_out = exp_temp.namelist
        codebase_to_use = IscaCodeBase
        
        #Make giant planet test case a lower resolution so that it runs in a finite time!
        nml_out['spectral_dynamics_nml']['num_fourier']=42
        nml_out['spectral_dynamics_nml']['num_spherical']=43
        nml_out['spectral_dynamics_nml']['lon_max']=128
        nml_out['spectral_dynamics_nml']['lat_max']=64
        nml_out['spectral_dynamics_nml']['cutoff_wn']=15

    if 'held_suarez' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/held_suarez/'))
        from held_suarez_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist
        codebase_to_use = IscaCodeBase
        
    if 'MiMA' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/MiMA/'))
        from MiMA_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist     
        codebase_to_use = IscaCodeBase
        
    if 'realistic_continents_fixed_sst' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/realistic_continents/'))
        from realistic_continents_fixed_sst_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles   
        nml_out = exp_temp.namelist        
        codebase_to_use = IscaCodeBase

    if 'realistic_continents_variable_qflux' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/realistic_continents/'))
        from realistic_continents_variable_qflux_test_case import exp as exp_temp
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist        
        codebase_to_use = IscaCodeBase

    if 'soc_realistic_continents_fixed_sst_with_linear_cld_scheme' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/simple_clouds/'))
        from socrates_aquaplanet import exp as exp_temp
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist       
        codebase_to_use=SocratesCodeBase

    if 'socrates_aquaplanet' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/socrates_test/'))
        from socrates_aquaplanet import exp as exp_temp
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist
        codebase_to_use=SocratesCodeBase
        
    if 'socrates_aquaplanet_cloud' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/socrates_test/'))
        from socrates_aquaplanet_cloud import exp as exp_temp
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist     
        codebase_to_use=SocratesCodeBase

    if 'top_down_test' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/top_down_test/'))
        from top_down_test_case import namelist as nml_out
        input_files = []
        codebase_to_use = IscaCodeBase

    if 'variable_co2_grey' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/variable_co2_concentration/'))
        from variable_co2_grey import exp as exp_temp
        input_files = exp_temp.inputfiles      
        nml_out = exp_temp.namelist           
        codebase_to_use = IscaCodeBase

    if 'variable_co2_rrtm' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/variable_co2_concentration/'))
        from variable_co2_rrtm import exp as exp_temp
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist              
        codebase_to_use = IscaCodeBase
                 
    if 'ape_aquaplanet' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/ape_aquaplanet/'))
        from socrates_ape_aquaplanet_T42 import exp as exp_temp
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist
        codebase_to_use=SocratesCodeBase
        
    if 'barotropic_vort_eq_stirring' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/barotropic_vorticity_equation/'))
        from barotropic_vor_eq_stirring_test_case import exp as exp_temp
        from isca import BarotropicCodeBase
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist
        codebase_to_use=BarotropicCodeBase

    if 'shallow_water_stirring' in test_case_name:
        sys.path.insert(0, os.path.join(GFDL_BASE, 'exp/test_cases/shallow_water/'))
        from shallow_water_stirring_test_case import exp as exp_temp
        from isca import ShallowCodeBase
        input_files = exp_temp.inputfiles
        nml_out = exp_temp.namelist
        codebase_to_use=ShallowCodeBase

    return nml_out, input_files, codebase_to_use

def list_all_test_cases_implemented_in_trip_test():

    #List of test cases to check
    exps_implemented = ['axisymmetric', 
                        'bucket_model', 
                        'frierson', 
                        'giant_planet', 
                        'held_suarez', 
                        'MiMA', 
                        'realistic_continents_fixed_sst', 
                        'realistic_continents_variable_qflux', 
                        #'simple_clouds', 
                        'socrates_aquaplanet', 
                        'socrates_aquaplanet_cloud',
                        'top_down_test', 
                        'variable_co2_grey', 
                        'variable_co2_rrtm', 
                        'ape_aquaplanet',
                        'barotropic_vort_eq_stirring',
                        'shallow_water_stirring',
                        'column_test',
                        'grey_mars',
                        'radiative_eq_mars',
                        #'socrates_mars', # requires Mars-specific Socrates spectral files not yet included in the repo - see exp/test_cases/socrates_mars/input/README.md
                        #'socrates_mars_dust', # requires Mars-dust-specific Socrates spectral files and dust climatology data not yet included in the repo - see exp/test_cases/socrates_mars_dust/input/README.md
                        #'frierson_dry_heating', # exercises the new local_heating feature, which the current ExeClim master doesn't have - included for opt-in testing, not run by default
                        ]

    return exps_implemented

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

def define_simple_diag_table_mars():
    """Defines a simple diag table for the Mars test cases, with the output cadence set to
    MARS_DAY_LENGTH_SECONDS (one Mars day) rather than the generic 1 Earth day used elsewhere -
    see the MARS_DAY_LENGTH_SECONDS comment above for why."""

    diag = DiagTable()
    diag.add_file('atmos_daily', MARS_DAY_LENGTH_SECONDS, 'seconds', time_units='days')

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

def define_simple_diag_table_2d(shallow_or_baro):
    """Defines a simple diag table for the 
    shallow water and barotropic vorticity test cases."""

    if shallow_or_baro=='shallow':
        diag_name = 'shallow_diagnostics'
    elif shallow_or_baro=='barotropic':
        diag_name = 'barotropic_diagnostics'
    else:
        raise NotImplementedError('incorrect option for 2d diag table')

    diag = DiagTable()
    diag.add_file('atmos_daily', 1, 'days', time_units='days')

    #Tell model which diagnostics to write
    diag.add_field(diag_name, 'ucomp', time_avg=True)
    diag.add_field(diag_name, 'vcomp', time_avg=True)
    diag.add_field(diag_name, 'vor', time_avg=True)

    return diag

def define_simple_diag_table_column():
    """Defines a simple diag table for the column model test case. The column
    model bypasses the dynamical core, so its diagnostics live under the
    'column' module rather than 'dynamics', and it has no vorticity/divergence."""

    diag = DiagTable()
    diag.add_file('atmos_daily', 1, 'days', time_units='days')

    diag.add_field('column', 'ps', time_avg=True)
    diag.add_field('column', 'bk')
    diag.add_field('column', 'pk')
    diag.add_field('column', 'ucomp', time_avg=True)
    diag.add_field('column', 'vcomp', time_avg=True)
    diag.add_field('column', 'temp', time_avg=True)

    return diag

def process_ids(base_commit_in, later_commit_in):


    if len(base_commit_in)==40:
        #Likely to be long-hash, rather than a tag
        base_commit_short = base_commit_in[0:7]
    else:
        base_commit_short = base_commit_in

    if len(later_commit_in)==40:
        #Likely to be long-hash, rather than a tag
        later_commit_short = later_commit_in[0:7]
    else:
        later_commit_short = later_commit_in

    return base_commit_short, later_commit_short

def conduct_comparison_on_test_case(base_commit, later_commit, test_case_name, repo_to_use='git@github.com:execlim/Isca', num_cores_to_use=4):
    """Process here is to checkout each commit in turn, compiles it if necessary, uses the appropriate nml for the test
    case under consideration, and runs the code with the two commits in turn. The output is then compared for all variables
    in the diag file. If there are any differences in the output variables then the test classed as a failure."""

    data_dir_dict = {}
    nml_use, input_files_use, codebase_obj  = get_nml_diag(test_case_name)
    
    if 'shallow_water' in test_case_name:
        diag_use = define_simple_diag_table_2d('shallow')
    elif 'barotropic_vort_eq' in test_case_name:
        diag_use = define_simple_diag_table_2d('barotropic')
    elif 'column_test' in test_case_name:
        diag_use = define_simple_diag_table_column()
    elif is_mars_test_case(test_case_name):
        diag_use = define_simple_diag_table_mars()
    else:
        diag_use = define_simple_diag_table()
        
    test_pass = True
    run_complete = True
    compile_successful=True

    #Do the run for each of the commits in turn
    for s in [base_commit, later_commit]:
        exp_name = test_case_name+'_trip_test_21_'+s
        cb = codebase_obj(repo=repo_to_use, commit=s)
        try:
            cb.compile()
            exp = Experiment(exp_name, codebase=cb)
            exp.namelist = nml_use.copy()
            exp.diag_table = diag_use
            exp.inputfiles = input_files_use

            #Only run for a short time to keep things short.
            if is_mars_test_case(test_case_name):
                #Override in whole Mars days (see MARS_DAY_LENGTH_SECONDS above) rather than
                #'days', since main_nml's 'days' always means Earth's fixed 86400s regardless
                #of planet, and 86400 is not a multiple of these test cases' dt_atmos=110s.
                exp.update_namelist({
                'main_nml': {
                'days': 0,
                'seconds': 3*MARS_DAY_LENGTH_SECONDS,
                }})
            else:
                exp.update_namelist({
                'main_nml': {
                'days': 3,
                }})
        except:
            run_complete = False
            test_pass = False      
            compile_successful=False                  
            continue            


        #The column model can currently only run on 1 core, regardless of -n.
        num_cores_for_test = 1 if 'column_test' in test_case_name else num_cores_to_use

        try:
            # run with a progress bar
            with exp_progress(exp, description=s) as pbar:
                exp.run(1, use_restart=False, num_cores=num_cores_for_test)
        except FailedRunError as e:
            #If run fails then test automatically fails
            run_complete = False
            test_pass = False
            continue



        data_dir_dict[s] = exp.datadir
    if run_complete:
        #For each of the diag files defined, compare the output
        for diag_file_entry in diag_use.files.keys():
            base_commit_dataset  = xar.open_dataset(data_dir_dict[base_commit] +'/run0001/'+diag_file_entry+'.nc', decode_times=False)
            later_commit_dataset = xar.open_dataset(data_dir_dict[later_commit]+'/run0001/'+diag_file_entry+'.nc', decode_times=False)

            diff = later_commit_dataset - base_commit_dataset

            #Check each of the output variables for differences
            for var in diff.data_vars.keys():
                maxval = np.abs(diff[var]).max()
                if maxval !=0.:
                    print('Test failed for '+var+' max diff value = '+str(maxval.values))
                    test_pass = False

            base_experiment_input_nml = f90nml.read(data_dir_dict[base_commit] +'/run0001/input.nml')
            later_commit_input_nml    = f90nml.read(data_dir_dict[later_commit] +'/run0001/input.nml')

            if base_experiment_input_nml!=later_commit_input_nml:
                raise AttributeError(f'The two experiments to be compared have been run using different input namelists, and so the results may be different because of this. This only happens when you have run the trip tests using one of the commit IDs before, and that you happen to have used a different version of the test cases on that previous occasion. Try removing both {data_dir_dict[base_commit]} and {data_dir_dict[later_commit]} and try again.')

        if test_pass:
            print('Test passed for '+test_case_name+'. Commit '+later_commit+' gives the same answer as commit '+base_commit)
            return_test_result = 'pass'
        else:
            print('Test failed for '+test_case_name+'. Commit '+later_commit+' gives a different answer to commit '+base_commit)
            return_test_result = 'fail'

    else:
        if compile_successful:
            #This means that the compiles were both successful, but at least one of the runs crashed.
            print('Test failed for '+test_case_name+' because the run crashed.')
        else:
            print('Test failed for '+test_case_name+' because at least one of the runs failed to compile.')

        return_test_result = 'fail'


    return return_test_result


def output_results_function(exp_outcome_dict, base_commit, later_commit):

    base_commit_short, later_commit_short = process_ids(base_commit, later_commit)

    #Decide if all tests passed or not
    overall_result = all([ k=='pass' for k in exp_outcome_dict.values() ])

    #Print results of each test case in turn, then overall results
    print('Results for all of the test cases ran comparing '+base_commit_short+' and '+later_commit_short+' are as follows...')
    for exp_key in exp_outcome_dict.keys():
        if exp_outcome_dict[exp_key]=='pass':
            print(exp_key+' : '+'\033[1;32m'+exp_outcome_dict[exp_key]+'\033[0;m')
        else:
            print(exp_key+' : '+'\033[1;31m'+exp_outcome_dict[exp_key]+'\033[0;m')

    if overall_result:
        print('Congratulations, all tests have passed')
    else:
        print('Nightmare, some tests have failed')

def run_all_tests(base_commit, later_commit, exps_to_check, repo_to_use='git@github.com:execlim/Isca', num_cores_to_use=4):

    exp_outcome_dict = {}

    #Run the test on each test case in turn
    for exp_name in exps_to_check:
        exp_outcome_dict[exp_name] = conduct_comparison_on_test_case(base_commit, later_commit, exp_name, repo_to_use = repo_to_use, num_cores_to_use=num_cores_to_use)

    output_results_function(exp_outcome_dict, base_commit, later_commit)

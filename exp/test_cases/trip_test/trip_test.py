# Run a parameter sweep of the Held Suarez model
# by varying the rotation rate from 1% to 1000% of Earth's rot rate
import numpy as np
from isca import Experiment, DryCodeBase, FailedRunError, GFDL_BASE, DiagTable
from isca.util import exp_progress
import xarray as xar
import pdb
import numpy as np
import sys

# sys.path.insert(0, GFDL_BASE+'/test_cases/frierson/')

def get_nml_diag(test_case_name):

    if 'held_suarez' in test_case_name:
        sys.path.insert(0, GFDL_BASE+'exp/test_cases/held_suarez/')
        from held_suarez_test_case import namelist as nml_out
       
    if 'bucket_model' in test_case_name:
        sys.path.insert(0, GFDL_BASE+'exp/test_cases/bucket_hydrology/')     
        from bucket_model_test_case import namelist as nml_out
        
    return nml_out      

def define_simple_diag_table():

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

def conduct_comparison_on_test_case(base_commit, later_commit, test_case_name, repo_to_use='git@github.com:execlim/Isca'):

    data_dir_dict = {}
    nml_use  = get_nml_diag(test_case_name)
    diag_use = define_simple_diag_table()

    for s in [base_commit, later_commit]:
        exp_name = test_case_name+'_trip_test_'+s
        cb = DryCodeBase(repo=repo_to_use, commit=s)
        cb.compile()
        exp = Experiment(exp_name, codebase=cb)        
        exp.namelist = nml_use.copy()
        exp.diag_table = diag_use

        exp.update_namelist({
        'main_nml': {
        'days': 3,
        }})

        try:
            # run with a progress bar with description showing omega
            with exp_progress(exp, description=s) as pbar:
                exp.run(1, use_restart=False, num_cores=4)

        except FailedRunError as e:
            # don't let a crash get in the way of good science
            # (we could try and reduce timestep here if we wanted to be smarter)
            continue
     
        data_dir_dict[s] = exp.datadir

    for diag_file_entry in diag_use.files.keys():
        base_commit_dataset  = xar.open_dataset(data_dir_dict[base_commit] +'/run0001/'+diag_file_entry+'.nc', decode_times=False)
        later_commit_dataset = xar.open_dataset(data_dir_dict[later_commit]+'/run0001/'+diag_file_entry+'.nc', decode_times=False)
    
        diff = later_commit_dataset - base_commit_dataset
    
        for var in diff.data_vars.keys():
            maxval = np.abs(diff[var]).max()
            if maxval !=0.:
                raise RuntimeError('Test failed for '+var+' '+str(maxval))
            
        print('Test passed for '+test_case_name+'. Commit '+later_commit+' gives the same answers as commit '+base_commit)
    
    
if __name__=="__main__":
    base_commit  = '155661f8c7945049cbac0dcf2019bb17fe7a6a8d'
    later_commit = 'HEAD'
    
    conduct_comparison_on_test_case(base_commit, later_commit, 'held_suarez_test_case')
    
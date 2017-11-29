# Run a parameter sweep of the Held Suarez model
# by varying the rotation rate from 1% to 1000% of Earth's rot rate
import numpy as np
from isca import Experiment, DryCodeBase, FailedRunError, GFDL_BASE
from isca.util import exp_progress
import xarray as xar
import pdb
import numpy as np

# import sys
# sys.path.insert(0, GFDL_BASE+'/test_cases/held_suarez/')

from held_suarez_test_case import namelist as hs_nml
from held_suarez_test_case import diag as  hs_diag

hs_nml['main_nml'] = {
        'dt_atmos': 600,
        'days': 3,
        'calendar': 'no_calendar'
}

base_commit  = '155661f8c7945049cbac0dcf2019bb17fe7a6a8d'
later_commit = 'HEAD'

data_dir_dict = {}


for s in [base_commit, later_commit]:
    exp_name = 'hs_trip_test_'+s
    cb = DryCodeBase(repo='git@github.com:execlim/Isca', commit=s)
    cb.compile()
    exp = Experiment(exp_name, codebase=cb)
    exp.namelist = hs_nml.copy()
    exp.diag_table = hs_diag

    try:
        # run with a progress bar with description showing omega
        with exp_progress(exp, description=s) as pbar:
            exp.run(1, use_restart=False, num_cores=16)

    except FailedRunError as e:
        # don't let a crash get in the way of good science
        # (we could try and reduce timestep here if we wanted to be smarter)
        continue
     
    data_dir_dict[s] = exp.datadir

for diag_file_entry in hs_diag.files.keys():
    base_commit_dataset  = xar.open_dataset(data_dir_dict[base_commit] +'/run0001/'+diag_file_entry+'.nc', decode_times=False)
    later_commit_dataset = xar.open_dataset(data_dir_dict[later_commit]+'/run0001/'+diag_file_entry+'.nc', decode_times=False)
    
    diff = later_commit_dataset - base_commit_dataset
    
    for var in diff.data_vars.keys():
        maxval = np.abs(diff[var]).max()
        if maxval !=0.:
            raise RuntimeError('Test failed for '+var+' '+str(maxval))
            
    print('Test passed for Held-Suarez. Commit '+later_commit+' gives the same answers as commit '+base_commit)
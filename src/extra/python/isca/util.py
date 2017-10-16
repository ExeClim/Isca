from contextlib import contextmanager
import json

import numpy as np
from tqdm import tqdm
import sh

@contextmanager
def exp_progress(exp, description='DAY {day}'):
    """Create a progress bar on the terminal output.

    DAY 39:  39%|====           | 39/100 [00:32<00:46,  1.31it/s, avgT=261, spd=40.6]

    Currently only works with a namelist with no_calendar and that
    specifies the run length in days.

    Use as a context manager.  e.g.
    with exp_progess(exp) as pbar:
        exp.run(1)
        # pbar.set_description('...')

    """
    exp.update_namelist({'spectral_dynamics_nml': {'json_logging': True}})
    with tqdm(total=exp.namelist['main_nml']['days']) as pbar:
        # handle the run output
        @exp.on('run:output')
        def parse_output(exp, line):
            try:
                data = json.loads(line)
                pbar.update(1)
                pbar.set_description(description.format(**data))
                pbar.set_postfix(spd=data['max_speed'], avgT=data['avg_T'])
            except ValueError as e:
                # wasn't valid JSON, just output as normal
                exp.log_output(line)

        yield pbar

        # after yield, we clean up.
        # we're done with logging so remove the temporary handler
        exp._events['run:output'].remove(parse_output)

def keep_only_certain_restart_files(exp, max_num_files, interval=12):
    try:
    #       sh.ls(sh.glob(P(exp.workdir,'restarts','res_*.cpio'))) #TODO get max_num_files calculated in line, rather than a variable to pass.

            #First defines a list of ALL the restart file numbers
        files_to_remove=range(0,max_num_files)

            #Then defines a list of the ones we want to KEEP
        files_to_keep  =range(0,max_num_files,interval)

            #Then we remove the files we want to keep from the list of all files, giving a list of those we wish to remove
        for x in files_to_keep:
               files_to_remove.remove(x)

            #Then we remove them.
        for entry in files_to_remove:
                sh.rm(P(exp.workdir,'restarts','res_'+str(entry)+'.cpio'))

    except sh.ErrorReturnCode_1:
        log.warning('Tried to remove some restart files, but the last one doesnt exist')


def runinterp(self, month, infile, outfile, var_names = '-a', p_levs = "EVEN", rm_input=False):
    """Interpolate data from sigma to pressure levels. Includes option to remove original file."""
    import subprocess
    pprocess = P(GFDL_BASE,'postprocessing/plevel_interpolation/scripts')
    interper = 'source '+pprocess+'/plevel.sh -i '
    inputfile = P(self.datadir, 'run%04d' % month, infile)
    outputfile = P(self.datadir, 'run%04d' % month, outfile)

    # Select from pre-chosen pressure levels, or input new ones in hPa in the format below.
    if p_levs == "MODEL":
        plev = ' -p "2 9 18 38 71 125 206 319 471 665 904 1193 1532 1925 2375 2886 3464 4115 4850 5679 6615 7675 8877 10244 11801 13577 15607 17928 20585 23630 27119 31121 35711 40976 47016 53946 61898 71022 81491 93503" '
    elif p_levs == "EVEN":
        plev = ' -p "100000 95000 90000 85000 80000 75000 70000 65000 60000 55000 50000 45000 40000 35000 30000 25000 20000 15000 10000 5000" '
    else:
        plev = p_levs
    command = interper + inputfile + ' -o ' + outputfile + plev + var_names
    subprocess.call([command], shell=True)
    if rm_input:
        sh.rm( inputfile)
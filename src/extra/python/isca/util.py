from contextlib import contextmanager
import json
from os.path import join as P

import numpy as np
from tqdm import tqdm
import xarray as xr
import sh

from isca import GFDL_BASE
from isca.create_alert import disk_space_alert

@contextmanager
def exp_progress(exp, description='DAY {n}'):
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
    main = exp.namelist['main_nml']
    total_days = int(main.get('seconds', 0) / 86400.0
        + main.get('days', 0)
        + main.get('months', 0)*30
        + main.get('years', 0)*360)
    with tqdm(total=total_days) as pbar:
        # handle the run output
        @exp.on('run:output')
        def parse_output(exp, line):
            try:
                data = json.loads(line)
            except ValueError as e:
                # wasn't valid JSON, just output as normal
                exp.log_output(line)
            else:
                pbar.update(1)
                data.update({'n': pbar.n})
                pbar.set_description(description.format(**data))
                if data.get('max_speed'):
                    pbar.set_postfix(spd=data['max_speed'], avgT=data['avg_T'])

        yield pbar

        # after yield, we clean up.
        # we're done with logging so remove the temporary handler
        exp._events['run:output'].remove(parse_output)


@contextmanager
def email_alerts(exp, email_address, limit=2000, cutoff=5):
    """A context manager for email alerts.
    e.g.

    with email_alerts(exp, 'myemail@example.com'):
        ...
        exp.run(...)
    """
    # add a handler to ready events
    @exp.on('run:ready')
    def check_disk_space(exp, month):
        dir = exp.datadir
        disk_space_alert(exp.datadir, exp.name, month, email_address, limit, cutoff)
    yield

    exp._events.remove(check_disk_space)

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

def clean_datadir(exp, run, keep_files=['input.nml', 'diag_table', 'field_table', 'git_hash_used.txt']):
    """Remove the `run` directory from output data, retaining only small
    configuration files."""
    outdir = exp.get_outputdir(run)
    for file in keep_files:
        filepath = P(outdir, 'run', file)
        if os.path.isfile(filepath):
            sh.cp(filepath, P(outdir, file))
            exp.log.info('Copied %s to %s' % (file, outdir))
    sh.rm('-r', P(outdir, 'run'))
    exp.log.info('Deleted %s directory' % P(outdir, 'run'))

def delete_all_restarts(exp, exceptions=None):
    """Remove the restart files for a given experiment except those given.

    e.g. remove_restarts(exp, [3,6,9,12])"""
    if exceptions:
        exceptions = [exp.restartfmt % i for i in exceptions]
    all_restarts = os.listdir(exp.restartdir)
    restarts_to_remove = [file for file in all_restarts if file not in exceptions]
    for file in restarts_to_remove:
        sh.rm(P(exp.restartdir, file))
        exp.log.info('Deleted restart file %s' % file)




def interpolate_output(infile, outfile, var_names=None, p_levs = "input"):
    """Interpolate data from sigma to pressure levels. Includes option to remove original file.

    This is a very thin wrapper around the plevel.sh script found in
    `postprocessing/plevel_interpolation/scripts/plevel.sh`.  Read the documentation
    in that script for more information.

    The interpolator must also be compiled before use.  See `postprocessing/plevel_interpolation/README`
    for instructions.

    infile: The path of a netcdf file to interpolate over.
    outfile: The path to save the output to.
    var_names: a list of fields to interpolate. if None, process all the fields in the input.
    p_levs: The list of pressure values, in pascals, to interpolate onto. Can be:
        * A list of integer pascal values
        * "input": Interpolate onto the pfull values in the input file
        * "even": Interpolate onto evenly spaced in Pa levels.
    Outputs to outfile.
    """
    interpolator = sh.Command(P(GFDL_BASE, 'postprocessing', 'plevel_interpolation', 'scripts', 'plevel.sh'))

    # Select from pre-chosen pressure levels, or input new ones in hPa in the format below.
    if isinstance(p_levs, str):
        if p_levs.upper() == "INPUT":
            with xr.open_dataset(infile, decode_times=False) as dat:
                levels = dat.pfull.data * 100
        elif p_levs.upper() == "MODEL":
        #     plev = ' -p "2 9 18 38 71 125 206 319 471 665 904 1193 1532 1925 2375 2886 3464 4115 4850 5679 6615 7675 8877 10244 11801 13577 15607 17928 20585 23630 27119 31121 35711 40976 47016 53946 61898 71022 81491 93503" '
            levels = [2, 9, 18, 38, 71, 125, 206, 319, 471, 665, 904, 1193, 1532, 1925, 2375, 2886, 3464, 4115, 4850, 5679, 6615, 7675, 8877, 10244, 11801, 13577, 15607, 17928, 20585, 23630, 27119, 31121, 35711, 40976, 47016, 53946, 61898, 71022, 81491, 93503]
        elif p_levs.upper() == "EVEN":
            #plev = ' -p "100000 95000 90000 85000 80000 75000 70000 65000 60000 55000 50000 45000 40000 35000 30000 25000 20000 15000 10000 5000" '
            levels = [100000, 95000, 90000, 85000, 80000, 75000, 70000, 65000, 60000, 55000, 50000, 45000, 40000, 35000, 30000, 25000, 20000, 15000, 10000, 5000]
    else:
        levels = p_levs

    plev = " ".join("{:.0f}".format(x) for x in reversed(sorted(levels)))
    if var_names is None:
        var_names = '-a'
    else:
        var_names = ' '.join(var_names)

    interpolator('-i', infile, '-o', outfile, '-p', plev, var_names)

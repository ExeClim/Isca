import os
import re

from f90nml import Namelist
from jinja2 import Environment, FileSystemLoader
import glob
import sh   
import pdb
import tarfile

# from gfdl import create_alert
# import getpass

from isca import GFDL_WORK, GFDL_DATA, GFDL_BASE, _module_directory, get_env_file, EventEmitter
from isca.diagtable import DiagTable
from isca.loghandler import Logger, clean_log_debug
from isca.helpers import destructive, useworkdir, mkdir

P = os.path.join

class CompilationError(Exception):
    pass

class FailedRunError(Exception): pass

class Experiment(Logger, EventEmitter):
    """A basic GFDL experiment"""

    RESOLUTIONS = {
        'T170': {
            'lon_max': 512,
            'lat_max': 256,
            'num_fourier': 170,
            'num_spherical': 171
        },

        'T85': {
            'lon_max': 256,
            'lat_max': 128,
            'num_fourier': 85,
            'num_spherical': 86
        },

        'T42': {
            'lon_max': 128,
            'lat_max': 64,
            'num_fourier': 42,
            'num_spherical': 43,
        },

        'T21': {
            'lon_max': 64,
            'lat_max': 32,
            'num_fourier': 21,
            'num_spherical': 22,
        },
    }

    runfmt = 'run%04d'
    restartfmt = 'res%04d.tar.gz'

    def __init__(self, name, codebase, safe_mode=False, workbase=GFDL_WORK, database=GFDL_DATA):
        super(Experiment, self).__init__()
        self.name = name
        self.codebase = codebase
        self.safe_mode = safe_mode

        # set the default locations of working directory,
        # executable directory, restart file storage, and
        # output data directory.
        self.workdir = P(workbase, 'experiment', self.name)
        self.rundir = P(self.workdir, 'run')          # temporary area an individual run will be performed
        self.datadir = P(database, self.name)        # where run data will be moved to upon completion
        self.restartdir = P(self.datadir, 'restarts') # where restarts will be stored
        self.template_dir = P(_module_directory, 'templates')

        self.env_source = get_env_file()

        self.templates = Environment(loader=FileSystemLoader(self.template_dir))

        self.diag_table = DiagTable()
        self.field_table_file = P(self.codebase.srcdir, 'extra', 'model', self.codebase.name, 'field_table')
        self.inputfiles = []

        self.namelist = Namelist()

    @destructive
    def rm_workdir(self):
        try:
            sh.rm(['-r', self.workdir])
        except sh.ErrorReturnCode:
            self.log.warning('Tried to remove working directory but it doesnt exist')

    @destructive
    def rm_datadir(self):
        try:
            sh.rm(['-r', self.datadir])
        except sh.ErrorReturnCode:
            self.log.warning('Tried to remove data directory but it doesnt exist')

    @destructive
    @useworkdir
    def clear_workdir(self):
        self.rm_workdir()
        mkdir(self.workdir)
        self.log.info('Emptied working directory %r' % self.workdir)

    @destructive
    @useworkdir
    def clear_rundir(self):
        #sh.cd(self.workdir)
        try:
            sh.rm(['-r', self.rundir])
        except sh.ErrorReturnCode:
            self.log.warning('Tried to remove run directory but it doesnt exist')
        mkdir(self.rundir)
        self.log.info('Emptied run directory %r' % self.rundir)

    def get_restart_file(self, i):
        return P(self.restartdir, self.restartfmt % i)

    def get_outputdir(self, run):
        return P(self.datadir, self.runfmt % run)

    def set_resolution(self, res, num_levels=None):
        """Set the resolution of the model, based on the standard triangular
        truncations of the spectral core.  For example,
            exp.set_resolution('T85', 25)
        creates a spectral core with enough modes to natively correspond to
        a 256x128 lon-lat resolution."""
        delta = self.RESOLUTIONS[res]
        if num_levels is not None:
            delta['num_levels'] = num_levels
        self.update_namelist({'spectral_dynamics_nml': delta})

    def update_namelist(self, new_vals):
        """Update the namelist sections, overwriting existing values."""
        for sec in new_vals:
            if sec not in self.namelist:
                self.namelist[sec] = {}
            nml = self.namelist[sec]
            nml.update(new_vals[sec])

    def write_namelist(self, outdir):
        namelist_file = P(outdir, 'input.nml')
        self.log.info('Writing namelist to %r' % namelist_file)
        self.namelist.write(namelist_file)

    def write_diag_table(self, outdir):
        outfile = P(outdir, 'diag_table')
        self.log.info('Writing diag_table to %r' % outfile)
        if self.diag_table.is_valid():
            if self.diag_table.calendar is None:
                # diagnose the calendar from the namelist
                cal = self.get_calendar()
                self.diag_table.calendar = cal
            self.diag_table.write(outfile)
        else:
            self.log.error("No output files defined in the DiagTable. Stopping.")
            raise ValueError()

    def write_field_table(self, outdir):
        self.log.info('Writing field_table to %r' % P(outdir, 'field_table'))
        sh.cp(self.field_table_file, P(outdir, 'field_table'))

    def log_output(self, outputstring):
        line = outputstring.strip()
        if 'warning' in line.lower():
            self.log.warn(line)
        else:
            self.log.debug(line)
        #return clean_log_debug(outputstring)

    def delete_restart(self, run):
        resfile = self.get_restart_file(run)
        if os.path.isfile(resfile):
            sh.rm(resfile)
            self.log.info('Deleted restart file %s' % resfile)

    def get_calendar(self):
        """Get the value of 'main_nml/calendar.
        Returns a string name of calendar, or None if not set in namelist.'"""
        if 'main_nml' in self.namelist:
            return self.namelist['main_nml'].get('calendar')
        else:
            return None

    def check_for_existing_output(self, i):
        outdir = P(self.datadir, self.runfmt % i)
        return os.path.isdir(outdir)

    @destructive
    @useworkdir
    def run(self, i, restart_file=None, use_restart=True, multi_node=False, num_cores=8, overwrite_data=False, save_run=False, run_idb=False, nice_score=0, mpirun_opts=''):
        """Run the model.0
            `num_cores`: Number of mpi cores to distribute over.
            `restart_file` (optional): A path to a valid restart archive.  If None and `use_restart=True`,
                                       restart file (i-1) will be used.
            `save_run`:  If True, copy the entire working directory over to GFDL_DATA
                         so that the run can rerun without the python script.
                         (This uses a lot of data storage!)

        """

        self.clear_rundir()

        indir =  P(self.rundir, 'INPUT')
        outdir = P(self.datadir, self.runfmt % i)
        resdir = P(self.rundir, 'RESTART')

        if self.check_for_existing_output(i):
            if overwrite_data:
                self.log.warning('Data for run %d already exists and overwrite_data is True. Overwriting.' % i)
                sh.rm('-r', outdir)
            else:
                self.log.warn('Data for run %d already exists but overwrite_data is False. Stopping.' % i)
                return False

        # make the output run folder and copy over the input files
        mkdir([indir, resdir, self.restartdir])

        self.codebase.write_source_control_status(P(self.rundir, 'git_hash_used.txt'))
        self.write_namelist(self.rundir)
        self.write_field_table(self.rundir)
        self.write_diag_table(self.rundir)

        for filename in self.inputfiles:
            sh.cp([filename, P(indir, os.path.split(filename)[1])])

        if multi_node:
            mpirun_opts += ' -bootstrap pbsdsh -f $PBS_NODEFILE'

        if use_restart and not restart_file and i == 1:
            # no restart file specified, but we are at first run number
            self.log.warn('use_restart=True, but restart_file not specified.  As this is run 1, assuming spin-up from namelist stated initial conditions so continuing.')
            use_restart = False

        if use_restart:
            if not restart_file:
                # get the restart from previous iteration
                restart_file = self.get_restart_file(i - 1)
            if not os.path.isfile(restart_file):
                self.log.error('Restart file not found, expecting file %r' % restart_file)
                raise IOError('Restart file not found, expecting file %r' % restart_file)
            else:
                self.log.info('Using restart file %r' % restart_file)

            self.extract_restart_archive(restart_file, indir)
        else:
            self.log.info('Running without restart file')
            restart_file = None

        vars = {
            'rundir': self.rundir,
            'execdir': self.codebase.builddir,
            'executable': self.codebase.executable_name,
            'env_source': self.env_source,
            'mpirun_opts': mpirun_opts,
            'num_cores': num_cores,
            'run_idb': run_idb,
            'nice_score': nice_score
        }

        runscript = self.templates.get_template('run.sh')

        # employ the template to create a runscript
        t = runscript.stream(**vars).dump(P(self.rundir, 'run.sh'))

        def _outhandler(line):
            handled = self.emit('run:output', self, line)
            if not handled: # only log the output when no event handler is used
                self.log_output(line)

        self.emit('run:ready', self, i)
        self.log.info("Beginning run %d" % i)
        try:
            #for line in sh.bash(P(self.rundir, 'run.sh'), _iter=True, _err_to_out=True):
            proc = sh.bash(P(self.rundir, 'run.sh'), _bg=True, _out=_outhandler, _err_to_out=True)
            self.log.info('process running as {}'.format(proc.process.pid))
            proc.wait()
            completed = True
        except KeyboardInterrupt as e:
            self.log.error("Manual interrupt, killing process.")
            proc.process.terminate()
            proc.wait()
            #log.info("Cleaning run directory.")
            #self.clear_rundir()
            raise e
        except sh.ErrorReturnCode as e:
            completed = False
            self.log.error("Run %d failed. See log for details." % i)
            self.log.error("Error: %r" % e)
            self.emit('run:failed', self)
            raise FailedRunError()

        self.emit('run:complete', self, i)
        self.log.info('Run %d complete' % i)
        mkdir(outdir)

        if num_cores > 1:
            # use postprocessing tool to combine the output from several cores
            codebase_combine_script = P(self.codebase.builddir, 'mppnccombine_run.sh')
            if not os.path.exists(codebase_combine_script):
                self.log.warning('combine script does not exist in the commit you are running Isca from.  Falling back to using $GFDL_BASE mppnccombine_run.sh script')
                sh.ln('-s',  P(GFDL_BASE, 'postprocessing', 'mppnccombine_run.sh'), codebase_combine_script)
            combinetool = sh.Command(codebase_combine_script)
            for file in self.diag_table.files:
                netcdf_file = '%s.nc' % file
                filebase = P(self.rundir, netcdf_file)
                combinetool(self.codebase.builddir, filebase)
                # copy the combined netcdf file into the data archive directory
                sh.cp(filebase, P(outdir, netcdf_file))
                # remove all netcdf fragments from the run directory
                sh.rm(glob.glob(filebase+'*'))
                self.log.debug('%s combined and copied to data directory' % netcdf_file)

            for restart in glob.glob(P(resdir, '*.res.nc.0000')):
                restartfile = restart.replace('.0000', '')
                combinetool(self.codebase.builddir, restartfile)
                sh.rm(glob.glob(restartfile+'.????'))
                self.log.debug("Restart file %s combined" % restartfile)

            self.emit('run:combined', self, i)
        else:
            for file in self.diag_table.files:
                netcdf_file = '%s.nc' % file
                filebase = P(self.rundir, netcdf_file)
                sh.cp(filebase, P(outdir, netcdf_file))
                sh.rm(glob.glob(filebase+'*'))
                self.log.debug('%s copied to data directory' % netcdf_file)

        # make the restart archive and delete the restart files
        self.make_restart_archive(self.get_restart_file(i), resdir)
        sh.rm('-r', resdir)

        if save_run:
            # copy the complete run directory to GFDL_DATA so that the run can
            # be recreated without the python script if required
            mkdir(resdir)
            sh.cp(['-a', self.rundir, outdir])
        else:
            # just save some useful diagnostic information
            self.write_namelist(outdir)
            self.write_field_table(outdir)
            self.write_diag_table(outdir)
            self.codebase.write_source_control_status(P(outdir, 'git_hash_used.txt'))

        self.clear_rundir()
        self.emit('run:finished', self, i)
        return True

    def make_restart_archive(self, archive_file, restart_directory):
        with tarfile.open(archive_file, 'w:gz') as tar:
            tar.add(restart_directory, arcname='.')
        self.log.info("Restart archive created at %s" % archive_file)

    def extract_restart_archive(self, archive_file, input_directory):
        with tarfile.open(archive_file, 'r:gz') as tar:
            tar.extractall(path=input_directory)
        self.log.info("Restart %s extracted to %s" % (archive_file, input_directory))

    def derive(self, new_experiment_name):
        """Derive a new experiment based on this one."""
        new_exp = Experiment(new_experiment_name, self.codebase)
        new_exp.namelist = self.namelist.copy()
        new_exp.diag_table = self.diag_table.copy()
        new_exp.inputfiles = self.inputfiles[:]

        return new_exp

    # TODO: replace this with util functionality
    # def run_parameter_sweep(self, parameter_values, runs=10, num_cores=16):
    #     # parameter_values should be a namelist fragment, with multiple values
    #     # for each study e.g. to vary obliquity:
    #     # exp.run_parameter_sweep({'astronomy_nml': {'obliq': [0.0, 5.0, 10.0, 15.0]}})
    #     # will run 4 independent studies and create data e.g.
    #     # <exp_name>/astronomy_nml_obliq_<0.0, 5.0 ...>/run[1-10]/daily.nc
    #     params = [(sec, name, values) for sec, parameters in parameter_values.items()
    #                     for name, values in parameters.items()]
    #     # make a list of lists of namelist section, parameter names and values
    #     params = [[(a,b,val) for val in values] for a,b,values in params]
    #     parameter_space = itertools.product(params)
    #     for combo in parameter_space:
    #         title = '_'.join(['%s_%s_%r' % (sec[:3], name[:5], val) for sec, name, val in combo])
    #         exp = self.derive(self.name + '_' + title)
    #         for sec, name, val in combo:
    #             exp.namelist[sec][name] = val
    #         exp.clear_rundir()
    #         exp.run(1, use_restart=False, num_cores=num_cores)
    #         for i in range(runs-1):
    #             exp.run(i+2)

# class RunSpec(Logger):
#     def __init__(self, exp):
#         self.exp = exp


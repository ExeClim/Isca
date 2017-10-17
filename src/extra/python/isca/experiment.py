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

from isca import GFDL_WORK, GFDL_DATA, _module_directory, get_env_file, EventEmitter
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
    }

    def __init__(self, name, codebase, safe_mode=False, workbase=GFDL_WORK, database=GFDL_DATA):
        super(Experiment, self).__init__()
        self.name = name
        self.codebase = codebase
        self.safe_mode = safe_mode

        # set the default locations of working directory,
        # executable directory, restart file storage, and
        # output data directory.
        self.workdir = P(workbase, 'experiment', self.name)
        self.restartdir = P(self.workdir, 'restarts') # where restarts will be stored
        self.rundir = P(self.workdir, 'run')          # temporary area an individual run will be performed
        self.datadir = P(database, self.name)        # where run data will be moved to upon completion
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
        sh.cd(self.workdir)
        try:
            sh.rm(['-r', self.rundir])
        except sh.ErrorReturnCode:
            self.log.warning('Tried to remove run directory but it doesnt exist')
        mkdir(self.rundir)
        self.log.info('Emptied run directory %r' % self.rundir)

    def get_restart_file(self, i):
        return P(self.restartdir, 'res%04d.tar.gz' % i)

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
        if len(self.diag_table.files):
            template = self.templates.get_template('diag_table')
            calendar = not self.namelist['main_nml']['calendar'].lower().startswith('no_calendar')
            vars = {'calendar': calendar, 'outputfiles': self.diag_table.files.values()}
            template.stream(**vars).dump(outfile)
        else:
            self.log.error("No output files defined in the DiagTable. Stopping.")
            raise ValueError()

    def write_field_table(self, outdir):
        self.log.info('Writing field_table to %r' % P(outdir, 'field_table'))
        sh.cp(self.field_table_file, P(outdir, 'field_table'))

    _day_re = re.compile('Integration completed through\s+([0-9]+) days')
    _month_re = re.compile('Integration completed through\s+([\w\s]+)\s([0-9]+:)')
    def log_output(self, outputstring):
        line = outputstring.strip()
        if 'warning' in line.lower():
            self.log.warn(line)
        else:
            self.log.debug(line)
        #return clean_log_debug(outputstring)

    @destructive
    @useworkdir
    def run(self, i, restart_file=None, use_restart=True, num_cores=8, overwrite_data=False, light=False, run_idb=False, experiment_restart=None, email_alerts=True, email_address_for_alerts=None, disk_space_limit=20, disk_space_cutoff_limit=5):

        indir =  P(self.rundir, 'INPUT')
        outdir = P(self.datadir, 'run%04d' % i)
        resdir = P(self.rundir, 'RESTART')

        self.codebase.write_source_control_status(P(self.rundir, 'git_hash_used.txt'))

        if os.path.isdir(outdir):
            if overwrite_data:
                self.log.warning('Data for run %d already exists and overwrite_data is True. Overwriting.' % i)
                sh.rm('-r', outdir)
            else:
                self.log.warn('Data for run %d already exists but overwrite_data is False. Stopping.' % i)
                return False

        # make the output run folder and copy over the input files
        mkdir([indir, resdir, self.restartdir])

        self.write_namelist(self.rundir)
        self.write_field_table(self.rundir)
        self.write_diag_table(self.rundir)

        for filename in self.inputfiles:
            sh.cp([filename, P(indir, os.path.split(filename)[1])])

        if use_restart:
            if not restart_file:
                # get the restart from previous iteration
                restart_file = self.get_restart_file(i - 1)
            if not os.path.isfile(restart_file):
                self.log.error('Restart file not found, expecting file %r' % restart_file)
                exit(2)
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
            'num_cores': num_cores,
            'run_idb': run_idb,
        }

        runscript = self.templates.get_template('run.sh')

        # employ the template to create a runscript
        t = runscript.stream(**vars).dump(P(self.rundir, 'run.sh'))

    # Check scratch space has enough disk space
        #if email_alerts:
        #    if email_address_for_alerts is None:
        #        email_address_for_alerts = getpass.getuser()+'@exeter.ac.uk'
        #    create_alert.run_alerts(self.execdir, GFDL_BASE, self.name, month, email_address_for_alerts, disk_space_limit)

        def _outhandler(line):
            handled = self.emit('run:output', self, line)
            if not handled: # only log the output when no event handler is used
                self.log_output(line)

        self.emit('run:ready', self)
        self.log.info("Beginning run %d" % i)
        try:
            #for line in sh.bash(P(self.rundir, 'run.sh'), _iter=True, _err_to_out=True):
            proc = sh.bash(P(self.rundir, 'run.sh'), _bg=True, _out=_outhandler, _err_to_out=True)
            proc.wait()
            completed = True
        except KeyboardInterrupt as e:
            self.log.error("Manual interrupt, killing process.")
            proc.process.kill()
            #log.info("Cleaning run directory.")
            #self.clear_rundir()
            raise e
        except sh.ErrorReturnCode as e:
            completed = False
            self.log.error("Run %d failed. See log for details." % i)
            self.log.error("Error: %r" % e)
            self.emit('run:failed', self)
            raise FailedRunError()

        self.emit('run:completed', self, i)
        self.log.info('Run %d complete' % i)
        mkdir(outdir)

        if num_cores > 1:
            combinetool = sh.Command(P(self.codebase.builddir, 'mppnccombine.x'))
            # use postprocessing tool to combine the output from several cores
            for file in self.diag_table.files:
                netcdf_file = '%s.nc' % file
                filebase = P(self.rundir, netcdf_file)
                combinetool(filebase)
                # copy the combined netcdf file into the data archive directory
                sh.cp(filebase, P(outdir, netcdf_file))
                # remove all netcdf fragments from the run directory
                sh.rm(glob.glob(filebase+'*'))
                self.log.debug('%s combined and copied to data directory' % netcdf_file)

            for restart in glob.glob(P(resdir, '*.res.nc.0000')):
                restartfile = restart.replace('.0000', '')
                combinetool(restartfile)
                sh.rm(glob.glob(restartfile+'.????'))
                self.log.debug("Restart file %s combined" % restartfile)

            self.emit('run:combined', self)

        # make the restart archive and delete the restart files
        self.make_restart_archive(self.get_restart_file(i), resdir)
        sh.rm('-r', resdir)

        sh.cp(['-a', self.rundir, outdir])
        self.clear_rundir()



        # TODO: replace this with util function
        # if light:
        #     os.system("cp -a "+self.rundir+"/*.nc "+outdir)
        #     sh.cp(['-a', P(self.restartdir, 'res_%d.cpio' % (month)), outdir])
        #     if month > 1:
        #         try:
        #             sh.rm( P(self.restartdir, 'res_%d.cpio' % (month-1)))
        #         except sh.ErrorReturnCode:
        #             log.warning('Previous months restart already removed')

        #         try:
        #             sh.rm( P(self.datadir, 'run03%d' % (month-1) , 'res_%d.cpio' % (month-1)))
        #         except sh.ErrorReturnCode:
        #             log.warning('Previous months restart already removed')

        # else:
        #     sh.cp(['-a', self.rundir+'/.', outdir])
        # self.clear_rundir()
        # sh.cd(self.rundir)
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
        new_exp.inputfiles = self.inputfiles.copy()
        # TODO: fix this
        # new_exp.commit_id = self.commit_id
        # new_exp.commit_id_base = self.commit_id_base
        # new_exp.git_status_output = self.git_status_output
        # new_exp.git_diff_output   = self.git_diff_output
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


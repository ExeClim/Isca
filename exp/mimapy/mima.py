#!/usr/bin/env python

import logging
import os

import f90nml
from jinja2 import Environment, FileSystemLoader
import sh

P = os.path.join

mkdir = sh.mkdir.bake('-p')



log = logging.getLogger('mima')
log.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
log.addHandler(ch)

def clean_log_info(s):
    if s.strip():
       log.info(s.strip())


try:
    GFDL_BASE = os.environ['GFDL_BASE']
    GFDL_WORK = os.environ['GFDL_WORK']
    GFDL_DATA = os.environ['GFDL_DATA']
except Exception, e:
    print('Environment variables GFDL_BASE, GFDL_WORK, GFDL_DATA must be set')
    exit(0)


class Experiment(object):
    """A basic MiMA experiment"""
    def __init__(self, name, commit=None, repo=None, overwrite_data=False):
        super(Experiment, self).__init__()
        self.name = name


        # set the default locations of working directory,
        # executable directory, restart file storage, and
        # output data directory.
        # These can be overridden e.g. if an experiment is to use
        # the same executable as another
        self.workdir = P(GFDL_WORK, self.name)
        self.execdir = P(self.workdir, 'exec')        # where executable will be compiled to / fectched from
        self.restartdir = P(self.workdir, 'restarts') # where restarts will be stored
        self.rundir = P(self.workdir, 'run')          # temporary area an individual run will be performed
        self.datadir = P(GFDL_DATA, self.name)        # where run data will be moved to upon completion


        if os.path.isdir(self.workdir):
            log.warning('Working directory for exp %r already exists' % self.name)
        else:
            log.debug('Making directory %r' % self.workdir)
            mkdir(self.workdir)

        # Running an experiment from a git repo
        self.repo = repo
        self.commit = commit
        if repo and commit:
            # do a checkout of the specific code and work from that source tree
            self.srcdir = P(self.workdir, 'source')
            self.clone_and_checkout()
        else:
            self.srcdir = GFDL_BASE

        self.mimapy_dir = P(self.srcdir, 'exp', 'mimapy')

        self.templates = Environment(loader=FileSystemLoader(P(self.mimapy_dir, 'templates')))

        self.diag_table_file = P(self.mimapy_dir, 'diag_table')
        self.field_table_file = P(self.mimapy_dir, 'field_table')

        self.path_names_file = P(self.mimapy_dir, 'path_names')
        self.path_names = self._get_default_path_names()

        self.namelist_files = [P(self.mimapy_dir, 'core.nml'), P(self.mimapy_dir, 'phys.nml')]
        self.namelist = self.rebuild_namelist()

        self.inputfiles = []

        self.overwrite_data = overwrite_data

    def use_template_file(self, filename, values):
        """Use a template file for compilation of the code"""
        new_pathname = P(self.rundir, filename)
        self.templates.get_template(filename).stream(values).dump(new_pathname)
        self.path_names.insert(0, new_pathname)

    def clear_workdir(self):
        sh.rm(['-r', self.workdir])
        mkdir(self.workdir)
        log.debug('emptied working directory %r' % self.workdir)

    def clear_rundir(self):
        sh.rm(['-r', self.rundir])
        mkdir(self.rundir)
        log.debug('emptied run directory %r' % self.rundir)

    def _get_default_path_names(self):
        with open(self.path_names_file) as pn:
            return [l.strip() for l in pn]

    def rebuild_namelist(self):
        namelist = f90nml.read(self.namelist_files[0])
        for nl in self.namelist_files[1:]:
            namelist.update(f90nml.read(nl))
        return namelist

    def clone_and_checkout(self, commit=None):
        c = commit or self.commit
        try:
            sh.cd(self.srcdir)
            sh.git.status()
        except Exception as e:
            log.debug('Repository not found at %r. Cloning.' % self.srcdir)
            log.debug(e.message)
            try:
                sh.git.clone(self.repo, self.srcdir)
            except Exception as e:
                log.error('Unable to clone repository %r' % self.repo)
                raise e
        try:
            log.debug('Checking out commit %r' % c)
            sh.cd(self.srcdir)
            sh.git.checkout(c)
        except Exception as e:
            log.error('Unable to checkout commit %r' % c)
            raise e

    def get_restart_file(self, month):
        return P(self.restartdir, 'res_%d.cpio' % (month))

    def write_path_names(self, outdir):
        log.debug('Writing path_names to %r' % P(outdir, 'path_names'))
        with open(P(outdir, 'path_names'), 'w') as pn:
            pn.writelines('\n'.join(self.path_names))

    def write_namelist(self, outdir):
        log.debug('Writing namelist to %r' % P(outdir, 'input.nml'))
        self.namelist.write(P(outdir, 'input.nml'))

    def write_diag_table(self, outdir):
        log.debug('Writing diag_table to %r' % P(outdir, 'diag_table'))
        sh.cp(self.diag_table_file, P(outdir, 'diag_table'))

    def write_field_table(self, outdir):
        log.debug('Writing field_table to %r' % P(outdir, 'field_table'))
        sh.cp(self.field_table_file, P(outdir, 'field_table'))

    def write_details_file(self, outdir):
        info = self.get_info()
        with open(P(outdir, 'details.txt'), 'w') as f:
            f.writelines(info)

    def get_info(self):
        git_commit = sh.git(['rev-parse', 'HEAD'])
        git_branch = sh.git(['symbolic-ref', 'HEAD'])
        git_desc   = sh.git(['describe', '--always'])
        git_show   = sh.git.show('--pretty')
        return {'git_desc': git_desc,
                'git_show': git_show}

    def compile(self):
        mkdir(self.execdir)

        vars = {
            'execdir': self.execdir,
            'mimapy_dir': self.mimapy_dir,
            'srcdir': self.srcdir,
            'workdir': self.workdir
        }

        self.write_path_names(self.workdir)

        self.templates.get_template('compile.sh').stream(**vars).dump(P(self.workdir, 'compile.sh'))
        log.debug('Running compiler')

        sh.bash(P(self.workdir, 'compile.sh'), _out=clean_log_info)
        log.debug('Compilation complete.')


    def runmonth(self, month, restart_file=None, use_restart=True, num_cores=8, overwrite_data=False):
        indir = P(self.rundir, 'INPUT')
        outdir = P(self.datadir, 'run%d' % month)

        if os.path.isdir(outdir):
            if self.overwrite_data or overwrite_data:
                log.debug('Data for month %d already exists and overwrite_data is True. Overwriting.' % month)
                sh.rm('-r', outdir)
            else:
                log.error('Data for month %d already exists but overwrite_data is False. Stopping.' % month)
                # exit(4)
                return False

        # make the output run folder and copy over the input files
        mkdir([indir, P(self.rundir, 'RESTART'), self.restartdir])

        self.write_namelist(self.rundir)
        self.write_field_table(self.rundir)
        self.write_diag_table(self.rundir)

        for filename in self.inputfiles:
            sh.cp([filename, P(indir, os.split(filename)[1])])


        if use_restart:
            if not restart_file:
                # get the restart from previous month
                restart_file = self.get_restart_file(month - 1)
            if not os.path.isfile(restart_file):
                log.error('Restart file not found, expecting file %r' % restart_file)
                exit(2)
            else:
                log.debug('Using restart file %r' % restart_file)
        else:
            log.debug('Running month %r without restart file' % month)
            restart_file = None

        vars = {
            'month': month,
            'datadir': outdir,
            'rundir': self.rundir,
            'execdir': self.execdir,
            'restart_file': restart_file,
            'num_cores': num_cores
        }

        runmonth = self.templates.get_template('runmonth.sh')

        # employ the template to create a runscript
        t = runmonth.stream(**vars).dump(P(self.rundir, 'runmonth.sh'))

        log.debug("Running GFDL for month %r" % month)
        sh.bash(P(self.rundir, 'runmonth.sh'), _out=clean_log_info)
        log.debug("Run for month %r complete" % month)

        mkdir(outdir)

        #restart_file = P(self.restartdir, 'res_%d.tar.gz' % month)
        #sh.tar('zcvf', restart_file, 'RESTART')
        restart_file = P(self.restartdir, 'res_%d.cpio' % month)
        sh.cd(P(self.rundir, 'RESTART'))
        state_files = sh.glob('*.res*')
        sh.cpio('-ov', _in='\n'.join(state_files), _out=restart_file)
        log.debug("Saved restart file %r" % restart_file)
        sh.rm('-r', P(self.rundir, 'RESTART'))

        sh.cp(['-a', self.rundir+'/.', outdir])
        self.clear_rundir()
        return True
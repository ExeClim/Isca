#!/usr/bin/env python

import logging
import os

import f90nml
from jinja2 import Environment, PackageLoader
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

mimapy_dir = P(GFDL_BASE, 'exp', 'mimapy')
mimapy_workdir = P(GFDL_WORK, 'mimapy')

mkdir(P(mimapy_workdir, 'exec'))
mkdir(P(mimapy_workdir, 'restarts'))

templates = Environment(loader=PackageLoader('mima', 'templates'))

variables = {
    'GFDL_BASE': GFDL_BASE,
    'GFDL_WORK': GFDL_WORK,
    'GFDL_DATA': GFDL_DATA
}


class Experiment(object):
    """A basic MiMA experiment"""
    def __init__(self, name, overwrite_data=False):
        super(Experiment, self).__init__()
        self.name = name

        # set the default locations of working directory,
        # executable directory, restart file storage and
        # data directory.
        # These can be overridden e.g. if an experiment is to use
        # the same executable as another
        self.workdir = P(mimapy_workdir, self.name)
        self.execdir = P(mimapy_workdir, 'exec', self.name)
        self.restartdir = P(mimapy_workdir, 'restarts', self.name)
        self.datadir = P(GFDL_DATA, self.name)


        if os.path.isdir(self.workdir):
            log.warning('Working directory for exp %r already exists' % self.name)
        else:
            log.debug('Making directory %r' % self.workdir)
            mkdir(self.workdir)

        self.variables = {
            'GFDL_BASE': GFDL_BASE,
            'GFDL_WORK': GFDL_WORK,
            'GFDL_DATA': GFDL_DATA,
            'workdir': self.workdir,
            'mimapy_dir': mimapy_dir
        }

        self.path_names = self._get_default_path_names()
        self.namelist = self._get_default_namelist()
        self.inputfiles = []

        self.overwrite_data = overwrite_data




    def clear_workdir(self):
        sh.rm(['-r', self.workdir])
        mkdir(self.workdir)
        log.debug('emptied working directory %r' % self.workdir)

    def _get_default_path_names(self):
        with open(P(mimapy_dir, 'path_names')) as pn:
            return [l.strip() for l in pn]

    def _get_default_namelist(self):
        namelist = f90nml.read(P(mimapy_dir, 'core.nml'))
        for nl in ['phys.nml']:
            filename = P(mimapy_dir, nl)
            namelist.update(f90nml.read(filename))
        return namelist

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
        sh.cp(P(mimapy_dir, 'diag_table'), P(outdir, 'diag_table'))

    def write_field_table(self, outdir):
        log.debug('Writing field_table to %r' % P(outdir, 'field_table'))
        sh.cp(P(mimapy_dir, 'field_table'), P(outdir, 'field_table'))

    def write_details_file(self, outdir):
        info = self.get_info()


    def get_info(self):
        git_commit = sh.git(['rev-parse', 'HEAD'])
        git_branch = sh.git(['symbolic-ref', 'HEAD'])
        git_desc   = sh.git(['describe', '--always'])
        git_show   = sh.git.show('--pretty')
        return {'git_desc': git_desc,
                'git_show': git_show}

    def compile(self):
        mkdir(self.execdir)

        vars = self.variables.copy()
        vars.update({
            'compile_dir': self.execdir,
            })

        self.write_path_names(self.workdir)

        templates.get_template('compile.sh').stream(**vars).dump(P(self.workdir, 'compile.sh'))
        log.debug('Running compiler')

        sh.bash(P(self.workdir, 'compile.sh'), _out=clean_log_info)
        log.debug('Compilation complete.')


    def runmonth(self, month, restart_file=None, use_restart=True, num_cores=8, overwrite_data=False):
        indir = P(self.workdir, 'INPUT')
        outdir = P(self.datadir, 'run%d' % month)

        if os.path.isdir(outdir):
            if self.overwrite_data or overwrite_data:
                log.debug('Data for month %d already exists and overwrite_data is True. Overwriting.' % month)
                sh.rm('-r', outdir)
            else:
                log.error('Data for month %d already exists but overwrite_data is False. Stopping.' % month)
                # exit(4)
                return False

        self.write_namelist(self.workdir)
        self.write_field_table(self.workdir)
        self.write_diag_table(self.workdir)

        # make the output run folder and copy over the input files
        mkdir([indir, outdir, P(self.workdir, 'RESTART'), self.restartdir])
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

        vars = self.variables.copy()
        vars.update({
            'month': month,
            'datadir': outdir,
            'workdir': self.workdir,
            'execdir': self.execdir,
            'restart_file': restart_file,
            'num_cores': num_cores
            })
        runmonth = templates.get_template('runmonth.sh')

        # employ the template to create a runscript
        t = runmonth.stream(**vars).dump(P(self.workdir, 'runmonth.sh'))

        log.debug("Running GFDL for month %r" % month)
        sh.bash(P(self.workdir, 'runmonth.sh'), _out=clean_log_info)
        log.debug("Run for month %r complete" % month)

        #restart_file = P(self.restartdir, 'res_%d.tar.gz' % month)
        #sh.tar('zcvf', restart_file, 'RESTART')
        restart_file = P(self.restartdir, 'res_%d.cpio' % month)
        sh.cd(P(self.workdir, 'RESTART'))
        state_files = sh.glob('*.res*')
        sh.cpio('-ov', _in='\n'.join(state_files), _out=restart_file)
        log.debug("Saved restart file %r" % restart_file)
        sh.rm('-r', P(self.workdir, 'RESTART'))

        sh.cp(['-a', self.workdir+'/.', outdir])
        self.clear_workdir()
        return True
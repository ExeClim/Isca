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

templates = Environment(loader=PackageLoader('mima', 'templates'))

variables = {
    'GFDL_BASE': GFDL_BASE,
    'GFDL_WORK': GFDL_WORK,
    'GFDL_DATA': GFDL_DATA
}


class Experiment(object):
    """A basic MiMA experiment"""
    def __init__(self, name):
        super(Experiment, self).__init__()
        self.name = name
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

    @property
    def workdir(self):
        "Base working directory"
        return P(mimapy_workdir, self.name)

    @property
    def execdir(self):
        "The directory where executables are compiled"
        return P(mimapy_workdir, 'exec', self.name)

    @property
    def datadir(self):
        "Directory where completed run data is stored."
        return P(GFDL_DATA, self.name)


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
        sh.bash(P(self.workdir, 'compile.sh'), _out=log.info)
        log.debug('Compilation complete.')


    def runmonth(self, month, restart=None, use_restart=True):
        indir = P(self.workdir, 'INPUT')
        outdir = P(self.datadir, 'run%d' % month)

        self.write_namelist(self.workdir)
        self.write_field_table(self.workdir)
        self.write_diag_table(self.workdir)

        # make the output run folder and copy over the input files
        mkdir([indir, outdir, P(self.workdir, 'RESTART')])
        for filename in self.inputfiles:
            sh.cp([filename, P(indir, os.split(filename)[1])])


        if use_restart:

        else:
            log.info('Running month %r without restart file' % month)
            restart_file = None

        vars = self.variables.copy()
        vars.update({
            'month': month,
            'datadir': outdir,
            'workdir': self.workdir,
            'execdir': self.execdir,
            'restart_file': restart_file,
            })
        runmonth = templates.get_template('runmonth.sh')

        # employ the template to create a runscript
        t = runmonth.stream(**vars).dump(P(self.workdir, 'runmonth.sh'))

        log.debug("Running GFDL for month %r" % month)
        sh.bash(P(self.workdir, 'runmonth.sh'), _out=log.info)
        log.debug("Run for month %r complete" % month)

        sh.cp(['-r', self.workdir, outdir])

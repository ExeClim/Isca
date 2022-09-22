from contextlib import contextmanager
import os
import socket

from jinja2 import Environment, FileSystemLoader
import sh

from isca import GFDL_WORK, GFDL_BASE, GFDL_SOC, _module_directory, get_env_file
from .loghandler import Logger
from .helpers import url_to_folder, destructive, useworkdir, mkdir, cd, git, P, git_run_in_directory

import pdb

class CodeBase(Logger):
    """The CodeBase.

    A CodeBase is a subset of the full FMS source code that is used to
    form a complete model. It can be used to compile a specific version from
    git source control, or compile from an existing directory.

    The CodeBase is a base class, use the derived models such as
    IscaCodeBase or ShallowCodeBase to compile codes.
    """
    # override these parameters in subclasses
    #templatedir = P(_module_directory, 'templates')
    path_names_file = None
    executable_name = None

    @classmethod
    def from_repo(cls, repo, commit=None, **kwargs):
        return cls(repo=repo, commit=commit, **kwargs)

    @classmethod
    def from_directory(cls, directory, **kwargs):
        return cls(directory=directory, **kwargs)

    def __init__(self, repo=None, commit=None, directory=None, storedir=P(GFDL_WORK, 'codebase'), safe_mode=False, socrates_version='1703'):
        """Create a new CodeBase object.

        A CodeBase can be created with either a git repository or a file directory as it's source.
        For example, to compile from a local directory where you are developing and changing the code:

            cb = gfdl.CodeBase(directory='/my/working/GFDLmoistModel')

        Alternatively code can be fetched from a git repository for a specfic commit hash, tag or branch:

            cb = gfdl.CodeBase(repo='git@github.com:execlim/GFDLmoistModel', commit='mytag0.2')

        Each directory or repo-commit can be compiled separately, allowing for multiple
        executables and prevent overwriting of known good states.
        The typical use of a CodeBase is to easily compile a model on a specfic
        computer/cluster, and then pass to an `Experiment` object for running.
        """

        if repo is None and directory is None:
            self.log.error('Not enough sources. Cannot create a CodeBase without either a source directory or source repository.')
            raise AttributeError('Either repo= or directory= required to create CodeBase.')
        if repo is not None and directory is not None:
            self.log.error('Too many sources. Cannot create a CodeBase with both a source directory and a source repository.')
            raise AttributeError('Either repo= or directory= required to create CodeBase.')


        self.safe_mode = safe_mode
        self.storedir = storedir

        if directory is not None:
            self.directory = directory
            self.repo = None
            self.commit = None
            workdir = url_to_folder(directory)
        else:
            self.repo = repo
            self.directory = None
            self.commit = 'HEAD' if commit is None else commit
            workdir = url_to_folder(self.repo) + '-' + self.commit

        # useful directory shortcuts
        self.workdir =  P(self.storedir, workdir)   # base for all codebase I/O actions
        self.codedir =  P(self.workdir, 'code')     # where code is checked out / symlinked to directory
        self.srcdir  =  P(self.codedir, 'src')      # ISCA_CODE/src
        self.builddir = P(self.workdir, 'build', self.executable_name.split('.')[0])
        self.templatedir = P(_module_directory, 'templates')  # templates are stored with the python isca module
        self.executable_fullpath = P(self.builddir, self.executable_name)

        # alias a version of git acting from within the code directory
        self.git = git_run_in_directory(GFDL_BASE, self.codedir)

        # check if the code is available.  If it's not, checkout the repo.
        if not self.code_is_available:
            if self.repo:
                self.log.info('Code not found. Checking out git repo.')
                self.checkout()
            else:
                self.link_source_to(directory)
        elif self.code_is_available and self.commit is not None:
            # problem is that if you try to checkout a specific commit, and it doesn't work, the next time you try it, the above code will only check if code exists, which it will, but it won't be at the correct commit. This will cause problems for e.g. the trip tests. Following code checks if the code that's checked out is the correct commit ID compared to what was asked for, and gives an error if they are different.         
            commit_at_HEAD_of_repo = self.git_commit.split('"')[1]
            commit_desired = self.commit
            if len(commit_desired)==len(commit_at_HEAD_of_repo):
                commit_to_compare_1 = commit_desired
                commit_to_compare_2 = commit_at_HEAD_of_repo
            elif len(commit_desired)>len(commit_at_HEAD_of_repo):
                commit_to_compare_1 = commit_desired[0:len(commit_at_HEAD_of_repo)]
                commit_to_compare_2 = commit_at_HEAD_of_repo
            else:
                commit_to_compare_1 = commit_desired
                commit_to_compare_2 = commit_at_HEAD_of_repo[0:len(commit_desired)]

            if commit_to_compare_1==commit_to_compare_2:
                self.log.info('commit requested successfully checked out')
            else:
                self.log.warn('commit requested is not the commit to be used')
                raise NotImplementedError("commit requested %s but commit supplied %s. This happens when you've previously tried to checkout a particular commit, but the commit was not found in the repo supplied. Try removing %s and trying again, making sure to select a repo that contains your desired commit." % (commit_to_compare_1, commit_to_compare_2, self.workdir ))


        #TODO 

        self.templates = Environment(loader=FileSystemLoader(self.templatedir))

        # read path names from the default file
        self.path_names = []
        self.extra_path_names = []
        self.compile_flags = []  # users can append to this to add additional compiler options

    @property
    def code_is_available(self):
        """Returns True if the repo has been checked out, or the directory
        points to a valid source directory.
        """
        # use the existence of the python directory as a simple test
        return os.path.isdir(P(self.srcdir, 'extra', 'python'))

    @property
    def is_clean(self):
        """Returns False if there are modified files or new files outside
        of git version control in the source directory."""
        raise NotImplementedError

    @property
    def git_commit(self):
        return self.git.log('-1', '--format="%H"').stdout.decode('utf8')

    # @property
    # def git_diff(self):
    #     """Returns the output of `git diff` run in the base directory, comparing to the existing commit."""
    #     return self.git.diff('--no-color')
        # if commit:
        #     commit_consistency_check = commit_id[0:len(commit)]==commit
        #     if not commit_consistency_check:
        #         raise ValueError('commit id specified and commit id actually used are not the same:' +commit+commit_id[0:len(commit)])

        # self.commit_id = commit_id

    def write_source_control_status(self, outfile):
        """Write the current state of the source code to a file."""

        gfdl_git = git_run_in_directory(GFDL_BASE, GFDL_BASE)

        with open(outfile, 'w') as file:
            # write out the git commit id of the compiled source code
            file.write("*---commit hash used for fortran code in workdir---*:\n")
            file.write(self.git_commit)

            # write out the git commit id of GFDL_BASE
            file.write("\n\n*---commit hash used for code in GFDL_BASE, including this python module---*:\n")
            file.write(gfdl_git.log('-1', '--format="%H"').stdout.decode('utf8'))

            # if there are any uncommited changes in the working directory,
            # add those to the file too
            source_status = self.git.status("-b", "--porcelain").stdout.decode('utf8')
            # filter the source status for changes in specific files
            filetypes = ('.f90', '.inc', '.c')
            source_status = [line for line in source_status.split('\n')
                    if any([suffix in line.lower() for suffix in filetypes])]

            # write the status and diff only when something is modified
            if source_status:
                file.write("\n#### Code compiled from dirty commit ####\n")
                file.write("*---git status output (only f90 and inc files)---*:\n")
                file.write('\n'.join(source_status))
                file.write('\n\n*---git diff output---*\n')
                source_diff = self.git.diff('--no-color').stdout.decode('utf8')
                file.write(source_diff)

    def read_path_names(self, path_names_file):
        with open(path_names_file) as pn:
            return [l.strip() for l in pn]

    @useworkdir
    @destructive
    def write_path_names(self, path_names):
        outfile = P(self.builddir, 'path_names')
        self.log.info('Writing path_names to %r' % outfile)
        with open(outfile, 'w') as pn:
            pn.writelines('\n'.join(path_names))

    @useworkdir
    @destructive
    def link_source_to(self, directory):
        # link workdir/code to the directory codebase for simplified paths
        if os.path.exists(self.codedir):
            self.log.info("Relinking %s to %s" % (self.codedir, directory))
            sh.rm(self.codedir)
        else:
            self.log.info("Linking %s to %s" % (self.codedir, directory))
        sh.ln('-s', directory, self.codedir)

    @useworkdir
    @destructive
    def checkout(self):
        if self.repo is None:
            self.log.warn('Cannot checkout a directory.  Use a CodeBase(repo="...") object instead.')
            return None

        try:
            self.git.status()
        except Exception as e:
            self.log.info('Repository not found at %r. Cloning.' % self.codedir)
            # self.log.debug(e.message)
            try:
                git.clone(self.repo, self.codedir)
            except Exception as e:
                self.log.error('Unable to clone repository %r' % self.repo)
                raise e

        if self.commit is not None:
            try:
                self.log.info('Checking out commit %r' % self.commit)
                self.git.checkout(self.commit)
            except Exception as e:
                self.log.error('Unable to checkout commit %r' % self.commit)
                raise e

    def _log_line(self, line):
        line = self.clean_log(line)
        if line is not None:
            if "warning" in line.lower():
                self.log.warn(line)
            else:
                self.log.info(line)

    @useworkdir
    @destructive
    def compile(self, debug=False, optimisation=None):
        env = get_env_file()
        mkdir(self.builddir)

        compile_flags = []
        # if debug:
        #     compile_flags.append('-g')
        #     compile_flags.append('-traceback')
        #     compile_flags.append('-debug all')

        # if optimisation is not None:
        #     compile_flags.append('-O%d' % optimisation)

        compile_flags.extend(self.compile_flags)
        compile_flags_str = ' '.join(compile_flags)

        # get path_names from the directory
        if not self.path_names:
            self.path_names = self.read_path_names(P(self.srcdir, 'extra', 'model', self.name, 'path_names')) + self.extra_path_names
        self.write_path_names(self.path_names)
        path_names_str = P(self.builddir, 'path_names')

        vars = {
            'execdir': self.builddir,
            'template_dir': self.templatedir,
            'srcdir': self.srcdir,
            'workdir': self.workdir,
            'compile_flags': compile_flags_str,
            'env_source': env,
            'path_names': path_names_str,
            'executable_name': self.executable_name,
            'run_idb': debug,
        }

        self.templates.get_template('compile.sh').stream(**vars).dump(P(self.builddir, 'compile.sh'))
        self.log.info('Running compiler')
        for line in sh.bash(P(self.builddir, 'compile.sh'), _iter=True, _err_to_out=True):
            self._log_line(line)

        self.log.info('Compilation complete.')



class IscaCodeBase(CodeBase):
    """The Full Isca Stack.
    This includes moist dynamics and the RRTM radiation scheme.
    """
    name = 'isca'
    executable_name = 'isca.x'

    def disable_soc(self):
        # add no compile flag
        self.compile_flags.append('-DSOC_NO_COMPILE')
        self.log.info('SOCRATES compilation disabled.')

    def __init__(self, *args, **kwargs):
        super(IscaCodeBase, self).__init__(*args, **kwargs)
        self.disable_soc()

class SocratesCodeBase(CodeBase):
    """Isca without RRTM but with the Met Office radiation scheme, Socrates.
    """
    #path_names_file = P(_module_directory, 'templates', 'moist_path_names')
    name = 'socrates'
    executable_name = 'soc_isca.x'

    def disable_rrtm(self):
        # add no compile flag
        self.compile_flags.append('-DRRTM_NO_COMPILE')
        self.log.info('RRTM compilation disabled.')

    def simlink_to_soc_code(self):
        #Make symlink to socrates source code if one doesn't already exist.
        socrates_desired_location = self.codedir+'/src/atmos_param/socrates/src/trunk'

        #First check if socrates is in correct place already
        if os.path.exists(socrates_desired_location):
            link_correct = os.path.exists(socrates_desired_location+'/src/')
            if link_correct:
                socrates_code_in_desired_location=True
            else:
                socrates_code_in_desired_location=False                
                if os.path.islink(socrates_desired_location):
                    self.log.info('Socrates source code symlink is in correct place, but is to incorrect location. Trying to correct.')
                    os.unlink(socrates_desired_location)
                else:
                    self.log.info('Socrates source code is in correct place, but folder structure is wrong. Contents of the folder '+socrates_desired_location+' should include a src folder.')
        else:
            socrates_code_in_desired_location=False
            self.log.info('Socrates source code symlink does not exist. Creating.')

        # If socrates is not in the right place already, then attempt to make symlink to location of code provided by GFDL_SOC
        if socrates_code_in_desired_location:
            self.log.info('Socrates source code already in correct place. Continuing.')
        else:
            if GFDL_SOC is not None:
                sh.ln('-s', GFDL_SOC, socrates_desired_location)
            elif GFDL_SOC is None:
                error_mesg = 'Socrates code is required for SocratesCodebase, but source code is not provided in location GFDL_SOC='+ str(GFDL_SOC)
                self.log.error(error_mesg)
                raise OSError(error_mesg)

    def read_version_specific_paths(self, socrates_version_to_use):
        self.extra_path_names = self.read_path_names(P(self.srcdir, 'extra', 'model', self.name, 'socrates_version_paths', socrates_version_to_use))        

    def __init__(self, *args, **kwargs):
        super(SocratesCodeBase, self).__init__(*args, **kwargs)
        self.disable_rrtm()
        self.simlink_to_soc_code()
        socrates_version_to_use = kwargs['socrates_version']
        self.read_version_specific_paths(socrates_version_to_use)

class GreyCodeBase(CodeBase):
    """The Frierson model.
    This is the closest to the Frierson model, with moist dynamics and a
    choice of grey radiation schemes.

    The Isca code can be configured to be run in exactly the same configuration
    as the Grey codebase, but doing so requires compilation of RRTM which
    can take a while during a development cycle.
    """
    #path_names_file = P(_module_directory, 'templates', 'moist_path_names')
    name = 'grey'
    executable_name = 'grey_isca.x'

    def disable_rrtm(self):
        # add no compile flag
        self.compile_flags.append('-DRRTM_NO_COMPILE')
        self.log.info('RRTM compilation disabled.')

    def disable_soc(self):
        # add no compile flag
        self.compile_flags.append('-DSOC_NO_COMPILE')
        self.log.info('SOCRATES compilation disabled.')

    def __init__(self, *args, **kwargs):
        super(GreyCodeBase, self).__init__(*args, **kwargs)
        self.disable_rrtm()
        self.disable_soc()

class DryCodeBase(GreyCodeBase):
    """The Held-Suarez model.

    Where the moist codebase uses a radiation scheme, incoming solar radiation
    and SSTs to force the model, the Dry model uses a prescribed 'equilibrium'
    temperature profile (Teq).  The model is relaxed towards Teq, generating
    a circulation in response.
    """
    #path_names_file = P(_module_directory, 'templates', 'dry_path_names')
    name = 'dry'
    executable_name = 'held_suarez.x'



# class ShallowCodeBase(CodeBase):
#     """The Shallow Water Equations.
#     """
#     name = 'shallow'
#     executable_name = 'shallow.x'

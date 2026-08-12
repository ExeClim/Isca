from contextlib import contextmanager
import os
import socket

from jinja2 import Environment, FileSystemLoader
import sh

from isca import GFDL_WORK, GFDL_BASE, GFDL_SOC, GFDL_SOC_DIR, _module_directory, get_env_file
from isca.socrates_paths import generate as generate_socrates_path_list
from .loghandler import Logger
from .helpers import url_to_folder, destructive, useworkdir, mkdir, git, P, git_run_in_directory, check_for_sh_stdout

import pdb

# The Socrates version vendored as a git submodule at
# src/atmos_param/socrates/src/<DEFAULT_SOCRATES_VERSION> and used when no
# socrates_version is given explicitly.
DEFAULT_SOCRATES_VERSION = '2026.07.1'

# Where SocratesCodeBase fetches Socrates source from when it isn't
# already available locally (neither vendored, nor under GFDL_SOC /
# GFDL_SOC_DIR). Socrates is BSD-3-Clause licensed and hosted at:
SOCRATES_REPO_URL = 'https://github.com/MetOffice/socrates.git'

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

    def __init__(self, repo=None, commit=None, directory=None, storedir=P(GFDL_WORK, 'codebase'), safe_mode=False):
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
        self.extra_path_names = []  # subclasses (e.g. SocratesCodeBase) can append version-specific files here
        self.compile_flags = []  # users can append to this to add additional compiler options
        # subclasses can set this where mkmf's automatic dependency detection
        # is known to miss a real build-order dependency (see
        # SocratesCodeBase for why) -- makes compile() do a best-effort
        # `make -k` pass before the normal one, so objects with no missing
        # dependency of their own get built regardless of what order mkmf
        # happened to schedule everything in
        self.make_retry_keep_going = False

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
        return check_for_sh_stdout(self.git.log('-1', '--format="%H"'))

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
            gfdl_git_out = check_for_sh_stdout(gfdl_git.log('-1', '--format="%H"'))
            file.write(gfdl_git_out)

            # if there are any uncommited changes in the working directory,
            # add those to the file too
            source_status = check_for_sh_stdout(self.git.status("-b", "--porcelain"))
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
                source_diff = check_for_sh_stdout(self.git.diff('--no-color'))
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
            self.path_names = self.read_path_names(
                P(self.srcdir, 'extra', 'model', self.name, 'path_names')) + self.extra_path_names
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
            'make_retry_keep_going': self.make_retry_keep_going,
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

    Socrates itself is not packaged inside the Isca repository. Isca
    depends on it via a git submodule pinned to DEFAULT_SOCRATES_VERSION
    (see .gitmodules), and this class is responsible for making a matching
    checkout of the requested `socrates_version` available at
    src/atmos_param/socrates/src/<version> and for compiling the right
    subset of it, in priority order:

    1. If that directory is already populated (e.g. the git submodule was
       initialised with `git submodule update --init`, or a previous run
       already fetched it), use it as-is.
    2. Else, if GFDL_SOC is set, symlink to it directly (legacy single-
       version override, kept for backwards compatibility).
    3. Else, if GFDL_SOC_DIR/<version> exists, symlink to that.
    4. Else, fetch <version> from Socrates' GitHub repository into
       GFDL_SOC_DIR (or a default cache under GFDL_WORK if GFDL_SOC_DIR
       isn't set) with a sparse partial clone -- only src/ and make/ are
       fetched (a few MB), not the full ~300MB repository, which is
       almost entirely example/data files Isca's build never touches.

    DEFAULT_SOCRATES_VERSION is a special case, exempt from steps 2-4
    above: that version's path is registered as a git submodule (see
    .gitmodules), and a plain symlink placed there instead -- even one
    pointing at a perfectly good Socrates checkout via GFDL_SOC/
    GFDL_SOC_DIR -- would make git treat the whole superproject repository
    as corrupt (`git status`/`git diff` fail outright with "expected
    submodule path ... not to be a symbolic link"), not just a problem
    local to this one path. So for that one version, unless it's already
    populated (step 1), Isca always runs `git submodule update --init`
    (with a blob filter and a post-init sparse-checkout, to keep it
    reasonably small) instead -- see _init_socrates_submodule. Request a
    different socrates_version to use a custom checkout via GFDL_SOC/
    GFDL_SOC_DIR.

    The list of Socrates source files to compile for a given version
    (`extra_path_names`) is read from
    src/extra/model/socrates/socrates_version_paths/<version> if that's
    been committed to Isca, or generated on the fly otherwise -- see
    isca.socrates_paths for how.
    """
    name = 'socrates'
    executable_name = None  # set per-instance in __init__, since it depends on socrates_version
    executable_prefix = 'soc_isca'  # overridden by SocColumnCodeBase

    def disable_rrtm(self):
        # add no compile flag
        self.compile_flags.append('-DRRTM_NO_COMPILE')
        self.log.info('RRTM compilation disabled.')

    def _socrates_source_dir(self):
        return self.codedir + '/src/atmos_param/socrates/src/' + self.socrates_version

    def _fetch_socrates_source(self, dest_dir):
        """Sparse partial clone of just src/, make/, and the spectral data
        files Isca's own test cases use (data/spectra/ga7, data/spectra/ga3_1
        -- a few MB, not the ~140MB of every spectral file Socrates ships)
        from Socrates' GitHub repository at the tag matching
        self.socrates_version."""
        self.log.info(
            'Socrates version %s not found locally. Fetching src/, make/ and '
            'the ga7/ga3_1 spectral data from %s (tag %s) into %s.' %
            (self.socrates_version, SOCRATES_REPO_URL, self.socrates_version, dest_dir))
        try:
            sh.git('clone', '--depth', '1', '--filter=blob:none', '--sparse',
                   '--branch', self.socrates_version, SOCRATES_REPO_URL, dest_dir)
            sh.git('-C', dest_dir, 'sparse-checkout', 'init', '--cone')
            sh.git('-C', dest_dir, 'sparse-checkout', 'set', 'src', 'make',
                   'data/spectra/ga7', 'data/spectra/ga3_1')
        except Exception as e:
            error_mesg = (
                'Failed to fetch Socrates version %s from %s. Check the '
                'version/tag exists at that repository, and that this '
                'machine has network access to GitHub. If it does not '
                '(e.g. on an HPC compute node), fetch it in advance '
                'somewhere with access and point GFDL_SOC_DIR/%s (or '
                'GFDL_SOC) at it instead. Underlying error: %s' %
                (self.socrates_version, SOCRATES_REPO_URL, self.socrates_version, e))
            self.log.error(error_mesg)
            raise OSError(error_mesg)

    def _init_socrates_submodule(self, dest_dir):
        """Populate dest_dir -- a path registered as a git submodule in
        .gitmodules, pinned to DEFAULT_SOCRATES_VERSION -- as a real,
        shallow, sparse submodule checkout, rather than the plain symlink
        ensure_socrates_source() uses for every other version/location.

        A plain symlink at a path git has registered as a submodule (via
        .gitmodules) makes git think the *superproject itself* is corrupt:
        `git status`/`git diff` fail outright, for the whole repository,
        with "error: expected submodule path ... not to be a symbolic
        link" -- not a failure localised to this one path, so it can't be
        worked around the way the ordinary fetch-and-symlink path can, and
        it hits every new user who follows the zero-config quick start
        (SocratesCodeBase.from_directory(...) with no GFDL_SOC/GFDL_SOC_DIR
        and no prior `git submodule update --init`), since the default
        version is exactly the one pinned as a submodule.
        """
        self.log.info(
            'Socrates version %s is pinned as a git submodule (see .gitmodules); '
            'initialising it as a real submodule checkout at %s (rather than the '
            'symlink used for other versions), so git does not treat a plain '
            'symlink at a submodule path as repository corruption.' %
            (self.socrates_version, dest_dir))
        try:
            # self.codedir can itself be a symlink into the real repo (see
            # CodeBase.link_source_to); resolve it to a real path first and
            # use a repo-root-relative pathspec, rather than handing git an
            # absolute path built through that symlink, since `git
            # submodule` matches pathspecs against its own canonicalised
            # repo root and a literal (unresolved) symlink prefix may not
            # match it. realpath() resolves the existing symlinked prefix
            # even though dest_dir itself doesn't exist on disk yet.
            repo_root = os.path.realpath(self.codedir)
            relative_submodule_path = os.path.relpath(os.path.realpath(dest_dir), repo_root)
            sh.git('-C', repo_root, 'submodule', 'update', '--init',
                   '--filter=blob:none', '--', relative_submodule_path)
            sh.git('-C', dest_dir, 'sparse-checkout', 'init', '--cone')
            sh.git('-C', dest_dir, 'sparse-checkout', 'set', 'src', 'make',
                   'data/spectra/ga7', 'data/spectra/ga3_1')
        except Exception as e:
            error_mesg = (
                'Failed to initialise the Socrates git submodule at %s. Try '
                'running `git submodule update --init` from the top of this '
                'repository yourself. Underlying error: %s' % (dest_dir, e))
            self.log.error(error_mesg)
            raise OSError(error_mesg)

    def ensure_socrates_source(self):
        """Make sure a Socrates checkout for self.socrates_version is
        available at src/atmos_param/socrates/src/<version> (symlinked in,
        except for the git-submodule-pinned default version -- see
        _init_socrates_submodule), fetching it if necessary. Returns the
        resolved source directory (containing src/ and make/)."""
        socrates_desired_location = self._socrates_source_dir()

        if os.path.exists(socrates_desired_location) and os.path.exists(socrates_desired_location + '/src/'):
            self.log.info('Socrates source for version %s already in place at %s.' %
                           (self.socrates_version, socrates_desired_location))
            return os.path.realpath(socrates_desired_location)

        is_empty_submodule_placeholder = (
            self.socrates_version == DEFAULT_SOCRATES_VERSION
            and os.path.isdir(socrates_desired_location)
            and not os.path.islink(socrates_desired_location)
            and not os.listdir(socrates_desired_location)
        )

        if os.path.islink(socrates_desired_location):
            # broken symlink (e.g. left over from a different GFDL_SOC_DIR)
            self.log.info('Socrates source symlink for version %s points somewhere invalid. Recreating.' %
                           self.socrates_version)
            os.unlink(socrates_desired_location)
        elif is_empty_submodule_placeholder:
            # A plain `git clone` of the superproject (without
            # --recurse-submodules, e.g. GitHub Actions' default checkout)
            # still materialises an empty directory at a submodule's path
            # even though its content was never fetched -- not the same as
            # the "missing entirely" case _init_socrates_submodule's `git
            # submodule update --init` otherwise runs against, but it
            # populates into an existing empty directory just as well, so
            # there's nothing to clean up here.
            self.log.info('Socrates source directory for version %s exists but is an empty, '
                           'uninitialised submodule placeholder. Initialising.' % self.socrates_version)
        elif os.path.exists(socrates_desired_location):
            error_mesg = ('%s exists but does not look like a Socrates checkout (no src/ '
                          'subdirectory found inside it).' % socrates_desired_location)
            self.log.error(error_mesg)
            raise OSError(error_mesg)

        if self.socrates_version == DEFAULT_SOCRATES_VERSION:
            # This path is registered as a git submodule (see .gitmodules),
            # so it must always end up a real git checkout, never a
            # symlink -- including when GFDL_SOC/GFDL_SOC_DIR are set,
            # since content from either (a hand-downloaded tarball
            # extraction, say) is not guaranteed to itself be a git
            # checkout that could safely be symlinked in here either. The
            # GFDL_SOC/GFDL_SOC_DIR overrides below remain available for
            # every *other* socrates_version.
            self._init_socrates_submodule(socrates_desired_location)
            return os.path.realpath(socrates_desired_location)

        if GFDL_SOC_DIR is None and GFDL_SOC is not None:
            # Legacy single-version override: whatever GFDL_SOC points at is used
            # regardless of the requested version.
            source_dir = GFDL_SOC
        else:
            cache_dir = GFDL_SOC_DIR if GFDL_SOC_DIR is not None else P(GFDL_WORK, 'socrates_src')
            source_dir = P(cache_dir, self.socrates_version)
            if not os.path.exists(source_dir):
                mkdir(cache_dir)
                self._fetch_socrates_source(source_dir)

        sh.ln('-s', source_dir, socrates_desired_location)
        return os.path.realpath(socrates_desired_location)

    def load_version_path_names(self, socrates_source_dir):
        committed_path_names = P(self.srcdir, 'extra', 'model', 'socrates',
                                  'socrates_version_paths', self.socrates_version)
        if os.path.exists(committed_path_names):
            self.extra_path_names = self.read_path_names(committed_path_names)
            return

        self.log.info(
            'No committed path list for Socrates version %s at %s. Generating one '
            'from the checked-out source (see isca.socrates_paths). Consider adding '
            'the generated file at that path to Isca so other users of this version '
            "don't need to regenerate it." % (self.socrates_version, committed_path_names))
        interface_dir = P(self.srcdir, 'atmos_param', 'socrates', 'interface')
        relpaths = generate_socrates_path_list(socrates_source_dir, interface_dir, warn=self.log.warn)
        self.extra_path_names = [
            'atmos_param/socrates/src/%s/src/%s' % (self.socrates_version, p) for p in relpaths]

    def __init__(self, *args, socrates_version=DEFAULT_SOCRATES_VERSION, **kwargs):
        self.socrates_version = socrates_version
        # self.executable_prefix resolves via the instance's actual class (so a
        # SocColumnCodeBase instance picks up its own override below), and must be
        # set before CodeBase.__init__ runs since that's where builddir and
        # executable_fullpath get computed from self.executable_name. Dots in the
        # version (e.g. '2026.07.1') are replaced with underscores here because
        # CodeBase.__init__ derives builddir from executable_name.split('.')[0],
        # i.e. everything before the *first* dot -- a version containing dots
        # would otherwise silently truncate and collide with other versions.
        self.executable_name = '%s_%s.x' % (self.executable_prefix, socrates_version.replace('.', '_'))
        super(SocratesCodeBase, self).__init__(*args, **kwargs)
        self.disable_rrtm()
        # mkmf's automatic dependency detection only scans each compiled
        # file's own literal USE/MODULE statements -- never the text of
        # files it pulls in via a plain Fortran INCLUDE (see bin/mkmf's
        # scanfile_for_keywords). Some Socrates source is structured that
        # way (e.g. aux/seaalbedo_driver.f both calls into and INCLUDEs the
        # whole of aux/seaalbedo.f -- see _drop_textually_included_duplicates
        # in isca.socrates_paths for why seaalbedo.f itself isn't compiled
        # independently), so a real module dependency inherited through the
        # INCLUDE (e.g. seaalbedo.f's own USE spline_evaluate_mod) is
        # invisible to mkmf, leaving that one file's build order relative to
        # its provider unconstrained in the generated Makefile. See
        # CodeBase.make_retry_keep_going for how this is worked around.
        self.make_retry_keep_going = True
        socrates_source_dir = self.ensure_socrates_source()
        self.load_version_path_names(socrates_source_dir)

class SocColumnCodeBase(SocratesCodeBase):
    """Isca without RRTM but with the Met Office radiation scheme, Socrates. THIS VERSION FOR SINGLE COLUMN USE.
    """
    name = 'socrates_column'
    executable_prefix = 'soc_column_isca'

    def column_model(self):
        self.compile_flags.append('-DCOLUMN_MODEL')
        self.log.info('USING SINGLE COLUMN MODEL')

    def __init__(self, *args, socrates_version=DEFAULT_SOCRATES_VERSION, **kwargs):
        super(SocColumnCodeBase, self).__init__(*args, socrates_version=socrates_version, **kwargs)
        self.column_model()

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

class ColumnCodeBase(CodeBase):
    """This contains code that will allow one to use all model physics in a single column configuration (i.e. without calling the dynamical core)
    """
    #path_names_file = P(_module_directory, 'templates', 'moist_path_names')
    name = 'column'
    executable_name = 'column_isca.x'

    def column_model(self):
        self.compile_flags.append('-DCOLUMN_MODEL')
        self.log.info('USING SINGLE COLUMN MODEL')

    def disable_soc(self):
        # add no compile flag
        self.compile_flags.append('-DSOC_NO_COMPILE')
        self.log.info('SOCRATES compilations diabled.') 

    def __init__(self, *args, **kwargs):
        super(ColumnCodeBase, self).__init__(*args, **kwargs)
        self.column_model()
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



class ShallowCodeBase(CodeBase):
    """The Shallow Water Equations.
    """
    name = 'shallow'
    executable_name = 'shallow.x'

class BarotropicCodeBase(CodeBase):
    """The Barotropic vorticity equations.
    """
    name = 'barotropic'
    executable_name = 'barotropic_isca.x'
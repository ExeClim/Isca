import os
from functools import wraps

import sh

mkdir = sh.mkdir.bake('-p')
# cd = sh.cd
git = sh.git.bake('--no-pager')

P = os.path.join

# DECORATORS
def destructive(fn):
    """functions decorated with `@destructive` are prevented from running
    in safe mode."""
    @wraps(fn)
    def _destructive(*args, **kwargs):
        self = args[0]
        if self.safe_mode:
            raise AttributeError('Cannot run destructive function %s in safe mode.' % fn.__name__)
        else:
            return fn(*args, **kwargs)
    return _destructive

def useworkdir(fn):
    """Decorator for functions that need to write to a working directory.
    Ensures before the function is run that the workdir exists."""
    @wraps(fn)
    def _useworkdir(*args, **kwargs):
        self = args[0]
        # TODO: check object has workdir attribute.
        if os.path.isdir(self.workdir):
            pass
            #self.log.warning('Working directory for exp %r already exists' % self.name)
        else:
            self.log.debug('Making directory %r' % self.workdir)
            mkdir(self.workdir)
        return fn(*args, **kwargs)
    return _useworkdir

def url_to_folder(url):
        """Convert a url to a valid folder name."""
        for sym in ('/:@'):
            url = url.replace(sym, '_')
        return url

def get_git_commit_id(directory):
    try:
        git_dir = P(directory, '.git')
        commit_id = git("--git-dir="+git_dir, "log", "--pretty=format:'%H'", "-n 1")
        commit_id = str(commit_id).split("'")[1]
    except:
        commit_id = None
    return commit_id

def git_diff(directory):
    git_diff_output = sh.git("-C", directory, "diff", "--no-color", ".")
    git_diff_output = str(git_diff_output).split("\n")
    return git_diff_output

def git_run_in_directory(GFDL_BASE_DIR, dir_in):
    """Function that bakes in git command to run in a specified directory.
       `Try except` is to catch different versions of this command
       in versions of git >1.8 (`git -C`) and older versions for which
       `git -C` is not a recognized command. """
       
    try:
        codedir_git = git.bake('-C', GFDL_BASE_DIR)        
        git_test = check_for_sh_stdout(codedir_git.log('-1', '--format="%H"'))
        baked_git_fn = git.bake('-C', dir_in)
    except:
        codedir_git = git.bake('--git-dir='+GFDL_BASE_DIR+'/.git', '--work-tree='+GFDL_BASE_DIR)
        git_test = check_for_sh_stdout(codedir_git.log('-1', '--format="%H"'))
        baked_git_fn = git.bake('--git-dir='+dir_in+'/.git', '--work-tree='+dir_in)        
        
    return baked_git_fn       

def check_for_sh_stdout(input_exp):
    """Versions of sh>2.* have started returning str types rather than a sh.RunningCommand type. To distinguish these possibilites, this function looks at the output of a sh expression, and asks for stdout and decodes it only if the type is a sh.RunningCommand."""

    if type(input_exp)==sh.RunningCommand:     
        input_exp=input_exp.stdout.decode('utf8')
    
    return input_exp
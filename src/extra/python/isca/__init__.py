import os
from collections import defaultdict

from isca.loghandler import log

_module_directory = os.path.dirname(os.path.realpath(__file__))

try:
    GFDL_BASE        = os.environ['GFDL_BASE']
    GFDL_WORK        = os.environ['GFDL_WORK']
    GFDL_DATA        = os.environ['GFDL_DATA']
except Exception as e:
    log.error('Environment variables GFDL_BASE, GFDL_WORK, GFDL_DATA must be set')
    raise ValueError('Environment variables GFDL_BASE, GFDL_WORK, GFDL_DATA must be set')

# GFDL_ENV: The environment on which the model is being run.
# Primarily, this determines which compilers and libraries are to be used
# to compile the code.
# A file with matching name in GFDL_BASE/src/extra/env is sourced before
# compilation or the model is run.  For example, if the user has `export GFDL_ENV=emps-gv`
# then before compilation or running the code, GFDL_BASE/src/extra/env/emps-gv
# is sourced to ensure that the same libraries will be used to run the code
# as those against which it was compiled.  To get the python scripts working
# on a new computer with different compilers, a new environment script will
# need to be developed and placed in this directory.
try:
    GFDL_ENV = os.environ['GFDL_ENV']
except:
    # if the user doesn't have the environment variable set, use the computer's
    # fully qualified domain name as GFDL_ENV
    import socket
    GFDL_ENV = socket.getfqdn()
    log.warning('Environment variable GFDL_ENV not set, using "%s".' % GFDL_ENV)

try:
    GFDL_SOC = os.environ['GFDL_SOC']
except:
    # if the user doesn't have the SOC variable set, then use None
    GFDL_SOC = None
    log.warning('Environment variable GFDL_SOC not set, but this is only required if using SocratesCodebase. Setting to '+str(GFDL_SOC))

def get_env_file(env=GFDL_ENV):
    filepath = os.path.join(GFDL_BASE, 'src', 'extra', 'env', env)
    if os.path.exists(filepath):
        return filepath
    else:
        log.error('Environment file %s not found' % filepath)
        raise IOError('Environment file %s not found' % filepath)

class EventEmitter(object):
    """A very simple event driven object to make it easier
    to tie custom functionality into the model."""
    def __init__(self):
        self._events = defaultdict(list)

    def on(self, event, fn=None):
        """Assign function `fn` to be called when `event` is triggered.

        Can be used in two ways:
        1. As a normal function:
            exp.on('start', fn)
        2. As a decorator:
            @exp.on('start')
            def fn(*args):
                pass
        """
        def _on(fn):
            self._events[event].append(fn)
            return fn

        if fn is None:
            return _on     # used as a decorator
        else:
            return _on(fn) # used as a normal function

    def emit(self, event, *args, **kwargs):
        """Trigger an event."""
        handled = False
        for callback in self._events[event]:
            callback(*args, **kwargs)
            handled = True
        return handled


from isca.experiment import Experiment, DiagTable, Namelist, FailedRunError
from isca.codebase import IscaCodeBase, SocratesCodeBase, DryCodeBase, GreyCodeBase #, ShallowCodeBase

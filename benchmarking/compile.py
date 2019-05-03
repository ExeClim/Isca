import numpy as np

from isca import IscaCodeBase, GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from isca.util import exp_progress
from util import get_core_list
import sh
import os
import namelist

# from ntfy import notify

NCORES = 16
CORE_LIST = get_core_list()

# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = GreyCodeBase.from_directory(GFDL_BASE)
cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

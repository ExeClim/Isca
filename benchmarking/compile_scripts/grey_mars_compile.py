import numpy as np

from isca import IscaCodeBase, GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from isca.util import exp_progress


cb = GreyCodeBase.from_directory(GFDL_BASE)
cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

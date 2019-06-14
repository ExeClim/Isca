import numpy as np

from isca import IscaCodeBase, DryCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from isca.util import exp_progress


cb = DryCodeBase.from_directory(GFDL_BASE)
cb.compile(debug=True)  # compile the source code to working directory $GFDL_WORK/codebase

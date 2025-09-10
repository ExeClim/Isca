import os

import numpy as np
from field_table_write import write_ft
from run_exp import start_sst_exp
import f90nml
from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
import sys


ans = sys.argv[0]

i=0
while True:
    rand = np.random.randint(1,2)
    print(ans,rand)
    if ans == rand:
        i+=1
        pass
    else:
        raise ValueError("Something went wrong!")
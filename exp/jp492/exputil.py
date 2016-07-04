import numpy as np
import gfdl.experiment

import namelists as _namelists
import diagtables as _diagtables

def get_tagged_source(tag):
    baseexp = gfdl.experiment.Experiment('base_%s' % tag,
        repo='git@github.com:jamesp/GFDLmoistModel.git',
        commit=tag)
    baseexp.disable_rrtm()
    baseexp.compile()
    return baseexp

def get_current_source():
    baseexp = gfdl.experiment.Experiment('_current_state')
    return baseexp


def get_namelist(namelist):
    return getattr(_namelists, namelist).copy()

def get_diagtable(diagtable):
    return getattr(_diagtables, diagtable).copy()

def new_experiment(name, tag, namelist, diagtable):
    baseexp = get_tagged_source(tag)
    baseexp.namelist = get_namelist(namelist)
    baseexp.diag_table = get_diagtable(diagtable)
    return baseexp.derive(name)
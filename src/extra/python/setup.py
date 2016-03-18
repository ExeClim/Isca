# This file provides the installation of the python library 'gfdl'
# Suggested installation procedure:
#       $ cd $GFDL_BASE/src/extra/python
#       $ pip install -e .
# This installs the package in *development mode* i.e. any changes you make to the python files
# or any additional files you add will be immediately available.
# In a new python console, from any directory, you can now access the execlim code:
#       >>> from gfdl import experiment
#       >>> exp = experiment.Experiment()
#       ...

from distutils.core import setup

setup(name='GFDL',
      version='0.1',
      description='GFDL utilities for running experiments and performing data analysis',
      author='James Penn',
      url='https://github.com/ExeClim/GFDLMoistModel',
      packages=['gfdl'],
      install_requires=[
        'sh',
        'jinja2',
        'f90nml',
        'numpy',
        'pandas',
        'xarray'
      ]
     )
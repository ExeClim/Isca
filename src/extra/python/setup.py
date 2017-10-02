# This file provides the installation of the python library 'isca'
# Suggested installation procedure:
#       $ cd $GFDL_BASE/src/extra/python
#       $ pip install -e .
# This installs the package in *development mode* i.e. any changes you make to the python files
# or any additional files you add will be immediately available.
# In a new python console, from any directory, you can now access the execlim code:
#       >>> from isca import experiment
#       >>> exp = experiment.Experiment()
#       ...

from distutils.core import setup

setup(name='Isca',
      version='0.2',
      description='Isca utilities for running experiments and performing data analysis',
      author='Isca Team',
      url='https://github.com/ExeClim/Isca',
      packages=['isca'],
      install_requires=[
        'sh',
        'jinja2',
        'f90nml',
        'numpy',
        'pandas',
        'xarray'
      ]
     )
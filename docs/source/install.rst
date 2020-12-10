.. _install:

Installation
============

Supported language:

- Python: 3.7

.. note::

   We highly recommend installing and using the free `Anaconda
   <https://www.anaconda.com/download/>`_ distribution of Python (or
   `Miniconda <https://conda.io/miniconda.html>`_, if you don't want
   all of the extra packages that come built-in with Anaconda), which
   works on Mac, Linux, and Windows, both on normal computers and
   institutional clusters and doesn't require root permissions.

Isca depends on the following libraries

- f90nml
- fortran-compiler
- jinja2
- libgfortran
- netcdf-fortran
- numpy
- openmpi
- pandas
- python=3.7
- pip
- pytest
- sh
- tqdm
- xarray

After the required packages are installed, Isca can be installed from source.


Clone from Github
-----------------

You can obtain it directly from the `Github repo <https://github.com/execlim/isca>`_ ::

  git clone https://www.github.com/execlim/isca.git
  cd  isca/src/extra/python
  pip install -e .

Verifying proper installation
-----------------------------

Once installed via any of these methods, you can run isca's suite of
tests using `py.test <http://doc.pytest.org/>`_.  From the top-level
directory of the isca installation ::

  conda install pytest  # if you don't have it already; or 'pip install pytest'
  py.test

If the pytest command results in any error messages or test failures,
something has gone wrong, and please refer to the Troubleshooting
information below.

Troubleshooting
---------------

Please search through the `Issues page`_ on Github if anybody else has had the same problem you're facing.
If none do, then please send open a new Issue.

.. _Issues page: https://github.com/execlim/isca/issues

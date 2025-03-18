Quick start guide for running Isca on NSF NCAR's Derecho HPC Machine
-- March 2025, Rory Kelly (rory@ucar.edu)
1. Download the Isca code
   ```bash
   git clone https://github.com/ExeClim/Isca
   ```

2. Set up the environment
   ```bash
   cd Isca

   module --force purge
   module load ncarenv/23.09 intel/2023.2.1 craype/2.7.31 cray-mpich/8.1.27 ncarcompilers/1.0.0 hdf5/1.12.2 netcdf/4.9.2 conda/latest

   mamba env create -f ci/environment-py3.12_ncar-derecho.yml
   conda activate isca_env

   cd src/extra/python
   pip install -e .
   ```
3. Add Isca settings to `~/.bashrc` (modify paths as desired). After adding these settings you'll
   need to source your ~/.bashrc file, or export the settings in the current shell as well.
   ```bash
   export GFDL_ENV=ncar-derecho
   export GFDL_MKMF_TEMPLATE=ncar-derecho-intel
   export GFDL_BASE=/glade/work/$USER/Isca
   export GFDL_WORK=/glade/derecho/scratch/$USER/isca_work
   export GFDL_DATA=/glade/derecho/scratch/$USER/isca_data
   ```
4. Build and Run a Test Case on a batch node
   ```bash
   cd $GFDL_BASE/exp/test_cases/held_suarez
   qsub isca-test.sh
   ```
   Example batch script `isca-test.sh` (you'll need to add you account code to the `-A` argument)
   ```bash
   #!/bin/bash
   #PBS -l walltime=00:30:00
   #PBS -l select=1:ncpus=128:mpiprocs=128
   #PBS -q main
   #PBS -A [ACCOUNT CODE]
   #PBS -N Isca-test
   #PBS -k eod
   #PBS -j oe
   #PBS -o Isca-test.out

   module --force purge
   module load ncarenv/23.09 intel/2023.2.1 craype/2.7.31 cray-mpich/8.1.27 ncarcompilers/1.0.0 hdf5/1.12.2 netcdf/4.9.2 conda/latest
   conda activate isca_env

   cd $GFDL_BASE/exp/test_cases/held_suarez
   python held_suarez_test_case.py
   ```

#!/usr/bin/env bash
# Compiles the mppncombine tool (tested on emps-gv2.ex.ac.uk) - copied from JP's compile.sh 


#-----------------------------------------------------------------------------------------------------
hostname=`hostname`
ppdir=./                 # path to directory containing the tool for combining distributed diagnostic output files
#-----------------------------------------------------------------------------------------------------



# 2. Load the necessary tools into the environment
module purge
source $GFDL_BASE/src/extra/env/$GFDL_ENV
module list
netcdf_flags=`nf-config --fflags --flibs`
#--------------------------------------------------------------------------------------------------------
# compile combine tool
#cd $ppdir
$CC -O -c `nf-config --cflags` mppnccombine.c
if [ $? != 0 ]; then
    echo "ERROR: could not compile combine tool"
    exit 1
fi
$CC -O -o mppnccombine.x `nf-config --fflags --flibs`  mppnccombine.o
if [ $? != 0 ]; then
    echo "ERROR: could not compile combine tool"
    exit 1
fi
#--------------------------------------------------------------------------------------------------------

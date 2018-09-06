#!/usr/bin/env bash
# Compiles the mppncombine tool (tested on emps-gv2.ex.ac.uk) - copied from JP's compile.sh 


#-----------------------------------------------------------------------------------------------------
hostname=`hostname`
ppdir=./                 # path to directory containing the tool for combining distributed diagnostic output files
#-----------------------------------------------------------------------------------------------------



# 2. Load the necessary tools into the environment
source $GFDL_BASE/src/extra/env/$GFDL_ENV
netcdf_flags=`nc-config --cflags --libs`
#--------------------------------------------------------------------------------------------------------
# compile combine tool
#cd $ppdir
$CC -O -c mppnccombine.c $netcdf_flags
if [ $? != 0 ]; then
    echo "ERROR: could not compile combine tool"
    exit 1
fi
$CC -O -o mppnccombine.x mppnccombine.o $netcdf_flags
if [ $? != 0 ]; then
    echo "ERROR: could not compile combine tool"
    exit 1
fi
#--------------------------------------------------------------------------------------------------------

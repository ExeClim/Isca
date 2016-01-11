#!/usr/bin/env bash
# Compiles the MiMa GFDL Moist Model (tested on emps-gv2.ex.ac.uk)


#-----------------------------------------------------------------------------------------------------
hostname=`hostname`
platform=ia64                                             # A unique identifier for your platform
template={{ GFDL_BASE }}/bin/mkmf.template.$platform      # path to template for your platform
mkmf={{ GFDL_BASE }}/bin/mkmf                             # path to executable mkmf
sourcedir={{ GFDL_BASE }}/src                             # path to directory containing model source code
pathnames={{ workdir }}/path_names                        # path to file containing list of source paths
ppdir={{ GFDL_BASE }}/postprocessing                      # path to directory containing the tool for combining distributed diagnostic output files
#-----------------------------------------------------------------------------------------------------
execdir={{ compile_dir }}        # where code is compiled and executable is created
executable=$execdir/fms_moist.x

netcdf_flags=`nf-config --fflags --flibs`

# 2. Load the necessary tools into the environment
source {{ GFDL_BASE }}/src/extra/loadmodule
module list
ulimit -s unlimited # Set stack size to unlimited

#--------------------------------------------------------------------------------------------------------
# compile combine tool
cd $ppdir
cc -O -c `nf-config --cflags` mppnccombine.c
if [ $? != 0 ]; then
    echo "ERROR: could not compile combine tool"
    exit 1
fi
cc -O -o mppnccombine.x `nf-config --libs`  mppnccombine.o
if [ $? != 0 ]; then
    echo "ERROR: could not compile combine tool"
    exit 1
fi
#--------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------
# setup directory structure
if [ ! -d $execdir ]; then
    echo "Creating exec directory $execdir"
    mkdir -p $execdir
fi
cd $execdir

# execute mkmf to create makefile
cppDefs="-Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DINTERNAL_FILE_NML -DOVERLOAD_C8"
$mkmf -a $sourcedir -t $template -p `basename $executable` -c "$cppDefs" $pathnames $sourcedir/shared/include $sourcedir/shared/mpp/include
if [ $? != 0 ]; then
   echo "ERROR: mkmf failed for fms_moist"
   exit 1
fi

# --- execute make ---
make `basename $executable`
if [ $? != 0 ]; then
    echo "ERROR: make failed for fms_moist"
    exit 1
fi

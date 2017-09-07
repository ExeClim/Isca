#!/usr/bin/env bash
# Compiles the MiMa GFDL Moist Model (tested on emps-gv2.ex.ac.uk)


#-----------------------------------------------------------------------------------------------------
hostname=`hostname`
template={{ template_dir }}/mkmf.template.ia64
mkmf={{ srcdir }}/bin/mkmf                             # path to executable mkmf
sourcedir={{ srcdir }}/src                             # path to directory containing model source code
pathnames={{ workdir }}/path_names                        # path to file containing list of source paths
ppdir={{ srcdir }}/postprocessing                      # path to directory containing the tool for combining distributed diagnostic output files
debug={{ run_idb }}                                     # logical to identify if running in debug mode or not
template_debug={{ template_dir }}/mkmf.template.debug
#-----------------------------------------------------------------------------------------------------
execdir={{ execdir }}        # where code is compiled and executable is created
executable=$execdir/fms_moist.x

netcdf_flags=`nc-config --fflags --flibs`

# 2. Load the necessary tools into the environment
module purge
source {{ env_source }}
module list
ulimit -s unlimited # Set stack size to unlimited
export MALLOC_CHECK_=0

if [[ ! -e {{ execdir }}/mppncombine.x ]]; then
  # --------------------------------------------------------------------------------------------------------
  # compile combine tool
  cd $ppdir
  cc -O -c `nc-config --cflags` mppnccombine.c
  if [ $? != 0 ]; then
      echo "ERROR: could not compile combine tool"
      exit 1
  fi
  cc -O -o mppnccombine.x `nc-config --libs`  mppnccombine.o
  if [ $? != 0 ]; then
      echo "ERROR: could not compile combine tool"
      exit 1
  fi
  #./compile_mppn.sh

  ln -s $ppdir/mppnccombine.x {{ execdir }}/mppnccombine.x
  #--------------------------------------------------------------------------------------------------------
fi

#--------------------------------------------------------------------------------------------------------

cd $execdir

if [ $debug == true ]; then

  echo "Compiling in debug mode"

  cppDefs="-Duse_netCDF  -Duse_LARGEFILE -DINTERNAL_FILE_NML -DOVERLOAD_C8 {{compile_flags}}"
  $mkmf  -a $sourcedir -t $template_debug -p `basename $executable` -c "$cppDefs" $pathnames $sourcedir/shared/include $sourcedir/shared/mpp/include

else
  # execute mkmf to create makefile
  cppDefs="-Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DINTERNAL_FILE_NML -DOVERLOAD_C8 {{compile_flags}}"
  $mkmf  -a $sourcedir -t $template -p `basename $executable` -c "$cppDefs" $pathnames $sourcedir/shared/include $sourcedir/shared/mpp/include
fi

make

# $mkmf $make_flags -a $source_dir  -p fms_moist.x -t   $template \
#     -c "-Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DINTERNAL_FILE_NML -DOVERLOAD_C8" $pathnames $sourcedir/shared/mpp/include $sourcedir/shared/constants $sourcedir/include
#     make

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

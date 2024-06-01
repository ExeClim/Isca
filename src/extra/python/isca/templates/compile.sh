#!/usr/bin/env bash
# Compiles the Isca Model

# 0. Source the environment file to load appropriate variables
source {{ env_source }}

# 1. Configuration
hostname=`hostname`
compiler=${GFDL_MKMF_TEMPLATE:-ia64}
echo compiler=$compiler
template={{ template_dir }}/mkmf.template.${compiler}
mkmf={{ srcdir }}/../bin/mkmf                             # path to executable mkmf
sourcedir={{ srcdir }}                             # path to directory containing model source code
pathnames={{ path_names }}                      # path to file containing list of source paths
ppdir={{ srcdir }}/../postprocessing                      # path to directory containing the tool for combining distributed diagnostic output files
debug={{ run_idb }}                                     # logical to identify if running in debug mode or not
template_debug={{ template_dir }}/mkmf.template.debug
#-----------------------------------------------------------------------------------------------------
execdir={{ execdir }}        # where code is compiled and executable is created
executable={{ executable_name }}

netcdf_flags=`nf-config --fflags --flibs`

ulimit -s unlimited # Set stack size to unlimited
export MALLOC_CHECK_=0

# 3. compile the mppncombine tool if it hasn't yet been done.
if [ ! -e "{{ execdir }}/mppnccombine.x" ]; then
  echo "Compiling postprocessing tools"
  cd $ppdir
  ./compile_mppn.sh
  # cc -O -c `nc-config --cflags` mppnccombine.c
  # if [ $? != 0 ]; then
  #     echo "ERROR: could not compile combine tool"
  #     exit 1
  # fi
  # cc -O -o mppnccombine.x `nc-config --libs`  mppnccombine.o
  # if [ $? != 0 ]; then
  #     echo "ERROR: could not compile combine tool"
  #     exit 1
  # fi

  ln -s $ppdir/mppnccombine.x {{ execdir }}/mppnccombine.x
  ln -s $ppdir/mppnccombine_run.sh {{ execdir }}/mppnccombine_run.sh
fi

cd $execdir

echo $pathnames


if [ $debug == True ]; then

 echo "Compiling in debug mode"

# execute mkmf to create makefile
cppDefs="-Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DINTERNAL_FILE_NML -DOVERLOAD_C8 {{compile_flags}}"
$mkmf  -a $sourcedir -t $template_debug -p $executable -c "$cppDefs" $pathnames $sourcedir/shared/include $sourcedir/shared/mpp/include

else

# execute mkmf to create makefile
cppDefs="-Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DINTERNAL_FILE_NML -DOVERLOAD_C8 {{compile_flags}}"
$mkmf  -a $sourcedir -t $template -p $executable -c "$cppDefs" $pathnames $sourcedir/shared/include $sourcedir/shared/mpp/include

fi

make

# $mkmf $make_flags -a $source_dir  -p fms_moist.x -t   $template \
#     -c "-Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DINTERNAL_FILE_NML -DOVERLOAD_C8" $pathnames $sourcedir/shared/mpp/include $sourcedir/shared/constants $sourcedir/include
#     make

if [ $? != 0 ]; then
   echo "ERROR: mkmf failed for $executable"
   exit 1
fi

# --- execute make ---
make $executable
if [ $? != 0 ]; then
    echo "ERROR: make failed for $executable"
    exit 1
fi

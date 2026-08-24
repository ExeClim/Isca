#!/usr/bin/env bash
# Compiles the Isca Model

# 0. Source the environment file to load appropriate variables
source {{ env_source }}

# 1. Configuration
hostname=`hostname`
compiler=${GFDL_MKMF_TEMPLATE:-ia64}
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

#netcdf_flags=`nf-config --fflags --flibs`

#ulimit -s unlimited # Set stack size to unlimited
export MALLOC_CHECK_=0

# Number of parallel make jobs. Defaults to the number of available cores,
# capped at 8 so we don't monopolise shared/multi-user machines by default.
# Override explicitly with e.g. ISCA_MAKE_JOBS=32 on a dedicated machine.
default_make_jobs=$(nproc 2>/dev/null || echo 1)
if [ "$default_make_jobs" -gt 8 ]; then
  default_make_jobs=8
fi
make_jobs=${ISCA_MAKE_JOBS:-$default_make_jobs}
echo "Using $make_jobs parallel make job(s) (set ISCA_MAKE_JOBS to override)"

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
cppDefs="-Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DINTERNAL_FILE_NML {{compile_flags}}"
$mkmf  -a $sourcedir -t $template_debug -p $executable -c "$cppDefs" $pathnames $sourcedir/shared/include $sourcedir/shared/mpp/include

else

# execute mkmf to create makefile
cppDefs="-Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DINTERNAL_FILE_NML ${CDEFS} {{compile_flags}}"
$mkmf  -a $sourcedir -t $template -p $executable -c "$cppDefs" $pathnames $sourcedir/shared/include $sourcedir/shared/mpp/include

fi

make -j"$make_jobs"

# $mkmf $make_flags -a $source_dir  -p fms_moist.x -t   $template \
#     -c "-Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DINTERNAL_FILE_NML -DOVERLOAD_C8" $pathnames $sourcedir/shared/mpp/include $sourcedir/shared/constants $sourcedir/include
#     make

if [ $? != 0 ]; then
   echo "ERROR: mkmf failed for $executable"
   exit 1
fi

# --- execute make ---
make -j"$make_jobs" $executable
if [ $? != 0 ]; then
    echo "ERROR: make failed for $executable"
    exit 1
fi

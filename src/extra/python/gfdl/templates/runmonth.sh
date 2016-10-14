#!/usr/bin/env bash
# Run a single month

module purge
source {{ srcdir }}/src/extra/loadmodule
module list

{% if experiment_restart %}
  echo 'Experiment restart is file: ', {{ experiment_restart }}
{% endif %}

ulimit -s unlimited

debug={{ run_idb }}                                     # logical to identify if running in debug mode or not

cd {{ rundir }}

# copy and extract the restart information
{% if restart_file %}
cd INPUT
cp {{ restart_file }} res
cpio -iv < res
cd {{ rundir }}
{% endif %}

export MALLOC_CHECK_=0

cp {{ execdir }}/fms_moist.x fms_moist.x

if [ $debug == True ]; then
   echo "Opening idb for debugging"
   idb -gdb  fms_moist.x
else
  mpirun  -np {{ num_cores }} fms_moist.x
fi

err_code=$?
if [[ $err_code -ne 0 ]]; then
	exit $err_code
fi

# combine output files
echo Month {{ month }} complete, combining nc files

if [ {{ num_cores }} > 1 ]; then
 for ncfile in `/bin/ls *.nc.0000`; do
    {{ execdir }}/mppnccombine.x $ncfile
    if [ $? == 0 ]; then
        rm -f "${ncfile%.*}".????
    fi
 done
fi



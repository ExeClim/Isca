#!/usr/bin/env bash
# Run a single month

rundir={{ rundir }}  # change this if you're rerunning from the output directory

source {{ env_source }}

ulimit -s unlimited

debug={{ run_idb }}                                     # logical to identify if running in debug mode or not

cd $rundir

export MALLOC_CHECK_=0

cp {{ execdir }}/{{ executable }} {{ executable }}

if [ $debug == True ]; then
   echo "Opening gdb for debugging"
   exec gdb  {{ executable}}
else
  # defaults to mpirun -> export EXECUTION_TYPE=APRUN for Cray systems 
  if [[ -z "${EXECUTION_TYPE}" ]] || [ "${EXECUTION_TYPE^^}" = "MPIRUN" ]; then
    exec nice -{{nice_score}} mpirun {{mpirun_opts}} -n {{ num_cores }} {{ execdir }}/{{ executable }}
  elif [ "${EXECUTION_TYPE^^}" = "APRUN" ]; then
    exec nice -{{nice_score}} aprun {{mpirun_opts}} -n {{ num_cores }} {{ execdir }}/{{ executable }}
  fi
fi

err_code=$?
if [[ $err_code -ne 0 ]]; then
	exit $err_code
fi

# # combine output files
# echo Month {{ month }} complete, combining nc files

# if [ {{ num_cores }} > 1 ]; then
#  for ncfile in `/bin/ls *.nc.0000`; do
#     {{ execdir }}/mppnccombine.x $ncfile
#     if [ $? == 0 ]; then
#         rm -f "${ncfile%.*}".????
#     fi
#  done
# fi



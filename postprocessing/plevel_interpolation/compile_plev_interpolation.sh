#!/usr/bin/env bash
#Script that compiles the plev.x executable

cd ./exec

source $GFDL_BASE/src/extra/env/$GFDL_ENV

../bin/mkmf -p plev.x -t $GFDL_BASE/src/extra/python/isca/templates/mkmf.template.ia64 -c "-Duse_netCDF" -a ../src ../src/path_names ../src/shared/mpp/include ../src/shared/include

make -f Makefile

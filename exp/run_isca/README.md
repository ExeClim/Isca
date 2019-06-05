## Running Isca from the command line

This example contains a script that will allow you to run Isca from the command line.  It uses the current source located at $GFDL_BASE and runs with an input.nml and diag_table you specify.  By default, these are looked for in the 'input' directory.  All other files in the input directory are also added the run e.g. if your run requires ozone or CO2 files, symlink or copy them into the input directory.

The only required argument for running from the command line is the experiment `name`, and all diagnostic output from the run will be stored at `$GFDL_DATA/<name>/run0001`.

Some examples:

* Run experiment "test1" with all the configuration files stored in `input/` directory, on 16 cores.

```$ ./isca test1```

* Run experiment "test2" as for "test1", but recompile the source before running.

```$ ./isca -c test2```

* Run experiment "test3" with an initial state from a previous run.

```$ ./isca test3 --restart restart/res0005.tar.gz```

* Run experiment "test4" with custom input.nml, diag_table and change the number of cores over which the model is distributed.

```$ ./isca --namelist my_input.nml --diag my_diag -n 16 test4```



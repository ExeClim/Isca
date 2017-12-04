# Running the Isca Trip tests

The idea of the `trip test` is to check whether any new code you have added to Isca changes the results of any of the test cases. This is to guard against new code having unexpected consequences.

The trip test code runs through each of the test cases in `exp/test_cases` and runs each of them with two different commit IDs, the first being the old code you knew worked, and the second being the new commit you want to test. 

## Basic examples

The trip tests are run using the command-line interface. The simplest example call is
```./trip_test_command_line 155661f ec29bf3```
This would run all of the test cases, comparing the results using commit ID `155661f` with commit ID `ec29bf3`.

To specify a subset of the test cases, use the `-e` option:
```./trip_test_command_line 155661f ec29bf3 -e 'axisymmetric' 'bucket_model' 'frierson'```

To change the number of cores Isca is run on, use the `-n` option:
```./trip_test_command_line 155661f ec29bf3 -e 'axisymmetric' 'bucket_model' 'frierson' -n 4```

To specify the github repo used for the tests (e.g. your fork rather than `Execlim/Isca`), use the `-r` option:
```./trip_test_command_line 155661f ec29bf3 -e 'axisymmetric' 'bucket_model' 'frierson' -r git@github.com:sit23/Isca```


## Example output

Running the command
```./trip_test_command_line 155661f8c7945049cbac0dcf2019bb17fe7a6a8d ec29bf389cf5ac53b50b23c363040479a6392e52```

Produces the following summary output:

```
Results for all of the test cases ran comparing 155661f and ec29bf3 are as follows...
variable_co2_grey : pass
realistic_continents_fixed_sst : pass
realistic_continents_variable_qflux : pass
held_suarez : pass
variable_co2_rrtm : pass
bucket_model : fail
top_down_test : pass
frierson : pass
axisymmetric : pass
MiMA : pass
giant_planet : pass
Nightmare, some tests have failed
```

The bucket test is the only one that fails in this instance. This is because the bucket model formulation was changed between these two commits, and so it is expected that the results with the two commits will differ. However, we note that the other tests have their results unchanged, meaning this is a safe modification to the code. Any unexpected failures should be investigated before submitting a pull-request.
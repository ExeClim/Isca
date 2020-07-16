Trip Tests
==============

Summary
-------
The trip test program takes two commit IDs as input, and compares the results of Isca's test cases using the two versions of Isca defined by the commit IDs. A comparison is then made between the output from each of the test cases. The output of the program shows which test cases produced bit-wise identical results and which did not.


Why do we ask you to run the trip tests?
-----------------------------------------

We are very keen for people to contribute their modifications and enhancements back into Isca's main version. However, before we accept your modifications, we require you to run the trip tests. By running these tests, we can check whether your modifications changes the results of any of the test cases.

We want to check this for three principal reasons:

    1. If you add a new feature, we want to make sure it is turned off by default. This is so that an existing Isca user who doesn't know about the feature isn't forced to use it when they upgrade to the newest Isca version. This prevents confusion where two versions of Isca give different answers for the same experimental configuration.
    2. If your modification contains a bug fix, we want to know which test cases are affected. We also want to know which P/R resulted in the bug fix.
    3. You may think your modifications won't have changed Isca's results, but we want to make sure that's the case. The trip tests are very good at highlighting unexpected consequences of your changes, and help us keep to a high-standard of reproducibility.


How do I run the trip tests?
---------------------------------

To run the trip tests, navigate to ``${GFDL_BASE}/exp/test_cases/trip_test/``. From this directory the trip tests are run using the command-line interface. The simplest example call is
::
    ./trip_test_command_line 155661f ec29bf3
This would run all of the test cases, comparing the results using commit ID ``155661f`` with commit ID ``ec29bf3``.

The optional arguments are
::
    -e 'all' - Runs all test experiments
    -n 4     - Uses 4 cores to run Isca
    -r 'git@github.com:execlim/Isca' - Uses the online Isca repo to checkout commits from

To specify a subset of the test cases, use the ``-e`` option
::
    ./trip_test_command_line 155661f ec29bf3 -e 'axisymmetric' 'bucket_model' 'frierson'

To change the number of cores Isca is run on, use the ``-n`` option
::
    ./trip_test_command_line 155661f ec29bf3 -e 'axisymmetric' 'bucket_model' 'frierson' -n 4

To specify the github repo used for the tests (e.g. your fork rather than ``Execlim/Isca``), use the ``-r`` option
::
    ./trip_test_command_line 155661f ec29bf3 -e 'axisymmetric' 'bucket_model' 'frierson' -r git@github.com:sit23/Isca

**Please note that by default the trip tests will try to connect to the ExeClim/Isca repository using SSH**. If you use HTTPs please use the option ``-r https://github.com/ExeClim/Isca.git``.

Example output
--------------

Running the command
::
    ./trip_test_command_line 155661f8c7945049cbac0dcf2019bb17fe7a6a8d ec29bf389cf5ac53b50b23c363040479a6392e52

Produces the following summary output
::
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

The bucket test is the only one that fails in this instance. This is because the bucket model formulation was changed between these two commits, and so it is expected that the results with the two commits will differ. However, we note that the other tests have their results unchanged, meaning this is a safe modification to the code. Any unexpected failures should be investigated before submitting a pull-request.

FAQs
----

* **What happens if the trip tests fail for my modifications?** If some of the trip tests fail, please investigate why you think this might be. If the tests fail on an isolated test case that you haven't used before, please post the results in your PR and we can discuss what is best to do. We will accept some PRs that fail the trip tests, but we will need to understand why first.

* **When I submit a pull request, which two commit IDs should I use?** Please use the commit ID of the current HEAD of Isca's master branch as one, and the HEAD of your pull-requested branch as the other. *Please note that in order for this to work, you will have to use the* `-r` *option to select your own fork, as described above.*

* **How long do the trip tests take to run?** The code is re-compiled several times during the trip tests, so the majority of time will be spent compiling. Each test case is only run for a few days, so once the code has compiled it shouldn't take more than 30 mins or so, depending on how many cores you're using. You can increase the number of cores by using the ``-n`` option.

* **I am getting an error message saying that git cannot find my commit ID - how do I fix this?** If you are specifying commit IDs that are present on your fork but not in the ExeClim version of Isca, you will need to specify your own fork as the repository. To do this, use the ``-r`` option, e.g. ``-r git@github.com:sit23/Isca``. Note that you can specify the repository location using the HTTPS option, or the SSH option depending on your preference.

* **Do I have to connect to GitHub to run the trip tests, or can I do it all locally?** It is possible to run the trip tests locally. Use the option ``-r PATH_TO_MY_LOCAL_REPO`` and it should work fine. 

* **Which of the test cases does the trip tests run**? A list of the available test cases can be found in the ``list_all_test_cases_implemented_in_trip_test`` function in ``trip_test_functions.py``.

* **I haven't setup Socrates - does that matter?** No - just post the results you can obtain.

* **If I add a new test case to the** ``test_cases`` **folder, how do I add it to the trip tests?** In the ``get_nml_diag`` function in ``trip_test_functions.py``, you will see several similar calls to import the namelist dictionary from each of the existing test cases. Just copy this syntax and edit it for your test case. You can then add the name of your test case to the list in the ``list_all_test_cases_implemented_in_trip_test`` function, and it should be run by default.

* **I have defined a new** ``codebase`` **object for my test case. How can I get the trip test to select it?** In the ``conduct_comparison_on_test_case`` the codebase object is selected. You can edit this section accordingly. 

* **Why don't you compare the results of the test cases to some standard output, rather than comparing the results from two commits?** We do this because different compilers and different hardware can produce slightly different results for a given test case. Therefore comparing to standard output is a much more difficult test for users to pass. We therefore compare the results of two different Isca versions *on the same hardware and with the same compilers* to make the tests easier to understand.

Authors
----------
This documentation was written by Stephen Thomson, peer reviewed by X, and quality controlled by Y.
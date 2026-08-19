Trip Tests
==============

Summary
-------
The primary purpose of the trip tests is to identify if proposed changes to the Isca master result in it giving different results for a range of test cases. 

The trip test program takes two commit IDs as input, and compares the results of Isca's test cases using the two versions of Isca that are defined by the commit IDs. A comparison is then made between the output from each of the test cases. The output of the program shows which test cases produced bit-wise identical results and which did not.

**Please note** that this guide assumes some knowledge of the `git` version control system. If you get stuck on this, there are lots of free guides online to using `git`, e.g. `GitHub's Git Handbook <https://guides.github.com/introduction/git-handbook/>`_.

Why do we ask you to run the trip tests?
-----------------------------------------

We are very keen for people to contribute their modifications and enhancements back into Isca's master - for instructions on how to contribute, see the :doc:`../contributing` . As part of this process, before we accept your modifications, we require you to run the trip tests. By running these tests, we can check whether your modifications changes the results of any of the test cases.

We want to check this for three principal reasons:

    1. If you add a new feature, we want to make sure it is turned off by default. This is so that an existing Isca user who doesn't know about the feature isn't forced to use it when they upgrade to the newest Isca version. This prevents confusion where two versions of Isca give different answers for the same experimental configuration.
    2. If your modification contains a bug fix, then this will likely alter Isca's output. If this is the case, we want to know which of test cases are affected by your bug fix.
    3. You may think your modifications won't have changed Isca's results, but we want to make sure that's the case. The trip tests are very good at highlighting unexpected consequences of your changes, and help us keep to a high-standard of reproducibility.


How do I run the trip tests?
---------------------------------
The first thing to do before you run the trip tests is to decide on which two commit IDs you would like to compare, and then obtain their commit IDs. When you submit a pull request, we ask you to use the current HEAD commit of Isca's master branch as one of the commits, and the HEAD of your pull-requested branch as the other. If you are using the trip tests as a testing tool, then you can choose any two commits you like.

To run the trip tests, navigate to ``${GFDL_BASE}/exp/test_cases/trip_test/``. From this directory the trip tests are run using the command-line interface. The simplest example call is
::
    ./trip_test_command_line 155661f ec29bf3

This would run all of the test cases, comparing the results using commit ID ``155661f`` with commit ID ``ec29bf3``. These commit IDs are for illustration purposes only. Please use two commit IDs relevant for your testing (i.e. the current HEAD of the Isca master and the HEAD of your branch).

The optional arguments and their defaults are
::
    -e 'all' - Runs all test experiments
    -n 4     - Uses 4 cores to run Isca
    -r 'git@github.com:ExeClim/Isca' - Uses the online Isca repo to checkout commits from

To specify a subset of the test cases, use the ``-e`` option
::
    ./trip_test_command_line 155661f ec29bf3 -e 'axisymmetric' 'bucket_model' 'frierson'

**Please note** that when you submit a pull request, you need to show results for *all* the test cases. Even if you think your modifications will not affect another test case, the trip tests are an easy way of making sure. For nearly all pull requests, we will require all the trip tests to pass. See :ref:`FAQs`. The option for running a subset of test cases is intended for subsequent testing if, for example, you find only one test fails, and you want to spend more focussed time investigating which commit caused the failure.

To change the number of cores Isca is run on, use the ``-n`` option.
::
    ./trip_test_command_line 155661f ec29bf3 -e 'axisymmetric' 'bucket_model' 'frierson' -n 4

To specify the github repo used for the tests (e.g. your fork rather than ``ExeClim/Isca``), use the ``-r`` option
::
    ./trip_test_command_line 155661f ec29bf3 -e 'axisymmetric' 'bucket_model' 'frierson' -r git@github.com:YOUR_GITHUB_USERNAME/Isca

**Please note that by default the trip tests will try to connect to the ExeClim/Isca repository using SSH**. If you use HTTPs please use the option ``-r https://github.com/ExeClim/Isca.git``.

How do I meet the trip test requirements for a pull request?
-----------------------------------------------------------------------

1. Identify the commit ID of the HEAD of Isca's master branch, and the HEAD of your pull-requested branch.
2. Run the trip test for all of the test cases (default) using these two commit IDs.
3. Post the summary output in your pull request (see below for example summary output).
4. Mark your pull request with the label ``trip tests passing`` if all of the tests pass. Mark your pull request with ``changes previous results`` if some of the trip tests fail.

If the tests all pass then your pull request should be accepted. If some of the tests fail then *please investigate why they have failed*. In general we will require all the tests to pass for a pull request to be accepted. If, however, you have fixed a bug and this causes a test case to fail, or similarly you think the results should change, then please make your case for this in the pull request.

**Please note** that you are required to run the trip tests when submitting a pull request that changes Isca's model code. Updates to the python front-end will be excluded from this requirement.


Example output
--------------

Running the command
::
    ./trip_test_command_line 155661f ec29bf3

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

The bucket test is the only one that fails in this instance. This is because the bucket model formulation was changed between these two commits, and so it is expected that the results with the two commits will differ. However, we note that the other tests have their results unchanged, meaning this is a safe modification to the code, and the pull request is likely to be approved. Any unexpected failures should be investigated before submitting a pull-request.

.. _FAQs:

FAQs
----

* **What happens if the trip tests fail for my modifications?** If some of the trip tests fail, please investigate why you think this might be. If the tests fail on an isolated test case that you haven't used before, please post the results in your pull request and we can discuss what is best to do. We will accept some pull requests that fail the trip tests, but we will need to understand why first.

* **When I submit a pull request, which two commit IDs should I use?** Please use the commit ID of the current HEAD of Isca's master branch as one, and the HEAD of your pull-requested branch as the other. *Please note that in order for this to work, you will have to use the* `-r` *option to select your own fork, as described above.*

* **How do I find a commit ID?** The commit ID is a 40-character string of letters and numbers associated with the commit, which you can find on GitHub under `Commits` or using `git log` on the command line. You require only the first 7 characters of the ID for these tests.

* **How long do the trip tests take to run?** The code is re-compiled several times during the trip tests, so the majority of time will be spent compiling. Each test case is only run for a few days, so once the code has compiled it shouldn't take more than 30 minutes or so, depending on how many cores you're using. You can increase the number of cores by using the ``-n`` option.

* **I am getting an error message saying that git cannot find my commit ID - how do I fix this?** If you are specifying commit IDs that are present on your fork but not in the ExeClim version of Isca, you will need to specify your own fork as the repository. To do this, use the ``-r`` option, e.g. ``-r git@github.com:YOUR_GITHUB_USERNAME/Isca``. Note that you can specify the repository location using the HTTPS option, or the SSH option depending on your preference.

* **Do I have to connect to GitHub to run the trip tests, or can I do it all locally?** It is possible to run the trip tests locally. Use the option ``-r PATH_TO_MY_LOCAL_REPO`` and it should work fine. 

* **Which of the test cases does the trip tests run**? A list of the available test cases can be found in the ``list_all_test_cases_implemented_in_trip_test`` function in ``${GFDL_BASE}/exp/test_cases/trip_test/trip_test_functions.py``.

* **I haven't setup Socrates - does that matter?** Socrates is an optional radiation scheme that Isca can be complied with, but is not by default (see :doc:`../modules/socrates`). So you are not required to run the Socrates test case if you have not set it up - just post the results you can obtain.

* **If I add a new test case to the** ``test_cases`` **folder, how do I add it to the trip tests?** In the ``get_nml_diag`` function in ``trip_test_functions.py``, you will see several similar calls to import the namelist dictionary from each of the existing test cases. Just copy this syntax and edit it for your test case. You can then add the name of your test case to the list in the ``list_all_test_cases_implemented_in_trip_test`` function, and it should be run by default. *Please note* that if you add a new trip test that uses new features, it will fail the trip test because the existing Isca master won't be able to run it. This is OK, and you should highlight this when submitting the pull request.

* **I have defined a new** ``codebase`` **object for my test case. How can I get the trip test to select it?** In the ``conduct_comparison_on_test_case`` the codebase object is selected. You can edit this section accordingly. 

* **Why don't you compare the results of the test cases to some standard output, rather than comparing the results from two commits?** We do this because different compilers and different hardware can produce slightly different results for a given test case. Therefore comparing to standard output is a much more difficult test for users to pass. We therefore compare the results of two different Isca versions *on the same hardware and with the same compilers* to make the tests easier to understand.

Authors
----------
This documentation was written by Stephen Thomson, peer reviewed by Penelope Maher, and quality controlled by Ross Castle.
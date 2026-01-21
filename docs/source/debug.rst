Ideas for ‘debugging’ problems in Isca
======================================

When Isca runs crash it can be difficult to diagnose the problem. We suggest trying these steps below in order to work out where the crash is happening. Not all may be applicable to your problem.

Write Statements
----------------
Add in some write commands into the Fortran files you suspect are causing trouble. These can be either to print out values to check if they seem sensible, or periodic statements (e.g. every 100 lines or after every major function) so see how far the code is running. A quick and dirty way of printing variables directly to the terminal is to use a command in the code like this (note this will print every timestep without further conditions):

.. code-block:: fortran

    write(6,*) "<Some description>: ", <variable>

Alternatively variables can be written to files, such as this `example <http://www.python.org/>`_.

Running on 1 core
-----------------
A lot of issues on Isca can stem from parallelisation across multiple computing cores. If you run a job on 1 core and it works, you know that the problem lies here. Quite often it can be the allocation of the start/end integers ``is, ie, js, je``

Column Model
------------
Using the column model can be useful to test that parameterisations like convection and radiation are working properly as you don’t have horizontal dynamics complicating things.

One Change At A Time
------------------
Generally when building a model, it is a good idea to create it iteratively. Start from somewhere you know works, for example one of our test cases, and make 1 or 2 changes at a time towards your ideal set up. Then you will know what change is causing the crash.

Diagnostics
-----------
Be sure to output diagnostics that can help you diagnose the problem. That said, sometimes turning off diagnostics can stop model crashes. If this is the case, your problem is likely data related. 

Recompile
---------
Ensure that you are recompiling the Fortran each time.

Timestep
--------
The usual quickest fix when you have a model run that is crashing part way through is to increase the timestep. Your timestep must obey the `CFL criterion <https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition>`_.

Check you’ve allocated large enough variable size
-------------------------------------------------
Particularly important when using an uncommon diagnostic. 


Increase the stack size
-----------------------
Ensure there is enough temporary memory to hold the run and data created.


Authors
-------
This documentation was put together by the Isca team from their experience with using Isca. 
Last updated 08/10/2021

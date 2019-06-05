#How to setup Isca to run with IDB

To run Isca with the Intel debugger, three changes should be made:

  1. Set NCORES=1
  2. set the 'debug=True' argument in the 'cb.compile'
  3. Add the 'run_idb=True' option to the first instance of exp.run

These changes have been made in the `frierson_debug.py` example in this folder.

# How to use IDB

Assuming you have managed to run the code with the ```run_idb=True``` flag and that your running in **serial** and not parallel (ie one core) . Here are some of the key commands you need to run idb.

1. The IDB GUI opens with the ```(idb)``` line ready for input.
2. The first thing you want to do is put some break points in (think of these as the points you want to stop the code). Typically I will add more than one.
 ```break filename.F90:line_number```
 3. Then run the code using ```run```
 4. The code will stop when it reached the first breakpoint.
 5. To print a variable use ```p var_name``` for a single value or ```p var_name(1,1)``` for the single value of a 2d array. Try not to print too many variables (IDB is smart enought to tell you no when relevant but this can cause you problems). To print a slice of a small-ish 2d array use ```p var_name(1,:)```
 5. Use ```n``` to move to the next line, ```c``` to move to the next breakpoint.
 6. To delete the first breakpoint use ```d 1```
 7. To see where your breakpoints are and how many times you have hit them use: ```info break```
 8. To stop the code use ```stop```
 9. To exit IDB use ```quit```
 
 The game them becomes how to place the break points in location that best suit your purpose.
 
 Have a look at the manual for more informaation https://software.intel.com/sites/default/files/m/8/4/c/5/7/6364-idb_manual.pdf

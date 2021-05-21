Isca Beginner's Guide
====================

Below is a list of reading and activites that will help you get comfortable using Isca. Be assured that your supervisor/tutor will not expect you to be fluent with this when you start using Isca, but it will help to have an idea of what to expect when you start running the model.

This document is essentially a suggestion of signposts. With this kind of work, self-study and initiative is very important. It is up to you to go and research the topics until you feel comfortable.

Note: users who are more familiar with scientific computing may find this document a little longwinded and may be better off looking at the `ReadMe <https://github.com/ExeClim/Isca#readme>`_ or this `guide on how to run Isca experiments <https://github.com/ExeClim/ictp-isca-workshop-2018/blob/master/experiments/isca_help_ictp.pdf>`_.

Essentials
-------------

Users who will be using Isca to run simple planetary models under close guidance of a supervisor/teacher should learn about the following topics. Do not worry if it doesn’t make sense right now, you will understand more as you go on. 

ssh/terminals
^^^^^^^^^^^^^
To run Isca you will need to be using quite a powerful computer, most laptops will not suffice, especially for high resolution runs. You’ll likely run Isca on a university owned *workstation* or *supercomputer*. Your supervisor will inform you of the name of the computer. However, these computers are *remote* – that is you do not sit in front of them and login to them as you might be used to. You have to login to them from another computer, e.g. your personal laptop. This is called *SSH*, or *Secure Shell*, a protocol that enables two computers to communicate and share data.

To do this you will need a *terminal* on your own personal computer. Mac and linux computers will have a factory installed app for this called *terminal*. Windows users will need to download a piece of software, the most common is `PuTTY <https://www.putty.org/>`_. These are *ssh clients*. This software is text based - everything is done using text commands from the *command line*, we talk about this more later. 

To login to the *unix server*, all you need to do is open the *terminal* and run something like:

``ssh USER@computername.ex.ac.uk``

Your supervisor will tell you the precise command. 

You may also need to be connected to a *VPN (Virtual Private Network)*. For Exeter, if you are not on campus you will certainly need to use `this <http://www.exeter.ac.uk/it/howdoi/vpn/>`_. 

There are some interesting videos on youtube on how the SSH protocol works: e.g. `here <https://www.youtube.com/watch?v=qWKK_PNHnnA>`_, but you don’t need to understand how it works in order to use it. 

Unix Servers
^^^^^^^^^^^^
Most scientific computing is done on computers/servers that use an *unix* operating system. These do not have a *graphical user interface (GUI)*, everything is done with text commands. There are good guides on the internet, like `here <https://ubuntu.com/tutorials/command-line-for-beginners#1-overview>`_. 

If you are an apple/linux user, you can just open the terminal app on your laptop and practise there. If you use windows, it maybe easiest to use this `emulator <https://cocalc.com/projects?session=default>`_, then select the drop down arrow to the right of `New` and select `_Terminal.term`. Note because you’re using an emulator it may not be possible to follow the steps of the guide exactly, but try and get a feel of the commands they suggest (`cd, cd .., ls, mkdir, rm, pwd, mv, cp, less`).

Python (namelists)
^^^^^^^^^^^^^^^^^^
Isca is configured using *python* scripts, although the actual model is coded in *FORTRAN* (see :ref:`FORTRAN` in Advanced). `Python <https://www.python.org>`_ is a powerful high-level programming language, similar to Matlab – but far better.

At this stage you don’t need to be able to code in python, as Isca has prewritten scripts called *test cases* which you can just edit in order to change the model set up. Editing will require changing one or more *namelists*. The namelists are actually part of the FORTRAN code but we use *python packages* to pass these variables between the python script and the FORTRAN model. 

For example say that in the model the value for the CO2 concentration was 300, and you wanted to make it 600, all you would need to do is change:

``'co2_conc' : 300.`` to ``'co2_conc' : 600.`` in your text editor (see :ref:`below<Text Editors>`).

If you haven’t used python before and want to become more familiar, there are hundreds of tutorials (e.g. `here <https://docs.python.org/3/tutorial/>`_) and videos. To practise you can use a *python notebook* on something like `Colab <https://colab.research.google.com>`_, or download software like `Anaconda <https://anaconda.org>`_. Python is comprised of the basic python *packages* and then additional *libraries* you have to install and *import*. In the future it may be useful to have python *environments* (see :ref:`Conda` in Intermediate). 

Text Editors
^^^^^^^^^^^^
You will need to be able to edit text based files, this includes code files like python and other scripts relevant to scientific computing (see :ref:`shell/bash<Shell Scripts>` in Advanced).
Every unix server will have *vim*, which is a powerful text editor, but it does take some getting used to. For example, unlike in 'MS Word' or editors with *GUIs*, you cannot click to place your cursor in vim.

You can try out this `tutorial <https://www.openvim.com>`_. Note: I have found that you can usually also use the `up/down/left/right` keys to navigate the text, not just `h/j/k/l` as is stated here. You can practise freely using terminal (Mac/Linux) or the unix `emulator <https://cocalc.com/projects?session=default>`_.

To open vim type on the command line:

``vim test.txt`` (if you’re lazy like me ``vi test.txt`` also works).

This opens a new text file in vim. You can edit an existing file in exactly the same way:

``vim alreadyexisted.txt``

vim is not the only option! *emacs* is a similar editor which is guaranteed to be installed. emacs is opened in the same way as vim.

Perhaps a better option is *gedit*. This is a simple text editor with a GUI which is usually installed on servers (it is on the GV machines). **This is what I would recommend using as a beginner, if available**. It’s a little clunky but more intuitive to use then the previous options. It is opened exactly the same as vim/emacs. In some cases you may need to set up *X11 forwarding* (see :ref:`X11 forwarding` in Intermediate).

As you get more comfortable with this scientific computing, you will likely find that you prefer a different text editor with a GUI which is far more user friendly. However, it will require a bit of setting up. Talk to your supervisor/research group about what they use and how they got it to work.

Intermediate
---------------

If the user will be running multiple experiments on their own and analysing the output, the following will likely be useful to them:

Isca Structure
^^^^^^^^^^^^^^
It may be useful for you to have a rough idea on how Isca works. The best way to do this is to look through the Isca `documentation <https://execlim.github.io/Isca/latest/html/>`_, especially the Isca structure page. You can also skim through the `source code <https://github.com/ExeClim/Isca/tree/master/src>`_, to get an idea of what files there are – there are lots, but you don’t need to worry about how they all work so do not be intimidated!

Conda
^^^^^
As mentioned earlier in the Python section, often Python libraries have to be installed, and you’ll need different libraries depending on what you’re doing. Python *environments* are very useful as loading them will load all the libraries you need for a given task. For example, there is an isca environment which is set up during the Isca installation, which has all the relevant python modules for running Isca. See `here <https://conda.io/projects/conda/en/latest/user-guide/concepts/environments.html>`_ for more details. 

Workstations
^^^^^^^^^^^^
Some terminology things to be aware of when running on servers/workstations:

- Workstations (for example the ‘GV machines’ at Exeter) have *cores* which are like groups of processors. So when running Isca you can run on a number of cores, generally the more cores the faster. Due to the way Isca works, you can only run on a number of cores that is a power of 2 (1, 2, 4, 8, 16, 32). We usually run at 8 or 16. 
- Unix has a feature called *screen* which allows you to leave something running and logout of a computer. When you’re logged in, simply type ``screen`` on the command line and a screen will start. You can then press ``CTRL+A+D`` to detach from the screen but leave your job running. Then you can log out of the computer. See `here <https://www.tecmint.com/screen-command-examples-to-manage-linux-terminals/>`_ for commands about reattaching, listing screens etc.
- Typing ``top`` on the command line will display a list of users/jobs that are happening at that time. This is useful to make sure you are not overloading the computer. For example, if you wanted a to run an 8 core job but the computer only had 4 cores free, you’d have to wait. 

X11 forwarding
^^^^^^^^^^^^^^
If you want to make plots and view them from a computer you have SSH’d into, you might need to set up some sort of *X11 forwarding*. It just allows images created in windows on another computer to appear as windows on your own computer.

Use software like `XQuartz <https://www.xquartz.org>`_ for macOS or Xming for Windows. You’ll also need to add the ``-Y`` or ``-X`` option to your ssh command (i.e. ``ssh –Y user@emps-gv1.ex.ac.uk``) . Getting it set up the first time may be a little tricky, but there is plenty of help available on google/your supervisor. 

netCDFs
^^^^^^^
Isca has to store the data it generates so that you can analyse it and make plots. The file type it uses is called a *netCDF* file which has a *.nc* suffix. For example, every month Isca can output a file called ``atmos_monthly.nc`` which contains all the variables asked for in the python run script (wind velocities, temperature, precipitation, etc). They are very useful for climate data because it allows variables to be stores on sets of *axis* like latitude, longitude, height* and time. This makes it easy to make plots and there are python libraries e.g. *netCDF4* which have many useful functions to make your life easier.

If you’re interested there is reams of documentation `here <https://www.unidata.ucar.edu/software/netcdf/docs/index.html>`_ but again, you don't need to understand it too much in order to use it.

*Note: In Isca’s case the ‘height’ axis is not measured in meters, but usually in `sigma pressure coordinates <https://glossary.ametsoc.org/wiki/Sigma_vertical_coordinate>`_.

Plotting/xarray
^^^^^^^^^^^^^^^
When Isca has finished it’s model run, you’ll want to look at the data created and analyse it and make plots. We have some scripts that will help get you started `here <https://github.com/ExeClim/ictp-isca-workshop-2018/tree/master/analysis>`_. These scripts are written using functions from python libraries called `xarray <http://xarray.pydata.org/en/stable/>`_, which is a very powerful way to work with datasets in python, and `matplotlib <https://matplotlib.org/2.0.2/api/pyplot_api.html>`_ which is a plotting library. You will need to install these libraries to a python environment to use them.

Transferring Files (SFTP/SCP)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Now you have made plots – or indeed any file you want to transfer between the computer you have SSH’d into and your own – you will need a way of transferring them. There are several ways of doing this.

*SFTP (SSH File Transfer Protocol)* is one, it will work on all operating systems and is the easiest for windows. One way of using SFTP is with an *SFTP client*, many are available. One of them is `Cyberduck <https://cyberduck.io>`_. It will require setting up but it is fairly straight forward. These clients tend to have a GUI so you can just drag and drop the files you want to transfer. It is also possible to view and transfer files using the native file browser if you're using Linux or macOS, using their built-in functions to connect via SFTP.

Other option is to use a command line function, for example ``scp``. This is a secure file copy protol, which uses SSH. The usage is simple, for example on the computer you want to transfer the file to, type:

``scp USER@COMPUTERNAME.ex.ac.uk:/path_to_file/file.png /path_to_destination/``

This uses the protocol to SSH into the computer with the file and copy it to the location specified on the RHS. Note to copy a directory you can use the ``-r`` (*recursion*) option. We also can use a ``.`` to copy to our current file location. 

``scp –r USER@COMPUTERNAME.ex.ac.uk:/path_to_directory/ ./``

See `here <https://www.ssh.com/ssh/scp/>`_ for more details.

Advanced
--------

Users who either intend to make changes to the Isca source code, or will use the model so often as to benefit from additional tools, should research the following:

Git
^^^
Git is a *version control software*, which allows you and every other user to have different copies of the Isca source code and modify it safely. Developers of Isca will have different *branches* on their own *fork*, which they can modify and improve. If the improvements are useful to everyone, the changes can be added to the `master copy <https://github.com/ExeClim/Isca>`_.

Here is a `video <https://www.youtube.com/watch?v=w3jLJU7DT5E>`_ about how git works. Here is a useful `cheat sheet <https://education.github.com/git-cheat-sheet-education.pdf>`_ on git commands.

Supercomputers
^^^^^^^^^^^^^^
You may be able to run Isca on a supercomputer, for example at Exeter we have *ISCA HPC (High Performance Computer)* - the same name get’s confusing. Your supervisor will help get you set up on this as they are a little more complicated, although usually faster. 

When you login to a supercomputer you are in fact logging in to a small *login node* which is not designed to run code. It is designed to allow you to submit your job to a *queue* which will then be run on the main computer (see :ref:`Slurm` below). Here is some `documentation for ISCA HPC <https://universityofexeteruk.sharepoint.com/sites/ExeterARC>`_, see the ISCA User Guide. 

Slurm
^^^^^
Submitting jobs to a queue requires you to use the supercomputers *workload manager*. ISCA HPC uses *Slurm*, but there is also *moab*. See here for a `slurm cheat sheet <http://www.physik.uni-leipzig.de/wiki/files/slurm_summary.pdf>`_. The important ones are ``sbatch`` and ``squeue``. 

FORTRAN
^^^^^^^
The actual Isca model is written in a coding language called *FORTRAN.90*. Therefor if you intend on modifying the source code, you’ll need to know a little FORTRAN. It is incredibly fast, but it has to be *compiled* before use (it is a *low level* language) and is slightly different from *high level* code. For example, you have to define variables before you can use them. There are plenty of FORTRAN tutorials around, e.g. `here <https://www.fortrantutorial.com>`_, however you will probably learn as you go by modifying the Isca code.

Shell Scripts
^^^^^^^^^^^^^
A *shell script* (``scriptname.sh``) is a useful tool if you have a series of command lines you have to write, especially if you do it often. For example, I have a shell script that transfers data from one server to another. The `example file <https://github.com/ExeClim/Isca/blob/master/exp/test_cases/isca_job.sh>`_ to submit a job to ISCA HPC is also a shell script. See `here <https://www.shellscript.sh>`_ for more details or google.

.bashrc Script (aliases)
^^^^^^^^^^^^^^^^^^^^^^^^
One particular shell script is your ``.bashrc`` script, see `here <https://www.journaldev.com/41479/bashrc-file-in-linux>`_. Your supervisor will set this up for you, as some Isca file locations need to be included in it. One very useful thing that you can set up in this script is *aliases*. This is where a text string is assigned to a command.

E.g. the line ``alias go_data='cd /scratch/USER/data_isca'`` will allow you to go to your data file location, just by typing ``go_data``.

Or the line ``alias i='source activate isca_env'`` will activate your ``isca`` python environment just by typing ``i``. 

Useful Links
------------

- `How to install isca <https://github.com/ExeClim/Isca/blob/master/ReadMe.md>`_
- Will Seviour's Scripts - useful code designed for `setting up Isca at Bristol university <https://github.com/wseviour/Bristol_Climate_Dynamics/blob/master/Isca_SOCRATES.md>`_ and `analysing data using a python notebook <https://github.com/wseviour/Bristol_Climate_Dynamics/blob/master/Anthropocene_Isca_analysis.ipynb>`_
- `ICPT workshop repo <https://github.com/ExeClim/ictp-isca-workshop-2018/tree/master/analysis>`_ (some lecture slides and analysis scripts).
- The 2018 `paper on Isca's release <https://gmd.copernicus.org/articles/11/843/2018/>`_
- The `Isca Website <https://execlim.github.io/IscaWebsite/index.html>`_

Authors
-------
This documentation was written by Ross Castle with input from the Isca team, notably Penny Maher, Denis Sergeev, Geoff Vallis and Will Seviour. It is hoped that this document will continue to be edited and improved, especially by masters and PhD students. 

Last updated 31/03/2021

Remote Access
==============

Summary
-------
This is a guide for how to edit remote files on a local text editor via port-
fowarding. By the end you should know: 

* how to quickly login to your remote server using an SSH config file
* how to edit remote files in a local text editor through ``rmate``
* how to edit remote python files in a local ``jupyter`` environment

Simplifying Logins and Port Fowarding
-------------------------------------
If you use a Unix-based operating system on your personal computer, you can make use
of an SSH config file to create shortcuts to your frequently used remote computers. We will take this one step further to simplify port-fowarding, a method which allows a user to redirect data from a specified remote host and port, through a secure tunnel, to a specified local port. Port-fowarding is helfpul because it will enable us to edit remote files locally.

As an example without any fancy tricks, let's set up an SSH tunnel that maps ``localhost`` port 3039 on my local machine to 8450 on my remote machine (the number is arbitrary as long as its between 1024 and 49150): ``$ ssh -l localhost:3039:host:8450 user@host``. You will then be required to enter in your password. This is cumbersome to repeat everytime we log in. Our goal will be to shorten the command to : ``$ ssh hostalias`` and without having to enter in your password. We give some instructions below:

**1.** In the home directory of your **local machine**, create a new directory called ``.ssh`` if it does not already exist. Navigate to this directory and create file called ``config``.
Put in the following contents (making sure to replace the text surrounded by ``**double astericks**`` with your own information)::

   Host **hostalias**
      Hostname **hostname**
      User **username**

By using an SSH config file, secure methods for copying (e.g. ``scp`` or ``sftp``) can now use the same host aliases. Let's say I want to copy a file from my local machine to my remove machine. This is now as simple as: ``$ scp localfile.txt hostalias:/path/to/directory``.

**2.** If your local machine is a Mac, you can eliminate the need to enter a password every time you want to log in by using an SSH key pair. To do this, on your **local machine** navigate to your ``~/.ssh`` directory and enter the following command to generate a set of RSA keys: ``$ ssh-keygen -t rsa``. You will then be prompted to supply a filename and a password. For the file name I recommend `id_rsa_hostalias` and for the password I recommend it to be the same as your remote machine's password. If you use a Mac, the operating system can use its internal keychain to remember your password, meaning you won't need to type it in every time you log in!

This command will generate 2 files:the one with the ``.pub`` extension is the "public key", the one without the ``.pub`` extension is the "private key". You keep your private 
key strictly on your local machine. You need to copy your public key to the remote machine you would like to use the key pair to log in to. First create the `~/.ssh` on your **remote machine**. Then from your **local machine** enter: ``$ scp ~/.ssh/id_rsa_hostalias.pub hostalias:~/.ssh/``. After copying the public key over to ``hostalias``, you need to create an ``authorized_keys`` file in your ``~/.ssh/`` on your **remote machine**.

**3.** Now from a terminal window on your **local machine**, you can try logging in to ``hostalias`` with the following command: ``$ ssh -i ~/.ssh/id_rsa_hostalias hostalias``. This will prompt you for the password you specified upon creating your key pair using ``ssh-keygen``. To always make use of your private key when logging in to ``hostalias``, add the following to your ``config`` file on your **local machine**::

   Hostname **hostname**
      User **username**
      LocalForward 8450 localhost:3039
      IdentityFile **~/.ssh/id_rsa_hostname**
      UseKeychain yes
      AddKeysToAgent Yes

Congratulations!


Edit Remote Files Locally
-------------------------
As an alternative to remote-based text editors such as ``vi`` and ``emacs``, we can
use port-fowarding to set up a local-based text editor like ``Atom`` which includes
features such as syntax highlighting and code completion.

Edit Remote python Files in a ``jupyter`` Environment
-----------------------------------------------------
The ``jupyter`` environment is a great environent for data exploration and integrating
your figures inline with your code.


References
----------
..
   Add relevant references. This is done in 2 steps:
   1. Add the reference itself to docs/source/references.rst
   2. Insert the citation key here, e.g. [Vallis2017]_
   
   See the Contributing guide for more info.

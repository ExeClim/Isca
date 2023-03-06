Remote Access
=============

Overview
--------
This is a guide for how to edit remote files on a local text editor via port-
fowarding. By the end you should know: 

* how to quickly login to your remote server using an SSH config file
* how to edit remote files in a local text editor through ``rmate``
* how to edit remote python files in a local ``jupyter`` environment

Simplifying logins and port-fowarding
-------------------------------------
If you use a Unix-based operating system on your personal computer, you can make use
of an SSH config file to create shortcuts to your frequently used remote computers. We will take this one step further to simplify port-fowarding, a method which allows a user to redirect data from a specified remote host and port, through a secure tunnel, to a specified local port. Port-fowarding is helfpul because it will enable us to edit remote files locally.

As an example without any fancy tricks, let's set up an SSH tunnel that maps ``localhost`` port 3039 on my local machine to 3039 on my remote machine (the number is arbitrary as long as its between 1024 and 49150): ``$ ssh -l localhost:3039:host:3039 user@host``. You will then be required to enter in your password. This is cumbersome to repeat everytime we log in. Our goal will be to shorten the command to : ``$ ssh hostalias`` and without having to enter in your password. We give some instructions below:

Set up an SSH config file
^^^^^^^^^^^^^^^^^^^^^^^^^
In the home directory of your **local machine**, create a new directory called ``.ssh`` if it does not already exist. Navigate to this directory            and create file called ``config``.
Put in the following contents (making sure to replace the text surrounded by ``**double asterisks**`` with your own information)::

   Host **hostalias**
      Hostname **hostname**
      User **username**

By using an SSH config file, secure methods for copying (e.g. ``scp`` or ``sftp``) can now use the same host aliases. Let's say I want to copy a file from my local machine to my remove machine. This is now as simple as: ``$ scp localfile.txt hostalias:/path/to/directory``.

Set up an SSH key pair
^^^^^^^^^^^^^^^^^^^^^^
If your local machine is a Mac, you can eliminate the need to enter a password every time you want to log in by using an SSH key pair. To do this, on your **local machine** navigate to your ``~/.ssh`` directory and enter the following command to generate a set of RSA keys: ``$ ssh-keygen -t rsa``. You will then be prompted to supply a filename and a password. For the file name I recommend ``id_rsa_hostalias`` and for the password I recommend it to be the same as your remote machine's password. If you use a Mac, the operating system can use its internal keychain to remember your password, meaning you won't need to type it in every time you log in!

This command will generate 2 files:the one with the ``.pub`` extension is the "public key", the one without the ``.pub`` extension is the "private key". You keep your private key strictly on your local machine. It is a good idea to change the permissions of the private key file. You need to copy your public key to the remote machine you would like to use the key pair to log in to. First create the ``~/.ssh`` on your **remote machine**. Then from your **local machine** enter: ``$ scp ~/.ssh/id_rsa_hostalias.pub hostalias:~/.ssh/``. After copying the public key over to ``hostalias``, you need to create an ``authorized_keys`` file in your ``~/.ssh/`` on your **remote machine**.

Now from a terminal window on your **local machine**, you can try logging in to ``hostalias`` with the following command: ``$ ssh -i ~/.ssh/id_rsa_hostalias hostalias``. This will prompt you for the password you specified upon creating your key pair using ``ssh-keygen``. To always make use of your private key when logging in to ``hostalias``, add the following to your ``config`` file on your **local machine**::

   Hostname **hostname**
      User **username**
      LocalForward 3039 localhost:3039
      IdentityFile **~/.ssh/id_rsa_hostname**
      UseKeychain yes
      AddKeysToAgent Yes

You should not be able to log in simply by typing ``$ ssh hostalias``. Congratulations!




Edit remote files in a local text editor using ``rmate``
--------------------------------------------------------
As an alternative to remote-based text editors such as ``vi`` and ``emacs``, we can
use port-fowarding to set up a `local-based text editor like Atom <https://atom.io>`_ which includes features such as syntax highlighting and code completion. For instructions to install ``rmate`` on your **local machine**, `click here <https://github.com/textmate/rmate>`_. Then to specifically use Atom to edit remote files, `click here <https://atom.io/packages/remote-atom>`_. You will need to add the following line to your ``~/.ssh/config`` file: ``RemoteForward 52698 localhost:52698``.


Edit remote python files in a ``jupyter`` environment
-----------------------------------------------------
The ``jupyter`` environment is a great environent for data exploration and integrating
your figures inline with your code. To open your first Jupyter notebook, log in to your **remote machine** and type: ``$ jupyter lab --no-browser --port=3039``. 

To make it even quicker, you can type the following on your **local machine**: ``$ssh remotehost "jupyter lab --no-browser --port=3039``. This should function because of all the work we put in during the port forwarding section. To shorten this command, add an alias to your ``~/.bashrc`` file on your **local machine**. I personally use the alias ``rjlab``.

Recap
-----
Your final ``~/.ssh/config`` file should look like this (making sure to replace the text surrounded by ``**double asterisks**`` with your own information)::

   Host **hostalias**
      Hostname **hostname**
      User **username**

   Hostname **hostname**
      User **username**
      LocalForward 3039 localhost:3039
      RemoteForward 52698 localhost:52698
      IdentityFile **~/.ssh/id_rsa_hostname**
      UseKeychain yes
      AddKeysToAgent Yes

**Final Notes:** 

* Do not use this with VPN, use ithome aka hashbang as proxy. If the connection is interrupted you can still reconnect, assuming the jupyter process is still running. But make sure not to leave zombie jupyter processes with open ports on remote hosts!

* Remember the port numbers chosen are arbitary. If you choose the same number as someone else on your network, their files may open up on your computer and vice versa!

Authors
-------
This documentation was written by Brett McKim, peer reviewed by Denis Sergeev, and quality controlled by Ross Castle.

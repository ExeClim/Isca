Contributor's Guide
===================

Isca takes contributions using pull requests on `GitHub <https://github.com/execlim/isca/pulls>`_.

Creating a Github fork
----------------------
1. Fork Isca on Github, so a copy of it appears in your Github profile and has a remote URL like `https://github.com/your_user_name/Isca <https://github.com/your_user_name/Isca>`_.
2. Git-clone it to your computer using the new URL (or change the remote of your existing folder).
3. It is better to create a new branch (e.g. :code:`git checkout -b your_new_branch`). If your contribution is relatively small, such as fixing one small bug or typo, you can skip this step.
3. Make changes, commit them and push to *your* remote. (Which again, should have your URL. You can check it by typing :code:`git remote -v`).
4. Go to your Isca page on Github and create a pull request (a prompt for this should show up on top of the page).


Contributing to the documentation
---------------------------------
We welcome contributions to improve the documentation of Isca.

1. Please make yourself familiar with the `reST formatting guide <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_.
2. Follow the steps 1-3 above to create a fork.
3. Make changes
   a. If you are adding a new section of documentation, make a copy of the template file (:code:`docs/source/_template.rst`), assign it an appropriate name and place it to the appropriate location within :code:`docs/source/`. Add a link to the newly created file to `index.rst` within the same folder.
   b. If you are changing an existing file, proceed.
4. Once you are done, remove all the comments from your .rst file if you are using the template.
5. Do not forget to add any relevant references to papers or textbooks.
   a. Create a new entry in :code:`docs/source/references.rst`, following the formatting style of the existing entries. For example,
   .. code-block::
      .. [VallisEtAl2018] Vallis, G. K. and Colyer, G. and Geen, R. and Gerber, E. and Jucker, M. and 
                 Maher, P. and Paterson, A. and Pietschnig, M. and Penn, J. and Thomson, S. I., 2018:
                 Isca, v1.0: a framework for the global modelling of the atmospheres of Earth and 
                 other planets at varying levels of complexity. *Geoscientific Model Development*,
                 **11(3)**, 843-859,
                 doi: `10.5194/gmd-11-843-2018 <https://doi.org/10.5194/gmd-11-843-2018>`_.

   b. Add the citation to your docs page by using the relevant citation key (note the underscore symbol at the end). For example:
   .. code-block::
      We use the Isca model ([VallisEtAl2018]_)
6. Create a pull request and wait for the Isca team to review it.

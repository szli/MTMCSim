This readme is for the developers.

Note that to generate figures, we need Graphviz and it needs to be in the PATH.

Do not inlucde cross-referenced source code in the output.

To publish, generate the docs and push to gh-pages branch. More specificially, go to the html directory, do the folling:
git clone https://github.com/szli/MTMCSim.git .
git checkout origin/gh-pages -b gh-pages  (the files from master will disappear)
git branch -d master

Then generate files, commit and push.

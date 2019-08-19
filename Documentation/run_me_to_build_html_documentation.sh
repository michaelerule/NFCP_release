#!/usr/bin/env bash

echo "This script should be run with its location as the working directory"

TITLE='Neural Field Cox Process package'
FILENAME='NFCP'
DIRNAME='_auto'

# Clean auto-generated python directory
rm -rf ./$DIRNAME
mkdir ./$DIRNAME

# Clean old .rst files 
mkdir ../._obsolete_rst
mv *.rst ../._obsolete_rst

# Extract matlab comments into a dummy python class for spinx

# old approach: everything in one file
#./extract_matlab_doc.py > ./$DIRNAME/$FILENAME.py

# new approach: one file per matlab source file
./extract_matlab_doc_separate_files.py ../NFCP/util/ ./$DIRNAME/util/
./extract_matlab_doc_separate_files.py ../NFCP/output/ ./$DIRNAME/output/
./extract_matlab_doc_separate_files.py ../NFCP/optimize/ ./$DIRNAME/optimize/
./extract_matlab_doc_separate_files.py ../NFCP/matrix/ ./$DIRNAME/matrix/
./extract_matlab_doc_separate_files.py ../NFCP/ ./$DIRNAME/
#touch ./$DIRNAME/util/__init__.py

# Patch
rm ./$DIRNAME/__init__.py

sphinx-apidoc -fe -o . ./$DIRNAME

# Patch modules to avoid orphan error
echo -e ':orphan:\n' > _modules.rst
cat modules.rst >> _modules.rst
mv _modules.rst modules.rst

# Patch to make functions top-level in TOC without being include twice
# ) :
cp index.template index.rst
cat modules.rst | sed -n '/maxdepth/ { s///; :a; n; p; ba; }' >> index.rst

# Run sphinx to generate HTML
make clean html
#make SPHINXOPTS='-W' clean html

# Post-Processing (VERY BRITTLE)
# Strip extra NFCP that appears in headers
sed -i -e "s/$FILENAME.//g" ./_build/html/$FILENAME.html
# Remove "function" label (everything is a function)
sed -i -e "s/\bfunction\b/$TITLE/g" ./_build/html/$FILENAME.html
# Remove module labels (matlab doesn't really have modules in the python sense)
sed -i -e "s/\bmodule\b//g" ./_build/html/*.html
# Remove the "_auto" text leftover from the auto-generation
sed -i -e "s/\b_auto\b/$FILENAME/g" ./_build/html/*.html
# Change "Other Parameter" to varargin options
# (TODO)
# Remove the wrapping lines to sneak matlab source past the 
# Python documentation builder
sed -i -e 's/<span class="k">pass<\/span><span class="c1">#SKIPME<\/span>//g' ./_build/html/_modules/*.html
sed -i -e 's/\&#39;\&#39;\&#39;#STARTCODE//g' ./_build/html/_modules/*.html
sed -i -e 's/<span class="sd">    \&#39;\&#39;\&#39;<\/span><span class="c1">#STOPCODE<\/span>//g' ./_build/html/_modules/*.html
sed -i -e 's/#!\/usr\/bin\/env python//g' ./_build/html/_modules/*.html
sed -i -e 's/# -\*- coding: UTF-8 -\*-//g' ./_build/html/_modules/*.html
sed -i -e 's/\bdef\s/function /g' ./_build/html/_modules/*.html


echo "If you got an error about fulltoc, try running"
echo "pip install --user sphinxcontrib-fulltoc"
echo "if that gives random errors, uninstall spinx and reinstall it?"
echo "you may also need this:"
echo "pip install --user sphinx_rtd_theme"

echo 
echo "If build succeeded, the documentation should be in _build/html"

#!/usr/bin/env bash

# Run code documentation generator
# This uses a kludge to render numpy style docstrings in Matlab
# through Sphinx.
cd ./Documentation
echo "Building documentation"
./run_me_to_build_html_documentation.sh
cd ../

echo "Generating static file browsing links"
# Generate github pages browsing links
# (stopgap until official documentation is prepared)
./Gittools/maketree.py

# KEEP THIS DISABLED TO AVOID DELETING PNG AND PDF OUTPUTS

echo "Deleting files specified in .gitignore"
# Do not track files that should be skipped by gitignore
# (in the case that these files get accidentally added, gitignore will not
#  automatically remove them)
cat .gitpurge | awk "/^[.\*]/" | sed 's/"/"\\""/g;s/.*/"&"/' |  xargs -E '' -I{} git rm -rf --cached {} 2>/dev/ null

# Clean up editor and temp files from the local directory (even if not 
# tracked by git)
echo "Deleting editor temporary files"
find . -name "*.pyc" -exec rm -rf {} \; 2>/dev/null
find . -name "*~" -exec rm -rf {} \;  2>/dev/null

# Add any new files, add all updates to all files
echo "Adding all changes"
git add --all . 
git add -u :/

# Check that there are no files over 100M in this commit
echo "Removing files over 100mb"
find -type f -size +100M -exec git rm --cached {} \;

# Commit using the message specified as first argument to this script
echo "Git commit"
git commit -m "… $1"

# Synchronize with master on github
echo "git pull"
git pull

echo "git push"
git push origin master

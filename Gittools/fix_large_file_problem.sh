#!/usr/bin/env bash

# Github will not store files larger than 100mb. If files larger than 100mb get 
# commited to a repository, they must be removed. Removing the file from the
# repo and from the current commit is not enough, the file must be removed from
# the entire commit history. This is nontrivial. Git "filter branch" is too
# complicated, difficult, and error prone. We use the BFG repo cleaner.
# ("bfg.jar")

#git gc
git repack

java -jar bfg.jar --strip-blobs-bigger-than 100M ../.git




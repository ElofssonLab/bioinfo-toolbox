#!/bin/bash -x

# Just a script that actually checks if everything is nicely packes in the way it should be


# First check if all tar.gz files really are compressed

find ./ -name "*.tar.gz" -exec gzip -t {} \; 2> notcompressed.lst
for i in `gawk '{print $2}' notcompressed.lst  | sed s/://g ` ; do j=`echo $i | sed s/.gz$//`;  mv $i $j ; gzip $j ; done


# Then we need to delete all directories that we already have as tar files.

dir=`pwd`


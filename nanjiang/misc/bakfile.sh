#!/bin/bash

tmpfile=$(mktemp /tmp/tmp.bakfile.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
trap 'rm -rf "$tmpfile"' INT TERM EXIT
find ./ -maxdepth 1  -type f  -name "*.py" -o -name "*.xsl" -o -name "*.gnu" -o -name "*.sh" -o -name "*.pl"  | grep -v "conflicted.*copy" > $tmpfile
backupfile.sh -l $tmpfile
rm -f $tmpfile

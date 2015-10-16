#!/usr/bin/env bash

tmpfile=$(mktemp /tmp/tmp.changepythonenv.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  

for file in $*; do
    sed 's/^\#\!\/usr\/bin\/python/#!\/usr\/bin\/env python/g' $file > $tmpfile
    swapfile "$file" "$tmpfile"
    chmod 755 "$file"
    echo "$file"
done

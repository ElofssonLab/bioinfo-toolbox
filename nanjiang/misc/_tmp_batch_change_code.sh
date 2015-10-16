#!/bin/bash

tmpfile=$(mktemp /tmp/tmp._tmp_batch_change_code.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  

for file in plot*.sh; do
    sed 's/epstopdf/my_epstopdf/g' $file > $tmpfile
    swapfile $file $tmpfile
    chmod 755 $file
    echo $file modified
done

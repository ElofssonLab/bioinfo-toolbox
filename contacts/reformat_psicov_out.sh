#!/bin/bash

IN=$1

cat ${IN} | sed 's/ 0 8//g' | sed 's/ /,/g'  > ${IN}.tmp
mv ${IN}.tmp ${IN}


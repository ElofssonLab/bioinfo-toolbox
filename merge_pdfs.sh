#!/bin/bash
#test
DIR=$1
echo "Merging ${DIR}/*.pdf into ${DIR}/merged.pdf". This may take a while...
gs -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${DIR}/merged.pdf `ls ${DIR}/*.pdf`
echo "Files successfully merged!"

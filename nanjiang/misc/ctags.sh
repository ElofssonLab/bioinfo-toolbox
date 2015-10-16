#!/bin/bash
osname=`uname -s`
exec_ctags=/usr/bin/ctags
case $osname in 
    Linux*)
        exec_ctags=/usr/bin/ctags;;
    Darwin*)
        exec_ctags=/usr/local/bin/ctags;;
    *)
        exec_ctags=/usr/local/bin/ctags;;
esac

$exec_ctags *.{h,py,c,cpp}  ../../../program/MyInclude/*.{c,cpp,h}

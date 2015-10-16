#!/bin/bash -f 

usage="
Usage: removehtmltag.sh  [-i infile] [-o outfile]

Examples:
   removehtmltag.sh < in.txt  > out.txt
   removehtmltag.sh -i in.txt -o out.txt
   cat in.txt | removehtmltag.sh > out.txt

Created 2013-02-11, updated 2013-02-11, Nanjiang Shu 
"
PrintHelp(){ #{{{
    echo "$usage"
}
#}}}

infile=
outfile=
isNonOptionArg=0
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then 
        echo Error! Wrong argument: $1 >&2; exit
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -i|--i|-infile|--infile)infile=$2;shift;;
            -o|--o) outfile=$2;shift;;
            -q|-quiet|--quiet) isQuiet=1;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else
        echo Error! Wrong argument: $1 >&2; exit
    fi
    shift
done

if [ "$outfile" != "" ]; then
    (
    if [ "$infile" != "" ]; then
        sed -n ' /^$/!{s/<[^>]*>//g;p;} ' < $infile
    else
        sed -n ' /^$/!{s/<[^>]*>//g;p;} '
    fi
    ) > $outfile
else
    if [ "$infile" != "" ]; then
        sed -n ' /^$/!{s/<[^>]*>//g;p;} ' < $infile
    else
        sed -n ' /^$/!{s/<[^>]*>//g;p;} '
    fi
fi


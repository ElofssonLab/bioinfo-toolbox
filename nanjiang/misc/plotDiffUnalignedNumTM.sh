#!/bin/bash
# draw the histogram given the file in the format
# #topology comparison of four states
# #
# format of the file to make the plot
# #Idx  NumDiffTM %# Count
# 0     1         50.0  130
# 1     2         ?    35
# 2     3         ?   33
# 3     5         ?   50
# 4     6         ?   60


usage="
Usage: $0 datafile [datafile ...]
Description: 
Options:
    -outstyle png|term|eps  set the output, default = terminal
    -outpath DIR            set the output path, default =./
    -h|--help               print this help message and exit

Created 2013-05-30, updated 2013-05-30, Nanjiang Shu 
"
PrintHelp() {
    echo "$usage"
}

MakePlot(){
    local basename=`basename "${dataFileList[0]}"`
    if [ $numFile -gt 1 ];then
        basename=${basename}.cmb${numFile}
    fi
    local outputSetting=
    case $outputStyle in
        term*)
            outputSetting=
            ;;
        png)
            outfile=$outpath/$basename.png
            outputSetting="
            set terminal png enhanced
            set output '$outfile' 
            "
            ;;
        eps)
            outfile=$outpath/$basename.eps
            pdffile=$outpath/$basename.pdf
            outputSetting="
            set term postscript eps color
            set output '$outfile' 
            "
            ;;
    esac
    local plotOption=
    local file=
    local colorindex=
    local title=
    local trailingcomma=
    local highestI=

    for ((i=0;i<$numFile;i++)); do
        file=${dataFileList[$i]}
        ((colorindex=i+1))
        title=
        trailingcomma=
        case $file in 
            *diffNumTM_Unaligned_Nterm*) title=N-terminal;;
            *diffNumTM_Unaligned_Cterm*) title=C-terminal;;
        esac
        ((highestI=numFile-1))
        if [ $highestI -eq $i ];then
            trailingcomma=
        else
            trailingcomma=,
        fi

        plotOption="$plotOption '$file' using 4:xtic(2) title '$title' ls $colorindex $trailingcomma"
    done


    /usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "green"
set style line 3 lt 1 lw 1 lc rgb "blue" 
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "cyan" 
set style line 6 lt 1 lw 1 lc rgb "orange"
set style line 7 lt 1 lw 1 lc rgb "grey40"
set style line 8 lt 1 lw 1 lc rgb "grey70"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"
#set key top center outside

set style data histogram 
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.8
set tics font 'Times-Roman,16'

plot $plotOption

EOF

    case $outputStyle in
        eps) my_epstopdf $outfile
            echo "Histogram image output to $pdffile"
            if [ -s "$pdffile" ]; then
                rm -f "$outfile"
            fi
            ;;
        *) echo "Histogram image output to $outfile" ;;
    esac

}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outpath=
dataFileList=()
isNonOptionArg=false
xlabel="Difference"
ylabel="Occurrences"

while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        dataFileList+=("$1")
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h|--help) PrintHelp; exit;;
            -outstyle|--outstyle) outputStyle=$2;shift;;
            -xlabel|--xlabel) xlabel="$2";shift;;
            -ylabel|--ylabel) ylabel="$2";shift;;
            -outpath|--outpath) outpath=$2;shift;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        dataFileList+=("$1")
    fi
    shift
done

numFile=${#dataFileList[@]}
if [ $numFile -eq 0  ]; then
    echo Input not set! Exit. >&2
    exit 1
fi


if [ "$outpath" == "" ]; then
    outpath=`dirname ${dataFileList[0]}`
elif [ ! -d "$outpath" ]; then
    mkdir -p $outpath
fi

MakePlot

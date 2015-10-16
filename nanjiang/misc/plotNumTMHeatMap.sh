#!/bin/bash
# draw numTM frequency heat map
# #
# format of the file to make the plot

usage="
Usage: $0 datafile
Description: 
Options:
    -o, -outstyle STR   set the output format,  (default: eps)
                        can be one of eps, png, term
    -outpath DIR        set the output path, default: \$dirname{datafile}
    -h,--help           print this help message and exit

Created 2012-09-25, updated 2013-01-23, Nanjiang Shu 
"
PrintHelp() {
    echo "$usage"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outpath=
dataFile=
isNonOptionArg=false
xlabel="Sequence identity"


font_label="Arial, 24"
font_tics="Arial, 24"
font_key="Arial, 24"


while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        dataFile=$1
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h|--help) PrintHelp; exit;;
            -outstyle|--outstyle|-o) outputStyle=$2;shift;;
            -xlabel|--xlabel) xlabel="$2";shift;;
            -ylabel|--ylabel) ylabel="$2";shift;;
            -outpath|--outpath) outpath=$2;shift;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        dataFile=$1
    fi
    shift
done


if [ ! -f "$dataFile" ]; then 
    echo "Error! dataFile = \"$dataFile\" does not exist. Exit..."
    exit
fi

if [ "$outpath" == "" ]; then
    outpath=`dirname $dataFile`
elif [ ! -d "$outpath" ]; then
    mkdir -p $outpath
fi

basename=`basename "$dataFile"` 

outputSetting=
case $outputStyle in
    term*)
    outputSetting=
    ;;
    png)
    outfile=$outpath/$basename.png
    outputSetting="
set terminal png color
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
    svg)
    outfile=$outpath/$basename.svg
    outputSetting="
set term svg
set output '$outfile'
    "
    ;;
esac

numLine=`awk 'BEGIN{cnt=0}/^[^#]/{cnt+=1}END{print cnt}' $dataFile`
((numLine--))
if [ $numLine -ge 14 ]; then
    numLine=14
fi
rangeMax=${numLine}.5

labelSetting=""
for ((i=1;i<=numLine;i++));do
    #((j=i-1))
    labelSetting="$labelSetting
    set label  $i \"$i\"  at  $i,$i  center font \"Symbol,24\" front tc rgb \"white\""
done

#echo "$labelSetting"

tmpdatafile=$(mktemp /tmp/tmp.plotNumTMHeatMap.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  

trap 'rm -f "$tmpdatafile"' INT TERM EXIT

awk '{for(i=1;i<=NF;i++){printf("%4g ", $i)} printf("\n")}' $dataFile > $tmpdatafile
#cat $tmpdatafile >/tmp/t1.txt


pdffile=$outpath/$basename.pdf
pngfile=$outpath/$basename.png
/usr/bin/gnuplot -persist<<EOF 
$outputSetting

#set style line 1 lt 1 lw 2 lc rgb "white"
set style line 11 lt 1 lw 2 lc rgb "white"


$labelSetting
unset key
#set palette rgbformulae 22,13,10
#color(gray) = gray**(1./gamma)
#set palette rgbformulae 30,31,32
#set palette rgbformulae  34,35,36
#set palette rgbformulae 3,3,3
#set palette rgbformulae 0,0,0
#set palette gray
#set palette functions sqrt(gray), sqrt(gray), sqrt(gray)
#set palette functions gray**2, gray**3, gray**3

#==================================
gamma = 0.35
#gamma = 1
#set logscale cb 2
#color(gray) = gray**(1./gamma)
#set palette model RGB functions color(gray), color(gray), color(gray)
#set palette negative
#==================================
set palette rgbformulae 3,3,3
set palette negative

#==================================
# color scheme
set logscale cb 2
set palette rgbformulae 22,13,10
set palette negative
#set palette rgbformulae 22,13,10
#==================================
#set palette model HSV functions gray, 1, 1


set xlabel "Number of TM helices" font "$font_label"
set ylabel "Number of TM helices" font "$font_label"
#set cbrange [0:100]
set cblabel "Normalized number of pairs" offset 2 font "$font_label"
set tic scale 1
set tics font "$font_tics"
#unset cbtics


#set xrange [0.5:$rangeMax]
#set yrange [0.5:$rangeMax]
set xrange [0.5:$rangeMax]
set yrange [0.5:$rangeMax]
#set xtic 1
#set ytic 1

plot '$tmpdatafile' matrix with image ti ''

#set pm3d interpolate 1,1
#unset surface; set pm3d map
#splot '$tmpdatafile' matrix with image
#set view map
#set pm3d map
#splot '$tmpdatafile' matrix

EOF

if [ -s "$outfile" ]; then
    case $outputStyle in
        eps) 
            my_epstopdf $outfile 
            convert -density 300 $pdffile $pngfile
            echo "Histogram image output to $pdffile $pngfile"
            ;;
        *)
            echo "Histogram image output to $outfile"
            ;;
    esac
fi

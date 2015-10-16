#!/bin/bash
# draw the histogram given the file in the format
# #topology comparison of four states
# #
# format of the file to make the plot
#  #Idx  SeqIDT         0         1         2   Maximum Occurrence
# 0       0-10    69.431    17.082    13.488    69.431     503959
# 1      10-20    82.700    10.024     7.276    82.700     484708
# 2      20-30    90.896     5.825     3.279    90.896     262787
# 3      30-40    94.371     3.555     2.074    94.371     303315
# 4      40-50    95.357     2.809     1.833    95.357     240253
# 5      50-60    95.549     2.710     1.741    95.549      86855
# 6      60-70    96.863     1.500     1.637    96.863      48070
# 7      70-80    97.592     1.157     1.250    97.592      29029
# 8      80-90    97.947     0.969     1.083    97.947      18464
# 9     90-100    98.337     0.725     0.937    98.337      16960
# #sum            85.119     8.585     6.296    85.119    1994400   
usage="
Usage: $0 datafile
Description: 
Options:
    -outstyle,-o STR    Set the output, (default: eps)
                        can be png, term, eps
    -outpath DIR        Set the output path, (default: the same as datafile)
    -xlabel STR         Set xlabel, (default: Sequence identity)
    -ylabel STR         Set Ylabel
    -plot1              Print merged bin, 0-10,10-20,20-30,30-40,40-100
    -h,--help           Print this help message and exit

Created 2012-09-11, updated 2013-05-03, Nanjiang Shu 
"
cmd="$*"
PrintHelp(){
    echo "$usage"
}
CreateDataFile1(){ #{{{
    local infile="$1"
    local outfile="$2"
awk -v method=$method -v idx=$idx '
BEGIN{
    totalsum = 0;
    COL=0;
    STARTINDEX = 4
    for (j=3;j<=9;j++){
        subsum[i] = 0;
    }
}

{
    if ($1>=STARTINDEX){
        COL=NF;
        totalsum += $NF;
        for (i=3;i<=NF-2;i++){
            cnt[i] = $i*$NF/100;
            subsum[i] +=  cnt[i]
#             print int(cnt[i]), int(subsum[i]);
        }
    }else if ($1 != "#sum"){
        print
    }
}
END{
    for (i=3;i<=COL-2;i++){
        freq[i] = subsum[i]/totalsum;
    }
    printf("%d ", STARTINDEX);
    printf("%d-100 ", STARTINDEX*10);
    maxf = -1.0
    for (i=3;i<=COL-2;i++){
        printf("%6.3f ", freq[i]*100);
        if (freq[i]*100 > maxf){
            maxf = freq[i]*100
        }
    }
    printf("%6.3f %d \n", maxf, totalsum);
}' $infile > $outfile
}
#}}}
CreateDataFile2(){ #{{{
    local infile="$1"
    local outfile="$2"
awk -v method=$method -v idx=$idx '
BEGIN{
#for merging  0-10 and 10-20
    totalsum0 = 0
    for (j=3;j<=9; j++){
        subsum0[i] = 0;
    }

    totalsum = 0;
    COL=-1;
    STARTINDEX = 4


    for (j=3;j<=9;j++){
        subsum[i] = 0;
    }
}

{
    if ($1 == 0 || $1 == 1){
        if(COL == -1){ COL=NF;}
        totalsum0 += $NF;
        for (i=3;i<=NF-2;i++){
            subsum0[i] +=  $i*$NF/100
        }
    }else if ($1 == 2){
        printf("%d ", 0);
        printf("0-20 ");
        maxf = -1.0
        for (i=3;i<=COL-2;i++){
            freq[i] = subsum0[i]/totalsum0;
            printf("%6.3f ", freq[i]*100);
            if (freq[i]*100 > maxf){
                maxf = freq[i]*100
            }
        }
        printf("%6.3f %d \n", maxf, totalsum0);

        printf("%d ", $1-1);
        printf("%s ", $2);
        for (i=3;i<=NF-2;i++){
            printf("%6.3f ", $i);
        }
        printf("%6.3f %d \n", $(NF-1), $NF);


    }else if ($1>=STARTINDEX){
        if(COL == -1){ COL=NF;}
        totalsum += $NF;
        for (i=3;i<=NF-2;i++){
            cnt[i] = $i*$NF/100;
            subsum[i] +=  cnt[i]
#             print int(cnt[i]), int(subsum[i]);
        }
    }else if (substr($1,1,1) != "#"){
        printf("%d ", $1-1);
        printf("%s ", $2);
        for (i=3;i<=NF-2;i++){
            printf("%6.3f ", $i);
        }
        printf("%6.3f %d \n", $(NF-1), $NF);
    }else if ($1 != "#sum"){
        print $0
    }
}
END{
    for (i=3;i<=COL-2;i++){
        freq[i] = subsum[i]/totalsum;
    }
    printf("%d ", STARTINDEX-1);
    printf("%d-100 ", STARTINDEX*10);
    maxf = -1.0
    for (i=3;i<=COL-2;i++){
        printf("%6.3f ", freq[i]*100);
        if (freq[i]*100 > maxf){
            maxf = freq[i]*100
        }
    }
    printf("%6.3f %d \n", maxf, totalsum);
}' $infile > $outfile
}
#}}}
MakePlot_shownumcase() {  #{{{
    local dataFile="$1"
    local outputSetting=
    local basename=`basename "$dataFile"` 
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

    /usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 lc rgb "orange"
set style line 2 lt 1 lw 1 lc rgb "grey20"
set style line 3 lt 1 lw 1 lc rgb "purple"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" font "Arial, 16"
set key top center outside
#set key bottom center outside
set key invert
set label 1 at 3.5, 110
set label 1 "Number of cases"

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set tics font 'Times-Roman,16'
set grid ytics
#set bmargin 5

plot \
    '$dataFile' using 3:xtic(2)   title 'TM2TM' ls 1 , \
    ''          using 4:xtic(2)   title 'TM2GAP' ls 2, \
    ''          using 5:xtic(2)   title 'TM2SEQ' ls 3, \
    ''          using (\$1):(\$3+\$4+\$5+3):7 title '' with labels
EOF
    local pdffile=$outpath/$basename.pdf
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
#}}}
MakePlot_no() {  #{{{
    local dataFile="$1"
    local basename=`basename "$dataFile"`
    local outfile=
    local pdffile=
    local outputSetting=""
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

    /usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 lc rgb "orange"
set style line 2 lt 1 lw 1 lc rgb "grey20"
set style line 3 lt 1 lw 1 lc rgb "purple"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" font "Arial, 16"
set key top center outside
#set key bottom center outside
set key invert
set yrange[0:100]

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set tics font 'Times-Roman,16'
set grid ytics

plot \
    '$dataFile' using 3:xtic(2)   title 'TM2TM' ls 1 , \
    ''          using 4:xtic(2)   title 'TM2GAP' ls 2, \
    ''          using 5:xtic(2)   title 'TM2SEQ' ls 3
EOF
    local pdffile=$outpath/$basename.pdf

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
#}}}
MakePlot_no_invert() {  #{{{
    local dataFile="$1"
    local basename=`basename "$dataFile"`
    local outfile=
    local pdffile=
    local outputSetting=""
    basename=$basename.revert
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

    /usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 lc rgb "orange"
set style line 2 lt 1 lw 1 lc rgb "grey20"
set style line 3 lt 1 lw 1 lc rgb "grey90"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" font "Arial, 16"
set key top center outside
#set key bottom center outside
set key invert
set yrange[0:6]

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set tics font 'Times-Roman,16'
set grid ytics

plot \
    '$dataFile' using 4:xtic(2)   title 'TM2GAP' ls 2, \
    ''          using 5:xtic(2)   title 'TM2SEQ' ls 3
EOF
    local pdffile=$outpath/$basename.pdf

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
#}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outpath=
dataFile=
xlabel="Sequence identity"
ylabel="Percentages of TM helix comparison classes"
isNonOptionArg=false
isShowNumCase=0
isPlot1=0

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
            -plot1|--plot1) isPlot1=1;;
            -showcase|--showcase) topt=$2; 
                if [ "${topt:0:1}" == "y" -o  "${topt:0:1}" == "Y" ]; then
                    isShowNumCase=1
                else
                    isShowNumCase=0
                fi
                shift
                ;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        dataFile=$1
    fi
    shift
done


if [ ! -f "$dataFile" ]; then 
    echo "Error! dataFile = \"$dataFile\" does not exist. Exit..."
    echo "$cmd"
    exit
fi

if [ "$outpath" == "" ]; then
    outpath=`dirname $dataFile`
elif [ ! -d "$outpath" ]; then
    mkdir -p $outpath
fi


stemname=${dataFile%.*}
dataFile1=$stemname.1.txt
dataFile2=$stemname.2.txt

if [ $isPlot1 -eq 1 ];then
    CreateDataFile1 "$dataFile" "$dataFile1"
    CreateDataFile2 "$dataFile" "$dataFile2"
fi


if [ $isShowNumCase -eq 1 ]; then
    MakePlot_shownumcase $dataFile
    if [ $isPlot1 -eq 1 ];then
        MakePlot_shownumcase $dataFile1
        MakePlot_shownumcase $dataFile2
    fi
else
    MakePlot_no $dataFile
    if [ $isPlot1 -eq 1 ];then
        for file in $dataFile1 $dataFile2; do
            MakePlot_no $file
            MakePlot_no_invert $file
        done
    fi
fi



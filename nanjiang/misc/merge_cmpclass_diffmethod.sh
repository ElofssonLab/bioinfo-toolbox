#!/bin/bash

# two methods, one is sampled average, another is numerical average
# default is numerical average

infile=$1
method=$2
idx=$3

avgmtd=0
if [ "$4" != "" ]; then 
    avgmtd=$4
fi
# avgmtd:
# 0: numerical average
# 1: sampled average
# default: 0

isGrep40_100=0
if [ "$5" == "40-100" ];then
    isGrep40_100=1  #temporary workaround
fi

if [ "$infile" == "" ]; then
    echo "usage: $0 infile"
    exit
fi

SampledAverage(){ #{{{
awk -v method=$method -v idx=$idx '
BEGIN{
    totalsum = 0;
    COL=0;
    for (i=3;i<=9;i++){
        subsum[i] = 0;
    }
}

{
    if ($1>=0){
        COL=NF;
        totalsum += $NF;
        for (i=3;i<=NF-2;i++){
            cnt[i] = $i*$NF/100;
            subsum[i] +=  cnt[i]
#             print int(cnt[i]), int(subsum[i]);
        }
    }
}
END{
    for (i=3;i<=COL-2;i++){
        freq[i] = subsum[i]/totalsum;
    }
    printf("%d ", idx);
    printf("%s ", method);
    for (i=3;i<=COL-2;i++){
        printf("%6.3f ", freq[i]*100);
    }
    printf("%d \n", totalsum);
}' $infile
}
#}}}
NumericalAverage(){ #{{{
awk -v method=$method -v idx=$idx '
BEGIN{
    totalsum = 0;
    COL = 0
    cntRow = 0;
    for (i=3;i<=9;i++){
        subsum[i] = 0;
    }
}

{
    if ($1>=0){
        COL = NF
        totalsum += $NF;
        cntRow += 1;
        for (i=3;i<=NF-2;i++){
            subsum[i] += $i;
        }
    }
}
END{
    for (i=3;i<=COL-2;i++){
        freq[i] = subsum[i]/cntRow/100;
    }
    printf("%d ", idx);
    printf("%s ", method);
    for (i=3;i<=COL-2;i++){
        printf("%6.3f ", freq[i]*100);
    }
    printf("%9d \n", totalsum);
}' $infile
}
#}}}

if [ $isGrep40_100 -eq 1 ];then
    awk -v method=$method -v idx=$idx '
    {
        if($2=="40-100"){
            printf("%d ", idx);
            printf("%s ", method);
            for (i=3;i<=NF-2;i++){
                printf("%s ", $i);
            }
            printf("%9d\n", $NF);
        }
    }
' $infile
    
else
    case $avgmtd in 
        0) NumericalAverage;;
        1) SampledAverage;;
    esac
fi

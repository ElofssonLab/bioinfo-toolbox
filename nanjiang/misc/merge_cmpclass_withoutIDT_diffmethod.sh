#!/bin/bash

infile=$1
method=$2
idx=$3

if [ "$infile" == "" ]; then
    echo "usage: ana_txt.sh infile"
    exit
fi

awk -v method=$method -v idx=$idx '
BEGIN{
    totalsum = 0;
    for (i=3;i<=7;i++){
        subsum[i] = 0;
    }
}

{
    if ($1>=3){
        totalsum += $8;
        for (i=3;i<=7;i++){
            cnt[i] = $i*$8/100;
            subsum[i] +=  cnt[i]
#             print int(cnt[i]), int(subsum[i]);
        }
    }
}
END{
    for (i=3;i<=7;i++){
        freq[i] = subsum[i]/totalsum;
    }
    printf("%d ", idx);
    printf("%s ", method);
    for (i=3;i<=7;i++){
        printf("%6.3f ", freq[i]*100);
    }
    printf("%d \n", totalsum);
}' $infile


#!/bin/bash
# merge_cmpclass_diffmethod_2class_non_accumulate.sh
infile=$1
method=$2
idx=$3
min_idx_to_count=$4

if [ "$infile" == "" ]; then
    echo "usage: $0 infile"
    exit
fi

awk -v method=$method -v idx=$idx -v min_idx_to_count=$min_idx_to_count '
BEGIN{
    totalsum = 0;
    for (i=3;i<=3;i++){
        subsum[i] = 0;
    }
}

{
    if ($1 == min_idx_to_count){
        totalsum += $NF;
        for (i=3;i<=3;i++){
            cnt[i] = $i*$NF/100;
            subsum[i] +=  cnt[i]
#             print int(cnt[i]), int(subsum[i]);
        }
    }
}
END{
    for (i=3;i<=3;i++){
        freq[i] = subsum[i]/totalsum;
    }
    printf("%-2d ", idx);
    printf("%-15s ", method);
    for (i=3;i<=3;i++){
        printf("%6.3f ", freq[i]*100);
    }
    i=3
    printf("%6.3f ", (1.0-freq[i])*100);
    printf("%d \n", totalsum);
}' $infile


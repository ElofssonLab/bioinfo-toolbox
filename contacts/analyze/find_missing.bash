#!/bin/bash



for dir in PF*/.
do
    #First delete all size 0 files
    find ${dir} -size 0 -exec rm {} \;
    for j in  $dir/*cm.tar.gz # $dir/conf*[04].tar.gz
    do
	i=`basename $j .tar.gz `
	if [ ! -s $dir/${i}_summary.csv ]
	then
	    k=`tar -ztvf $dir/${i}.tar.gz | grep -c ".pdb" `
	    if [ $k -lt 10  ]
	    then
		echo "TAR: "$dir/${i}.tar.gz $k
	    fi
	    if [ ! -s $dir/${i}_cns.out ]
	    then
		echo "CNS: "$dir/${i}_cns.out
	    fi
	    if [ ! -s $dir/${i}_proq3.tar.gz ]
	    then
		echo "ProQ3: "$dir/${i}_proq3.tar.gz
	    fi
	    if [ ! -s $dir/${i}.raw ]
	    then
		echo "Pcons: "$dir/${i}.raw
	    fi
	    if [ ! -s $dir/${i}_TM.out ]
	    then
		echo "TM: "$dir/${i}_TM.out
	    else
		k=`grep -c TM-score $dir/${i}_TM.out`
		if [ $k -lt 20  ]
		then
		    echo "TM: "$dir/${i}_TM.out $k
		fi
	    fi
	fi
    done
done

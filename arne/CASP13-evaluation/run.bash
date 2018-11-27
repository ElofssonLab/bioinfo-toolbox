#!/bin/bash -x
for i in `cat targets.txt` ; do cat $i.txt ; done | gawk '{print $2, $4}'    | grep -v Model | grep TS | sed "s/-D1//g" | sort > GDT_TS.tsv
for i in `cat targets.txt` ; do cat $i.txt ; done | gawk '{print $2, $11}'   | grep -v Model | grep TS | sed "s/-D1//g" | sort > GDT_HA.tsv
for i in `cat targets.txt` ; do cat $i.txt ; done | gawk '{print $2, $35}'   | grep -v Model | grep TS | sed "s/-D1//g" | sort > lDDT.tsv
for i in `cat targets.txt` ; do cat $i.txt ; done | gawk '{print $2, $37}'   | grep -v Model | grep TS | sed "s/-D1//g" | sort > CAD.tsv
for i in `cat targets.txt` ; do cat $i.txt ; done | gawk '{print $2, $46}'   | grep -v Model | grep TS | sed "s/-D1//g" | sort > TM.tsv
# sed s/TS.*//g TM.tsv | sort -u > targets.txt  # THis was edited manually due to domains
for i in ../../qa/QA*tsv ; do j=`basename $i .tsv` ; for k in GDT_*.tsv lDDT.tsv CAD.tsv TM.tsv ; do l=`basename $k .tsv` ; join $i $k > $j-$l.tsv ; done ; done


ls QA*TM.tsv  | sed s/\-TM.tsv//g | sort -u   > methods.txt


for j in TM CAD GDT_TS GDT_HA lDDT
do
    for i in QA*-$j.tsv
    do
	echo -n $i " "
	~/bin/corr 2 3 $i
    done | sort -rnk 8 | sed s/QA267_2/ProQ3D_TM/ | sed s/QA457_2/Pcomb/ | sed s/QA139_2/ProQ3D/ | sed s/QA044_2/ProQ2/ | sed s/QA022_2/Pcons/ | sed s/QA187_2/ProQ3/  | sed s/QA360_2/ProQ3D_lDDT/ | sed s/QA198_2/ProQ3D_CAD/  | sed s/QA440_2/ProQ4/ | gawk '{print $1,$3,$8}' | sed "s/\-/ /g" | sed s/.tsv//g | sed "s/\s+/ /g"> corr-$j.tsv
done


for k in TM CAD GDT_TS GDT_HA lDDT
do
    for j in `cat methods.txt `
    do
	echo -n $j " "
	for i in `cat targets.txt | sed "s/\-D1//g"`
	do
	    sort -nk 2 $j-$k.tsv | grep $i  | tail -1
	done | gawk '{i++;sum+=$3};END{print i,sum,sum/i}'
    done | sort -rnk 4 | sed s/QA267_2/ProQ3D_TM/ | sed s/QA457_2/Pcomb/ | sed s/QA139_2/ProQ3D/ | sed s/QA044_2/ProQ2/ | sed s/QA022_2/Pcons/ | sed s/QA187_2/ProQ3/  | sed s/QA360_2/ProQ3D_lDDT/ | sed s/QA198_2/ProQ3D_CAD/  | sed s/QA440_2/ProQ4/| sed "s/\s+/ /g" > sum-$k.tsv
done


for i in QA*tsv
do
    j=`echo $i | sed "s/QA/QQ/g" ` 
    sed "s/TS/ TS/g" $i > $j
done

for k in TM CAD GDT_TS GDT_HA lDDT
do
    for l in `cat methods.txt `
    do
	j=`echo $l | sed "s/QA/QQ/g" `
	echo -n $j " " 
	for m in `cat targets.txt | sed "s/\-D1//g" | sort -u`
	do
	    grep -w $m ${j}-$k.tsv  |  gawk -v v1=3 -v v2=4 -f $HOME/bin/correlation.awk | gawk '{print $6}'
	done   | gawk '{i+=1;sum+=$1};END{ if (i>0){print i,sum,sum/i}}'
    done | sort -rnk 4 | sed s/QQ267_2/ProQ3D_TM/ | sed s/QQ457_2/Pcomb/ | sed s/QQ139_2/ProQ3D/ | sed s/QQ044_2/ProQ2/ | sed s/QQ022_2/Pcons/ | sed s/QQ187_2/ProQ3/  | sed s/QQ360_2/ProQ3D_lDDT/ | sed s/QQ198_2/ProQ3D_CAD/  | sed s/QQ440_2/ProQ4/ | sed s/QQ/QA/g | sed "s/\s+/ /g" > pertarget-$k.tsv
done
    
grep -ni p  sum-*tsv  > summary.txt
grep -ni p  corr-*tsv  > summary-corr.txt
grep -ni p  pertarget-*tsv  > summary-pertarget.txt

for i in sum-*.tsv corr-*tsv pertarget*tsv
do
    j=`basename $i .tsv`
    Rscript bin/Zscore.R $i Z-$j.tsv
done

for k in TM CAD GDT_TS GDT_HA lDDT
do
    Rscript bin/merge.R Z-sum-$k.tsv Z-corr-$k.tsv Z-pertarget-$k.tsv average-$k.tsv
done

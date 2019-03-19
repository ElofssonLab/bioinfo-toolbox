#!/bin/bash -x
for i in T*txt
do
    echo -n $i " "
    grep -v \# $i | gawk '{i++;sum+=$4;if($4>best){best=$4}};END{print i,sum,sum/i,best}'
done | sort -nrk 4 | sed s/.txt//g > targets-difficulty.txt
gawk '{print $1}' targets-difficulty.txt | sort > alltargets.txt

# Manally edited to create targets.txt

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
    done | sort -rnk 8 | sed s/QA267_2/ProQ3D_TM/ | sed s/QA457_2/Pcomb/ | sed s/QA139_2/ProQ3D/ | sed s/QA044_2/ProQ2/ | sed s/QA022_2/Pcons/ | sed s/QA187_2/ProQ3/  | sed s/QA360_2/ProQ3D_lDDT/ | sed s/QA198_2/ProQ3D_CAD/  | sed s/QA440_2/ProQ4/ | gawk '{print $1,$3,$8}' | sed "s/\-/ /g" | sed s/.tsv//g | sed s/QA359_2/3DCNN/g |     sed s/QA014_2/Bhattacharya-ClustQ/g |      sed s/QA102_2/Bhattacharya-Server/g |      sed s/QA170_2/Bhattacharya-SingQ/g |      sed s/QA471_2/CPClab/g |      sed s/QA349_2/Davis-EMAconsensus/g |      sed s/QA413_2/FALCON-QA/g | sed s/QA027_2/FaeNNz/g |      sed s/QA196_2/Grudinin/g |      sed s/QA065_2/Jagodzinski-Cao-QA/g | sed s/QA344_2/Kiharalab/g |      sed s/QA067_2/LamoureuxLab/g |      sed s/QA146_2/MASS1/g |      sed s/QA415_2/MASS2/g |      sed s/QA197_2/MESHI/g |      sed s/QA289_2/MESHI-enrich-server/g |      sed s/QA347_2/MESHI-server/g |      sed s/QA113_2/MUFold-QA/g | sed s/QA312_2/MUFold_server/g |      sed s/QA243_2/MULTICOM-CONSTRUCT/g |      sed s/QA023_2/MULTICOM-NOVEL/g |      sed s/QA058_2/MULTICOM_CLUSTER/g |      sed s/QA107_2/MUfold-QA2/g | sed s/QA211_2/MUfold-QA-T/g | sed s/QA275_2/ModFOLD7/g |      sed s/QA213_2/ModFOLD7_cor/g |      sed s/QA272_2/ModFOLD7_rank/g |      sed s/QA373_2/ModFOLDclust2/g |      sed s/QA209_2/PLU-Angular-QA/g | sed s/QA134_2/PLU-Top-QA/g | sed s/QA083_2/Pcomb/g |      sed s/QA022_2/Pcons/g |      sed s/QA044_2/ProQ2/g |      sed s/QA187_2/ProQ3/g |      sed s/QA139_2/ProQ3D/g |      sed s/QA198_2/ProQ3D-CAD/g |      sed s/QA267_2/ProQ3D-TM/g |      sed s/QA360_2/ProQ3D-lDDT/g |      sed s/QA440_2/ProQ4/g |      sed s/QA334_2/RaptorX-Deep-QA/g | sed s/QA220_2/SASHAN/g |      sed s/QA135_2/SBROD/g |      sed s/QA207_2/SBROD-plus/g |      sed s/QA364_2/SBROD-server/g |      sed s/QA194_2/UOSHAN/g |      sed s/QA339_2/VoroM-QA-A/g | sed s/QA030_2/VoroM-QA-B/g | sed s/QA457_2/Wallner/g |       sed "s/\s+/ /g"> corr-$j.tsv
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
    done | sort -rnk 4 | sed s/QA267_2/ProQ3D_TM/ | sed s/QA457_2/Pcomb/ | sed s/QA139_2/ProQ3D/ | sed s/QA044_2/ProQ2/ | sed s/QA022_2/Pcons/ | sed s/QA187_2/ProQ3/  | sed s/QA360_2/ProQ3D_lDDT/ | sed s/QA198_2/ProQ3D_CAD/ | sed s/QA359_2/3DCNN/g |     sed s/QA014_2/Bhattacharya-ClustQ/g |      sed s/QA102_2/Bhattacharya-Server/g |      sed s/QA170_2/Bhattacharya-SingQ/g |      sed s/QA471_2/CPClab/g |      sed s/QA349_2/Davis-EMAconsensus/g |      sed s/QA413_2/FALCON-QA/g | sed s/QA027_2/FaeNNz/g |      sed s/QA196_2/Grudinin/g |      sed s/QA065_2/Jagodzinski-Cao-QA/g | sed s/QA344_2/Kiharalab/g |      sed s/QA067_2/LamoureuxLab/g |      sed s/QA146_2/MASS1/g |      sed s/QA415_2/MASS2/g |      sed s/QA197_2/MESHI/g |      sed s/QA289_2/MESHI-enrich-server/g |      sed s/QA347_2/MESHI-server/g |      sed s/QA113_2/MUFold-QA/g | sed s/QA312_2/MUFold_server/g |      sed s/QA243_2/MULTICOM-CONSTRUCT/g |      sed s/QA023_2/MULTICOM-NOVEL/g |      sed s/QA058_2/MULTICOM_CLUSTER/g |      sed s/QA107_2/MUfold-QA2/g | sed s/QA211_2/MUfold-QA-T/g | sed s/QA275_2/ModFOLD7/g |      sed s/QA213_2/ModFOLD7_cor/g |      sed s/QA272_2/ModFOLD7_rank/g |      sed s/QA373_2/ModFOLDclust2/g |      sed s/QA209_2/PLU-Angular-QA/g | sed s/QA134_2/PLU-Top-QA/g | sed s/QA083_2/Pcomb/g |      sed s/QA022_2/Pcons/g |      sed s/QA044_2/ProQ2/g |      sed s/QA187_2/ProQ3/g |      sed s/QA139_2/ProQ3D/g |      sed s/QA198_2/ProQ3D-CAD/g |      sed s/QA267_2/ProQ3D-TM/g |      sed s/QA360_2/ProQ3D-lDDT/g |      sed s/QA440_2/ProQ4/g |      sed s/QA334_2/RaptorX-Deep-QA/g | sed s/QA220_2/SASHAN/g |      sed s/QA135_2/SBROD/g |      sed s/QA207_2/SBROD-plus/g |      sed s/QA364_2/SBROD-server/g |      sed s/QA194_2/UOSHAN/g |      sed s/QA339_2/VoroM-QA-A/g | sed s/QA030_2/VoroM-QA-B/g | sed s/QA457_2/Wallner/g |       sed s/QA440_2/ProQ4/| sed "s/\s+/ /g" > sum-$k.tsv
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
    done | sort -rnk 4 | sed s/QQ267_2/ProQ3D_TM/ | sed s/QQ457_2/Pcomb/ | sed s/QQ139_2/ProQ3D/ | sed s/QQ044_2/ProQ2/ | sed s/QQ022_2/Pcons/ | sed s/QQ187_2/ProQ3/  | sed s/QQ360_2/ProQ3D_lDDT/ | sed s/QQ198_2/ProQ3D_CAD/  | sed s/QQ440_2/ProQ4/  | sed s/QQ359_2/3DCNN/g |     sed s/QQ014_2/Bhattacharya-ClustQ/g |      sed s/QQ102_2/Bhattacharya-Server/g |      sed s/QQ170_2/Bhattacharya-SingQ/g |      sed s/QQ471_2/CPClab/g |      sed s/QQ349_2/Davis-EMAconsensus/g |      sed s/QQ413_2/FALCON-QA/g | sed s/QQ027_2/FaeNNz/g |      sed s/QQ196_2/Grudinin/g |      sed s/QQ065_2/Jagodzinski-Cao-QA/g | sed s/QQ344_2/Kiharalab/g |      sed s/QQ067_2/LamoureuxLab/g |      sed s/QQ146_2/MASS1/g |      sed s/QQ415_2/MASS2/g |      sed s/QQ197_2/MESHI/g |      sed s/QQ289_2/MESHI-enrich-server/g |      sed s/QQ347_2/MESHI-server/g |      sed s/QQ113_2/MUFold-QA/g | sed s/QQ312_2/MUFold_server/g |      sed s/QQ243_2/MULTICOM-CONSTRUCT/g |      sed s/QQ023_2/MULTICOM-NOVEL/g |      sed s/QQ058_2/MULTICOM_CLUSTER/g |      sed s/QQ107_2/MUfold-QA2/g | sed s/QQ211_2/MUfold-QA-T/g | sed s/QQ275_2/ModFOLD7/g |      sed s/QQ213_2/ModFOLD7_cor/g |      sed s/QQ272_2/ModFOLD7_rank/g |      sed s/QQ373_2/ModFOLDclust2/g |      sed s/QQ209_2/PLU-Angular-QA/g | sed s/QQ134_2/PLU-Top-QA/g | sed s/QQ083_2/Pcomb/g |      sed s/QQ022_2/Pcons/g |      sed s/QQ044_2/ProQ2/g |      sed s/QQ187_2/ProQ3/g |      sed s/QQ139_2/ProQ3D/g |      sed s/QQ198_2/ProQ3D-CAD/g |      sed s/QQ267_2/ProQ3D-TM/g |      sed s/QQ360_2/ProQ3D-lDDT/g |      sed s/QQ440_2/ProQ4/g |      sed s/QQ334_2/RaptorX-Deep-QA/g | sed s/QQ220_2/SASHAN/g |      sed s/QQ135_2/SBROD/g |      sed s/QQ207_2/SBROD-plus/g |      sed s/QQ364_2/SBROD-server/g |      sed s/QQ194_2/UOSHAN/g |      sed s/QQ339_2/VoroM-QA-A/g | sed s/QQ030_2/VoroM-QA-B/g | sed s/QQ457_2/Wallner/g |       sed "s/\s+/ /g" > pertarget-$k.tsv
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
    Rscript bin/merge.R Z-sum-$k.tsv Z-corr-$k.tsv Z-pertarget-$k.tsv foo-$k.tsv
done

for k in TM CAD GDT_TS GDT_HA lDDT
do
    sed s/QA267_2/ProQ3D_TM/ foo-$k.tsv | sed s/QA457_2/Pcomb/ | sed s/QA139_2/ProQ3D/ | sed s/QA044_2/ProQ2/ | sed s/QA022_2/Pcons/ | sed s/QA187_2/ProQ3/  | sed s/QA360_2/ProQ3D_lDDT/ | sed s/QA198_2/ProQ3D_CAD/  | sed s/QA440_2/ProQ4/ |  sed s/QA359_2/3DCNN/g |     sed s/QA014_2/Bhattacharya-ClustQ/g |      sed s/QA102_2/Bhattacharya-Server/g |      sed s/QA170_2/Bhattacharya-SingQ/g |      sed s/QA471_2/CPClab/g |      sed s/QA349_2/Davis-EMAconsensus/g |      sed s/QA413_2/FALCON-QA/g | sed s/QA027_2/FaeNNz/g |      sed s/QA196_2/Grudinin/g |      sed s/QA065_2/Jagodzinski-Cao-QA/g | sed s/QA344_2/Kiharalab/g |      sed s/QA067_2/LamoureuxLab/g |      sed s/QA146_2/MASS1/g |      sed s/QA415_2/MASS2/g |      sed s/QA197_2/MESHI/g |      sed s/QA289_2/MESHI-enrich-server/g |      sed s/QA347_2/MESHI-server/g |      sed s/QA113_2/MUFold-QA/g | sed s/QA312_2/MUFold_server/g |      sed s/QA243_2/MULTICOM-CONSTRUCT/g |      sed s/QA023_2/MULTICOM-NOVEL/g |      sed s/QA058_2/MULTICOM_CLUSTER/g |      sed s/QA107_2/MUfold-QA2/g | sed s/QA211_2/MUfold-QA-T/g | sed s/QA275_2/ModFOLD7/g |      sed s/QA213_2/ModFOLD7_cor/g |      sed s/QA272_2/ModFOLD7_rank/g |      sed s/QA373_2/ModFOLDclust2/g |      sed s/QA209_2/PLU-Angular-QA/g | sed s/QA134_2/PLU-Top-QA/g | sed s/QA083_2/Pcomb/g |      sed s/QA022_2/Pcons/g |      sed s/QA044_2/ProQ2/g |      sed s/QA187_2/ProQ3/g |      sed s/QA139_2/ProQ3D/g |      sed s/QA198_2/ProQ3D-CAD/g |      sed s/QA267_2/ProQ3D-TM/g |      sed s/QA360_2/ProQ3D-lDDT/g |      sed s/QA440_2/ProQ4/g |      sed s/QA334_2/RaptorX-Deep-QA/g | sed s/QA220_2/SASHAN/g |      sed s/QA135_2/SBROD/g |      sed s/QA207_2/SBROD-plus/g |      sed s/QA364_2/SBROD-server/g |      sed s/QA194_2/UOSHAN/g |      sed s/QA339_2/VoroM-QA-A/g | sed s/QA030_2/VoroM-QA-B/g | sed s/QA457_2/Wallner/g |       sed "s/\s+/ /g" > average-$k.tsv
done


# best results


# We also need to find best method and best possible individual method


grep _1 QA237_2-CAD.tsv | sed "s/.*TS/TS/g" | gawk '{print $1}' > servers.txt

# best possible

for k in TM CAD GDT_TS GDT_HA lDDT
do
    #echo -n $k " " 
    for i in  `cat servers.txt`
    do
	echo -n $i " " 
	grep $i QA237_2-$k.tsv | gawk '{i++;sum+=$3};END{print i,sum,sum/i}' 
    done  | sort -rn -k 3 > servers-$k.tsv
done 


    
# best possible

for k in TM CAD GDT_TS GDT_HA lDDT
do
    echo -n $k " " 
    for i in  `cat targets.txt | sed "s/\-D1//g" | sort -u`
    do
	grep $i QA237_2-$k.tsv | sort -g -k 3 | tail -1 
    done  | gawk '{i++;sum+=$3};END{print i,sum,sum/i}' 
done > best.tsv





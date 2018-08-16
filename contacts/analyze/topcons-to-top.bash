

# just a reminder ...

wget http://topcons.net/static/result/rst_0DRMXv/rst_0DRMXv/query.result.txt
~/git/bioinfo-toolbox/arne/misc/topcons2top.pl query.result.txt > list.top
mkdir tmp
~/git/bioinfo-toolbox/arne/PfamProteomes/extract-seq.py list.top tmp/
for i in tmp/*fa ; do j=`basename $i .fa` ; k=`echo $j | sed "s/.[0-9]..._[a-zA-Z]//g" ` ; mv $i $k/$j.top ; done
rmdit tmp

 for i in *.pdb ; do j=`basename $i .pdb`  ; gawk '{print $2, $4, $35}' $j.txt ; done | grep -E "^T" | sort  | sed s/-D1//g > data.txt
for i in *.pdb ; do j=`basename $i .pdb`  ;  cat "$j/${j}QA339_2"  | grep -E "^T[0-9]" | gawk '{print $1,$2}' ;  done | sort > VoroMQA.txt
for i in *.pdb ; do j=`basename $i .pdb`  ;  cat "$j/${j}QA139_2"  | grep -E "^T[0-9]" | gawk '{print $1,$2}' ;  done | sort > ProQ3D.txt
for i in *.pdb ; do j=`basename $i .pdb`  ;  cat "$j/${j}QA440_2"  | grep -E "^T[0-9]" | gawk '{print $1,$2}' ;  done | sort > ProQ4.txt
 for i in *.pdb ; do j=`basename $i .pdb`  ;  cat "$j/${j}QA196_2"  | grep -E "^T[0-9]" | gawk '{print $1,$2}' ;  done | sort > Ornate.txt
for i in *.pdb ; do j=`basename $i .pdb`  ;  cat "$j/${j}QA359_2"  | grep -E "^T[0-9]" | gawk '{print $1,$2}' ;  done | sort > 3DCNN.txt

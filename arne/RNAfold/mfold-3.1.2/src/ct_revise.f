      read(7,6060,end=999) n,ctlabel
 6060 format(i5,a50)
      do i = 1,n
         basepr(i) = 0
      enddo
      do i = 1,n
         read(7,6070,end=998,err=997) record
 6070    format(a80)
         read(record,6080,err=997) k,seq(i)
 6080    format(i5,1x,a1)
         read(record(8:80),*,err=997) itmp,itmp,basepr(i),hstnum(i)
         if (i.ne.k) go to 997
         if (basepr(i).gt.0) basepr(basepr(i)) = i
      enddo

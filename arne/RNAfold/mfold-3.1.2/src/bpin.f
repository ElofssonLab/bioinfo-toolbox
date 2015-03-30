      subroutine bpin(n,l,bp,indeof)
c
c     (FS382)
c
      character*70 l
      character*80 ctrec
      integer*4 bp(1000)
c
c     READ BASE PAIR INFORMATION FOR ONE STRUCTURE
c
      indeof = 0
      read (2,101,end=20) n,l
      do i=1,n
         read (2,102,end=20) ctrec
         read (ctrec(1:5),*,end=90) k
         read (ctrec(8:80),*,end=90) itmp,itmp,bp(k)
      enddo   
      return
c
c     EXPECTED END OF FILE
   20 indeof = 1
      return
c
c     UNEXPECTED END OF FILE
   90 write (6,190)
      stop
c
 101  format (i5,1x,a70)
 102  format (a80)
 190  format (/' Unexpected end of file: no data following label')
      end

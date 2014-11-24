      implicit integer (a-z)

      include 'maxn2.inc'
      parameter (maxtloops=100)
      parameter (maxtriloops=50)

      integer vst(maxn*maxn)
      integer nsave(2),asint3(6,6,5,5,5),asint5(6,6,5,5,5,5,5),bulge(30),
     .        dangle(5,5,5,2),eparam(9),hairpin(30),inter(30),poppen(4),
     .        sint2(6,6,5,5),sint4(6,6,5,5,5,5),sint6(6,6,25,5,5,5,5),
     .        stack(5,5,5,5),tloop(maxtloops,2),triloop(maxtriloops,2),
     .        tstkh(5,5,5,5),tstki(5,5,5,5),maxbp

      character*1  lorc
      character*50 argline,sfile,seqlab,pnumfile
      real r1,r2,prec

      common /block/ n,nsave,vst
      data infinity/999999/

c     Determine minimum helix length in suboptimal foldings

      if (iargc().lt.1) then
         lmin = 1
      else
         call getarg(1,argline)
         read(argline,*) lmin
         if (lmin.lt.1) lmin = 1
      endif

c     Prompt for SAVE file name.

      write (6,3)
3     format(' Enter save file name (default fold.sav)')
      read (5,15,end=999) sfile
      if (sfile.eq.'         ') sfile= 'fold.sav'
      open(30,err=999,file=sfile,status='old',form='UNFORMATTED')

c     Read save file
 
      read(30,err=999) n,break,nsave,vmin,listsz,seqlab,lorc,maxbp
      read(30,err=999) asint3,asint5,bulge,dangle,eparam,hairpin,
     . inter,sint2,sint4,sint6,stack,tloop,triloop,tstkh,tstki,prec
      read(30,err=999) (vst(i),i=1,n*n)
      close(30)

c     Get info for plot file

      prec = 100.0
      r1 = float(vmin)/prec

      write (6,10)
 10   format(' Enter p-num file name '/' >')
      read  (5,15,err=999,end=999) pnumfile
 15   format(a50)

      open(unit=7,file=pnumfile,status='unknown')

 31   write (6,35) r1
 35   format('vmin = ',f7.1,'. Enter energy increment >')
      read  (5,fmt=*,err=999,end=999) r2
      vinc = r2*prec

c     Filter plot file - suboptimal helices must have length lmin or more
      do l = 2,2*n-2
         if (l.le.n) then
            i = 1
            j = l
         else
            i = l + 1 - n
            j = n
         endif
c        Stop at diagonal
         sum = 0
         do while (j-i.gt.0)
            dg = vst((n-1)*(j-1) + i) + vst((n-1)*(i-1) + j+n)
            if (dg.le.vmin+vinc) then
               sum = sum + 1
            else
               if (sum.gt.0.and.sum.lt.lmin) then
                  do isum = 1,sum
                     ip = i - isum
                     jp = j + isum
                     dgp = vst((n-1)*(jp-1)+ip) + vst((n-1)*(ip-1)+jp+n)
                     if (dgp.gt.vmin) then
                        vst((n-1)*(jp-1)+ip) = infinity
                        vst((n-1)*(ip-1)+jp+n) = infinity
                     endif
                  enddo
                  sum = 0
               endif
            endif
            i = i + 1
            j = j - 1
         enddo
      enddo

      do i = 1,n
         sum = 0
         do j = 1,n
            if (j.lt.i) then
               if (vst((n-1)*(j-1)+i)+vst((n-1)*(i-1)+j+n).le.vmin+vinc) sum = sum + 1
            elseif (j.gt.i) then
               if (vst((n-1)*(i-1)+j)+vst((n-1)*(j-1)+i+n).le.vmin+vinc) sum = sum + 1
            endif
         enddo
         write (7,*) i,sum
      enddo

999   stop
      end


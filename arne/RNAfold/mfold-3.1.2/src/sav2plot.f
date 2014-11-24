      implicit integer (a-z)

      include 'maxn.inc'
      parameter (maxtloops=100)
      parameter (maxtriloops=50)
      integer*2 vst(maxn*maxn)
      integer nsave(2),asint3(6,6,5,5,5),asint5(6,6,5,5,5,5,5),bulge(30),
     .        dangle(5,5,5,2),eparam(9),hairpin(30),inter(30),poppen(4),
     .        sint2(6,6,5,5),sint4(6,6,5,5,5,5),sint6(6,6,25,5,5,5,5),
     .        stack(5,5,5,5),tloop(maxtloops,2),triloop(maxtriloops,2),
     .        tstkh(5,5,5,5),tstki(5,5,5,5),maxbp

      character*1 lorc
      character*50 sfile,seqlab,plotfile
      real r1,r2,prec

      common /block/ n,nsave,vst

c     Prompt for SAVE file name.

      write (6,3)
3     format(' Enter save file name (default fold.sav)')
      read (5,15,end=99) sfile
      if (sfile.eq.'         ') sfile= 'fold.sav'
      open(30,err=999,file=sfile,status='old',form='unformatted')

c     Read save file
 
      read(30,err=999) n,break,nsave,vmin,listsz,seqlab,lorc,maxbp
      read(30,err=999) asint3,asint5,bulge,dangle,eparam,hairpin,
     . inter,sint2,sint4,sint6,stack,tloop,triloop,tstkh,tstki
      read(30,err=999) (vst(i),i=1,n*n)
      close(30)
 
c     Get info for plot file

      prec = 10.0
      r1 = float(vmin)/prec

      write (6,10)
 10   format(' Enter helix file name '/' >')
      read  (5,15,err=99,end=99) plotfile
 15   format(a50)

      open(unit=7,file=plotfile,status='unknown')

 31   write (6,35) r1
 35   format('vmin = ',f7.1,'. Enter energy increment >')
      read  (5,fmt=*,err=99,end=99) r2
      vinc = r2*prec

      if (vinc.lt.0) then
         go to 31
      elseif (vinc.gt.0) then
         write (6,40)
 40      format('Enter number of levels >')
         read  (5,fmt=*,err=99,end=99) levels
      else
         write (6,42)
 42      format('Energy increment = 0. Number of levels is set to 1.')
         levels = 1
      endif
c
      write(7,43)
43    format('   level  length istart jstart energy')
      if(levels.gt.1) icrit = vinc/(levels-1)
      do diag = 1,2*n-1
         flagid = 0
         i = (diag+1)/2
         j = (diag+2)/2
         do while (i.ge.1.and.j.le.n)
            level = 0
            check = v(i,j)+v(j,n+i)
            if(levels.gt.1) then
               k = (check + icrit - 1 - vmin)/icrit + 1
               if(k.le.levels) then
                  level = k
               else
                  if(check.le.vmin+vinc) level = levels
               endif
            else
               if(check.le.vmin+vinc) level = 1
            end if
            if (level.gt.0) then
               if (flagid.ne.level) then
                  if (flagid.gt.0) call plotout(i+1,j-1,
     .               istart-i,flagid,v(i+1,j-1)+v(j-1,n+i+1),0)
                  flagid = level
                  istart = i
                  jstart = j
               endif
            elseif (flagid.gt.0) then
               call plotout(i+1,j-1,istart-i,flagid,
     .                    v(i+1,j-1)+v(j-1,n+i+1),0)
               flagid = 0
            endif
            if (i.eq.1.or.j.eq.n) then
               if (flagid.gt.0) call plotout(i,j,
     .                          istart-i+1,flagid,v(i,j)+v(j,n+i),0)
            endif
            i = i - 1
            j = j + 1
         enddo
      enddo
      call plotout(0,0,0,0,0,1)
      close(7)
99    stop
999   write(6,*) 'STOP: Premature end of save file. '
      call exit(1)
      end

c     Used to recall values of V which are actually stored in VST.
      function v(i,j)
      implicit integer (a-z)

      include 'maxn.inc'
      integer nsave(2)
      integer*2 vst(maxn*maxn)
      common /block/ n,nsave,vst

        v = vst((n-1)*(i-1)+j)

      return
      end

      subroutine plotout(istart,jstart,length,level,energy,dump)
      implicit integer (a-z)

      include 'maxn.inc'
      integer nsave(2)
      integer*2 vst(maxn*maxn)
      common /block/ n,nsave,vst

      integer energy,dump,stack(5,5),k
      data k/0/

      if (dump.eq.0) then
         k = k + 1
         stack(1,k) = level
         stack(2,k) = length
         stack(3,k) = nsave(1) + istart - 1
         stack(4,k) = nsave(1) + jstart - 1
         stack(5,k) = energy
      endif

      if (k.eq.5.or.(k.gt.0.and.dump.eq.1)) then
         write(7,fmt='(5i7)') ((stack(i,j),i=1,5),j=1,k)
         k = 0
      endif

      return
      end

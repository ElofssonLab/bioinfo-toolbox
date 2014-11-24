      SUBROUTINE fileplot2
      INCLUDE 'rna.inc'    
      CHARACTER*30 plotfile
      REAL r1,r2

      r1 = float(vmin)/10.

      write (6,10)
 10   format(' Enter helix file name '/' >')
      read  (5,15,err=99,end=99) plotfile
 15   format(a30)

      open(unit=7,file=plotfile,status='unknown')

 31   write (6,35) r1
 35   format('vmin = ',f7.1,'. Enter energy increment >')
      read  (5,fmt=*,err=99,end=99) r2
      vinc = r2*10.

      if (vinc.lt.0) then
         go to 31
      elseif (vinc.gt.0) then
         write (6,40)
 40      format('Enter number of levels >')
         read  (5,fmt=*,err=99,end=99) levels
      else
         write (6,42)
 42      format('Energy increment = 0. Number of levels is set to 1.')
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
               if(k.le.levels) level = k
            else
               if(check.le.vmin+vinc) level = 1
            end if
            if (level.gt.0) then
               if (flagid.ne.level) then
                  if (flagid.gt.0) call plotout(hstnum(i+1),hstnum(j-1),
     .               istart-i,flagid,v(i+1,j-1)+v(j-1,n+i+1),0)
                  flagid = level
                  istart = i
                  jstart = j
               endif
            elseif (flagid.gt.0) then
               call plotout(hstnum(i+1),hstnum(j-1),istart-i,flagid,
     .                    v(i+1,j-1)+v(j-1,n+i+1),0)
               flagid = 0
            endif
            if (i.eq.1.or.j.eq.n) then
               if (flagid.gt.0) call plotout(hstnum(i),hstnum(j),
     .                          istart-i+1,flagid,v(i,j)+v(j,n+i),0)
            endif
            i = i - 1
            j = j + 1
         enddo
      enddo
      call plotout(0,0,0,0,0,1)
      close(7)
      return
c     error return.
 99   print *,' Read error or end of file...'
c     call ringbe
      return
      end
      subroutine plotout(istart,jstart,length,level,energy,dump)
      integer energy,dump,stack(5,5),k
      data k/0/

      if (dump.eq.0) then
         k = k + 1
         stack(1,k) = level
         stack(2,k) = length
         stack(3,k) = istart
         stack(4,k) = jstart
         stack(5,k) = energy
      endif

      if (k.eq.5.or.(k.gt.0.and.dump.eq.1)) then
         write(7,fmt='(5i7)') ((stack(i,j),i=1,5),j=1,k)
         k = 0
      endif

      return
      end

c     Program efn - RNA folding energy function for linear or
c                   circular RNA
c
c     Input - a CT file containing one or more structures
c           - energy files : stack.???, tstack.??? etc.
c     Output - the energy of each of the foldings according
c              to the input datasets
c             - optional text output and thermodynamic DETails output
c     Objective - this program serves as a check on energy
c                 calculations performed by nafold
c               - energies can be calculated according to different
c                 rules (e.g. temperature change)
c
c     force() = 4 option not yet implemented. Thus WXYZ are recognized as
c     ACGU (respectively) and are allowed to pair anywhere.
c     M. Zuker - 11/15/97
c
      include 'efn.inc'
      real energy
      character*1 ans,lorc
      character*80 ctrec
      logical flag
      data lorc/'c'/
 
c     Linear or circular RNA
 
      write(6,1005)
1005  format(' Linear or circular RNA ? (L [default] or C) : ')
      read(5,1007) ans
1007  format(a1)
      if (ans.eq.'C'.or.ans.eq.'c')  then
         ans = 'c'
      else
         ans = 'l'
      endif
 
c     Initial setup for run : get CT file name and open file.

5     write(6,1010)
1010  format(' Enter CT file name (default=fold.ct) : ')
      read(5,1020,end=999) ctnam
      if (ctnam.eq.'     ') ctnam = 'fold.ct'
1020  format(a50)
      open(unit=7,file=ctnam,status='old',err=5)
      iter = 0
c
c     Determine output specifications.
c
      call outputs
c
c     Read energy information
c
      call enefiles
      call ergread
c
c     Read foldings one by one.
c
 10   read(7,*,end=999) n,ctlabel
      do i = 1,n
         basepr(i) = 0
      enddo
      do i = 1,n
         read(7,1035,end=998,err=997) ctrec
 1035    format(a80)
         read(ctrec,*,err=997) k,seq(i)
         kb = index(ctrec,seq(i))
         read(ctrec(kb+1:80),*,err=997) itmp,itmp,basepr(i),hstnum(i)
         if (i.ne.k) go to 997
         if (basepr(i).gt.0) basepr(basepr(i)) = i
      enddo
      iter = iter + 1
      if (ans.eq.lorc) then
         break = 3*n
      else
         break = n
      endif
c
c     Check for knots
c
      do i = 1,n
         if (basepr(i).ne.0) then
            j = max(i,basepr(i))
            ip = min(i,basepr(i))
            if (basepr(ip).ne.j.or.basepr(j).ne.ip) then
               write(6,1041) hstnum(ip),hstnum(j)
               call exit(1)
            endif
            do k = ip+1,j-1
               if (basepr(k).ne.0) then
                  l = basepr(k)
                  if (l.le.ip.or.l.ge.j) then
                     write(6,1042) hstnum(ip),hstnum(j),hstnum(k),hstnum(l)
                     call exit(1)
                  endif
               endif
            enddo
         endif
      enddo
1041  format(' Base pair ',i5,'.',i5,39H has 5' or 3' end in another base pair.)
1042  format(' Base pair ',i5,'.',i5,' conflicts with ',i5,'.',i5)
c
c     Process sequence to find base types.
c
      call process
c
c     Double up if circular
c
      if (ans.eq.lorc) then
         do i = 1,n
            basepr(i+n) = 0
            seq(i+n) = seq(i)
            force(i+n) = force(i)
            hstnum(i+n) = hstnum(i)
            numseq(i+n) = numseq(i)
         enddo
      endif
c
c     Compute energy of folding
c
      if (ans.eq.lorc) then
c
c     Find a base pair that closes a hairpin loop.
c     First find a base pair.
c
         i = 1
         do while (basepr(i).eq.0.and.i.le.n)
            i = i + 1
         enddo
         if (i.gt.n) then
c
c     No structure !
c
            write(6,1050) iter,ctlabel
1050        format(' ',i5,'.  ',a50/' No structure found.'/)
            go to 10
         endif
c
         flag = .false.
         do while (.not.flag)
            j = i + 1
            do while (basepr(j).eq.0.and.j.lt.basepr(i))
               j = j + 1
            enddo
            if(j.eq.basepr(i)) then
               flag = .true.
            else
               i = j
            endif
         enddo
c
c     i.basepr(i) close a hairpin loop
c
         energy = float(erg4(i,basepr(i)))
c
c     Circular permutation of structure : start at basepr(i)
c
         do j = 1,i
            if (basepr(j).gt.0) then
               if (basepr(j).lt.i) then
                  basepr(n+j) = basepr(j) + n
                  basepr(n+basepr(j)) = n + j
               else
                  basepr(basepr(j)) = n + j
                  basepr(n+j) = basepr(j)
               endif
            endif
         enddo
         call efn(nergy,basepr(i),n+i,1)
         energy = ( energy + float(nergy) )/100.0
c
c     Restore folding
c
         do j = n + 1,n + i
            if(basepr(j).le.n) basepr(basepr(j)) = j - n
         enddo
 
      else
c
c     linear molecule
c
         call efn(nergy,1,n,1)
         energy = float(nergy)/100.0
c
      endif
c
      write(6,2010) iter,ctlabel,energy
2010  format(' ',i5,'.  ',a50/' Computed energy = ',f8.1/)
c
c     Text output.
c
      if(cntrl(2).eq.1.or.cntrl(2).eq.3) call linout(1,n,energy,0,0,err)
c
c     Read next structure
c
      go to 10
c
c     Error in CT file.
c
997   write(6,9080)
9080  format('STOP: Error in CT file.')
      call exit(1)
c
c     Error - incomplete CT file.
c
998   write(6,9090)
9090  format('STOP:  Premature end of CT file.')
      call exit(1)
c
999   stop
      end
 
      subroutine efn(e,ii,ji,ext_open)
      include 'efn.inc'
 
      e = 0
      error = 0
      eold = 0

      call initst
      call push(ii,ji,ext_open,0)
 
100   stz = pull(i,j,open,null)
      if (stz.ne.0) return
      write(22,102) hstnum(i),hstnum(j)
 102  format('New segment: (',i5,',',i5,').')
 
c     Do I and J base-pair with one another?
      if (basepr(i).eq.j) goto 300
 
      if (open.eq.0) then
 
         do while (basepr(i).eq.0.and.basepr(i+1).eq.0)
c            Whittle away from the 5' end.
             i = i + 1
             e = e + eparam(6)
             if (i.ge.j-1) goto 100
         enddo
         do while (basepr(j).eq.0.and.basepr(j-1).eq.0)
c            Whittle away from the 3' end.
             j = j - 1
             e = e + eparam(6)
             if (i.ge.j-1) goto 100
         enddo
 
         if (basepr(i).eq.0.and.basepr(i+1).gt.i+1) then
c            I dangles over I+1,basepr(I+1).
             e = e + min0(0,erg6(basepr(i+1),i+1,i,2)) + eparam(6)
             i = i + 1
         endif
         if (basepr(j).eq.0.and.basepr(j-1).ne.0.and.basepr(j-1).lt.j-1) then
c            J dangles over basepr(J-1),J-1.
             e = e + min0(0,erg6(j-1,basepr(j-1),j,1)) + eparam(6)
             j = j - 1
         endif
 
      else
 
         do while (basepr(i).eq.0.and.basepr(i+1).eq.0)
c            Whittle away from the 5' end.
             i = i + 1
             if (i.ge.j-1) goto 100
         enddo
         do while (basepr(j).eq.0.and.basepr(j-1).eq.0)
c            Whittle away from the 3' end.
             j = j - 1
             if (i.ge.j-1) goto 100
         enddo
 
         if (basepr(i).eq.0.and.basepr(i+1).gt.i+1) then
c            I dangles over I+1,basepr(I+1).
             e = e + min0(0,erg6(basepr(i+1),i+1,i,2)) 
             i = i + 1
         endif
         if (basepr(j).eq.0.and.basepr(j-1).ne.0.and.basepr(j-1).lt.j-1) then
c            J dangles over basepr(J-1),J-1.
             e = e + min0(0,erg6(j-1,basepr(j-1),j,1)) 
             j = j - 1
         endif
      endif
 
      if (basepr(i).ne.j) then
c        Cannot chop away at the ends any more and still the ends do not
c        base-pair with one another. Structure MUST bifucate (OPEN).
         k = basepr(i)
         kp = basepr(j)
         if (k .ge. kp) then
            error = 50
            write (6,1010) hstnum(i),hstnum(k),hstnum(kp),hstnum(j)
1010        format(' Knot : ',i5,'.',i5,' conflicts with ',i5,'.',i5)
            e = infinity
            return
         endif
         if (basepr(k+1).ne.0) then
            call push(i,k,open,0)
            call push(k+1,j,open,0)
            write(22,1012) hstnum(i),hstnum(k),hstnum(k+1),hstnum(j)
 1012       format('Bifurcation: (',i5,',',i5,') and (',i5,',',i5,').')
         elseif (basepr(k+2).eq.0) then
            call push(i,k+1,open,0)
            call push(k+2,j,open,0)
            write(22,1012) hstnum(i),hstnum(k+1),hstnum(k+2),hstnum(j)
         elseif (erg6(k,i,k+1,1).le.erg6(basepr(k+2),k+2,k+1,2)) then
            call push(i,k+1,open,0)
            call push(k+2,j,open,0)
            write(22,1012) hstnum(i),hstnum(k+1),hstnum(k+2),hstnum(j)
         else
            call push(i,k,open,0)
            call push(k+1,j,open,0)
            write(22,1012) hstnum(i),hstnum(k),hstnum(k+1),hstnum(j)
         endif
         goto 100
      endif
c
c     Add penalty for accessible base pair in a multi-loop.
c
 300  if(open.eq.0) e = e + eparam(9)
      e = e + au_pen(i,j)
 310  open = 0
c     Write energy of loop (or stack) closed by (i,j) and the PREVIOUS
c     base pair (if none, value = 0)

      write(22,2020) seq(i),seq(j),i,j,e,e - eold
 2020 format('Base pair ',2a1,' (',i4,',',i4,'), Energy = ',i7,' E-inc = ',i6)
      eold = e

c     Perhaps I,J stacks over I+1,J-1?
      if (basepr(i+1).eq.j-1) then
         e = e + erg2(i,j)
         i = i + 1
         j = j - 1
         goto 310
      endif
 
      sum = 0
      k = i + 1
      do while (k.lt.j)
         if (basepr(k).gt.k) then
            sum = sum + 1
            ip  = k
            k   = basepr(k) + 1
            jp  = k - 1
            if (k.gt.j) then
               write (6,*) 'ERROR'
               error = 51
               e = infinity
               return
            endif
         elseif (basepr(k).eq.0) then
            k = k + 1
         endif
      enddo
      if (sum.eq.0) then
         e = e + erg4(i,j)
         write(22,104) hstnum(i),hstnum(j),erg4(i,j)
 104     format('Hairpin loop closed by: (',i5,',',i5,') E = ',i6)
         goto 100
      elseif (sum.eq.1) then
         e = e + erg3(i,j,ip,jp)
         i = ip
         j = jp
         goto 310
      else
c     Closed bifurcation - multi-branch loop
c     The free energy of this loop is em
         em = eparam(5) + eparam(9) + au_pen(i,j)
         k = i + 1
         do while (k.lt.j)
            if (basepr(k).eq.0) then
               predangle = 0
               if (basepr(k-1).gt.0) predangle = min0(0,erg6(k-1,basepr(k-1),k,1))
               postdangle = 0
               if (basepr(k+1).gt.0) postdangle = min0(0,erg6(basepr(k+1),k+1,k,2))
               em = em + eparam(6) + min0(predangle,postdangle)
            else
               em = em + eparam(9) + au_pen(k,basepr(k))
               k = basepr(k)
            endif
         k = k + 1
         enddo
         write(22,105) hstnum(i),hstnum(j),em
 105     format('Multi-branch loop closed by: (',i5,',',i5,') E = ',i6)
         is = i+1
         js = j-1
         e = e + eparam(5) + eparam(9) + au_pen(i,j)
         if (basepr(i+1).eq.0.and.basepr(i+2).ne.0) then
            if (erg6(i,j,i+1,1).le.erg6(basepr(i+2),i+2,i+1,2)) then
               is = i+2
               e = e + min0(0,erg6(i,j,i+1,1)) + eparam(6) 
            endif
         endif
         if (basepr(i+1).eq.0.and.basepr(i+2).eq.0) then
            is = i+2
            e = e + min0(0,erg6(i,j,i+1,1)) + eparam(6)
         endif
         if (basepr(j-1).eq.0.and.basepr(j-2).ne.0) then
            if (erg6(i,j,j-1,2).le.erg6(j-2,basepr(j-2),j-1,1)) then
               js = j-2
               e = e + min0(0,erg6(i,j,j-1,2)) + eparam(6)
            endif
         endif
         if (basepr(j-1).eq.0.and.basepr(j-2).eq.0) then
            js = j-2
            e = e + min0(0,erg6(i,j,j-1,2)) + eparam(6)
         endif
         call push (is,js,0,0)
         goto 100
      endif
      end

c     Energy funtions.
c     ergk (k=2,3,4,5) is the energy of a loop closed by i,j (new numbering).
c     ip,jp is the other closing base-pair when MODE = 2 or 3.
c     The ends of the sequence cannot be contained in a hairpin, bulge
c     or interior loop. By convention, the ends of the sequence are
c     put into a special kind of multi-loop. This can be called an
c     exterior loop or an open multi-loop.
c
      function erg2(i,j)
      include 'efn.inc'
      integer e(4)
      integer*2 tlink,tlptr,itemp
 
      erg2 = 0
  
c     Molecule is not circular. n is not covalently bonded to n+1.
      if (i.eq.break.or.j.eq.break+1) then
         erg2 = infinity
         return
      endif
c     Stacking energy.
      erg2 = erg2 + stack(numseq(i),numseq(j),numseq(i+1),numseq(j-1))
     .          + eparam(1)
      return
      end

      function erg3(i,j,ip,jp)
      include 'efn.inc'
      dimension e(4)
      integer*2 tlink,tlptr,itemp

      erg3 = 0
      if ((i.le.break.and.ip.gt.break).or.(jp.le.break.and.j.gt.break)) then
c          Loop is not allowed to contain the ends of the sequence.
           erg3 = infinity
           return
      endif
c
      size1 = ip - i - 1
      size2 = j - jp - 1
      if (size1.eq.0.or.size2.eq.0) then
          size = size1+size2
c         Bulge loop energy.
          if (size.gt.1) erg3 = erg3 + au_pen(i,j) + au_pen(ip,jp)
          if (size.eq.1) then
             erg3 = erg3 + stack(numseq(i),numseq(j),numseq(ip),numseq(jp))
     .              + bulge(size) + eparam(2)
          elseif (size.gt.30) then
             loginc = nint(prelog*log((float(size)/30.0)))
             erg3 = erg3 + bulge(30) + loginc + eparam(2)
          else
             erg3 = erg3 + bulge(size) + eparam(2)
          endif
          return
      else
          size = size1+size2
          lopsid = abs((size1-size2))
c         Interior loop.
          if (size.gt.30) then
             loginc = nint(prelog*log((float(size)/30.0)))
             if ((size1.eq.1.or.size2.eq.1).and.eparam(16).eq.1) then
c     Apply GAIL Rule (Grossly Asymmetric Interior Loop Rule)
                erg3 = erg3 + tstki(numseq(i),numseq(j),1,1)
     .           + tstki(numseq(jp),numseq(ip),1,1)
     .           + inter(30) + loginc + eparam(3)
     .           + min0(maxpen,(lopsid*poppen(min0(4,size1,size2))))
             else
                erg3 = erg3 + tstki(numseq(i),numseq(j),numseq(i+1),numseq(j-1))
     .           + tstki(numseq(jp),numseq(ip),numseq(jp+1),numseq(ip-1))
     .           + inter(30) + loginc + eparam(3)
     .           + min0(maxpen,(lopsid*poppen(min0(4,size1,size2))))
             endif

c   Kevin's changes start here:
c#        elseif(lopsid.eq.1.and.(size.eq.3.or.size.eq.5))then
          elseif(lopsid.eq.1.and.(size.eq.3))then
c   Asymmetric interior loop with size1 < size2
             if(size1.lt.size2)then
                if((numseq(i)+numseq(j)).eq.5)then
                   a=numseq(i)
                elseif(numseq(i).eq.3.and.numseq(j).eq.4)then
                   a=5
                elseif(numseq(i).eq.4.and.numseq(j).eq.3)then 
                   a=6
                endif

                if((numseq(ip)+numseq(jp)).eq.5)then
                   b=numseq(ip)
                elseif(numseq(ip).eq.3.and.numseq(jp).eq.4)then
                   b=5
                elseif(numseq(ip).eq.4.and.numseq(jp).eq.3)then
                   b=6
                endif

c       Size = 3
                if(size.eq.3)then
                   erg3 = erg3 + eparam(3)+asint3(a,b,numseq(i+1),numseq(j-1),
     .             numseq(jp+1))
                   return

c       Size = 5
c#              elseif(size.eq.5)then
c#                 erg3 = erg3 + eparam(3)+asint5(a,b,numseq(i+1),numseq(ip-1),
c#   .             numseq(j-1),numseq(j-2),numseq(jp+1))
c#                 return
                endif

c   Asymmetric interior loop with size1 > size2
             else
                if((numseq(jp)+numseq(ip)).eq.5)then
                   a=numseq(jp)
                elseif(numseq(jp).eq.3.and.numseq(ip).eq.4)then
                   a=5
                elseif(numseq(jp).eq.4.and.numseq(ip).eq.3)then
                   a=6
                endif

                if((numseq(j)+numseq(i)).eq.5)then
                   b=numseq(j)
                elseif(numseq(j).eq.3.and.numseq(i).eq.4)then
                   b=5
                elseif(numseq(j).eq.4.and.numseq(i).eq.3)then
                   b=6
                endif

c       Size = 3
                if(size.eq.3)then
                   erg3 = erg3 + eparam(3)+asint3(a,b,numseq(jp+1),
     .             numseq(ip-1),numseq(i+1))
                   return

c       Size = 5
c#              else
c#                 erg3 = erg3 + eparam(3)+asint5(a,b,numseq(jp+1),
c#   .             numseq(j-1),numseq(ip-1),numseq(ip-2),numseq(i+1))
c#                 return
                endif
             endif

c#        elseif(lopsid.ne.0.or.size.gt.6)then
          elseif(lopsid.ne.0.or.size.gt.4)then

             if ((size1.eq.1.or.size2.eq.1).and.eparam(16).eq.1) then
c     Apply GAIL Rule (Grossly Asymmetric Interior Loop Rule)
                erg3 = erg3 + tstki(numseq(i),numseq(j),1,1)
     .           + tstki(numseq(jp),numseq(ip),1,1) + inter(size) + eparam(3)
     .           + min0(maxpen,(lopsid*poppen(min0(4,size1,size2))))
             else 
                erg3 = erg3 + tstki(numseq(i),numseq(j),numseq(i+1),numseq(j-1))
     .           + tstki(numseq(jp),numseq(ip),numseq(jp+1),numseq(ip-1))
     .           + inter(size) + eparam(3)
     .           + min0(maxpen,(lopsid*poppen(min0(4,size1,size2))))
             endif

          else
             test1=numseq(i)+numseq(j)
          
             if(test1.eq.5)then
                lft=numseq(i)
             elseif(test1.eq.7.and.(numseq(i).eq.4.or.numseq(i).eq.3)) then
                lft=numseq(i)+2
             else
                erg3 = infinity
                return
             endif

             test2=numseq(ip)+numseq(jp)

             if(test2.eq.5)then
                rt=numseq(ip)
             elseif(test2.eq.7.and.(numseq(ip).eq.4.or.numseq(ip).eq.3)) then
                rt=numseq(ip)+2
             else
                erg3 = infinity
                return
             endif

c   Symmetric interior loop of size 2
             if(size.eq.2) then
                erg3 = erg3 + sint2(lft,rt,numseq(i+1),numseq(j-1)) + eparam(3)
                return

c   Symmetric interior loop of size 4
             elseif(size.eq.4) then
                erg3 = erg3 + sint4(lft,rt,numseq(i+1),numseq(j-1),
     .          numseq(ip-1),numseq(jp+1)) + eparam(3)
                return
c   Symmetric interior loop of size 6
c#           elseif(size.eq.6) then
c#              erg3 = erg3 + sint6(lft,rt,(numseq(i+1)-1)*5+numseq(j-1),
c#   .          numseq(i+2),numseq(j-2),numseq(ip-1),numseq(jp+1)) + eparam(3)
c#              return
             endif
          endif
         return
      endif
c   Kevin's modifications end --------------------
      end

      function erg4(i,j)
      include 'efn.inc'
      integer e(4),strand,d3,d5
      integer*2 tlink,tlptr,itemp

      erg4 = 0

      if (i.le.break.and.j.gt.break) then
c         Hairpin loop must not contain the ends of the sequence.
          erg4 = infinity
          return
      endif
      size = j-i-1
      if ((size.eq.3).and.seq(hstnum(i+1)).eq.' ') then
c        Closed excision
         erg4 = 0
         return
      endif

c      if ((size.eq.3).and.seq(hstnum(i+1)).eq.'L') then
cc        Treat as a dimer - done for John SantaLucia
c         erg4 = eparam(15) + au_pen(i,j)
c         return
c      endif
      strand = 0
      do k=i+1,j-1
         if (seq(hstnum(k)).eq.'L') strand = strand + 1
      enddo
      if (strand.eq.3) then
c        Treat as a dimer - done for John SantaLucia (and hybridization server)
         erg4 = eparam(15) + au_pen(i,j)
         d3 = dangle(numseq(i),numseq(j),numseq(i+1),1)
         if (seq(hstnum(i+1)).ne.'L'.and.d3.lt.0) erg4 = erg4 + d3
         d5 = dangle(numseq(i),numseq(j),numseq(j-1),2)
         if (seq(hstnum(j-1)).ne.'L'.and.d5.lt.0) erg4 = erg4 + d5
         return
      endif

c     Check for poly-C loop
      k = i + 1
      do while (numseq(k).eq.2.and.k.lt.j)
         k = k + 1
      enddo
      if (k.eq.j) then
         if (size.eq.3) then
            erg4 = eparam(14)
         else
            erg4 = eparam(13) + size*eparam(12)
         endif
      endif

c     Check for GGG hairpin loop and worry about circular RNA/DNA
      if ((i.gt.2.and.j.le.n).or.i-break.gt.2) then
         if (numseq(i).eq.3.and.numseq(i-1).eq.3.and.numseq(i-2).eq.3.
     .        and.numseq(j).eq.4) erg4 = erg4 + eparam(11)
      else
         if (ans.eq.lorc) then
            if (i.eq.1) then
               if (numseq(1).eq.3.and.numseq(n).eq.3.and.numseq(n-1).eq.3.
     .            and.numseq(j).eq.4) erg4 = erg4 + eparam(11)
            elseif (i.eq.2) then
               if (numseq(2).eq.3.and.numseq(1).eq.3.and.numseq(n).eq.3.
     .        and.numseq(j).eq.4) erg4 = erg4 + eparam(11)
            else
               if (numseq(i).eq.3.and.numseq(i-1).eq.3.and.numseq(i-2).eq.3.
     .        and.numseq(j).eq.4) erg4 = erg4 + eparam(11)
            endif
         endif
      endif

      if (size.gt.30) then
         loginc = nint(prelog*log((float(size)/30.0)))
         erg4 = erg4 + tstkh(numseq(i),numseq(j),numseq(i+1),numseq(j-1))
     .          + hairpin(30) + loginc + eparam(4)
      elseif (size .lt. 4) then
c
c* Special case for hairpin of 3
c
c------------------Kevin's modifications begin---------------------------
         tlink=0
         if(size.eq.3)then
            key=(((numseq(i+4)*8+numseq(i+3))*8+numseq(i+2))*8+
     +       numseq(i+1))*8+numseq(i)

            itemp=1
            dowhile(itemp.le.numtriloops.and.triloop(itemp,1).ne.key)
               itemp=itemp+1
            enddo   
            if(triloop(itemp,1).eq.key)then
               tlink=triloop(itemp,2)
            endif

         endif

         erg4 = erg4 + hairpin(size) + eparam(4) + au_pen(i,j) + tlink

         return
      else
c
         tlink=0
         if (size.eq.4) then
            key=((((numseq(i+5)*8+numseq(i+4))*8+numseq(i+3))*8+
     +       numseq(i+2))*8+numseq(i+1))*8+numseq(i)
            tlptr=1
            do while (tloop(tlptr,1).ne.key.and.tlptr.le.numtloops)
               tlptr=tlptr+1
            enddo

            if (tloop(tlptr,1).eq.key) tlink=tloop(tlptr,2)
         endif
            erg4 = erg4 + tstkh(numseq(i),numseq(j),numseq(i+1),numseq(j-1))
     .          + hairpin(size) + eparam(4) + tlink

         return
      endif
      return
      end
 
      function erg6(i,j,ip,jp)
      include 'efn.inc'
      integer e(4)
      integer*2 tlink,tlptr,itemp

c     Dangling base stacking energy. ip dangles over the i,j base-pair.
c     3' or 5' dangle if jp = 1 or 2 respectively.
600   erg6 = dangle(numseq(i),numseq(j),numseq(ip),jp)
      return
      end

c     Used to penalize non CG GC closings of helices in multi-branch
c     and external loops

      function au_pen(i,j)
      include 'efn.inc'
      integer inc2(5,5)
      data inc2/0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0/

      au_pen = inc2(numseq(i),numseq(j))*eparam(10)
      return
      end

      subroutine process
c     Process RNA sequence.
      include 'efn.inc'
c     Selected fragment is from hstnum(1) to hstnum(n) in historical
c     numbering. ( 1 to n in internal numbering )
 
      do k = 1,n
c           A - type 1
c           B - an A accessible to nuclease cleavage
c           W - a modified A that can pair only at a helix end
c           C - type 2
c           Z - a C accessible to nuclease cleavage
c           X - a modified C that can pair only at a helix end
c           G - type 3
c           H - a G accessible to nuclease cleavage
c           Y - a modified G that can pair only at a helix end
c           U/T - type 4
c           V/W - a U/T accessible to nuclease cleavage
c           Z - a modified U/T that can pair only at a helix end
c           anything else - type 5
c           HSTNUM stores historical numbering
c           NUMSEQ stores nucleotide type.
c           AUX stores auxiliary force or prohibit information
c           This information may be used to find invalid constraints
c           or to display the folded fragment annotated with constraints.
            numseq(k) = 5
            force(k) = 0
            if (seq(k).eq.'A'.or.seq(k).eq.'a') numseq(k) = 1
            if (seq(k).eq.'B'.or.seq(k).eq.'b') then
               numseq(k) = 1
               force(k) = 3
            endif
            if (seq(k).eq.'W'.or.seq(k).eq.'w') then
               numseq(k) = 1
               force(k) = 4
            endif
            if (seq(k).eq.'C'.or.seq(k).eq.'c') numseq(k) = 2
            if (seq(k).eq.'Z'.or.seq(k).eq.'z') then
               numseq(k) = 2
               force(k) = 3
            endif
            if (seq(k).eq.'X'.or.seq(k).eq.'x') then
               numseq(k) = 2
               force(k) = 4
            endif
            if (seq(k).eq.'G'.or.seq(k).eq.'g') numseq(k) = 3
            if (seq(k).eq.'H'.or.seq(k).eq.'h') then
               numseq(k) = 3
               force(k) = 3
            endif
            if (seq(k).eq.'Y'.or.seq(k).eq.'y') then
               numseq(k) = 3
               force(k) = 4
            endif
            if (seq(k) .eq. 'U'.or.seq(k).eq.'u') numseq(k) = 4
            if (seq(k) .eq. 'V'.or.seq(k).eq.'v') then
               numseq(k) = 4
               force(k) = 3
            endif
            if (seq(k).eq.'Z'.or.seq(k).eq.'z') then
               numseq(k) = 4
               force(k) = 4
            endif
            if (seq(k) .eq. 'T'.or.seq(k).eq.'t') numseq(k) = 4
            if (seq(k) .eq. 'W'.or.seq(k).eq.'w') then
               numseq(k) = 4
               force(k) = 3
            endif


      enddo
      return
      end

c     Used in reading the energy files.
      function convt(str)
      implicit integer (a-z)
      character*6 str
      logical neg
 
      neg = .false.
      place = 0
      convt = 0
 
      do i = 6,1,-1
        if (str(i:i).eq.'-') then
          neg = .true.
        else
          if (str(i:i).ge.'0'.and.str(i:i).le.'9') then
             convt = convt + 10**place * (ichar(str(i:i)) - ichar('0'))
             place = place+1
          endif
        endif
      enddo
      if (neg) convt = -1 * convt
      return
      end
 
c     Reads energy file names and open the files for reading.
      subroutine enefiles
      character*80 filen,path
 
      call getenv('MFOLDLIB',path)
      in = index(path,' ')
      if (path.eq.'     ') then
         call getenv('MFOLD',path)
         in = index(path,' ')
         path(in:in+4) = '/dat/'
      else
         path(in:in) = '/'
      endif
 3    write(6,*) 'Enter asymmetric interior loop of size 3 energy file name'
      write(6,*) '(default asint1x2.dat)'
      read (5,100,end=1) filen
      if (filen.eq.'         ') filen = 'asint1x2.dat'
      open(8,file=filen,status='old',err=4)
      goto 5
 4    open(8,file=path(1:index(path,' ')-1)//filen,status='old',err=3)

 5    write(6,*) 'Enter asymmetric interior loop of size 5 energy file name'
      write(6,*) '(default asint2x3.dat)'
      read (5,100,end=1) filen
      if (filen.eq.'         ') filen = 'asint2x3.dat'
c      open(9,file=filen,status='old',err=6)
      goto 10
c 6    open(9,file=path(1:index(path,' ')-1)//filen,status='old',err=5)

10    write(6,*) 'Enter dangle energy file name (default dangle.dat)'
      read (5,100,end=1) filen
      if (filen.eq.'         ') filen = 'dangle.dat'
      open(10,file=filen,status='old',err=11)
      goto 20
11    open(10,file=path(1:index(path,' ')-1)//filen,status='old',err=10)
 
20    write(6,*) 'Enter loop energy file name (default loop.dat)'
      read (5,100,end=1) filen
      if (filen.eq.'         ') filen = 'loop.dat'
      open(11,file=filen,status='old',err=21)
      goto 25
21    open(11,file=path(1:index(path,' ')-1)//filen,status='old',err=20)
 
25    write(6,*) 'Enter misc. loop energy file name (default miscloop.dat)'
      read (5,100,end=1) filen
      if (filen.eq.'         ') filen = 'miscloop.dat'
      open(32,file=filen,status='old',err=26)
      goto 30
26    open(32,file=path(1:index(path,' ')-1)//filen,status='old',err=25)
 
30    write(6,*)'Enter symmetric interior loop of size 2 energy file'
      write(6,*)' name (default sint2.dat)'
      read(5,100,end=1)filen
      if(filen.eq.'         ')filen='sint2.dat'
      open(33,file=filen,status='old',err=31)
      goto 35
31    open(33,file=path(1:index(path,' ')-1)//filen,status='old',err=30)

35    write(6,*)'Enter symmetric interior loop of size 4 energy file'
      write(6,*)' name (default sint4.dat)'
      read(5,100,end=1)filen
      if(filen.eq.'         ')filen='sint4.dat'
      open(34,file=filen,status='old',err=36) 
      goto 37
36    open(34,file=path(1:index(path,' ')-1)//filen,status='old',err=35)

37    write(6,*)'Enter symmetric interior loop of size 6 energy file'
      write(6,*)' name (default sint6.dat)'
      read(5,100,end=1)filen
      if(filen.eq.'         ')filen='sint6.dat'
c      open(35,file=filen,status='old',err=38)
      goto 40
c38    open(35,file=path(1:index(path,' ')-1)//filen,status='old',err=37)
 
40    write(6,*) 'Enter stack energy file name (default stack.dat)'
      read (5,100,end=1) filen
      if (filen.eq.'         ') filen = 'stack.dat'
      open(12,file=filen,status='old',err=41)
      goto 45
41    open(12,file=path(1:index(path,' ')-1)//filen,status='old',err=40)

45    write(6,*) 'Enter tetraloop energy file name (default tloop.dat)'
      read (5,100,end=1) filen
      if (filen.eq.'         ') filen = 'tloop.dat'
      open(29,file=filen,status='old',err=46)
      goto 47
46    open(29,file=path(1:index(path,' ')-1)//filen,status='old',err=45)

47    write(6,*) 'Enter triloop energy file name (default triloop.dat)'
      read(5,100,end=1) filen
      if (filen.eq.'         ') filen = 'triloop.dat'
      open(39,file=filen,status='old',err=48)
      goto 50
48    open(39,file=path(1:index(path,' ')-1)//filen,status='old',err=47)
 
50    write(6,*)'Enter tstackh energy file name (default tstackh.dat)'
      read(5,100,end=1)filen
      if(filen.eq.'         ') filen = 'tstackh.dat'
      open(13,file=filen,status='old',err=51)
      goto 55
51    open(13,file=path(1:index(path,' ')-1)//filen,status='old',err=50)
 
55    write(6,*) 'Enter tstacki energy file name (default tstacki.dat)'
      read (5,100,end=1) filen
      if (filen.eq.'         ') filen = 'tstacki.dat'
      open(14,file=filen,status='old',err=56)
      goto 2
56    open(14,file=path(1:index(path,' ')-1)//filen,status='old',err=55)
 
100   format(a40)
      goto 2
1     stop
2     return
      end
 

c     Initialize the stack.
      subroutine initst
      implicit integer (a-z)
      dimension stk(500,4)
      common /stk/ stk,sp
 
      sp = 0
      return
      end
c     Add A,B,C,D to the bottom of the stack.
      subroutine push(a,b,c,d)
      implicit integer (a-z)
      dimension stk(500,4)
      common /stk/ stk,sp
 
      sp = sp + 1
      if (sp.gt.500) then
         write (6,*) 'STOP: STACK OVERFLOW'
         call exit(1)
      endif
      stk(sp,1) = a
      stk(sp,2) = b
      stk(sp,3) = c
      stk(sp,4) = d
      return
      end
c     Retrieve A,B,C,D from the bottom of the stack and decrease the
c     stack size by one.
      function pull(a,b,c,d)
      implicit integer (a-z)
      dimension stk(500,4)
      common /stk/ stk,sp
 
      if (sp.eq.0) then
         pull = 1
         return
      endif
      a = stk(sp,1)
      b = stk(sp,2)
      c = stk(sp,3)
      d = stk(sp,4)
      sp = sp - 1
      pull = 0
      return
      end

c     Text output of a secondary structure.
      subroutine linout(n1,n2,energy,iret,jret,error)
c
c      include 'rna.inc'
      include 'efn.inc'
      character*1 array(6,5000)
      real energy
      integer unit
c
      data amax/5000/
c
c     Write sequence label and computed energy 
c
      unit = cntrl(4)
      hstn1 = hstnum(n1)
      hstn2 = hstnum(n2)
      write(unit,103) hstn1,hstn2,ctlabel,energy
c      write(unit,103) hstn1,hstn2,seqlab,energy
 103  format('Folding bases ',i5,' to ',i5,' of ',a50,/
     .       'Computed energy  =  ',f8.1/)
c     .       ' dG = ',f8.1/)
c
c     Initialize traceback
c
      call initst
      call push(n1,n2,0,0)
c
c     Fill in output array
c
      do while (pull(i,j,countr,xx).eq.0)
c
c     Look for dangling ends
c
         ip = i
         jp = j
c        Chew off dangling 5' bases
         do while (basepr(ip).eq.0.and.ip.lt.j)
            ip = ip+1
         enddo
         if (ip.eq.jp) then
c        Hairpin loop found - dump array
            size = j - i + 1
            hsize = (size-1)/2
            if (hsize.gt.0) then
               do k = 1,hsize
                  array(1,countr+k) = ' '
                  array(2,countr+k) = seq(hstnum(i+k-1))
                  if(10*(hstnum(i+k-1)/10).eq.hstnum(i+k-1))
     .                 call digit(1,countr+k,hstnum(i+k-1),amax,array)
                  array(3,countr+k) = ' '
                  array(4,countr+k) = ' '
                  array(5,countr+k) = seq(hstnum(j-k+1))
                  array(6,countr+k) = ' '
                  if (10*(hstnum(j-k+1)/10).eq.hstnum(j-k+1))
     .                 call digit(6,countr+k,hstnum(j-k+1),amax,array)
               enddo
               countr = countr + hsize + 1
               if (2*hsize.eq.size-1) then
                  array(1,countr) = ' '
                  array(2,countr) = ' '
                  array(3,countr) = '\\'
                  array(4,countr) = seq(hstnum(i+hsize))
                  array(5,countr) = ' '
                  array(6,countr) = ' '
               else
                  array(1,countr) = ' '
                  array(2,countr) = ' '
                  array(3,countr) = seq(hstnum(i+hsize))
                  array(5,countr) = ' '
                  array(6,countr) = ' '
                  if (10*(hstnum(i+hsize)/10).eq.hstnum(i+hsize))
     .                 call digit(1,countr,hstnum(i+hsize),amax,array)
                  array(4,countr) = seq(hstnum(i+hsize+1))
                  if (10*(hstnum(i+hsize+1)/10).eq.hstnum(i+hsize+1))
     .                 call digit(6,countr,hstnum(i+hsize+1),amax,array)
               endif
            endif
            call dump_array(amax,array,countr,unit)
            do k = 1,countr
               do k1 = 1,6
                  array(k1,k) = ' '
               enddo
            enddo
         else
c           Chew off dangling 3' bases
            do while (basepr(jp).eq.0)
               jp = jp-1
            enddo
c
c           Test for bifurcation
c
            if (basepr(ip).lt.basepr(jp)) then
               do k1 = 1,6
                  do k2 = 1,2
                     array(k1,countr+k2) = ' '
                  enddo
               enddo
               if (basepr(i).gt.0) then
                  array(3,countr+1) = '-'
                  array(3,countr+2) = '-'
               else
                  array(2,countr+1) = '.'
                  array(2,countr+2) = '-'
               endif
               array(5,countr+1) = '\\'
               countr = countr + 2
               call push(basepr(ip)+1,j,countr,0) 
               call push(i,basepr(ip),countr,0)
            else
               k = max0(ip-i,j-jp)
               if (k.gt.0) then
c              Unpaired bases are in an interior or bulge loop.
                  do kk = 1,k
                     array(1,countr+kk) = ' '
                     array(3,countr+kk) = ' '
                     array(4,countr+kk) = ' '
                     array(6,countr+kk) = ' '
                     if (i+kk-1.lt.ip) then
                        array(2,countr+kk) = seq(hstnum(i+kk-1))
                        if (10*(hstnum(i+kk-1)/10).eq.hstnum(i+kk-1))
     .                       call digit(1,countr+kk,hstnum(i+kk-1),amax,array)
                     else
                        array(2,countr+kk) = '-'
                     endif
                     if (j-kk+1.gt.jp) then
                        array(5,countr+kk) = seq(hstnum(j-kk+1))
                        if (10*(hstnum(j-kk+1)/10).eq.hstnum(j-kk+1))
     .                       call digit(6,countr+kk,hstnum(j-kk+1),amax,array)
                     else
                        array(5,countr+kk) = '-'
                     endif
                  enddo
                  countr = countr + k
               endif
c      
c              Stacking must occur
c      
               i = ip
               j = jp
               do while (basepr(i).eq.j)
c              Base pair case
                  countr = countr + 1
                  array(1,countr) = ' '
                  array(2,countr) = ' '
                  array(5,countr) = ' '
                  array(6,countr) = ' '
                  array(3,countr) = seq(hstnum(i))
                  if (10*(hstnum(i)/10).eq.hstnum(i))
     .                 call digit(1,countr,hstnum(i),amax,array)
                  array(4,countr) = seq(hstnum(j))
                  if (10*(hstnum(j)/10).eq.hstnum(j))
     .                 call digit(6,countr,hstnum(j),amax,array)
                  if(i.eq.iret.and.j.eq.jret) then
                     array(2,countr) = '|'
                     array(5,countr) = '^'
                  endif
                  i = i + 1
                  j = j - 1
               enddo
               call push(i,j,countr,0)
            endif
         endif
      enddo
c     Normal exit.
      return
      end
 
      subroutine dump_array(amax,array,countr,unit)
      integer amax,countr,unit,k1,k2
      character*1 array(6,amax)
      do k1 = 1,6
         write(unit,100) (array(k1,k2),k2=1,countr)
 100     format(5000a1)
      enddo
      write(unit,101)
 101  format(/)
      return
      end

c     Puts the number base_num in row 'row' and column 'col' of the array.
c     The least significant digit ends up in column 'col'. If the number
c     is too large to fit, a period is put in column 'col'.
      subroutine digit(row,col,base_num,amax,array)
      implicit integer (a-z)
      integer amax,row,col,base_num,p,q,r
      character*1 array(6,amax)

      position = col
      p = base_num/10
      q = base_num - 10*p
      do while (q.gt.0.or.p.gt.0)
         write(array(row,position),100) q
 100     format(i1)
         if (position.eq.1.and.p.gt.0) then
c     Number overflows to the left. Put a dot in column 'col' and return.
            do k = position,col
               array(row,k) = ' '
            enddo
            array(row,col) = '.'
            return
         endif
         position = position - 1
         r = p
         p = r/10
         q = r - 10*p
      enddo
      return
      end

c     Set up output units and files for RNA folding.
      subroutine outputs
      include 'efn.inc'
      character*1 in
      character*40 str,dstr
      data dstr/'                                        '/

c     Examine CT file name to get default names for output files.
 
      k = 1
      do while ((ctnam(k:k).lt.'A'.or.ctnam(k:k).gt.'Z').and.
     .          (ctnam(k:k).lt.'a'.or.ctnam(k:k).gt.'z'))
         k = k + 1
      enddo
      slen = min0(30,25+k)
      do while (ctnam(slen:slen).eq.' ')
         slen = slen - 1
      enddo
      j = 1
      do i = k,slen
         if ((ctnam(i:i).ge.'A'.and.ctnam(i:i).le.'Z').or.
     .       (ctnam(i:i).ge.'a'.and.ctnam(i:i).le.'z').or.
     .       (ctnam(i:i).ge.'0'.and.ctnam(i:i).le.'9')) then
 
            dstr(j:j) = ctnam(i:i)
         else
            dstr(j:j) = '_'
         endif
         j = j + 1
      enddo
      slen = j - 1
 
c     Line printer output. Get name and open file for write.
      cntrl(2) = 0
      write (6,5010)
      read (5,5000,end=1) in
      if (in.ne.'N'.and.in.ne.'n') then
          cntrl(2) = 1
          write (6,5011)
          read (5,5000,end=1) in
          if (in.eq.'N'.or.in.eq.'n') then
             dstr(slen+1:slen+4) = '.out'
51           write (6,5012) dstr
             read (5,5001) str
             if (str.eq.'     ') str = dstr
             cntrl(4) = 20
             open(20,file=str,status='unknown',err=51)
          else
             cntrl(4) = 6
          endif
      endif
 
c     Thermodynamic DETails output. Get name and open file for write.
      write (6,5030)
      read (5,5000,end=1) in
      if (in.eq.'Y'.or.in.eq.'y') then
         cntrl(2) = cntrl(2) + 2
         dstr(slen+1:slen+4) = '.det'
53       write (6,5031) dstr
         read (5,5001) str
         if (str.eq.'     ') str = dstr
         open(22,file=str,status='unknown')
      endif
      write (6,*) ' '
      return
1     stop
5000  format(a1)
5001  format(a40)
5010  format(' Do you want printer output? (Y,n) ')
5011  format(' Output to terminal? (Y,n) ')
5012  format(' Enter output file name (default ',a35,')')
5030  format(' Do you want the thermodynamic details table? (y,N) ')
5031  format(' Enter details table file name (default ',a35,')')
      end
 
c     Reads energy files.
      subroutine ergread
 
      include 'efn.inc'
      logical  endfile,find
      character*96 inrec
      character*6 temp
      real a,b,c,d
 
c     TLoop INFORMATION IN
      call gettloops
 
c     TriLoop INFORMATION IN
      call gettri

c     Get misc loop info
      if(find(32,3,'-->')) then
         write(6,*) 'STOP: Premature end of miscloop file'
         call exit(1)
      endif
      read (32,*) prelog
      prelog = nint(prelog*100)
c     asymmetric internal loops: the ninio equation
      if(find(32,3,'-->')) then
         write(6,*) 'STOP: Premature end of miscloop file'
         call exit(1)
      endif
      read (32,*) a
      maxpen = nint(a*100.0)
      if(find(32,3,'-->')) then
         write(6,*) 'STOP: Premature end of miscloop file'
         call exit(1)
      endif
      read (32,*) a,b,c,d
      poppen(1) = nint(a*100.0)
      poppen(2) = nint(b*100.0)
      poppen(3) = nint(c*100.0)
      poppen(4) = nint(d*100.0)
c     Set default values of eparam.
      eparam(1) = 0
      eparam(2) = 0
      eparam(3) = 0
      eparam(4) = 0
      eparam(7) = 30
      eparam(8) = 30
c      eparam(9) = -500 Bonus energies are no longer used!
c     multibranched loops
      if(find(32,3,'-->')) then
         write(6,*) 'STOP: Premature end of miscloop file'
         call exit(1)
      endif
      read (32,*) a,b,c
      eparam(5) = nint(a*100.0)
      eparam(6) = nint(b*100.0)
      eparam(9) = nint(c*100.0)
c     efn2 multibranched loops
      if(find(32,3,'-->')) then
c     Version 2.3 rules and DNA rules do not need extra parameters
         write(6,*) 'End of miscloop file. Parameters 10 -> end set to 0.'
         do k= 10,16
            eparam(k) = 0
         enddo
      else
         read (32,*) a,b,c
c        Don't need these parameters yet in nafold
c        terminal AU penalty
         if(find(32,3,'-->')) then
            write(6,*) 'STOP: Premature end of miscloop file'
            call exit(1)
         endif
         read (32,*) a
         eparam(10) = nint(a*100.0)
c        bonus for GGG hairpin
         if(find(32,3,'-->')) then
            write(6,*) 'STOP: Premature end of miscloop file'
            call exit(1)
         endif
         read (32,*) a
         eparam(11) = nint(a*100.0)
c        c hairpin slope
         if(find(32,3,'-->')) then
            write(6,*) 'STOP: Premature end of miscloop file'
            call exit(1)
         endif
         read (32,*) a
         eparam(12) = nint(a*100.0)
c        c hairpin intercept
         if(find(32,3,'-->')) then
            write(6,*) 'STOP: Premature end of miscloop file'
            call exit(1)
         endif
         read (32,*) a
         eparam(13) = nint(a*100.0)
c        c hairpin of 3
         if(find(32,3,'-->')) then
            write(6,*) 'STOP: Premature end of miscloop file'
            call exit(1)
         endif
         read (32,*) a
         eparam(14) = nint(a*100.0)
c        Intermolecular initiation free energy
         if(find(32,3,'-->')) then
            write(6,*) 'STOP: Premature end of miscloop file'
            call exit(1)
         endif
         read (32,*) a
         eparam(15) = nint(a*100.0)
c        GAIL (Grossly Asymmetric Interior Loop) Rule (on/off <-> 1/0)
         if(find(32,3,'-->')) then
            write(6,*) 'STOP: Premature end of miscloop file'
            call exit(1)
         endif
         read (32,*) eparam(16)
      endif

c     DANGLE IN
 
      do a = 1,5
        do b = 1,5
          do c = 1,5
            do d = 1,2
              dangle(a,b,c,d) = 0
            enddo
          enddo
        enddo
      enddo
      endfile = find(10,3,'<--')
      if (.not.endfile) then
         do var4 = 1,2
            do var1 = 1,4
               if (endfile) goto 150
               read(10,100,end=150,err=91) inrec
               do var2 = 1,4
                  do var3 = 1,4
                     j = 0
                     tstart = (var2-1)*24 + (var3-1)*6 + 1
                     temp = inrec(tstart:tstart+5)
                     do i = 2,5
                        if (temp(i-1:i+1).eq.' . ') j = infinity
                     enddo
                     if (temp(1:1).eq.'.'.or.temp(6:6).eq.'.') j = infinity
                     if (j.eq.0) j = convt(temp)
                     if (j.ne.infinity) dangle(var1,var2,var3,var4) = j
                  enddo
               enddo
               endfile = find(10,3,'<--')
            enddo
         enddo
      else
         write(6,*) 'STOP: DANGLE ENERGY FILE NOT FOUND'
         call exit(1)
      endif
 
100   format(a96)
      goto 200
 
150   write(6,*) 'STOP: PREMATURE END OF DANGLE ENERGY FILE'
      call exit(1)
 
 
c    INTERNAL,BULGE AND HAIRPIN IN
 
200   endfile = find(11,5,'-----')
      i = 1
201   read(11,100,end=300) inrec
      j = -1
      do ii = 1,3
         j = j + 6
         do while (inrec(j:j).eq.' ')
            j = j + 1
         enddo
         temp = inrec(j:j+5)
         k = 0
         do jj = 2,4
            if (temp(jj-1:jj+1).eq.' . ') k = infinity
         enddo
         if (temp(1:1).eq.'.'.or.temp(6:6).eq.'.') k = infinity
         if (k.eq.0) k = convt(temp)
         if (ii.eq.1) inter(i) = k
         if (ii.eq.2) bulge(i) = k
         if (ii.eq.3) hairpin(i) = k
      enddo
      i = i + 1
      if (i.le.30) goto 201
 
c     STACK IN
 
300   do a = 1,5
        do b = 1,5
          do c = 1,5
            do d = 1,5
              stack(a,b,c,d) = infinity
            enddo
          enddo
        enddo
      enddo
      endfile = find(12,3,'<--')
      if (.not.endfile) then
         do var1 = 1,4
            do var3 = 1,4
               if (endfile) goto 350
               read(12,100,end=350,err=92) inrec
               do var2 = 1,4
                  do var4 = 1,4
                     j = 0
                     tstart = (var2-1)*24 + (var4-1)*6 + 1
                     temp = inrec(tstart:tstart+5)
                     do i = 2,5
                        if (temp(i-1:i+1).eq.' . ') j = infinity
                     enddo
                     if (temp(1:1).eq.'.'.or.temp(6:6).eq.'.') j = infinity
                     if (j.eq.0) j = convt(temp)
                     stack(var1,var2,var3,var4) = j
                  enddo
               enddo
            enddo
            endfile = find(12,3,'<--')
         enddo
      else
         write(6,*) 'STOP: STACK ENERGY FILE NOT FOUND'
         call exit(1)
      endif
      call stest(stack,'STACK ')
 
      goto 400

 350   write(6,*) 'STOP: PREMATURE END OF STACK ENERGY FILE'
       call exit(1)

400   do a = 1,5
        do b = 1,5
          do c = 1,5
            do d = 1,5
              tstkh(a,b,c,d) = 0
            enddo
          enddo
        enddo
      enddo
      endfile = find(13,3,'<--')
      if (.not.endfile) then
         do var1 = 1,4
            do var3 = 1,4
               if (endfile) goto 350
               read(13,100,end=450,err=93) inrec
               do var2 = 1,4
                  do var4 = 1,4
                     j = 0
                     tstart = (var2-1)*24 + (var4-1)*6 + 1
                     temp = inrec(tstart:tstart+5)
                     do i = 2,5
                        if (temp(i-1:i+1).eq.' . ') j = infinity
                     enddo
                     if (temp(1:1).eq.'.'.or.temp(6:6).eq.'.') j = infinity
                     if (j.eq.0) j = convt(temp)
                     tstkh(var1,var2,var3,var4) = j
                  enddo
               enddo
            enddo
            endfile = find(13,3,'<--')
         enddo
      else
         write(6,*) 'STOP: TSTACKH ENERGY FILE NOT FOUND'
         call exit(1)
      endif

c**   CALL STEST(TSTK,'TSTACK')
      goto 500

 450  write(6,*) 'STOP: PREMATURE END OF TSTACKH ENERGY FILE'
      call exit(1)
 500  do a = 1,5
        do b = 1,5
          do c = 1,5
            do d = 1,5
              tstki(a,b,c,d) = 0
            enddo
          enddo
        enddo
      enddo

      endfile = find(14,3,'<--')
      if (.not.endfile) then
         do var1 = 1,4
            do var3 = 1,4
               if (endfile) goto 450
               read(14,100,end=550,err=94) inrec
               do var2 = 1,4
                  do var4 = 1,4
                     j = 0
                     tstart = (var2-1)*24 + (var4-1)*6 + 1
                     temp = inrec(tstart:tstart+5)
                     do i = 2,5
                        if (temp(i-1:i+1).eq.' . ') j = infinity
                     enddo
                     if (temp(1:1).eq.'.'.or.temp(6:6).eq.'.') j = infinity
                     if (j.eq.0) j = convt(temp)
                     tstki(var1,var2,var3,var4) = j
                  enddo
               enddo
            enddo
            endfile = find(14,3,'<--')
         enddo
      else
         write(6,*) 'STOP: TSTACKI ENERGY FILE NOT FOUND'
         call exit(1)
      endif

c**   call stest(tstki,'TSTACKI')
      goto 600

550   write(6,*) 'STOP: PREMATURE END OF TSTACK ENERGY FILE'
      call exit(1)

c  Read in symmetric interior loop energies
600   call symint

c  Read in asymmetric interior loop energies
      call asmint

      call symtest

      close(8)
      close(9)
      close(10)
      close(11)
      close(12)
      close(13)
      close(14)
      close(33)
      close(34)
      close(35)
      close(39)

      return
 91   write(6,*) 'STOP: ERROR reading dangle energy file'
      call exit(1)

 92   write(6,*) 'STOP: ERROR reading stacking energy file'
      call exit(1)

 93   write(6,*) 'STOP: ERROR reading tstackh energy file'
      call exit(1)

 94   write(6,*) 'STOP: ERROR reading tstacki energy file'
      call exit(1)
      end

c     Symmetry test on stacking and terminal stacking energies.
c     For all i,j,k,l between 1 and 4, stack(i,j,k,l) MUST equal
c     stack(l,k,j,i). If this fails at some i,j,k,l; these numbers
c     are printed out and the programs grinds to an abrupt halt!
      subroutine stest(stack,sname)
      integer stack(5,5,5,5),a,b,c,d
      character*6 sname
 
      do a = 1,4
        do b = 1,4
          do c = 1,4
            do d = 1,4
              if (stack(a,b,c,d).ne.stack(d,c,b,a)) then
                 write(6,*) 'SYMMETRY ERROR'
                 write(6,101) sname,a,b,c,d,stack(a,b,c,d)
                 write(6,101) sname,d,c,b,a,stack(d,c,b,a)
                 stop
              endif
            enddo
          enddo
        enddo
      enddo
      return
101   format(5x,a6,'(',3(i1,','),i1,') = ',i10)
      end
 
c     Used in reading the energy files.
c     Locates markers in the energy files so that data can be read
c     properly.
      function find(unit,len,str)
      implicit integer (a-z)
      logical find,flag
      character*20 str
      character*96 inrec
 
      find = .false.
      flag = .false.
      do  while(.not.flag)
         read(unit,100,end=200) inrec
         count = 1
         do 101 i = 1,80-len+1
           if (inrec(i:i).eq.str(count:count)) then
             count = count + 1
             if (count.gt.len) flag = .true.
             if (inrec(i+1:i+1).ne.str(count:count)) count = 1
           endif
101      continue
      enddo
 
      return
100   format(a96)
200   find = .true.
      return
      end
 
      subroutine gettloops
c------------Kevin added below-----------------------------------------
      include 'efn.inc'
      integer flag,i,j,numseqn(6)
      real energy
      character row*80,seqn*6

 2020 format(a80)
 2021 format(1x,a6,1x,f6.2)

      flag=0
      do while(flag.eq.0)
         read(29,2020,err=91)row
         if(index(row,'---').gt.0)then
            flag=1
         endif
      enddo    

      flag=0
      j=0
      do while(flag.eq.0)
         j=j+1
         read(29,2021,end=99,err=91)seqn,energy

         if(index(seqn,'    ').gt.0.or.j.eq.maxtloops)then
            flag=1
            numtloops=j-1
         else
            do i=1,6
               call letr2num(seqn(i:i),numseqn(i))
            enddo

            tloop(j,1)=((((numseqn(6)*8+numseqn(5))*8+numseqn(4))*8+
     .       numseqn(3))*8+numseqn(2))*8+numseqn(1)
            tloop(j,2)=nint(100.0*energy)
         endif
      enddo   

      close(29,status='KEEP')

      return

 91   write(6,*) 'STOP: ERROR reading tetraloop file'
      call exit(1)

 99   close(29,status='KEEP')
      numtloops=j-1

      return

      end

c---------------Kevin added above------------------------------------

c   Reads sequence and energy for triloop.d__ and stores info in array triloop
      subroutine gettri

      include 'efn.inc'
      integer i,j,flag,numseqn(5)
      real energy
      character row*80,seqn*5

 3030 format(a80)
 3031 format(1x,a5,2x,f6.2)

      flag=0
      do while(flag.eq.0)
         read(39,3030,err=91)row
         if(index(row,'------').gt.0)then
            flag=1
         endif
      enddo
      flag=0
      i=0
      do while(flag.eq.0)
         i=i+1

         if(i.eq.maxtriloops)then
            numtriloops=i
            flag=1
         endif

         read(39,3031,end=98,err=91)seqn,energy

         if(index(seqn,'   ').gt.0)then
            numtriloops=i-1
            flag=1
         else
            do j=1,5
               call letr2num(seqn(j:j),numseqn(j))
            enddo
 
            triloop(i,1)=(((numseqn(5)*8+numseqn(4))*8+numseqn(3))*8+
     .       numseqn(2))*8+numseqn(1)
            triloop(i,2)=nint(100.0*energy)
         endif
      enddo

      return

 98   continue
      numtriloops=i-1

      return

 91   write(6,*) 'STOP: ERROR reading triloop file'
      call exit(1)
      end
c----------------------------------------------------------------
         subroutine letr2num(letr,num)
         
         integer num
         character*1 letr

         if(letr.eq.'A')then
            num=1
         elseif(letr.eq.'C')then
            num=2
         elseif(letr.eq.'G')then
            num=3
         elseif(letr.eq.'U'.or.letr.eq.'T')then
            num=4
         endif

      return
1     format(/)
2     format (a)
5     format (1x,'Too many characters in numeric field of this line of',/,
     1        1x,'tloop.dat file: ',a)
      end
c---------------------------------------------------------------------- 
      subroutine symint

      include 'efn.inc'
      integer test1,test2,i,j,test3,test4,i1,j1,i2,j2,x1,y1,ip,jp,
     +i3,j3,i4,test5,flag
      real rsint2(6,6,5,5),rsint4(6,6,5,5,5,5),
     +rsint6(6,6,25,5,5,5,5)
      character rowcha*144,row2*96

c  Explanation of sint6(a,b,c,e,f,g,h); variable changes as corresponding
c  pair changes.
c  a  c   e   g  b
c     U   W   Y
c  A             C    Symmetric interior loop of size 6
c  B             D
c     V   X   Z
c         f   h 

c  Reads in energies for symmetric interior loop of size 2

      i1=0
      flag=0

 3030 format(a144)
 3031 format(24f6.2)
                     
      do while(flag.lt.4)

         read(33,3030,end=93,err=91)rowcha

         test1=index(rowcha(25:144),'<--')
         test2=index(rowcha(1:24),'                        ')

         if(test2.gt.0)then
            flag=flag+1
         elseif(test1.gt.0)then
            i1=i1+1

            do i2=1,4
               read(33,3031,err=91)((rsint2(i1,j1,i2,j2),j2=1,4),j1=1,6)
            enddo
            flag=0
         else
            flag=0
         endif
      enddo

 93   continue

      do i=1,6
         do j=1,6
            worst = 0
            do ip=1,4
               do jp=1,4
                  sint2(i,j,ip,jp) = nint(100*rsint2(i,j,ip,jp))
                  worst = max0(worst,nint(100*rsint2(i,j,ip,jp)))
               enddo
            enddo
            do ip = 1,5
               sint2(i,j,ip,5) = worst
               sint2(i,j,5,ip) = worst
            enddo
         enddo
      enddo

c   Reads in energies for symmetric interior loop of size 4

 3032 format(a96)
 3033 format(16f6.2)

      test3=0

      do while(test3.eq.0)
         read(34,3032) row2
         test3=index(row2,'<-----')
      enddo

      do x1=1,6
         do y1=1,6

            test4=0

            do while(test4.eq.0)
               read(34,3032,err=92) row2
               test4=index(row2,'<------')
            enddo

            do i2=1,4
               do j2=1,4
                  read(34,3032,err=92) row2
                  read(row2,3033,err=92)((rsint4(x1,y1,i2,j2,i3,j3),j3=1,4),
     .            i3=1,4)
               enddo
            enddo
         enddo
      enddo

      do i=1,6
         do j=1,6
            worst = 0
            do ip=1,4
               do jp=1,4
                  do i1=1,4
                     do j1=1,4
                        sint4(i,j,ip,jp,i1,j1) = nint(100*rsint4(i,j,ip,jp,i1,j1))
                        worst = max0(worst,nint(100*rsint4(i,j,ip,jp,i1,j1)))
                     enddo
                  enddo
               enddo
            enddo
c            do ip = 1,5
c               do jp = 1,5
c                  do i1 = 1,5
c                     sint4(i,j,ip,jp,i1,5) = worst
c                     sint4(i,j,ip,jp,5,i1) = worst
c                     sint4(i,j,ip,5,jp,i1) = worst
c                     sint4(i,j,5,ip,jp,i1) = worst
c                  enddo
c               enddo
c            enddo
         enddo
      enddo

c   Reads in energies for symmetric interior loop size 6
c
c   KEY
c      c  e  g  
c   a           b    Capital letters:  RNA/DNA bases
c      U--W--Y       Lower case letters:  variables assigned to 
c   A_/       \_C                         specific bases
c   B \       / D
c      V--X--Z       sint6(a,b,c,e,f,g,h)
c      
c         f  h

c3034 format(24F6.2)

c     test1=0

c     do while(test1.eq.0)
c        read(35,3030,err=94)rowcha
c        test1=index(rowcha(41:120),'-->')
c     enddo
c     
c     do i=1,6
c        do i1=1,19
c            if((i1.ne.5.and.i1.ne.10).and.i1.ne.15)then
c               do i3=1,4
c                 do j3=1,4
c                     test5=0
c                     do while(test5.eq.0)
c                       read(35,3030,err=94)rowcha
c                       test5=index(rowcha(41:120),'<--')
c                    enddo
c              
c                    do i2=1,4
c                       read(35,3034,err=94)((rsint6(i,j,i1,i2,j2,
c    .                  i3,j3),j2=1,4),j=1,6)
c                    enddo
c                 enddo                     
c              enddo
c           endif
c        enddo
c     enddo
c      do i=1,6
c        do j=1,6
c           do i1=1,25
c              do i2=1,5
c                 do j2=1,5
c                    do i3=1,5
c                       do j3=1,5
c                          if((nint((real(i1)-.5)/5.)+1).eq.5)then
c                             sint6(i,j,i1,i2,j2,i3,j3)=infinity
c                          elseif((i1-nint((real(i1)-.5)/5.)*5).eq.5)then
c                             sint6(i,j,i1,i2,j2,i3,j3)=infinity
c                          elseif(i2.eq.5.or.j2.eq.5)then
c                             sint6(i,j,i1,i2,j2,i3,j3)=infinity
c                          elseif(i3.eq.5.or.j3.eq.5)then
c                             sint6(i,j,i1,i2,j2,i3,j3)=infinity
c                          else
c                             sint6(i,j,i1,i2,j2,i3,j3)=
c    .                        nint(rsint6(i,j,i1,i2,j2,i3,j3)*100.0)
c                          endif
c                       enddo
c                    enddo
c                 enddo
c              enddo
c           enddo
c        enddo
c     enddo
c  
      return

 91   write(6,*) 'STOP: ERROR reading symmetric interior of size 2 file'
      call exit(1)

 92   write(6,*) 'STOP: ERROR reading symmetric interior of size 4 file'
      call exit(1)

c94   write(6,*) 'STOP: ERROR reading symmetric interior of size 6 file'
      call exit(1)
      end


c-----------------------------------------------------------------
      subroutine asmint

      include 'efn.inc'
      integer a,b,u,w,x,y,z,test
      character row*120
      real raint3(6,6,4,4,4),raint5(6,6,4,4,4,4,4)

 11   format(a144)
 12   format(24f6.2)

c   Reads in asymmetric interior loop of size 3 energies

      do a=1,6
         do b=1,6
            do x=1,5
               do y=1,5
                  do z=1,5
                     asint3(a,b,x,y,z)=infinity
                     
                     do u=1,5
                        do w=1,5
                           asint5(a,b,x,w,y,u,z)=infinity
                        enddo
                     enddo

                  enddo
               enddo
            enddo
         enddo
      enddo

      do a=1,6
         do z=1,4
            test=0
            do while(test.eq.0)
               read(8,11,err=91)row
               test=index(row(41:120),'<--')
            enddo

            do x=1,4
               read(8,12,err=91)((raint3(a,b,x,y,z),y=1,4),b=1,6)

               do b=1,6
                  do y=1,4
                     asint3(a,b,x,y,z)=nint(100.0*raint3(a,b,x,y,z))
                  enddo
               enddo

            enddo
         enddo
      enddo

c   Reads in asymmetric interior loop of size 5 energies
c----------------------------------------------------------------
c   Diagram and explanation of asint5 dimensions:
c
c     a  x w    b
c        A-C--          asint5(a,b,x,w,y,u,z)
c  ===U_/     \_A===
c  ===G \     / U===
c        G-C-A
c        y u z           *  raint5 has same dimensions
c-----------------------------------------------------------------
cTEMP do a=1,6
cTEMP    do x=1,4
cTEMP       do y=1,4
cTEMP          do u=1,4
cTEMP             test=0
cTEMP             do while(test.eq.0)
cTEMP                read(9,11,err=92)row
cTEMP                test=index(row(41:120),'<--')
cTEMP             enddo
cTEMP             
cTEMP             do w=1,4
cTEMP                read(9,12,err=92)((raint5(a,b,x,w,y,u,z),z=1,4),b=1,6)
cTEMP                 do b=1,6
cTEMP                   do z=1,4
cTEMP                      asint5(a,b,x,w,y,u,z)=nint(100.0*
cTEMP.                     raint5(a,b,x,w,y,u,z))
cTEMP                   enddo
cTEMP                enddo
cTEMP              enddo
cTEMP          enddo
cTEMP       enddo
cTEMP    enddo
cTEMP  enddo

      return

 91   write(6,*) 'STOP: ERROR reading asymmetric interior 1 x 2 file'
      call exit(1)

 92   write(6,*) 'STOP: ERROR reading asymmetric interior 2 x 3 file'
      call exit(1)
      end

c------------------------------------------------------------------
      subroutine symtest
      
      include 'efn.inc'
      integer diff,i,j,ii,jj,ip,jp,i1,i1p,j1p,i2,j2,i2p

c   Symmetry test for symmetric interior loops of size 2
      do i=1,6
         do j=1,6
            if(i.gt.4.and.i.le.6)then
               ip=11-i
            else
               ip=5-i
            endif
            if(j.gt.4.and.j.le.6)then
               jp=11-j
            else
               jp=5-j
            endif

            do ii=1,4
               do jj=1,4
                  diff=sint2(i,j,ii,jj)-sint2(jp,ip,jj,ii)
                  if(diff.ne.0)then
                     write(6,*) 'STOP: sint2 FAILED SYMMETRY TEST'
                     write(6,93) i,j,ii,jj,sint2(i,j,ii,jj)
 93                  format('sint2(',4i2,')= ',i6,' BUT ')
                     write(6,94) jp,ip,jj,ii,sint2(jp,ip,jj,ii)
 94                  format('sint2(',4i2,')= ',i6,'!')
                     call exit(1)
                  endif
               enddo
            enddo
         enddo
      enddo
      
c   Symmetry test for symmetric interior loops of size 4
      do i=1,6
         do j=1,6
            if(i.gt.4.and.i.le.6)then
               ip=11-i
            else
               ip=5-i
            endif
            if(j.gt.4.and.j.le.6)then
               jp=11-j
            else
               jp=5-j
            endif

            do ii=1,4
               do jj=1,4
                  do i2=1,4
                     do j2=1,4
                        diff=sint4(i,j,ii,jj,i2,j2)-
     .                  sint4(jp,ip,j2,i2,jj,ii)
                        if(diff.ne.0)then
                           write(6,*) 'STOP: sint4 FAILED SYMMETRY TEST'
                           write(6,*) i,j,ii,jj,i2,j2,sint4(i,j,ii,jj,i2,j2)
                           write(6,*) jp,ip,j2,i2,jj,ii,sint4(jp,ip,j2,i2,jj,ii)
                           call exit(1)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

c   Symmetry test for symmetric interior loops of size 6
      do i=1,6
         do j=1,6
            if(i.gt.4.and.i.le.6)then
               ip=11-i
            else
               ip=5-i
            endif
            if(j.gt.4.and.j.le.6)then
               jp=11-j
            else
               jp=5-j
            endif
            
            do i1=1,19
               if((i1.ne.5.and.i1.ne.10).and.i1.ne.15)then
                  i1p=nint((real(i1)-.5)/5.)+1
                  j1p=i1-nint((real(i1)-.5)/5.)*5
                  do i2=1,4
                     do j2=1,4
                        i2p=(j2-1)*5+i2
                        do ii=1,4
                           do jj=1,4
                              diff=sint6(i,j,i1,ii,jj,i2,j2)-
     .                        sint6(jp,ip,i2p,jj,ii,j1p,i1p)
                              if(diff.ne.0)then
                                 write(6,*) 'STOP: sint6 FAILED SYMMETRY TEST'
                                 call exit(1)
                              endif
                           enddo
                        enddo
                     enddo
                  enddo
               endif
            enddo
         enddo
      enddo

      return
      end

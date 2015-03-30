c     Energy funtions.
c     ergk (k=2,3,4,5) is the energy of a loop closed by i,j (new numbering).
c     ip,jp is the other closing base-pair when MODE = 2 or 3.
c     The ends of the sequence cannot be contained in a hairpin, bulge
c     or interior loop. By convention, the ends of the sequence are
c     put into a special kind of multi-loop. This can be called an
c     exterior loop or an open multi-loop.
c
      function erg2(i,j)
      include 'rna.inc'
      integer e(4)
      integer*2 tlink,tlptr,itemp
      logical fce
 
      erg2 = 0
c     Do not allow prohibited bases to pair.
      if (force(i).eq.1.or.force(j).eq.1) then
         erg2 = infinity
         return
      endif
  
c     Nucleotide accessibility option.
      if (force(i).eq.3.or.force(j-1).eq.3) then
         erg2 = infinity
         return
      endif
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
      include 'rna.inc'
      dimension e(4)
      integer*2 tlink,tlptr,itemp
      logical fce

      erg3 = 0
c     Do not allow prohibited bases to pair.
      if (force(i).eq.1.or.force(j).eq.1) then
         erg3 = infinity
         return
      endif
      if ((i.le.break.and.ip.gt.break).or.(jp.le.break.and.j.gt.break)) then
c          Loop is not allowed to contain the ends of the sequence.
           erg3 = infinity
           return
      endif
c
      size1 = ip - i - 1
      size2 = j - jp - 1
c     Check for single stranded bases that MUST be paired
      if ((size1.gt.0.and.fce(i,ip)).or.(size2.gt.0.and.fce(jp,j))) then
         erg3 = infinity
         return
      endif
      if (size1.eq.0.or.size2.eq.0) then
c         Check for nucleotide accessibility.
          if (size1.eq.0.and.force(i).eq.3) then
             erg3 = infinity
             return
          endif
          if (size2.eq.0.and.force(jp).eq.3) then
             erg3 = infinity
             return
          endif
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
                else
                   erg3 = infinity
                   return
                endif

                if((numseq(ip)+numseq(jp)).eq.5)then
                   b=numseq(ip)
                elseif(numseq(ip).eq.3.and.numseq(jp).eq.4)then
                   b=5
                elseif(numseq(ip).eq.4.and.numseq(jp).eq.3)then
                   b=6
                else
                   erg3 = infinity
                   return
                endif

c       Size = 3
                if(size.eq.3) then
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
      include 'rna.inc'
      integer e(4)
      integer*2 tlink,tlptr,itemp
      logical fce

      erg4 = 0
c     Do not allow prohibited bases to pair.
      if (force(i).eq.1.or.force(j).eq.1) then
         erg4 = infinity
         return
      endif
 
      if (i.le.break.and.j.gt.break) then
c         Hairpin loop must not contain the ends of the sequence.
          erg4 = infinity
          return
      endif
c     Hairpin loop must not contain bases that are forced to be double
c     stranded
      if (fce(i,j)) then
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

       if (strand(i).eq.1.and.strand(j).eq.2) then
c      Treat as a dimer - done for John SantaLucia (and hybridization server)
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
         if (lorc.eq.'c') then
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
 
      function erg5(i,j)
      include 'rna.inc'
      integer e(4),loop
      integer*2 tlink,tlptr,itemp
      logical fce
      data loop/3/

      erg5 = 0
c     Do not allow prohibited bases to pair.
      if (force(i).eq.1.or.force(j).eq.1) then
         erg5 = infinity
         return
      endif

c     Multi-branch  (or multi-) loop closed by i,j.
      do ii = 1,4
         e(ii) = infinity
      enddo
c     Open bifurcation is possible when j > n
      if (j.gt.break) then
         e(1) = w3(i+1) + w5(j-n-1)
         if (i.ne.break) e(2) = dangle(numseq(i),numseq(j),numseq(i+1),1) 
     .                      + w3(i+2) + w5(j-n-1)
         if (j.ne.break+1) e(3) = dangle(numseq(i),numseq(j),numseq(j-1),2)
     .                        + w3(i+1) + w5(j-n-2)
         if (i.ne.break.and.j.ne.break+1) then
            e(4) = dangle(numseq(i),numseq(j),numseq(i+1),1) + w3(i+2) +
     .             dangle(numseq(i),numseq(j),numseq(j-1),2) + w5(j-n-2)
         endif
      endif
#     write(6,*) 'FIRST',i,j,e
      if (j-i.gt.2*loop+4) then
c        Cases:
c        No dangling ends next to the i,j base-pair.
         if (i.lt.break) e(1) = 
     .            min0(e(1),wmb(i+1,mod(j-1,3)) + eparam(5) + eparam(9))
c        i+1 dangles on the i,j base-pair.
         if (i+1.lt.break) e(2) = 
     .            min0(e(2),dangle(numseq(i),numseq(j),numseq(i+1),1)
     .          + wmb(i+2,mod(j-1,3)) + eparam(5) + eparam(6) + eparam(9))
c        j-1 dangles on the i,j base-pair.
         if (i.lt.break.and.(j.lt.break.or.j.gt.break+1)) e(3) =
     .          min0(e(3),dangle(numseq(i),numseq(j),numseq(j-1),2)
     .          + wmb(i+1,mod(j-2,3)) + eparam(5) + eparam(6) + eparam(9))
c        Both i+1 and j-1 dangle on the i,j base-pair.
         if (i+1.lt.break.and.(j.lt.break.or.j.gt.break+1)) e(4) = 
     .          min0(e(4),dangle(numseq(i),numseq(j),numseq(i+1),1)
     .          + dangle(numseq(i),numseq(j),numseq(j-1),2) + 
     .          wmb(i+2,mod(j-2,3)) + eparam(5) + 2*eparam(6) + eparam(9))
      endif
#     write(6,*) 'SECOND',i,j,e
c     Check for bases that must pair.
      if (i.ne.break.and.force(i+1).eq.2) then
         e(2) = infinity
         e(4) = infinity
      endif
      if (j.ne.break+1.and.force(j-1).eq.2) then
         e(3) = infinity
         e(4) = infinity
      endif
      erg5 = erg5 + au_pen(i,j) + min0(e(1),e(2),e(3),e(4))
#     write(6,*) 'THIRD',i,j,erg5
      return
      end

      function erg6(i,j,ip,jp)
      include 'rna.inc'
      integer e(4)
      integer*2 tlink,tlptr,itemp
      logical fce

      erg6 = 0
c     Do not allow prohibited bases to pair.
      if (force(i).eq.1.or.force(j).eq.1) then
         erg6 = infinity
         return
      endif

c     Dangling base stacking energy. ip dangles over the i,j base-pair.
c     3' or 5' dangle if jp = 1 or 2 respectively.
600   erg6 = erg6 + dangle(numseq(i),numseq(j),numseq(ip),jp)
      if (force(ip).eq.2) erg6 = infinity
      return
      end
 
      subroutine fill(nofold)
c     This subroutine computes the arrays of optimal energies.
      include 'rna.inc'
      integer inc(5,5),e(5),vtemp
      data loop/3/,inc/0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0
     ./
 
c     Sweep along anti-diagonals and mark non base pairs and 
c     isolated base pairs. 
c     Also mark ANY base pair involving a modified base that COULD
c     be an interior base in a helix! This is overkill, but the 
c     penalty/bonus method does not appear to work.

      do l = 2,2*n-2
         if (l.le.n) then
            i = 1
            j = l
         else
            i = l + 1 - n
            j = n
         endif
c        Stop at minimum hairpin loop size
         test1 = 0
         test2 = inc(numseq(i),numseq(j))
         if (j-i.gt.loop+2) then
            test3 = inc(numseq(i+1),numseq(j-1))
         else
            test3 = 0
         endif
         do while (j-i.gt.loop)
c           test1: i-1.j+1, if it exists
c           test2: Current base pair, i.j
c           test3: i+1.j-1, if it can exist
            if(test2.eq.0.or.(test2.eq.1.and.(test1+test3.eq.0))) then
cnofilter   if(test2.eq.0) then
               vst((n-1)*(i-1)+j) = infinity
               vst((n-1)*(j-1)+i+n) = infinity
            endif
            if(test2.eq.1.and.(force(i).eq.4.or.force(j).eq.4).
     .           and.test1+test3.eq.2) then
               vst((n-1)*(i-1)+j) = infinity
               vst((n-1)*(j-1)+i+n) = infinity
            endif
            i = i + 1
            j = j - 1
            test1 = test2
            test2 = test3
            if (j-i.gt.loop+2) then
               test3 = inc(numseq(i+1),numseq(j-1))
            else
               test3 = 0
            endif
         enddo
      enddo

      vmin = infinity
      do j = 1,2*n-1

c        How far along is the computation?
         if(j.eq.20*(j/20)) then
            open(19,file='mfold.log',status='unknown')
            write(19,*) j/2,','
            close(19)
            if (usage.eq.'text') then
               call system('cat mfold.log | tr -cd ",[0-9]" > /dev/tty')
            endif
         endif

         do i = 1,2*n
            wmb(i,mod(j,3)) = infinity
         enddo
         do i = min0(j,n),max0(1,j-n+1),-1
           vij =  infinity
           wij = infinity
           wmbij = infinity
           if (j.le.break) then
              if (j-i.le.loop) goto 300
           else
              if (i.eq.break.or.j.eq.break+1) goto 100
           endif

c          Use maxbp constraint, if applicable. 
           if (lorc.eq.'l') then
              if (j.le.n) then
                 if (j-i.gt.maxbp) goto 200
              else
                 if (n+i-j.gt.maxbp) goto 200
              endif
           else
              if ((j-i.gt.maxbp).and.(n+i-j.gt.maxbp)) goto 200
           endif

#          write(6,2091) 'VST&INC',i,j,vst((n-1)*(i-1)+j),
#    .  inc(numseq(i),numseq(j)),'SEQIJ',numseq(i),numseq(j)
#2091      format(/,a7,1x,4i6,1x,a5,2i2)

c          Test for a prohibited base-pair or a pair which cannot form.
           if (vst((n-1)*(i-1)+j).gt.0.or.inc(numseq(i),numseq(j)).eq.0) goto 200
c          Compute vij, the minimum energy of the fragment from i to j
c          inclusive where i and j base-pair with one another.
c          Perhaps i,j closes a hairpin loop.
           vij = min0(vij,erg4(i,j))
#          write(6,*) 'HLOOP',i,j,erg4(i,j),'VIJ=',vij
           if (j-i-1.ge.loop+2.or.j.gt.break) then
c             Perhaps i,j stacks over i+1,j-1.
              vij = min0(vij,erg2(i,j)+v(i+1,j-1))
           endif
c          Search for the best bulge or interior loop closed by i,j.
           if (j-i-1.ge.loop+3.or.j.gt.break) then
              do d = j-i-3,1,-1
                 do ip = i+1,j-1-d
                    jp = d+ip
                    if (j-i-2-d.gt.eparam(7)) goto 100
                    if (abs(ip-i+jp-j).le.eparam(8).and.
     .                  inc(numseq(ip),numseq(jp)).eq.1) then
                       if (ip.gt.n) then
                          vij = min0(vij,erg3(i,j,ip,jp)+vst((n-1)*(ip-n-1)+jp-n))
                       else
                          vij = min0(vij,erg3(i,j,ip,jp)+vst((n-1)*(ip-1)+jp))
                       endif
                    endif
                 enddo
              enddo
           endif
#          write(6,*) 'VBI VIJ=',vij
c          Search for the best multi-loop closed by i,j.

100        if (j-i-1.ge.2*loop+4.or.j.gt.break) then
              if (vst((n-1)*(i-1)+j).gt.0.or.inc(numseq(i),numseq(j)).eq.0) goto 200
              vij=min0(vij,erg5(i,j))
#             write(6,*) 'MLOOP',i,j,erg5(i,j)
           endif

c          Compute wij.
c          A multi-loop containing n and 1 (ie. m+1) as single-stranded
c          bases is called an exterior loop. wij is the minimum folding
c          energy of a non-empty folding on i to j inclusive where an
c          exterior loop is given an energy of eparam(5) plus eparam(6)
c          per single-stranded exterior base plus eparam(9) per
c          double-stranded exterior base-pair in addition to
c          possible dangling base energies. 
c          Version 3.0 adds a penalty of eparam(10) to each exterior
c          base pair that is AU, UA, GU or UG. This is also included
c          in w5 and w3 below, whereas eparam(9) is not.

200        do ii = 1,5
             e(ii) = infinity
           enddo
           if (i.ne.break) then
c            Add single-stranded i to an optimal structure containing
c            the base-pair i+1,j.
             e(1) = v(i+1,j) + eparam(9) + eparam(6) + erg6(j,i+1,i,2) + au_pen(i+1,j)
             e(4) = w(i+1,j) + eparam(6)
c            Check if i must base pair
             if (force(i).eq.2) then
                e(1) = infinity
                e(4) = infinity
             endif
           endif
           if (j.ne.break+1) then
c            Add single-stranded j to an optimal structure containing
c            the base-pair i,j-1.
             e(2) = v(i,j-1) + eparam(9) + eparam(6) + erg6(j-1,i,j,1) + au_pen(i,j-1)
             e(5) = w(i,j-1) + eparam(6)
c            Check if j must base pair
             if (force(j).eq.2) then
                e(2) = infinity
                e(5) = infinity
             endif
           endif
           if (i.ne.break.and.j.ne.break+1) then
c            Add single-stranded i and j to an optimal structure containing
c            the base-pair i+1,j-1.
             e(3) = v(i+1,j-1) + eparam(9) + 2*eparam(6) + au_pen(i+1,j-1) +
     .                 erg6(j-1,i+1,i,2) + erg6(j-1,i+1,j,1)
c            Check if i or j must base pair
             if (force(i).eq.2.or.force(j).eq.2) then
                e(3) = infinity
             endif
           endif
 
           wij = min0(eparam(9)+au_pen(i,j)+vij,e(1),e(2),e(3),e(4),e(5))
 
           if (j-i-1.gt.2*loop+2.or.j.gt.break) then
              index = (n-1)*(i-1)
c             Search for an open bifurcation.
              do k = i,j-1
                 if (k.ne.break) then
                    wmbij = min0(wmbij,wst(index+k)+work(k+1,mod(j,3)))
                 endif
              enddo
           endif
           wij = min0(wij,wmbij)
 
c          Store vij and wij. They can be regarded as elements
c          v(i,j) and w(i,j) in two dimensional arrays. They
c          are actually stored in one dimensional arrays VST and WST
c          in position (n-1)*(i-1) + j.
c          Columns j,j-1 and j-2 of w are also stored in the work array,
c          work. This is done to reduce virtual memory swaps.
c          Columns j and j-1 of wmb (best closed multi-branched loop
c          energies) and stored in the temporary array, wmb.

300        vst((n-1)*(i-1)+j) = vij
           wst((n-1)*(i-1)+j) = wij
           work(i,mod(j,3)) = wij
           wmb(i,mod(j,3)) = wmbij
           if (i.eq.1) w5(j) = infinity
           if (j.eq.n) w3(i) = infinity
#          write(6,*) 'END VIJ=',vij
           if (j.gt.n) then
c             vmin is the minimum folding energy of the entire sequence.
              vtemp = vst((n-1)*(i-1)+j) + vst((n-1)*(j-n-1)+i)
              vmin = min0(vmin,vtemp)
           endif
         enddo

c        Compute the best folding energies for the fragments 1-->i
c        (stored in w5(i)) and i-->n (stored in w3(i)).
         if (j.eq.break) then
            w5(-1) = 0
            w5(0) = 0
            do i = 1,loop+1
               w5(i) = w5(i-1)
               if(force(i).eq.2) w5(i) = infinity
            enddo
            do i = loop+2,n
               w5(i) = w5(i-1)
               do k = 1,4
                  e(k) = infinity
               enddo
               do k = 0,i-4
                  e(1) = min0(e(1),w5(k)+vst((n-1)*k+i)+au_pen(k+1,i))
                  e(2) = min0(e(2),w5(k)+erg6(i,k+2,k+1,2)+vst((n-1)*(k+1)+i)+au_pen(k+2,i))
                  e(3) = min0(e(3),w5(k)+erg6(i-1,k+1,i,1)+vst((n-1)*k+i-1)+au_pen(k+1,i-1))
                  e(4) = min0(e(4),w5(k)+erg6(i-1,k+2,k+1,2)+au_pen(k+2,i-1)+
     .                        erg6(i-1,k+2,i,1)+vst((n-1)*(k+1)+i-1))
               enddo
c              Check if i must base pair.
               if (force(i).eq.2) then
                  w5(i) = infinity
                  e(3) = infinity
                  e(4) = infinity
               endif
               w5(i) = min0(w5(i),e(1),e(2),e(3),e(4))
            enddo
            w3(n+1) = 0
            w3(n+2) = 0
            do i = n,n-loop,-1
               w3(i) = w3(i+1)
               if (i.le.n.and.force(i).eq.2) w3(i) = infinity
            enddo
            do i = n-loop-1,1,-1
               w3(i) = w3(i+1)
               do k = 1,4
                  e(k) = infinity
               enddo
               do k = n+1,i+4,-1
                  e(1) = min0(e(1),vst((n-1)*(i-1)+k-1)+w3(k)+au_pen(i,k-1))
                  e(2) = min0(e(2),vst((n-1)*i+k-1)+erg6(k-1,i+1,i,2)+w3(k)+au_pen(i+1,k-1))
                  e(3) = min0(e(3),vst((n-1)*(i-1)+k-2)+erg6(k-2,i,k-1,1)+
     .                        au_pen(i,k-2)+w3(k))
                  e(4) = min0(e(4),vst((n-1)*i+k-2)+erg6(k-2,i+1,i,2) +
     .                        erg6(k-2,i+1,k-1,1)+au_pen(i+1,k-2)+w3(k))
               enddo
c              Check if i must base pair.
               if (force(i).eq.2) then
                  w3(i) = infinity
                  e(2) = infinity
                  e(4) = infinity
               endif
               w3(i) = min0(w3(i),e(1),e(2),e(3),e(4))
            enddo
         endif

c        Fill in some work array values before beginning on
c        the next column.
         if (j.ge.n) then
            do k = j+1,n+1,-1
               work(k,mod(j+1,3)) = wst((k-n-1)*(n-1)+j+1-n)
            enddo
         endif
      enddo

c     Dump of i,j base pairs and energies.
c     Useful only for debugging.

c      write(46,*) '      i          j        v(i,j)       w(i,j)'
c      do diag = 1,n-1
c         do i = 1,n
c            j = diag + i
c            write(46,*) i,j,v(i,j),w(i,j)
c         enddo
c      enddo
c      write(46,*) '      j        w5(j)       w3(j)'
c      do j = 1,n
c         write(46,*) j,w5(j),w3(j)
c      enddo
c      write(46,*) vmin
       if (usage.eq.'text') then
          call system('echo " " > /dev/tty')
          call system('echo "End of Fill" > /dev/tty')
       endif

c     Check for existence of a reasonable folding

      if (vmin.gt.5000) then
         write(6,*) 'No folding possible in this segment.'
         nofold = 1
      endif
      return
      end
c     Used to recall values of V which are actually stored in VST.
      function v(i,j)
      include 'rna.inc'
      if(j.lt.i) then
         v = 0
         return
      endif
      if (i.gt.n) then
        v = vst((n-1)*(i-n-1)+j-n)
      else
        v = vst((n-1)*(i-1)+j)
      endif
      return
      end
 
c     Used to recall values of W which are actually stored in WST.
      function w(i,j)
      include 'rna.inc'
 
      if (i.gt.n) then
        w = wst((n-1)*(i-n-1)+j-n)
      else
        w = wst((n-1)*(i-1)+j)
      endif
      return
      end
 
c     Used to penalize non CG GC closings of helices in multi-branch
c     and external loops

      function au_pen(i,j)
      include 'rna.inc'
      integer inc2(5,5)
      data inc2/0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0/

      au_pen = inc2(numseq(i),numseq(j))*eparam(10)
      return
      end

c     Computes an optimal structure on the subsequence from ii to ji
c     where ii and ji must base-pair with each other. ERROR = 0
c     indicates a normal termination.
c     Base-pair information is stored in the array BASEPR.
      subroutine trace(ii,ji,error)
      include 'rna.inc'

      error = 0
c     Zero the appropriate region of BASEPR.
      if (ji.le.n) then
        do k=ii,ji
          basepr(k) = 0
        enddo
      else
        do k=1,ji-n
          basepr(k) = 0
        enddo
        do k = ii,n
          basepr(k) = 0
        enddo
      endif

c     Initialize the stack of outstanding base-pairs and push ii,ji,
c     v(ii,ji) and 0 on to the stack.
      call initst
      i = ii
      j = ji
      e = v(i,j)
      open = 0
      go to 300
 
100   i = j
c     Pull a fragment ( i to j ) and its expected energy ( e ) from
c     the stack. open = 1 indicates that the free bases are part of
c     an exterior loop. open = 0 (ie. closed) indicates that the
c     free bases are part of a multi-loop.
      stz = pull(i,j,e,open)
#     write(6,*) 'Pull from stack: i,j,e,open = ',i,j,e,open
      if (stz.ne.0) return

      if (open.eq.0) then
c     Multibranch loop
             do while (e.eq.w(i+1,j)+eparam(6).and.force(i).ne.2)
c            Whittle away from the 5' end.
             i = i + 1
             e = w(i,j)
#            write(6,*) 'Reduce 5'' end (mb): i,j,e = ',i,j,e
             if (i.ge.j) goto 100
         enddo
         do while (e.eq.w(i,j-1)+eparam(6).and.force(j).ne.2)
c            Whittle away from the 3' end.
             j = j - 1
             e = w(i,j)
#            write(6,*) 'Reduce 3'' end: (mb) i,j,e = ',i,j,e
             if (i.ge.j) goto 100
         enddo
 
         if (e.eq.v(i+1,j)+au_pen(i+1,j)+eparam(9)+eparam(6)+erg6(j,i+1,i,2)) then
c            i dangles over i+1,j.
             i = i + 1
             e = v(i,j)
         elseif (e.eq.v(i,j-1)+au_pen(i,j-1)+eparam(9)+eparam(6)+erg6(j-1,i,j,1)) then
c            j dangles over i,j-1.
             j = j - 1
             e = v(i,j)
         elseif (e.eq.v(i+1,j-1) + au_pen(i+1,j-1) + eparam(9) + 2*eparam(6) +
     .                erg6(j-1,i+1,i,2) + erg6(j-1,i+1,j,1) ) then
c            Both i and j dangle over i+1,j-1.
             i = i + 1
             j = j - 1
             e = v(i,j)
         endif
c        Check for stem closing a multi-loop.
         if (e.eq.v(i,j)+au_pen(i,j)+eparam(9)) e = v(i,j)
 
      else
 
         do while (j.eq.break.and.e.eq.w3(i+1))
c            Whittle away at the 5' end.
             i = i + 1
#            write(6,*) 'Reduce 5'' end (free): i,j,e = ',i,j,e
             if (i.ge.j) goto 100
         enddo
         do while (i.eq.1+break-n.and.e.eq.w5(j-1))
c            Whittle away at the 3' end.
             j = j - 1
#            write(6,*) 'Reduce 3'' end (free): i,j,e = ',i,j,e
             if (i.ge.j) goto 100
         enddo
 
         if (e.eq.v(i+1,j) + au_pen(i+1,j) + erg6(j,i+1,i,2)) then
c            i dangles over i+1,j.
             i = i + 1
             e = v(i,j)
         elseif (e.eq.v(i,j-1) + au_pen(i,j-1) + erg6(j-1,i,j,1)) then
c            j dangles over i,j-1.
             j = j - 1
             e = v(i,j)
         elseif (e.eq.v(i+1,j-1)+au_pen(i+1,j-1)+erg6(j-1,i+1,i,2)+erg6(j-1,i+1,j,1)) then
c            Both i and j dangle over i+1,j-1.
             i = i + 1
             j = j - 1
             e = v(i,j)
         endif
c        Check for stem closing an exterior loop.
         if (e.eq.v(i,j)+au_pen(i,j)) e = v(i,j)
 
      endif

c     Cannot chop away at the ends any more and still the ends do not
c     base-pair with one another. Structure MUST bifucate (open).
      if (e.ne.v(i,j)) then

         k = i
         do while (k.lt.j-2)
            if (open.eq.0.and.e.eq.w(i,k) + w(k+1,j)) then
c              Best structure on i,j splits into best structures on i,k
c              and k+1,j. Push these fragments on to the stack. (open = 0)
               call push(i,k,w(i,k),0)
               call push(k+1,j,w(k+1,j),0)
               goto 100
            elseif (open.eq.1.and.i.eq.1+break-n) then
               if (e.eq.w5(k) + au_pen(k+1,j) + v(k+1,j)) then
c              Best structure on 1,j splits into best structure on 1,k
c              and the best structure in which k+1 pairs with j.
                  call push(1,k,w5(k),1)
                  i = k + 1
                  goto 300
               elseif (e.eq.w5(k) + erg6(j,k+2,k+1,2) + au_pen(k+2,j) + v(k+2,j)) then
c              Best structure on 1,j splits into best structure on 1,k
c              and the best structure in which k+2 pairs with j and
c              k+1 stacks on this base pair.
                  call push(1,k,w5(k),1)
                  i = k + 2
                  goto 300
               elseif (e.eq.w5(k) + erg6(j-1,k+1,j,1) + au_pen(k+1,j-1) + v(k+1,j-1)) then
c              Best structure on 1,j splits into best structure on 1,k
c              and the best structure in which k+1 pairs with j-1 and
c              j stacks on this base pair.
                  call push(1,k,w5(k),1)
                  i = k + 1
                  j = j - 1
                  goto 300
               elseif (e.eq.w5(k) + erg6(j-1,k+2,k+1,2) + erg6(j-1,k+2,j,1)
     .                  + au_pen(k+2,j-1) + v(k+2,j-1)) then
c              Best structure on 1,j splits into best structure on 1,k
c              and the best structure in which k+2 pairs with j-1 and
c              k+1 and j both stack on this base pair.
                  call push(1,k,w5(k),1)
                  i = k + 2
                  j = j - 1
                  goto 300
               endif
            elseif (open.eq.1.and.j.eq.break) then
               if (e.eq.v(i,k+2)+au_pen(i,k+2)+w3(k+3)) then
c              Best structure on i,n splits into best structure in which
c              i pairs with k+2 and the best structure on k+3 to n.
                  call push(k+3,n,w3(k+3),1)
                  j = k + 2
                  goto 300
               elseif (e.eq.v(i+1,k+2)+au_pen(i+1,k+2)+erg6(k+2,i+1,i,2)+w3(k+3)) then
c              Best structure on i,n splits into best structure in which
c              i+1 pairs with k+2 with i stacked on this base pair and the 
c              best structure on k+3 to n.
                  call push(k+3,n,w3(k+3),1)
                  i = i + 1
                  j = k + 2
                  goto 300
               elseif (e.eq.v(i,k+1)+au_pen(i,k+1)+erg6(k+1,i,k+2,1)+w3(k+3)) then
c              Best structure on i,n splits into best structure in which
c              i pairs with k+1 with k+2 stacked on this base pair and the 
c              best structure on k+3 to n.
                  call push(k+3,n,w3(k+3),1)
                  j = k + 1
                  goto 300
               elseif (e.eq.v(i+1,k+1)+au_pen(i+1,k+1)+erg6(k+1,i+1,i,2)+
     .                      erg6(k+1,i+1,k+2,1)+w3(k+3)) then
c              Best structure on i,n splits into best structure in which
c              i+1 pairs with k+1 with both i and k+2 stacked on this
c              base pair and the best structure on k+3 to n.
                  call push(k+3,n,w3(k+3),1)
                  i = i + 1
                  j = k + 1
                  goto 300
               endif

            endif
            k = k + 1
         enddo
c        Structure will not split. Error
         ii = hstnum(i)
         ji = hstnum(j)
         error = 10
         return

      endif
 
c     Base-pair found. All base-pairs are stored in the range 1 <= i < j <= n.
c     If i and j form a base-pair, then basepr(i) = j and basepr(j) = i.
300   e = v(i,j)
      if (j.le.n) then
         basepr(i) = j
         basepr(j) = i
      elseif (i.gt.n) then
         basepr(i-n) = j-n
         basepr(j-n) = i-n
         i = i - n
         j = j - n
      else
         basepr(j-n) = i
         basepr(i) = j-n
      endif
#     write(6,*) 'Base pair: ',i,j
      open = 0

c     Perhaps i,j stacks over i+1,j-1?
      if (i.ne.break.and.j.ne.break+1) then
         if (e.eq.erg2(i,j) + v(i+1,j-1)) then
            i = i + 1
            j = j - 1
            e = v(i,j)
            goto 300
         endif
      endif
 
c     Perhaps i,j closes a hairpin loop?
      if (e.eq.erg4(i,j)) goto 100
 
c     Perhaps i,j closes a multi-loop?
      k = i
      do while (k.le.j)
         if (k.ne.break.and.k.ge.i+2.and.k.le.j-3) then
            if (e.eq.w(i+1,k)+w(k+1,j-1)+au_pen(i,j)+eparam(9)+eparam(5)) then
c              Multi-loop. No dangling ends on i,j.
               call push(i+1,k,w(i+1,k),0)
               call push(k+1,j-1,w(k+1,j-1),0)
               goto 100
            elseif (e.eq.erg6(i,j,i+1,1)+w(i+2,k)+w(k+1,j-1)+au_pen(i,j)+
     .                eparam(9) + eparam(6) + eparam(5)) then
c              Multi-loop. i+1 dangles over the i,j base-pair.
               call push(i+2,k,w(i+2,k),0)
               call push(k+1,j-1,w(k+1,j-1),0)
               goto 100
            elseif (e.eq.erg6(i,j,j-1,2)+w(i+1,k)+w(k+1,j-2)+au_pen(i,j)+
     .                eparam(9) + eparam(6) + eparam(5)) then
c              Multi-loop. j-1 dangles over the i,j base-pair.
               call push(i+1,k,w(i+1,k),0)
               call push(k+1,j-2,w(k+1,j-2),0)
               goto 100
            elseif (e.eq.erg6(i,j,i+1,1)+erg6(i,j,j-1,2)+w(i+2,k)
     .         +w(k+1,j-2)+au_pen(i,j)+eparam(9)+2*eparam(6)+eparam(5)) then
c              Multi-loop. Both i+1 and j-1 dangle over the i,j base-pair.
               call push(i+2,k,w(i+2,k),0)
               call push(k+1,j-2,w(k+1,j-2),0)
               goto 100
            endif
         elseif (k.eq.break.and.j.gt.n) then
            if (e.eq.w3(i+1) + w5(j-n-1)+au_pen(i,j)) then
c              Exterior loop. No ends dangling on i,j.
               if (i+1.lt.n) call push(i+1,n,w3(i+1),1)
               if (1.lt.j-n-1) call push(1,j-n-1,w5(j-n-1),1)
               goto 100
            elseif  (e.eq.erg6(i,j,i+1,1)+w3(i+2)+w5(j-n-1)+au_pen(i,j)) then
c              Exterior loop. i+1 dangles over the i,j base-pair.
               if (i+2.lt.n) call push(i+2,n,w3(i+2),1)
               if (1.lt.j-n-1) call push(1,j-n-1,w5(j-n-1),1)
               goto 100
            elseif (e.eq.erg6(i,j,j-1,2)+w3(i+1)+w5(j-n-2)+au_pen(i,j)) then
c              Exterior loop. j-1 dangles over the i,j base-pair.
               if (i+1.lt.n) call push(i+1,n,w3(i+1),1)
               if (1.lt.j-n-2) call push(1,j-n-2,w5(j-n-2),1)
               goto 100
            elseif (e.eq.erg6(i,j,i+1,1)+erg6(i,j,j-1,2)+w3(i+2)+
     .               w5(j-n-2)+au_pen(i,j)) then
c              Exterior loop. Both i+1 and j-1 dangle over the i,j base-pair.
               if (i+2.lt.n) call push(i+2,n,w3(i+2),1)
               if (1.lt.j-n-2) call push(1,j-n-2,w5(j-n-2),1)
               goto 100
            endif
         endif
         k = k + 1
      enddo
  
c     None of the above work. i,j MUST close a bulge or interior loop.
500   do d = j-i-3,1,-1
        do ip = i+1,j-1-d
           jp = d+ip
           if (j-i-2-d.gt.eparam(7)) then
c             Error, bulge or interior loop not found.
              ii = hstnum(i)
              ji = hstnum(j)
              error = 11
              return
           endif
           if (abs(ip-i+jp-j).le.eparam(8)) then
              if (e.eq.erg3(i,j,ip,jp)+v(ip,jp)) then
                 i = ip
                 j = jp
                 e = v(i,j)
                 goto 300
              endif
           endif
        enddo
      enddo
c    Error, bulge or interior loop not found.
      ii = hstnum(i)
      ji = hstnum(j)
      error = 11
      return
      end
 
c     Store results of a SAVE run for a continuation run.
      subroutine putcont
      include 'rna.inc'
 
      write(30) n,break,nsave,vmin,listsz,seqlab,lorc,maxbp
      write(30) asint3,asint5,bulge,dangle,eparam,hairpin,inter,
     .sint2,sint4,sint6,stack,tloop,triloop,tstkh,tstki
      write(30) (vst(i),i=1,n*n)
      write(30) (wst(i),i=1,n*n)
      write(30) (w5(i),i=0,n)
      write(30) (w3(i),i=1,n+1)
      write(30) (seq(i),i=nsave(1),nsave(2))
      write(30) ((list(i,j),i=1,listsz),j=1,4)
      write(30) tloop,numtloops,triloop,numtriloops
      write(30) (poppen(i),i=1,4),maxpen,prelog
      return
      end
c     Read results from a SAVE run for a CONTINUATION run.
      subroutine getcont
      include 'rna.inc'
 
      read(30,err=10) n,break,nsave,vmin,listsz,seqlab,lorc,maxbp
      read(30,err=10)asint3,asint5,bulge,dangle,eparam,hairpin,
     .inter,sint2,sint4,sint6,stack,tloop,triloop,tstkh,tstki
      read(30,err=10) (vst(i),i=1,n*n)
      read(30,err=10) (wst(i),i=1,n*n)
      read(30,err=10) (w5(i),i=0,n)
      read(30,err=10) (w3(i),i=1,n+1)
      read(30,err=10) (seq(i),i=nsave(1),nsave(2))
      read(30,err=10) ((list(i,j),i=1,listsz),j=1,4)
      read(30,err=10) tloop,numtloops,triloop,numtriloops
      read(30,err=10) (poppen(i),i=1,4),maxpen,prelog
      goto 11
 
10    call errmsg(40,0,0)
 
11    return
      end

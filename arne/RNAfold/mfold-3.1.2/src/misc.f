      subroutine process
c     Process RNA sequence to be folded.
      include 'rna.inc'
c     Selected fragment is from nsave(1) to nsave(2) in historical
c     numbering.
      do i = 1,n
         newnum(i) = 0
      enddo
c     list contains information on excisions, and on forced or prohibited
c     base-pairs.
      ptr = 0
      maxbp = infinity
      do while (ptr.lt.listsz)
         ptr = ptr + 1
         if (list(ptr,1).eq.4) then
c     Closed excision beween list(ptr,2) & list(ptr,3) (historical numbering).
            do i = list(ptr,2)+4,list(ptr,3)-1
               newnum(i) = 1
            enddo
         else if (list(ptr,1).eq.5) then
c     Open excision beween list(ptr,2) & list(ptr,3) (historical numbering).
            do i = list(ptr,2),list(ptr,3)
               newnum(i) = 1
            enddo
         else if (list(ptr,1).eq.9) then
c     Set maxbp, the maximum distance between paired bases.
            maxbp = list(ptr,2)
         endif
      enddo

      n = 0
      do k = nsave(1),nsave(2)
c     Generate new numbering of fragment ( 1 to n ).
         if (newnum(k).eq.0) then
            n = n+1
            newnum(k) = n
         else
c           An excised base gets a new numbering of 0.
            newnum(k) = 0
         endif
      enddo
 
c     When folding circular RNA/DNA, the "break" point is
c     set to 3*n, which is effectively infinite. Thus no
c     special rules pertaining to loops containing the sequence
c     ends are invoked.
      if (lorc.eq.'c') then
         break = 3*n
      else
         break = n
      endif

      if (n*2.gt.fldmax) goto 700
 
c     Initialize the force and vst arrays.
c     'strand' and 'scheck' were added on July 22, 2002 to quickly determine
c     if 2 strands are present using the 'LLL' trick.
      scheck = 0
      do i = 1,n
        force(i) = 0
        strand(i) = 1
        if (cntrl(1).ne.2) then
           do j = i,i+n-1
             vst((n-1)*(i-1)+j) = 0
           enddo
        endif
      enddo

      do k = nsave(1),nsave(2)
         i = newnum(k)
         if (i.gt.0) then
c           Non-excised bases are examined to determine their type.
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
            hstnum(i) = k
            numseq(i) = 5
            aux(i) = ' '
            if (seq(k).eq.'A'.or.seq(k).eq.'a') then
               numseq(i) = 1
               if (scheck.eq.-1) strand(k) = 2
            elseif (seq(k).eq.'B'.or.seq(k).eq.'b') then
               numseq(i) = 1
               force(i) = 3
               if (scheck.eq.-1) strand(k) = 2
            elseif (seq(k).eq.'W'.or.seq(k).eq.'w') then
               numseq(i) = 1
               force(i) = 4
               if (scheck.eq.-1) strand(k) = 2
            elseif (seq(k).eq.'C'.or.seq(k).eq.'c') then
               numseq(i) = 2
               if (scheck.eq.-1) strand(k) = 2
            elseif (seq(k).eq.'Z'.or.seq(k).eq.'z') then
               numseq(i) = 2
               force(i) = 3
               if (scheck.eq.-1) strand(k) = 2
            elseif (seq(k).eq.'X'.or.seq(k).eq.'x') then
               numseq(i) = 2
               force(i) = 4
               if (scheck.eq.-1) strand(k) = 2
            elseif (seq(k).eq.'G'.or.seq(k).eq.'g') then
               numseq(i) = 3
               if (scheck.eq.-1) strand(k) = 2
            elseif (seq(k).eq.'H'.or.seq(k).eq.'h') then
               numseq(i) = 3
               force(i) = 3
               if (scheck.eq.-1) strand(k) = 2
            elseif (seq(k).eq.'Y'.or.seq(k).eq.'y') then
               numseq(i) = 3
               force(i) = 4
               if (scheck.eq.-1) strand(k) = 2
            elseif (seq(k) .eq. 'U'.or.seq(k).eq.'u') then
               numseq(i) = 4
               if (scheck.eq.-1) strand(k) = 2
            elseif (seq(k) .eq. 'V'.or.seq(k).eq.'v') then
               numseq(i) = 4
               force(i) = 3
               if (scheck.eq.-1) strand(k) = 2
            elseif (seq(k).eq.'Z'.or.seq(k).eq.'z') then
               numseq(i) = 4
               force(i) = 4
               if (scheck.eq.-1) strand(k) = 2
            elseif (seq(k) .eq. 'T'.or.seq(k).eq.'t') then
               numseq(i) = 4
               if (scheck.eq.-1) strand(k) = 2
            elseif (seq(k) .eq. 'W'.or.seq(k).eq.'w') then
               numseq(i) = 4
               force(i) = 3
               if (scheck.eq.-1) strand(k) = 2
             elseif (seq(k) .eq. 'L'.or.seq(k).eq.'l') then
c               seq(k) = ' '
                seq(k) = 'L'
                if (scheck.ge.0) scheck = scheck + 1
                if (scheck.eq.3) then
                   if (seq(k-2).eq.'L'.and.seq(k-1).eq.'L') scheck = -1
                endif
            endif
         endif
      enddo
 
      ptr = 0
500   if (ptr.eq.listsz) goto 600
      ptr = ptr + 1
      if (list(ptr,1).lt.0) goto 580
      if (list(ptr,1).gt.7) goto 500
      i = newnum(list(ptr,2))
      j = newnum(list(ptr,3))
      k = list(ptr,4)
      if (list(ptr,1).eq.2.or.list(ptr,1).eq.6) k = list(ptr,3)
      goto (500,520,530,540,500,560,570),list(ptr,1)

c     Version 3.0 marks a new era in constrained base pairs.
c     Prohibited base pairs will be handled as before.
c
c     When a base, I, is forced to pair, it will be marked as '2'
c     in the 'force' array. This is as before.
c
c     When a base pair I.J is forced, both I and J are marked as '2'
c     in 'force' AND all base pairs I'.J' (!=I.J) are prohibited as in 
c     570 below. Thus I and J can only pair with each other, and marking 
c     them with '2' in 'force' reduces the problem to forcing single
c     bases to be paired.
c     The single bit 'fce' array will be used differently. Instead of 
c     fce(i,j)=true designating a forced base pair, fce(i,j)=true will 
c     now mean that the interval from i to j (NON-INCLUSIVE) contains at
c     least one base that must be paired. Why is this done? 
c     For example, if fce(i,j)=true, then i,j CANNOT close a hairpin loop.
c     If i < i' < j' < j, and fce(i,i')=true, then i.j and i'.j' cannot
c     close an interior or bulge loop. Large penalty energies will
c     be used to enforce these rules.

c     Also new is the Prohibit region option. The user inputs 4 numbers, 
c     i,j,k and l satisfying 1 <= i < j < k < l <= n (n is 3' base).
c     These numbers define 2 distinct segments; the 5' segment,
c     5'-i- ... -j-3' and the 3' segment, 5'-k- ... -l-3'.
c     The constraint is that base pairs between the 5' and 3' segments
c     are forbidden.

c     Force bases I to I+K-1 to be double-stranded.
520   do x = i,i+k-1
         force(x) = 2
         aux(x) = 'F'
      enddo
      goto 500
c     Force base-pairs I.J , I+1.J-1 , ... I+K-1.J-K+1.
530   do x = 0,k-1
        force(i+x) = 2
        force(j-x) = 2
        aux(j-x) = ')'
        aux(i+x) = '('
         do l = 1,n
            if (l.lt.i+x) then 
               vst((n-1)*(l-1) + i+x) = infinity
               vst((n-1)*(i+x-1) + l+n) = infinity
               vst((n-1)*(l-1) + j-x) = infinity
               vst((n-1)*(j-x-1) + l+n) = infinity
            else if (l.gt.i+x.and.l.lt.j-x) then
               vst((n-1)*(i+x-1) + l) = infinity
               vst((n-1)*(l-1) + i+x+n) = infinity
               vst((n-1)*(l-1) + j-x) = infinity
               vst((n-1)*(j-x-1) + l+n) = infinity
            else if (l.gt.j-x) then
               vst((n-1)*(i+x-1) + l) = infinity
               vst((n-1)*(l-1) + i+x+n) = infinity
               vst((n-1)*(j-x-1) + l) = infinity
               vst((n-1)*(l-1) + j-x+n) = infinity
            endif
         enddo
      enddo
c     Eliminate potential pseudoknots to the helix: 
c     Prohibit ALL base pairs l1.l2
c     where l1 < i, i+k-1 < l2 or l1 < j-k+1, j < l2
      if (i.gt.1) then
         do l1 = 1,i-1
            do l2 = i+k,j-k
               vst((n-1)*(l1-1) + l2) = infinity
               vst((n-1)*(l2-1) + l1+n) = infinity
            enddo
         enddo
      endif
      if (j.lt.n) then
         do l1 = i+k,j-k
            do l2 = j+1,n
               vst((n-1)*(l1-1) + l2) = infinity
               vst((n-1)*(l2-1) + l1+n) = infinity
               aux(j-x) = '}'
               aux(i+x) = '{'
            enddo
         enddo
      endif
      goto 500
c     Force the ends of a closed excision to base-pair.
 540  do ii = i+1,i+3
        seq(ii) = ' '
      enddo
      force(i) = 2
      force(j) = 2
      aux(j) = ')'
      aux(i) = '('
      do l = 1,n
         if (l.lt.i) then 
            vst((n-1)*(l-1) + i) = infinity
            vst((n-1)*(i-1) + l+n) = infinity
         else if (l.gt.i) then
            vst((n-1)*(i-1) + l) = infinity
            vst((n-1)*(l-1) + i+n) = infinity
         endif
         if (l.lt.j) then 
            vst((n-1)*(l-1) + j) = infinity
            vst((n-1)*(j-1) + l+n) = infinity
         else if (l.gt.j) then
            vst((n-1)*(j-1) + l) = infinity
            vst((n-1)*(l-1) + j+n) = infinity
         endif
      enddo
      goto 500
c     Prohibit bases I to I+K-1 from base-pairing.
560   do ii = i,i+k-1
          force(ii) = 1
          aux(ii) = 'P'
      enddo
      goto 500
c     Prohibit the base-pairs I.J , I+1,J-1 , ... I+K-1.J-K+1.
570   if (cntrl(1).ne.2) then
         do x = 0,k-1
           vst((n-1)*(i+x-1)+j-x) = infinity
           vst((n-1)*(j-x-1)+i+x+n) = infinity
         enddo
      endif
      goto 500
c     Prohibit base pairs between 5'-i- ... -j-3' and 5'-k- ... -l-3'
580   if (cntrl(1).ne.2) then
         i = newnum(-list(ptr,1))
         j = newnum(list(ptr,2))
         k = newnum(list(ptr,3))
         l = newnum(list(ptr,4))
         if (-list(ptr,1).lt.nsave(1)) then
            i = 1
            if (list(ptr,2).lt.nsave(1)) then
               j = 0
            elseif (list(ptr,2).gt.nsave(2)) then
               j = n
            endif
         elseif (-list(ptr,1).gt.nsave(2)) then
            i = 1
            j = 0
         else
            if (list(ptr,2).gt.nsave(2)) j = n
         endif
         if (list(ptr,3).lt.nsave(1)) then
            k = 1
            if (list(ptr,4).lt.nsave(1)) then
               l = 0
            elseif (list(ptr,4).gt.nsave(2)) then
               l = n
            endif
         elseif (list(ptr,3).gt.nsave(2)) then
            k = 1
            l = 0
         else
            if (list(ptr,4).gt.nsave(2)) l = n
         endif
         do ii = i,j
            do jj = k,l
               if (ii.lt.jj) then
                  vst((n-1)*(ii-1)+jj) = infinity
                  vst((n-1)*(jj-1)+ii+n) = infinity
               elseif (ii.gt.jj) then
                  vst((n-1)*(jj-1)+ii) = infinity
                  vst((n-1)*(ii-1)+jj+n) = infinity
               endif
            enddo
         enddo
      endif
      goto 500
c     Double up the sequence.
600   do i = 1,n
        hstnum(i+n) = hstnum(i)
        force(i+n) = force(i)
        numseq(i+n) = numseq(i)
      enddo
c     Now that all the forced bases and base pairs have been
c     catalogued, it is time to fill the fce array.
c     fce(i,j) = true iff force(k) = 2 for some i < k < j.
c     fce will be filled recursively.
      do j = 3,n
         i = j-2
         do while (force(i+1).ne.2.and.i.gt.0)
            i = i - 1
         enddo
         if (i.gt.0) then
            do ip = i,1,-1
               call sfce(ip,j)
               call sfce(j,ip+n)
            enddo
         endif
      enddo
700   return
      end
 
c     Error message subroutine.
      subroutine errmsg(err,i,j)
      include 'rna.inc'
 
      if (err.eq.10) then
         write(6,10) i,j
         call exit(1)
      endif
      if (err.eq.11) then
         write(6,11) i,j
         call exit(1)
      endif
      if (err.eq.12) then
         write(6,12) i,j
         call exit(1)
      endif
      if (err.eq.20) then
         write(6,20) i,j
         call exit(1)
      endif
      if (err.eq.21) then
         write(6,21)
         err = 0
      endif
      if (err.eq.30) write(6,30) i
      if (err.eq.31) then
         write(6,31) sortmax,i,j
         err = 0
      endif
      if (err.eq.40) then
         write(6,40)
         call exit(1)
      endif
      return
 
10    format('STOP: Open bifurcation not found between ',i6,' and ',i6)
11    format('STOP: Bulge or interior loop closed by (',i6,',',i6,') not found')
12    format('STOP: Closed bifurcation not found between ',i6,' and ',i6)
20    format('STOP: Base pair between ',i6,' and ',i5,' conflicts with ',
     . 'at least one other pair')
21    format(' Buffer overflow in lineout')
30    format(' End reached at traceback ',i5)
31    format(' More than ',i6,' basepairs in sort at (',i6,',',i6,')')
40    format('STOP: Premature end of save file')
      end

c     Initialize the stack.
      subroutine initst
      implicit integer (a-z)
      integer stk(500,4),sp
      common /stk/ stk,sp
 
      sp = 0
      return
      end
c     Add A,B,C,D to the bottom of the stack.
      subroutine push(a,b,c,d)
      implicit integer (a-z)
      integer stk(500,4),sp
      common /stk/ stk,sp
 
      sp = sp + 1
      if (sp.gt.500) then
         write(6,*) 'STOP: STACK OVERFLOW'
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
      integer stk(500,4),sp
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
      include 'rna.inc'
c      include 'efn.inc'
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
c      write(unit,103) hstn1,hstn2,ctlabel,energy
      write(unit,103) hstn1,hstn2,seqlab,energy
 103  format('Folding bases ',i6,' to ',i6,' of ',a50,/,
     .       ' dG = ',f8.1/)
c     .       'Computed energy  =  ',f8.1/)
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

c     Generates a CT file.
      subroutine ct(r)
      include 'rna.inc'
      real r
 
      write(21,100) n,r,seqlab
      do k = 1,n
         k1 = k - 1
         k2 = k + 1
         if (k.eq.1) then
            if (lorc.eq.'c') k1 = n
         endif
         if (k.eq.n) then
            if (lorc.eq.'c') then
               k2 = 1
            else
               k2 = 0
            endif
         endif
         write(21,200) k,seq(hstnum(k)),k1,k2,basepr(k),hstnum(k)
      enddo
      return
 
100   format(i5,1x,'dG = ',f7.1,4x,a50)
200   format(i5,1x,a1,3i6,i7)
      end

c     This subroutine performs a detailed calculation of the energy 
c     contributions of every loop and stack in an RNA folding.
c     The folding is assumed to be on seq(nsave(1))---seq(nsave(2))
c     and the base pairs are stored in basepr(i), i=1, 2, ... ,n
c
      subroutine erg_det(deltag,iret,jret)
      include 'rna.inc'
      logical instack
 
c     Initialize the stack of outstanding base-pairs and push 1,n
c     0 and 1 on to the stack.
c     The 4th argument in the stack is the 'open' argument.
c     open = 1 : exterior loop; ends do not necessarily pair with each other
c     For circular nucleic acids, the ends cannot pair with each other.
c     open = 0 : ends pair

      call initst
      call push(1,n,0,1)
      write(22,1005) seqlab,float(deltag)/prec
 1005 format(a50/' dG = ',f8.2/)

c     Pull a fragment ( i to j ) and +-v(i,j) or +-v(j,i) from the stack.
c     This energy is 0 when open = 1.
      instack = .false.
      do while (pull(i,j,ddg,open).eq.0)
         single = 0
         if (open.eq.1) then
            double = 0
            ip = i
         else
            ip = i + 1
            double = 1
         endif
         do while (ip.le.j)
            if (basepr(ip).eq.0) then
               single = single + 1
            elseif (basepr(ip).gt.ip) then
               double = double + 1
               if (ip.lt.iret.and.jret.lt.basepr(ip)) then
                  ddg = ddg + v(basepr(ip),ip+n)
                  call push(ip,basepr(ip),-v(basepr(ip),ip+n),0)
               elseif (ip.eq.iret.and.jret.eq.basepr(ip)) then
                  ddg = ddg + v(basepr(ip),ip+n)
                  call push(ip,basepr(ip),v(ip,basepr(ip)),0)
               else
                  ddg = ddg - v(ip,basepr(ip))
                  call push(ip,basepr(ip),v(ip,basepr(ip)),0)
               endif
               ip = basepr(ip)
            endif
            ip = ip + 1
         enddo
         if (open.eq.1.and.lorc.eq.'l') then
            write(22,1010) float(ddg)/prec,single,double
 1010       format('External loop:	ddG =',f7.2,1x,i3,' ss bases & ',
     .              i2,' closing helices.')
         elseif (open.eq.1.and.lorc.eq.'c') then
            if (double.eq.1) then
               write(22,1012) float(ddg)/prec
 1012          format('Sequence ends are in a hairpin loop. ddg =',f7.2)
            elseif (double.eq.2) then
               if (single.eq.0) then
                  write(22,1014) float(ddg)/prec
 1014             format('Sequence ends are in a stack. ddg =',f7.2)
                  dghelix = ddg
                  instack = .true.
               elseif ((basepr(1).ne.0.and.basepr(n).ne.0).or.
     .                 basepr(basepr(ip-1)-1).gt.0) then
                  write(22,1016) float(ddg)/prec
 1016             format('Sequence ends are in a bulge loop. ddg =',f7.2)
               else
                  write(22,1018) float(ddg)/prec
 1018             format('Sequence ends are in an interior loop. ddg =',f7.2)
               endif
            else
               write(22,1019) float(ddg)/prec
 1019          format('Sequence ends are in a multi-loop. ddg =',f7.2)
            endif
         elseif (double.gt.2) then
            if (instack) then
               write(22,1024) float(dghelix)/prec,hsize
 1024          format('Helix:		ddG =',f7.2,1x,i3,' base pairs.')
               instack = .false.
            endif
            write(22,1020) float(ddg)/prec,seq(hstnum(i)),hstnum(i),
     .                    seq(hstnum(j)),hstnum(j),single,double
 1020       format('Multi-loop:	ddG =',f7.2,' External closing pair is '
     .           ,a1,'(',i6,')-',a1,'(',i6,')',/,'			',
     .           i3,' ss bases & ',i2,' closing helices.')
         elseif (double.eq.2) then
            if (single.eq.0) then
               write(22,1030) float(ddg)/prec,seq(hstnum(i)),hstnum(i),
     .                       seq(hstnum(j)),hstnum(j)
 1030          format('Stack:		ddG =',f7.2,
     .          ' External closing pair is ',a1,'(',i6,')-',a1,'(',i6,')')
               if (instack) then
                  dghelix = dghelix + ddg
                  hsize = hsize + 1
               else
                  dghelix = ddg
                  instack = .true.
                  hsize = 2
               endif
            elseif (basepr(i+1)+basepr(j-1).gt.0) then
               if (instack) then
                  write(22,1024) float(dghelix)/prec,hsize
                  instack = .false.
               endif
               write(22,1040) float(ddg)/prec,seq(hstnum(i)),hstnum(i),
     .                       seq(hstnum(j)),hstnum(j)
 1040          format('Bulge loop:	ddG =',f7.2,
     .             ' External closing pair is ',a1,'(',i6,')-',a1,'(',i6,')')
            else
               if (instack) then
                  write(22,1024) float(dghelix)/prec,hsize
                  instack = .false.
               endif
               write(22,1050) float(ddg)/prec,seq(hstnum(i)),hstnum(i),
     .                       seq(hstnum(j)),hstnum(j)
 1050          format('Interior loop:	ddG =',f7.2,
     .              ' External closing pair is ',a1,'(',i6,')-',a1,'(',i6,')')
            endif
         elseif (double.eq.1) then
            if (instack) then
               write(22,1024) float(dghelix)/prec,hsize
               instack = .false.
            endif
            write(22,1060) float(ddg)/prec,seq(hstnum(i)),hstnum(i),
     .                    seq(hstnum(j)),hstnum(j)
 1060       format('Hairpin loop:	ddG =',f7.2,
     .            '          Closing pair is ',a1,'(',i6,')-',a1,'(',i6,')')
         else
            write(6,*) 'STOP: Big mess in erg_det function'
            call exit(1)
         endif
      enddo
      return
      end
c     This subroutine performs a detailed calculation of the energy 
c     contributions of every loop and stack in an RNA folding.
c     The folding is assumed to be on seq(nsave(1))---seq(nsave(2))
c     and the base pairs are stored in basepr(i), i=1, 2, ... ,n
c
      subroutine erg_det_html(deltag,iret,jret,strnum)
      include 'rna.inc'
      character*80 path
      character*29 bp_label,tmp_lab
      logical instack
 
c     Initialize the stack of outstanding base-pairs and push 1,n
c     0 and 1 on to the stack.
c     The 4th argument in the stack is the 'open' argument.
c     open = 1 : exterior loop; ends do not necessarily pair with each other
c     For circular nucleic acids, the ends cannot pair with each other.
c     open = 0 : ends pair

      call initst
      call push(1,n,0,1)

c     Get path for ddG image
      call getenv('MFOLDLIB',path)
      in = index(path,' ')
      if (path.eq.'     ') then
         call getenv('MFOLD',path)
         in = index(path,' ')
         path(in:in+12) = '/dat/ddG.gif"'
      else
         path(in:in+8) = '/ddG.gif"'
      endif

      write(22,1005) seqlab,float(deltag)/prec,path
 1005 format('<BR>',a50,'<BR> dG = ',f8.2,/'<P>',
     .  '<TABLE BORDER=3 CELLPADDING=4><TR><TH>Structural element</TH><TH>'
     .  ,'<IMG'/'SRC="',a80/' ALT="ddG"></TH><TH>Information</TH></TR>')

c     Pull a fragment ( i to j ) and +-v(i,j) or +-v(j,i) from the stack.
c     This energy is 0 when open = 1.
      instack = .false.
      do while (pull(i,j,ddg,open).eq.0)
         single = 0
         if (open.eq.1) then
            double = 0
            ip = i
         else
            ip = i + 1
            double = 1
         endif
         do while (ip.le.j)
            if (basepr(ip).eq.0) then
               single = single + 1
            elseif (basepr(ip).gt.ip) then
               double = double + 1
               if (ip.lt.iret.and.jret.lt.basepr(ip)) then
                  ddg = ddg + v(basepr(ip),ip+n)
                  call push(ip,basepr(ip),-v(basepr(ip),ip+n),0)
               elseif (ip.eq.iret.and.jret.eq.basepr(ip)) then
                  ddg = ddg + v(basepr(ip),ip+n)
                  call push(ip,basepr(ip),v(ip,basepr(ip)),0)
               else
                  ddg = ddg - v(ip,basepr(ip))
                  call push(ip,basepr(ip),v(ip,basepr(ip)),0)
               endif
               ip = basepr(ip)
            endif
            ip = ip + 1
         enddo
         tmp_lab = bp_label(strnum,i,j)
         if (open.eq.1.and.lorc.eq.'l') then
            write(22,1010) float(ddg)/prec,single,double
 1010       format('<TR><TD>External loop</TD><TD ALIGN=RIGHT>',f7.2,
     .         '</TD><TD>',i3,' ss bases & ',i2,' closing helices.</TD></TR>')
         elseif (open.eq.1.and.lorc.eq.'c') then
            if (double.eq.1) then
               write(22,1012) float(ddg)/prec
 1012          format('<TR><TD>Sequence ends are<BR>in a hairpin loop.',
     .         '</TD><TD ALIGN=RIGHT>',f7.2,'</TD><TD>&nbsp;</TD></TR>')
            elseif (double.eq.2) then
               if (single.eq.0) then
                  write(22,1014) float(ddg)/prec
 1014             format('<TR><TD>Sequence ends are<BR>in a stack.</TD>',
     .                 '<TD ALIGN=RIGHT>',f7.2,'</TD><TD>&nbsp;</TD></TR>')
                  dghelix = ddg
                  instack = .true.
               elseif ((basepr(1).ne.0.and.basepr(n).ne.0).or.
     .                 basepr(basepr(ip-1)-1).gt.0) then
                  write(22,1016) float(ddg)/prec
 1016             format('<TR><TD>Sequence ends are<BR>in a bulge loop.'
     .            ,'</TD><TD ALIGN=RIGHT>',f7.2,'</TD><TD>&nbsp;</TD></TR>')
               else
                  write(22,1018) float(ddg)/prec
 1018             format('<TR><TD>Sequence ends are<BR>in an interior loop.'
     .              ,'</TD><TD ALIGN=RIGHT>',f7.2,'</TD><TD>&nbsp;</TD></TR>')
               endif
            else
               write(22,1019) float(ddg)/prec
 1019          format('<TR><TD>Sequence ends are<BR>in a multi-loop.',
     .            '</TD><TD ALIGN=RIGHT>',f7.2,'</TD><TD>&nbsp;</TD></TR>')
            endif
         elseif (double.gt.2) then
            if (instack) then
               write(22,1024) float(dghelix)/prec,hsize
 1024          format('<TR><TH>Helix</TH><TD ALIGN=RIGHT>',f7.2,
     .              '</TD><TD>',i3,' base pairs.</TD></TR>')
               instack = .false.
            endif
            write(22,1020) float(ddg)/prec,tmp_lab,seq(hstnum(i)),hstnum(i),
     .                    seq(hstnum(j)),hstnum(j),single,double
 1020       format('<TR><TD>Multi-loop</TD>','<TD ALIGN=RIGHT>',f7.2,
     .           '</TD><TD><A NAME=',a29,'>External closing pair is</A> '
     .           ,a1,'<SUP>',i6,'</SUP>-',a1,'<SUP>',i6,'</SUP><BR>',i3,
     .           ' ss bases & ',i2,' closing helices.</TD></TR>')
         elseif (double.eq.2) then
            if (single.eq.0) then
               write(22,1030) tmp_lab,float(ddg)/prec,
     .              seq(hstnum(i)),hstnum(i),seq(hstnum(j)),hstnum(j)
 1030          format('<TR><TD><A NAME=',a29,'>Stack</A></TD><TD ALIGN=RIGHT>'
     .                ,f7.2,'</TD><TD>','External closing pair is ',a1,
     .                '<SUP>',i6,'</SUP>-',a1,'<SUP>',i6,'</SUP></TD></TR>')
               if (instack) then
                  dghelix = dghelix + ddg
                  hsize = hsize + 1
               else
                  dghelix = ddg
                  instack = .true.
                  hsize = 2
               endif
            elseif (basepr(i+1)+basepr(j-1).gt.0) then
               if (instack) then
                  write(22,1024) float(dghelix)/prec,hsize
                  instack = .false.
               endif
               write(22,1040) tmp_lab,float(ddg)/prec,seq(hstnum(i)),hstnum(i),
     .                       seq(hstnum(j)),hstnum(j)
 1040          format('<TR><TD><A NAME=',a29,'>Bulge loop</A></TD>',
     .                '<TD ALIGN=RIGHT>',f7.2,'</TD>',
     .                '<TD>External closing pair is ',a1,'<SUP>',i6,
     .                '</SUP>-',a1,'<SUP>',i6,'</SUP></TD></TR>')
            else
               if (instack) then
                  write(22,1024) float(dghelix)/prec,hsize
                  instack = .false.
               endif
               write(22,1050) tmp_lab,float(ddg)/prec,seq(hstnum(i)),
     .                        hstnum(i),seq(hstnum(j)),hstnum(j)
 1050          format('<TR><TD><A NAME=',a29,'>Interior loop</A></TD>',
     .                '<TD ALIGN=RIGHT>',f7.2,'</TD>'
     .               ,'<TD>External closing pair is ',a1,'<SUP>',i6,
     .                '</SUP>-',a1,'<SUP>',i6,'</SUP></TD></TR>')
            endif
         elseif (double.eq.1) then
            if (instack) then
               write(22,1024) float(dghelix)/prec,hsize
               instack = .false.
            endif
            write(22,1060) tmp_lab,float(ddg)/prec,seq(hstnum(i)),hstnum(i),
     .                    seq(hstnum(j)),hstnum(j)
 1060       format('<TR><TD><A NAME=',a29,'>Hairpin loop</A></TD>',
     .       '<TD ALIGN=RIGHT>',f7.2,'</TD><TD>Closing pair is ',a1,'<SUP>',
     .       i6,'</SUP>-',a1,'<SUP>',i6,'</SUP></TD></TR>')
         else
            write(6,*) 'STOP: Big mess in erg_det function'
            call exit(1)
         endif
      enddo
      write(22,*) '</TABLE><P>'
      return
      end

c     Function to create anchor label without blank spaces.
      character*29 function bp_label(strnum,i,j)
      integer strnum,i,j,k,l
      bp_label = '                             '
      write(bp_label,10) strnum,i,j
 10   format('Structure_',i5,'_',i6,'_',i6)
      k = 11
      l = 12
      do while (l.le.29)
         if (bp_label(k:k).eq.' ') then
            if (bp_label(l:l).eq.' ') then
               l = l + 1
            else
               bp_label(k:k) = bp_label(l:l)
               bp_label(l:l) = ' '
               k = k + 1
               l = l + 1
            endif   
         else
            k = k + 1
            l = l + 1
         endif   
      enddo
      return
      end

c     Menu subroutine for RNA folding program.
c     Allows the user to set energy parameters and to
c     add auxiliary information.
      subroutine menu
 
      include 'rna.inc'
      data listsz/0/
 
10    if (listsz.ge.3000) goto 800
      write(6,900)
50    write(6,901)
      read  (5,*,end=1,err=1) choice
      if (choice.lt.1.or.choice.gt.12) goto 50
      goto (100,200,300,400,400,200,300,350,370,800,60,70),choice
 
60    call listout(6)
      goto 10
 
70    listsz = 0
      goto 10
1     stop
 
 
100   write(6,1000) (eparam(i),i=1,16)
101   write(6,1001)
      read (5,1002,end=10,err=10) parm
      if (parm.lt.1.or.parm.gt.10) goto 10
      write(6,1003)
      read (5,*,end=10,err=10) val
      eparam(parm) = val
      goto 100
1000  format(/,
     .  10x,'   Energy Parameters (10ths kcal/mole)',//,
     .  10x,' 1 Extra stack energy                        [',i5,']',/,
     .  10x,' 2 Extra bulge energy                        [',i5,']',/,
     .  10x,' 3 Extra loop energy (interior)              [',i5,']',/,
     .  10x,' 4 Extra loop energy (hairpin)               [',i5,']',/,
     .  10x,' 5 Extra loop energy (multi)                 [',i5,']',/,
     .  10x,' 6 Multi loop energy/single-stranded base    [',i5,']',/,
     .  10x,' 7 Maximum size of interior loop             [',i5,']',/,
     .  10x,' 8 Maximum lopsidedness of an interior loop  [',i5,']',/,
     .  10x,' 9 Multi loop energy/closing base pair       [',i5,']',/,
     .  10x,'10 Helix penalty per non-GC closing pair     [',i5,']',/,
     .  10x,'11 GGG hairpin bonus                         [',i5,']',/,
     .  10x,'12 poly-C loop penalty: slope                [',i5,']',/,
     .  10x,'13 poly-C loop penalty: intercept            [',i5,']',/,
     .  10x,'14 poly-C loop penalty: hairpin of size 3    [',i5,']',/,
     .  10x,'15 Intermolecular initiation free energy     [',i5,']',/,
     .  10x,'16 GAIL Rule - 1 for version 3 rules only    [',i5,']',//)
1001   format(' Enter Parameter to be changed (<return> for main menu)  ')
1002   format(i6)
1003   format(' Enter new value  ')
 
200   write(6,2001)
      read (5,*,end=10,err=10) i,k
      listsz = listsz + 1
      list(listsz,1) = choice
      list(listsz,2) = i
      list(listsz,3) = k
      list(listsz,4) = -1
      goto 10
2001  format(' Enter base and length  ')
 
300   write(6,3001)
      read (5,*,end=10,err=10) i,j,k
      listsz = listsz + 1
      list(listsz,1) = choice
      list(listsz,2) = i
      list(listsz,3) = j
      list(listsz,4) = k
      goto 10
3001  format(' Enter base pair and length    ')

350   write(6,3501)
      read (5,*,end=10,err=10) i,j,k,l
      listsz = listsz + 1
      list(listsz,1) = -i
      list(listsz,2) = j
      list(listsz,3) = k
      list(listsz,4) = l
      goto 10
3501  format(' Enter the 4 ends defining the region  ')

370   write(6,3701)
      read (5,*,end=10,err=10) i
      listsz = listsz + 1
      list(listsz,1) = choice
      list(listsz,2) = i
      goto 10
3701  format(' Enter the maximum permitted distance between paired bases ')
 
400   write(6,4001)
      read (5,*,end=10,err=10) i,j
      listsz = listsz + 1
      list(listsz,1) = choice
      list(listsz,2) = i
      list(listsz,3) = j
      list(listsz,4) = -1
      goto 10
4001  format(' Enter begining and end    ')
 
800   return
 
900   format(/,
     .  10x,'1  Energy Parameter           7  Double Prohibit ',/,
     .  10x,'2  Single Force               8  Prohibit Range  ',/,
     .  10x,'3  Double Force               9  Maximum distance',/,
     .  10x,'4  Closed Excision           10  Begin folding   ',/,
     .  10x,'5  Open Excision             11  Show current    ',/,
     .  10x,'6  Single Prohibit           12  Clear current   ',/)
901   format(' Enter Choice   ')
      end
 
 
      subroutine listout(u)
c     This subroutine lists current choices on excisions and on
c     forced or prohibited base-pairs.
      include 'rna.inc'
      integer u,i
      character*20 choices(9)
      data choices/'Energy Parameter   ','Single Force       ',
     .'Double Force       ','Closed Excision    ','Open Excision      ',
     .'Single Prohibit    ','Double Prohibit    ','Prohibit range     ',
     .'Maximum distance   '/
 
      if (listsz.eq.0) then
         write(u,*) ' No choices currently defined'
      else
         write(u,*) ' '
         write(u,*) ' Current Choices'
         do 100 i = 1,listsz
             if (list(i,1).eq.3.or.list(i,1).eq.7) then
                write(u,1000) choices(list(i,1)),(list(i,k),k = 2,4)
             else if (list(i,1).lt.0) then
                write(u,1002) choices(-list(i,1)),-list(i,1),(list(i,k),k = 2,4)
             else if (list(i,1).eq.8) then
                write(u,1000) choices(list(i,1)),list(i,2)
             else
                write(u,1001) choices(list(i,1)),(list(i,k),k = 2,3)
             endif
 100      continue
         write(u,*) ' '
      endif
      return
 1000 format(10x,a20,': (',i5,',',i5,') ',i5)
 1001 format(10x,a20,':  ',i5,',',i5)
 1002 format(10x,a20,':  ',i5,'-',i5,'   ',i5,'-',i5)
      end
 
c     Control subroutine for RNA folding.
      subroutine device
      include 'rna.inc'
      character*40 sfile,str
 
c     What kind of run is this? ( regular, save or continuation )
      write(6,2000)
      read (5,2001,end=1) cntrl(1)
      write(6,*) ' '
      if (cntrl(1).lt.0.or.cntrl(1).gt.2) cntrl(1) = 0
 
      if (cntrl(1).eq.1) then
         cntrl(7) = 1
      else
c        What mode is the program to be run in?
c        dot plot, automatic sorted tracebacks of one sequence
c        fragment or suboptimal foldings of every complete
c        sequence in a file.
9        write(6,1002)
         read (5,2001,end=1) cntrl(7)
         if (cntrl(7).lt.0.or.cntrl(7).gt.2) cntrl(7) = 0
         if (cntrl(1).eq.2.and.cntrl(7).eq.2) then
            write(6,*) 'Combination of continuation run and multiple foldings disallowed'
            write(6,*) ' '
            goto 9
         endif
         write(6,*) ' '
      endif
c     Folding multiple sequences is treated as a sort run.
c     Find total number of sequences to be folded in a multiple
c     sequence run.
 
         if (cntrl(7).eq.2) then
             cntrl(5) = 0
             call mseq(cntrl(5))
         endif
 
         if (cntrl(7).eq.0) then
            write(6,1001)
         elseif (cntrl(1).ne.1) then
c           Prompt for controls on sort.
            write(6,1004)
            read (5,2001,end=1) cntrl(8)
            write(6,1003)
         endif
         if (cntrl(1).ne.1) then
            read (5,2001,end=1) cntrl(6)
            if (cntrl(6).lt.1) cntrl(6) = 1
 
            write(6,1005)
            read (5,2001,end=1) cntrl(9)
            write(6,*) ' '
         endif
c        Prompt for SAVE file name for a save/continuation run.
         if (cntrl(1).ne.0) then
4           write(6,3000)
            read (5,3001,end=1) sfile
 5          if (cntrl(1).eq.1) then
                str = 'unknown'
            else
                str = 'old'
            endif
            if (sfile.eq.'         ') sfile= 'fold.sav'
            open(30,err=2,file=sfile,status=str,form='unformatted')
            goto 3
2           if (cntrl(1).eq.2) goto 4
3           if (cntrl(1).eq.2) call getcont
         endif
c        Obtain sequence. Original length is N.
c        A fragment from nsave(1) to nsave(2) is selected.
c        After process, n becomes the length of the processed sequence
c        to be folded.
         if (cntrl(1).ne.2.and.cntrl(7).ne.2) then
            call formid(seqlab,seq,n,maxsiz)
            write(6,4000) seqlab,n
            write(6,*) 'Enter start of fragment (default 1)'
            read (5,4001,end=1) nsave(1)
            if (nsave(1).le.0) nsave(1) = 1
            write(6,4002) n
            read (5,4001,end=1) nsave(2)
            if (nsave(2).le.0) nsave(2) = n
         endif
 
1001  format(/,' Enter minimum vector size for plot (default 1)  ')
1002  format(1x,'Enter run mode',/,5x,'0  Sub-optimal plot (default)',
     .   /,5x,'1  N-best',/,5x,'2  Multiple Molecules')
1003  format(/,' Enter number of tracebacks (default 1)  ')
1004  format(/,' Enter percentage for sort (default 0)  ')
1104  format('A negative number will be interpreted as',
     .          ' -[energy increment] in 10ths of a kcal/mole')
1005  format(/,' Enter window size (default 0)  ')
2000  format(1x,'Enter run type',/,5x,'0  Regular run (default)',
     .   /,5x,'1  Save run',/,5x,'2  Continuation run')
2001  format(i6)
3000  format(' Enter save file name (default fold.sav)')
3001  format(a30)
4000  format(/,' ',a50,5x,i6,' nucleotides',/)
4001  format(i10)
4002  format(1x,'Enter end of fragment (default ',i6,')')
      return
1     stop
      end
 
c     Obtain multiple sequences from a sequence file using MULTID.
      subroutine mseq(i)
      include 'rna.inc'
 
      if (i.eq.0) then
         call multid(seqlab,seq,n,maxsiz,i)
         write(6,*) ' '
      else
         call multid(seqlab,seq,n,maxsiz,i)
         write(6,4000) seqlab,n
      endif
      nsave(1) = 1
      nsave(2) = n
      return
4000  format(/,' ',a50,5x,i6,' nucleotides',/)
      end
 
c     Set up output units and files for RNA folding.
      subroutine outputs
      include 'rna.inc'
      character*1 in
      character*60 str,dstr
      data dstr/'                                                            '/
 
c     Examine sequence label to get default names for output files.
 
      do i = 1,50
         if ((seqlab(i:i).ge.'A'.and.seqlab(i:i).le.'Z').or.
     .       (seqlab(i:i).ge.'a'.and.seqlab(i:i).le.'z').or.
     .       (seqlab(i:i).ge.'0'.and.seqlab(i:i).le.'9')) then
 
            dstr(i:i) = seqlab(i:i)
         else
            dstr(i:i) = '_'
         endif
      enddo
      i = 50
      do while (dstr(i:i).eq.'_'.and.i.gt.1)
         i = i - 1
      enddo
      slen = i
 
c     Line printer output. Get name and open file for write.
      cntrl(2) = 0
      write(6,5010)
      read (5,5000,end=1) in
      if (in.ne.'N'.and.in.ne.'n') then
          cntrl(2) = 1
          write(6,5011)
          read (5,5000,end=1) in
          if (in.eq.'N'.or.in.eq.'n') then
             dstr(slen+1:slen+4) = '.out'
51           write(6,5012) dstr
             read (5,5001) str
             if (str.eq.'         ') str = dstr
             cntrl(4) = 20
             open(20,file=str,status='unknown',err=51)
          else
             cntrl(4) = 6
          endif
      endif
 
c     CT file output. Get name and open file for write.
      write(6,5020)
      read (5,5000,end=1) in
      if (in.eq.'Y'.or.in.eq.'y') then
         cntrl(2) = 2 + 2*cntrl(2)
         dstr(slen+1:slen+4) = '.ct '
52       write(6,5021) dstr
         read (5,5001) str
         if (str.eq.'         ') str = dstr
         open(21,file=str,status='unknown',err=52)
      endif
 
c     Region table output. Get name and open file for write.
      write(6,5030)
      read (5,5000,end=1) in
      if (in.eq.'Y'.or.in.eq.'y') then
         if (cntrl(2).eq.1.or.cntrl(2).eq.2) cntrl(2) = cntrl(2) + 1
         cntrl(2) = cntrl(2) + 3
         dstr(slen+1:slen+4) = '.det'
53       write(6,5031) dstr
         read (5,5001) str
         if (str.eq.'         ') str = dstr
         open(22,file=str,status='unknown',err=53)
      endif
      write(6,*) ' '
      return 
1     stop
 
5000  format(a1)
5001  format(a60)
5010  format('Do you want printer output? (Y,n) ')
5011  format('Output to terminal? (Y,n) ')
5012  format('Enter output file name (default ',a60,')')
5020  format('Do you want ct file? (y,N) ')
5021  format('Enter ct file name (default ',a60,')')
5030  format('Do you want a detailed ddG table? (y,N) ')
5031  format('Enter detailed ddG table file name (default ',a60,')')
 
      end
 
      subroutine cdump
      include 'rna.inc'
      character*40 name
      character yn
 
      write(6,*) 'Enter file name for continuation dump (return for terminal)'
      read (5,100,end=1) name
      if (name.eq.' ') then
         u = 6
      else
         u = 31
         open(31,status='unknown',file=name)
      endif
      call listout(u)
      write(u,101) 'Energy Parameters'
      write(u,1000) eparam
      write(6,*) 'Listing of energy files? (y/N)'
      read(5,102) yn
      if (yn.eq.'Y'.or.yn.eq.'y') then
         call out(u)
      endif
      return
1     stop
100   format(a30)
101   format(a20,/)
102   format(a1)
1000  format(/,
     .  10x,'   Energy Parameters (10ths kcal/mole)',//,
     .  10x,' 1 Extra stack energy                        [',i5,']',/,
     .  10x,' 2 Extra bulge energy                        [',i5,']',/,
     .  10x,' 3 Extra loop energy (interior)              [',i5,']',/,
     .  10x,' 4 Extra loop energy (hairpin)               [',i5,']',/,
     .  10x,' 5 Extra loop energy (multi)                 [',i5,']',/,
     .  10x,' 6 Multi loop energy/single-stranded base    [',i5,']',/,
     .  10x,' 7 Maximum size of interior loop             [',i5,']',/,
     .  10x,' 8 Maximum lopsidedness of an interior loop  [',i5,']',/,
     .  10x,' 9 Multi loop energy/closing base pair       [',i5,']',/,
     .  10x,'10 Helix penalty per non-GC closing pair     [',i5,']',/,
     .  10x,'11 GGG hairpin bonus                         [',i5,']',/,
     .  10x,'12 poly-C loop penalty: slope                [',i5,']',/,
     .  10x,'13 poly-C loop penalty: intercept            [',i5,']',/,
     .  10x,'14 poly-C loop penalty: hairpin of size 3    [',i5,']',/,
     .  10x,'15 Intermolecular initiation free energy     [',i5,']',/,
     .  10x,'16 GAIL Rule - 1 for version 3 rules only    [',i5,']',//)
      end

c     Program ct_compare : RNA folding comparison (scoring)
c
c     Input1 - a CT file containing one structure (i.e. the correct one)
c     Input2 - a CT file containing one or more structures
c     Output - i. The number of helices of the Input1 structure that are
c                 contained in each Input2 structure. A helix is defined
c                 as in Jaeger et al. (PNAS 86, 7706-7710, 1989)
c             ii. The number of base pairs of the Input1 structure that are
c                 contained in each Input2 structure.
c            iii. The first and the overall best structures from Input2 that 
c                 matches the Input1 structure (output file for iii. is ctout,
c                 where ctout is 'Input1 name before the dot'.compare).
c
      implicit integer (a-z)
 
      parameter (maxn=5000)
      character*50 ctlabel1,ctnam1,ctlabel2,ctnam2,ctout
      character*80 ctrec
      dimension basepr1(maxn),basepr2(maxn)
c 
c     Initial setup for run : get first CT file name and open file.
c 
      case = 0
      if (iargc().eq.0) then
         write(6,1010)
 1010    format(' Enter first CT file name : ')
         read(5,1020,end=1000) ctnam1
 1020    format(a50)
      else
         call getarg(1,ctnam1)
      endif
      if (ctnam1.eq.'   ') stop
      open(unit=7,file=ctnam1,status='old',err=1000)
      case = case + 1
      iter = 0
c
c     Open summary output file.
c
      if (case.eq.1) then
         ctout = '                                                  '
         i = 1
         do while (i.lt.45.and.ctnam1(i:i).ne.' '.and.ctnam1(i:i+2).ne.'.ct')
            ctout(i:i) = ctnam1(i:i)
            i = i + 1
         enddo
         ctout(i:i+7) = '.compare'
         open(unit=11,file=ctout,status='unknown')
      endif
c
c     Read first folding.
c
      base_count1 = 0
      read(7,1030,end=999) n1,ctlabel1
 1030 format(i5,a50)
      do i = 1,n1
         basepr1(i) = 0
      enddo
      do i = 1,n1
         read(7,1035,end=998,err=997) ctrec
 1035    format(a80)
         read(ctrec(1:6),*,err=997) k
         read(ctrec(8:80),*,err=997) itmp,itmp,basepr1(i)
         if (i.ne.k) go to 997
         if (basepr1(i).gt.0) then
             basepr1(basepr1(i)) = i
             if (basepr1(i).gt.i) base_count1 = base_count1 + 1
         endif
      enddo
      write (6,1021) ctlabel1,base_count1
 1021 format(' ',a50/' # of base pairs in reference structure is ',i4/)
c 
c     Get second CT file name and open file.
c 
      close (7)
      if (iargc().lt.2) then
         write(6,1012)
 1012    format(' Enter second CT file name (default=fold.ct) : ')
         read(5,1020,end=999) ctnam2
      else
         call getarg(2,ctnam2)
      endif
      if (ctnam2.eq.' ') ctnam2 = 'fold.ct'
      open(unit=7,file=ctnam2,status='old',err=1000)
c
c     Read foldings one by one.
c
      iter = 0
      score_max = 0
      helix_max = 0
      iter_max = 0
      base_count_max = 0
 10   read(7,1030,end=999) n2,ctlabel2
      do i = 1,n2
         basepr2(i) = 0
      enddo
      do i = 1,n2
         read(7,1035,end=998,err=997) ctrec
         read(ctrec(1:6),*,end=998,err=997) k
         read(ctrec(8:80),*,end=998,err=997) itmp,itmp,basepr2(i)
         if (i.ne.k) go to 997
         if (basepr2(i).gt.0) basepr2(basepr2(i)) = i
      enddo
      iter = iter + 1
      if (n1.ne.n2) then
         write(6,1070) n1,n2
 1070    format(' Sequences of unequal length : ',i5,' and ',i5/)
         stop
      endif
      if (iter.eq.1) then
         write (11,1022) ctlabel1,ctlabel2
 1022    format(/a50/'versus'/a50/)
      endif
c
c     Score structures.
c
c     Counts base pairs
      base_count2 = 0
      do i = 1,n1
         if (basepr1(i).gt.i) then
            if (basepr2(i).eq.basepr1(i)) base_count2 = base_count2 + 1
         endif
      enddo
c
      score = 0
      helix_size = 1
      helix_count = 0
      all_but = 0
      i = 1
      do while (i.lt.n1-2)
         if (basepr1(i).gt.i) then
            if (basepr2(i).ne.basepr1(i)) all_but = all_but + 1
            if (basepr1(i+1).eq.basepr1(i)-1) then
               i = i + 1
               helix_size = helix_size + 1
            elseif (basepr1(i+1).eq.basepr1(i)-2) then
               i = i + 1
               helix_size = helix_size + 1
            elseif (basepr1(i+1).eq.basepr1(i)-3) then
               i = i + 1
               helix_size = helix_size + 1
            elseif (basepr1(i+2).eq.basepr1(i)-1) then
               i = i + 2
               helix_size = helix_size + 1
            elseif (basepr1(i+3).eq.basepr1(i)-1) then
               i = i + 3
               helix_size = helix_size + 1
            elseif (basepr1(i+2).eq.basepr1(i)-2) then
               i = i + 2
               helix_size = helix_size + 1
            else
               if (helix_size.gt.2) then
                  helix_count = helix_count + 1  
                  if (all_but.le.2) score = score + 1
               endif
               all_but = 0
               helix_size = 1
               i = i + 1
            endif
         else
             i = i + 1
         endif
      enddo
      if (iter.eq.1) then
          score_1 = score
          iter_1 = iter
          helix_1 = helix_count
          base_count_1 = base_count2
      endif
      if (score.gt.score_max) then
          score_max = score
          iter_max = iter
          helix_max = helix_count
          base_count_max = base_count2
      elseif (score.eq.score_max.and.base_count2.gt.base_count_max) then
          score_max = score
          iter_max = iter
          helix_max = helix_count
          base_count_max = base_count2
      endif
      write (6,1080) iter,score,helix_count,base_count2,base_count1
 1080 format('Structure ',i4,' : ',i4,' out of ',i4,' helices;',i4,
     .       ' out of ',i4,' base pairs.')
c
c     Read next structure
c
      go to 10
c
c     Error in CT file.
c
 997  write(6,9080)
 9080 format(' Error in CT file.')
      stop
c
c     Error - incomplete CT file.
c
 998  write(6,9090)
 9090 format(' Premature end of CT file.')
      stop
c
c     Summary of run : First and best folding identified.
c
 999  write (11,1082) iter_1,iter,score_1,helix_1
      write (11,1083) base_count_1,base_count1
      write (11,*) ' '
      write (11,1082) iter_max,iter,score_max,helix_max
      write (11,1083) base_count_max,base_count1
 1082 format('Structure ',i4,' of ',i4,' : ',i4,' out of ',i3,' helices.')
 1083 format('                         ',i4,' out of ',i4,' base pairs.')
c
 1000 stop
      end

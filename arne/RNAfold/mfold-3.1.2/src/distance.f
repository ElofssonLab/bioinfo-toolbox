c     PROGRAM FP493 - DISTANCE
c
c     COMPUTES DISTANCES BETWEEN PAIRS OF RNA STRUCTURES
c
      character*80 ctfile,outfil
      character*70 label
      character*4 term,termlc
      integer*4 d,bp1(1000),bp2(1000),dm(15,15)
c
      data term/'TERM'/, termlc/'term'/, maxs/15/
c
      ns1 = 1
      nstr = 0
      iout = 6
c
c     PRINT PROGRAM ID AND REQUEST NAME OF INPUT FILE, THEN OPEN FILE
c
      write (6,101)
   10 write (6,105)
      read (5,106) ctfile
      write (6,107) ctfile
      open(unit=2,file=ctfile,status='old',err=10)
c
c     GET INSTRUCTIONS FROM USER FOR COMPARISONS TO BE DONE AND FOR OUTPUT
c
c        INS = 1: COMPARE 1st STRUCTURE WITH EVERY OTHER ONE IN THE FILE
c        INS = 2: COMPARE ALL POSSIBLE PAIRS OF STRUCTURES IN THE FILE
c
      write (6,110)
      read (5,111) ins
      write (6,112)
      read (5,106) outfil
      write (6,113)
      if (outfil.eq.term.or.outfil.eq.termlc) go to 15
c
c     OPEN FILE FOR OUTPUT AND WRITE INITIAL INFORMATION
      iout = 3
      open(unit=3,file=outfil,status='NEW')
      write (iout,101)
      write (iout,107) ctfile
c
c     OPEN FILE FOR PROGRAM VERIFICATION OUTPUT
c     (INDICES OF BASE PAIRS THAT GIVE DISTANCE SOLUTION)
   15 continue
!     OPEN(UNIT=4,FILE='DISTANCE.INDX',STATUS='NEW')
!     WRITE (4,107) CTFILE
c
c     GET BASE PAIR INFORMATION FOR ONE STRUCTURE
c     (DUMMY READ THROUGH STRUCTURES ALREADY COMPARED IF DOING MULTIPLE
c                                                           COMPARISONS)
   20 nstr = nstr + 1
      call bpin(nbases,label,bp1,indeof)
      if (indeof.ne.0) go to 50
      if (nstr.lt.ns1) go to 20
      if (ns1.eq.1) write (iout,115) nstr,nbases,label
c
c     GET BASE PAIR INFORMATION FOR ANOTHER STRUCTURE
c
   30 nstr = nstr + 1
      call bpin(nbs2,label,bp2,indeof)
      if (indeof.ne.0) go to 50
      if (nstr.gt.maxs) go to 92
      if (ns1.eq.1) write (iout,115) nstr,nbs2,label
c
c     COMPUTE DISTANCE BETWEEN A PAIR OF STRUCTURES AND PRINT RESULTS
c     (THIS REQUIRES TWO SEARCH LOOPS, ONE FOR EACH STRUCTURE)
c
      d = 0
      call bpsrch(nbases,bp1,nbs2,bp2,d)
      call bpsrch(nbs2,bp2,nbases,bp1,d)
      dm(nstr,ns1) = d
c
c     GO ON TO NEXT STRUCTURE (IF ANY)
!     WRITE (4,111)
      go to 30
c
c     END OF FILE FOUND - DECIDE WHAT TO DO NEXT
c
   50 nstr = nstr - 1
c
c     CHECK FOR INSUFFICIENT INPUT
      if (nstr.eq.0) go to 90
      if (nstr.eq.1) go to 91
c
c     CHECK FOR COMPLETION OF COMPARISONS
      if (ins.eq.1) go to 60
      if (ns1.eq.(nstr-1)) go to 60
c
c     SET CONTROLS FOR NEXT SET OF COMPARISONS AND REWIND INPUT UNIT
      ns1 = ns1 + 1
      nstr = 0
      rewind 2
      go to 20
c
c     OUTPUT
c
   60 write (iout,120) (j,j = 1,nstr)
      write (iout,111)
      if (ins.eq.2) go to 65
      write (iout,121) (dm(j,1),j = 1,nstr)
      stop
c
   65 do j = 1,nstr
          dm(i,i) = 0
          write (iout,125) j,(dm(j,i),i = 1,j)
          write (iout,111)
      end do
      stop
c
c     ERROR CONDITIONS
c
c     NOT ENOUGH DATA IN INPUT FILE TO DO DIFFERENCE CALCULATION
   90 write (6,190)
      stop
   91 write (6,191)
      stop
c
c     GREATER THAN MAX ALLOWABLE STRUCTURES
c     WRITE MESSAGE (1st ENCOUNTER ONLY) AND PROCEED AS FOR END OF INPUT
   92 if (ns1.eq.1) write (iout,192)
      go to 50
c
  101 format (//' Program DISTANCE computes the distance between pairs ',
     1 'of RNA structures')
  105 format (/' Enter name of input data file: ')
  106 format (a80)
  107 format (/' File: ',a70/)
  110 format (' For distance between 1st structure and every other ',
     1 'one in the file, enter 1' / ' For distance between all ',
     2 'possible pairs of structures in the file,  enter 2'/' ? ')
  111 format (i1)
  112 format (/'  Enter name of output file'/' (enter  term  if you ',
     1 'wish output to terminal): ')
  113 format (/)
  115 format ( 1x,i2,i6,' bases     ',a60)
  120 format (//'     Distances'// 5x,15i5)
  121 format ('  1  ',15i5)
  125 format (1x,i2,2x,15i5)
c
  190 format (/' No data in input file')
  191 format (/' Only one structure in input file')
  192 format (/ ' Too many structures in input file...'/ ' Only the ',
     1 'first 15 will be processed')
      end

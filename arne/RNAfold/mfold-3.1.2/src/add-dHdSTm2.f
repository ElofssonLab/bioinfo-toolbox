c     Read .out, .out.html, .det or .det.html file and a corresponding .dh file
c     For each "dG =" line in the first input file, read the dG
c     Read the corresponding dH in the .dh file (stop for discrepancies)
c     Compute dS and Tm based on simple melt.
c     Add dH, dS  and Tm to "dG = " line
c     Output is same as first input, but with extra information on 
c     the "dG =" lines.
c     This version also reads the total nucleic acid concentration
c     and outputs the correct Tm (still with a 2-state model) using
c     this concentration. Rln2 is added to the entropy for homo-dimers.

      implicit real (a-h,o-z)
      integer amax
      character*3 homodimer
      character*4 mode
      character*60 infile1,infile2,outfile
      character*500 record1,record2
      data amax/500/,rgas/0.0019872/

      if (iargc().lt.5) then
         write(6,*) '2 input files, temperature, mode, nucleic acid 
     .concentration and "YES" or "NO" for homodimer are required'
         call exit(1)
      endif
      call getarg(1,infile1)
      call getarg(2,infile2)
      call getarg(3,record1)
      read(record1,*) t
      t = t + 273.15
      call getarg(4,mode)
      call getarg(5,record1)
      read(record1,*) conc
      call getarg(6,homodimer)
      labstart = index(infile2,'.dh')
      outfile = infile2
      outfile(labstart:labstart+6) = '.dHdSTm'
      open(2,file=infile1,status='old',err=91)
      open(3,file=infile2,status='old',err=92)
      open(4,file=outfile,status='unknown')
      
      read(2,100,end=93) record1
 100  format(a500)
      do while (1.eq.1)
         do while (index(record1,' dG =').eq.0.and.index(record1,' dF =').eq.0)
            call dump(amax,record1)
            read(2,100,end=99) record1
         enddo
         read(record1(index(record1,'=')+1:amax),*) dg
c     Now look for dh in second input file.
         read(3,100,end=94) record2
         do while (index(record2,'Computed').eq.0)
            read(3,100,end=94) record2
         enddo
         read(record2(index(record2,'=')+1:amax),*) dh

c     Now compute. Round dg to 1 decimal.
         idg = nint(10.0*dg)
         dg = float(idg)/10.0
         ds = (dh - dg)/t
         if (homodimer.eq.'yes'.or.homodimer.eq.'YES') then
            factor = 1.0
         else
            factor = 4.0
         endif
         tm = dh/(ds + rgas*alog(conc/factor))
         tm = tm - 273.15
         istart = index(record1,'dG =') + 13
         if (mode.eq.'html') then
            write(record1(istart:amax),105) dh,1000.0*ds,tm
 105        format(' &nbsp; dH = ',f8.1,' &nbsp; dS = ',f8.1,' &nbsp; ',
     .             'T<SUB>m</SUB> = ',f5.1)
         else
            write(record1(istart:amax),107) dh,1000.0*ds,tm
 107        format('  dH = ',f8.1,'  dS = ',f8.1,'  Tm = ',f6.1)
         endif
         call dump(amax,record1)
         record1 = '                                         '
      enddo

 91   stop 'Error opening first input file'
 92   stop 'Error opening second input file'
 93   stop 'First input file is empty'
 94   stop 'Missing information in dh file'

 99   call exit(0)
      end
      subroutine dump(amax,record)
      implicit real (a-h,o-z)
      integer amax
      character*1 record(amax)
      m = amax
      do while (record(m).eq.' '.and.m.gt.1)
         m = m - 1
      enddo
      write(4,100) (record(l),l=1,m)
 100  format(500a1)
      return
      end

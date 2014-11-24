c   Program to get free energies at 37 degrees celsius and enthalpies from
c   input files and calculate the free energies at a different temperature, 
c   which is entered by the user.  User also tells program whether to deal 
c   with RNA or DNA.  
c   example command line:   newtemp rna 50
c                              or
c                           newtemp RNA 50
c--------------------------------------------------------------------------
      program newtemp

      character acid*3,naconc*10,mgconc*10,olipoly*7
      real rnaconc,rmgconc,ion

      character*3 chtemp,gsuffix,hsuffix
      real salt_corr1,salt_corr2
      integer temper
      common /block1/chtemp,gsuffix,hsuffix
      common /block2/salt_corr1,salt_corr2
      common /block3/temper

      call testdate
      call getarg(1,acid)
      call getarg(2,chtemp)
      call getarg(3,naconc)
      call getarg(4,mgconc)
      call getarg(5,olipoly)
         
 900  format(a3)
 905  format(a10)
 910  format(i3)

c   If type of nucleic acid was not entered in command line, ask user for it
      if((index(acid,'DNA').eq.0.and.index(acid,'dna').eq.0).and.
     . (index(acid,'RNA').eq.0.and.index(acid,'rna').eq.0))then
         write(6,*) 'Do you want to use DNA or RNA files? (default = RNA)'
         read(5,900) acid
         if(index(acid,'  ').gt.0)then
            acid = 'RNA'
         endif
      endif

c   If the temperature was not entered in the command line, ask user for it
      if(index(chtemp,'   ').gt.0)then
         write(6,*) 'What is the new temperature? (integer in degrees celsius)'
         read(5,900) chtemp
         if(index(chtemp,'   ').gt.0)then
            write(6,*) 'What is the new temperature? (integer in degrees celsius)'
            read(5,900) chtemp
            if(index(chtemp,'   ').gt.0)then
               write(6,*) 'Temperature must be entered.'
               write(6,*) 'Rerun program and enter correct temperature.'
               call exit(1)
            endif
         endif
      endif

      read(chtemp,910) temper
      
c     If Na+ or Mg++ concentration not entered in the command line, ask for it.
      if (naconc.eq.'          ') then
         write(6,*) 'Enter Na+ concentration (molar).'
         read(5,905) naconc
      endif
      if (mgconc.eq.'          ') then
         write(6,*) 'Enter Mg++ concentration (molar).'
         read(5,905) mgconc
      endif
      read(naconc,*) rnaconc
      read(mgconc,*) rmgconc
c     Test for input that specifies mM.
      if (index(naconc,'mM').gt.0) rnaconc = rnaconc/1000.0
      if (index(mgconc,'mM').gt.0) rmgconc = rmgconc/1000.0
      write(6,810) acid,rnaconc,rmgconc
 810  format(a3,' free energies. [Na+] = ',f6.4,' M [Mg++] = ',f6.4,' M')
      if (olipoly.eq.'       '.or.(olipoly.ne.'polymer'.and.
     .          olipoly.ne.'POLYMER')) then
         olipoly = 'oligo  '
         salt_corr1 = ion(acid,olipoly,'      ',float(temper),rnaconc,rmgconc)
         salt_corr2 = salt_corr1
      else
         salt_corr1 = ion(acid,olipoly,'affine',float(temper),rnaconc,rmgconc)
         salt_corr2 = ion(acid,olipoly,'linear',float(temper),rnaconc,rmgconc)
      endif

c   Determine what the filename suffix should be based on the nucleic acid
      if(index(acid,'RNA').gt.0.or.index(acid,'rna').gt.0)then
         gsuffix = 'dg '
         hsuffix = 'dh '
      else
         gsuffix = 'dgd'
         hsuffix = 'dhd'
      endif

c    Asymmetric interior 1x2 and then 2x3
      call asint1x2
c      call asint2x3  Zuker disables this for now

c    Dangles
      call dangle

c    Loop
      call loop

c   Miscellaneous Loops
      call miscloop

c   Symmetric interior loops of sizes 2, 4, and 6
      call sint2
      call sint4
c      call sint6  Zuker disables this for now

c   Stacking energies
      call stack

c   Tetra-loops
      call tloop

c   Triloops
      call triloop

c   Terminal stacking energies (hairpin loops and internal loops)
      call tstack

      stop
      end
c---------------------------------------------------------------------------
      subroutine asint1x2

      character as1gfile*12,as1hfile*12,path*80,
     +a1header(50)*70,a1tables(24,0:30)*144,rowdown*20,a1outfile*12,row*20
      integer itemp,jtemp,iptemp
      real ra1garray(24,4,24),ra1harray(24,4,24),ras1out(24,4,24)

      character*3 chtemp,gsuffix,hsuffix
      real salt_corr1,salt_corr2
      integer temper
      common /block1/chtemp,gsuffix,hsuffix
      common /block2/salt_corr1,salt_corr2
      common /block3/temper

 900  format(a3)
 920  format(a9)
 930  format(a20)
 940  format(24f6.2)
 950  format(a1)
 960  format(a9)
 970  format(a70)
 980  format(a144)

c   Assigns the filenames of the input files
      as1gfile(1:9) = 'asint1x2.'
      as1hfile(1:9) = 'asint1x2.'
      as1gfile(10:12) = gsuffix
      as1hfile(10:12) = hsuffix

c     Get path for energy files
      call getenv('MFOLDLIB',path)
      in = index(path,' ')
      if (path.eq.'     ') then
         call getenv('MFOLD',path)
         in = index(path,' ')
         path(in:in+4) = '/dat/'
      else
         path(in:in) = '/'
      endif

c   Opens input files
      open(10,file=path(1:index(path,' ')-1)//as1gfile,status='old',err=91)
      open(20,file=path(1:index(path,' ')-1)//as1hfile,status='old',err=91)

c   Reads headers and data from asint1x2.dg(d)
      itemp=0
      do while(index(a1header(itemp),'(U)').eq.0.and.index(a1header(itemp),'(T)').eq.0)
         itemp=itemp+1
         read(10,970)a1header(itemp)
      enddo

c   Reads in header from each table and stores it in a1tables array
      do itemp=1,24
         jtemp=0
         do while(index(a1tables(itemp,jtemp),'<--').eq.0)
            jtemp=jtemp+1
            read(10,980)a1tables(itemp,jtemp)
         enddo

c   Reads in real, numerical, free energy data and stores it in ra1garray
         do jtemp=1,4
            read(10,940)(ra1garray(itemp,jtemp,iptemp),iptemp=1,24)
         enddo
      enddo

c   Reads data from asint1x2.dh(d)

c   Skips main header
      do while(index(row,'(U)').eq.0.and.index(row,'(T)').eq.0)
         read(20,930)row
      enddo

      do itemp=1,24
         
         rowdown = '                    '

c   Skips down to beginning of numerical data
         do while(index(rowdown,'<--').eq.0)
            read(20,930)rowdown
         enddo

c   Reads in real, numerical, enthalpy data and stores it in ra1harray
         do jtemp=1,4
            read(20,940)(ra1harray(itemp,jtemp,iptemp),iptemp=1,24)
         enddo
      enddo

c   Calculates free energy at specified temperature
      do itemp=1,24
         do jtemp=1,4
            do iptemp=1,24
               ras1out(itemp,jtemp,iptemp) = (1./310.15)*(((37.0 -
     .           real(temper))*ra1harray(itemp,jtemp,iptemp)) + ((273.15+
     .           real(temper))*ra1garray(itemp,jtemp,iptemp))) + salt_corr1
     .           + 1.5*salt_corr2
            enddo
         enddo
      enddo

c   Opens output file (asint1x2."temperature")
      a1outfile(1:9) = 'asint1x2.'
      read(chtemp,900)a1outfile(10:12)

      open(30,file=a1outfile,status='unknown')

c   Writes headers and output data to asint1x2."temperature"
      read(chtemp,900)a1header(2)(18:20)
      
c   Writes main header to output file
      itemp=0
      do while(index(a1header(itemp),'(U)').eq.0.and.index(a1header(itemp),'(T)').eq.0)
         itemp=itemp+1
         write(30,970)a1header(itemp)
      enddo
      
      do itemp=1,24
c   Writes header of each table to the output file
         jtemp=0
         do while(index(a1tables(itemp,jtemp),'<--').eq.0)
            jtemp=jtemp+1
            write(30,980)a1tables(itemp,jtemp)
         enddo
c   Writes the calculated, numerical output to the output file
         do jtemp=1,4
            write(30,940)(ras1out(itemp,jtemp,iptemp),iptemp=1,24)
         enddo
      enddo

      write(6,*) 'Output successfully written to ',a1outfile
      
      return

 91   write(6,*) 'Error opening asint1x2.__ input data file(s)'
      call exit(1)
      end

c----------------------------------------------------------------------------
c   Asymmetric interior 2 x 3
      subroutine asint2x3

      character as2gfile*12,as2hfile*12,path*80,
     +a2header(50)*70,a2tables(384,0:30)*120,rowdown*20,a2outfile*12,row*20
      integer itemp,jtemp,iptemp
      real ra2garray(384,4,24),ra2harray(384,4,24),ras2out(384,4,24)

      character*3 chtemp,gsuffix,hsuffix
      real salt_corr1,salt_corr2
      integer temper
      common /block1/chtemp,gsuffix,hsuffix
      common /block2/salt_corr1,salt_corr2
      common /block3/temper

 900  format(a3)
 920  format(a9)
 930  format(a20)
 940  format(24f6.2)
 950  format(a1)
 960  format(a9)
 970  format(a70)
 980  format(a120)

      as2gfile(1:9) = 'asint2x3.'
      as2hfile(1:9) = 'asint2x3.'
      as2gfile(10:12) = gsuffix
      as2hfile(10:12) = hsuffix

c     Get path for energy data 
      call getenv('MFOLDLIB',path)
      in = index(path,' ')
      if (path.eq.'     ') then
         call getenv('MFOLD',path)
         in = index(path,' ')
         path(in:in+4) = '/dat/'
      else
         path(in:in) = '/'
      endif
      open(30,file=path(1:index(path,' ')-1)//as2gfile,status='old',err=92)
      open(40,file=path(1:index(path,' ')-1)//as2hfile,status='old',err=92)

c   Read in data from asint2x3.dg(d) -- free energies
      itemp=0
      do while(index(a2header(itemp),'(U)').eq.0.and.index(a2header(itemp),'(T)').eq.0)
         itemp=itemp+1
         read(30,970)a2header(itemp)
      enddo

      do itemp=1,384
         jtemp=0
         do while(index(a2tables(itemp,jtemp),'<--').eq.0)
            jtemp=jtemp+1
            read(30,980)a2tables(itemp,jtemp)
         enddo

         do jtemp=1,4
            read(30,940)(ra2garray(itemp,jtemp,iptemp),iptemp=1,24)
         enddo
      enddo

c   Read in data from asint2x3.dh(d) -- enthalpies
      do while(index(row,'(U)').eq.0.and.index(row,'(T)').eq.0)
         read(40,930)row
      enddo

      do itemp=1,384
         if(itemp.ge.2)then
            rowdown = '                    '
         endif

         do while(index(rowdown,'<--').eq.0)
            read(40,930)rowdown
         enddo

         do jtemp=1,4
            read(40,940)(ra2harray(itemp,jtemp,iptemp),iptemp=1,24)
         enddo
      enddo

c   Calculates free energy at specified temperature
      do itemp=1,384
         do jtemp=1,4
            do iptemp=1,24
               ras2out(itemp,jtemp,iptemp)=(1./310.15)*(((37.0 -
     .          real(temper))*ra2harray(itemp,jtemp,iptemp))+((273.15 +
     .           real(temper))*ra2garray(itemp,jtemp,iptemp))) + salt_corr1
     .            + 2.5*salt_corr2
            enddo
         enddo
      enddo

      a2outfile(1:9) = 'asint2x3.'
      read(chtemp,900)a2outfile(10:12)

      open(50,file=a2outfile,status='unknown')

      read(chtemp,900)a2header(2)(18:20)
      
c   Writes new free energies to output file (asint2x3."temperature")
      itemp=0
      do while(index(a2header(itemp),'(U)').eq.0.and.index(a2header(itemp),'(T)').eq.0)
         itemp=itemp+1
         write(50,970)a2header(itemp)
      enddo
      
      do itemp=1,384
         jtemp=0
         do while(index(a2tables(itemp,jtemp),'<--').eq.0)
            jtemp=jtemp+1
            write(50,980)a2tables(itemp,jtemp)
         enddo
         do jtemp=1,4
            write(50,940)(ras2out(itemp,jtemp,iptemp),iptemp=1,24)
         enddo
      enddo

      write(6,*) '              ''''               ',a2outfile
      
      close(30)
      close(40)
      close(50)

      return

 92   write(6,*) 'Error opening asint2x3.__ input data file(s)'
      call exit(1)
      end
      
c-----------------------------------------------------------------------------
c   Dangles
      subroutine dangle

      character freein*10,enthpyin*10,path*80,rowheader(8,0:15)*96,rowskip*20,
     .          outfile*10,freeen(8,16)*6,enthalpies(8,16)*6,output(8,16)*6
      integer itemp,jtemp
      real fenum,enthnum,outnum

      character*3 chtemp,gsuffix,hsuffix
      real salt_corr1,salt_corr2
      integer temper
      common /block1/chtemp,gsuffix,hsuffix
      common /block2/salt_corr1,salt_corr2
      common /block3/temper

 1010 format(a7)
 1020 format(a3)
 1030 format(a96)
 1040 format(16a6)
 1050 format(a20)
 1060 format(f6.2)

      freein(1:7) = 'dangle.'
      read(gsuffix(1:3),1020) freein(8:10)
      enthpyin(1:7) = 'dangle.'
      read(hsuffix(1:3),1020) enthpyin(8:10)
      
c     Get path for energy data 
      call getenv('MFOLDLIB',path)
      in = index(path,' ')
      if (path.eq.'     ') then
         call getenv('MFOLD',path)
         in = index(path,' ')
         path(in:in+4) = '/dat/'
      else
         path(in:in) = '/'
      endif
      open(60,file=path(1:index(path,' ')-1)//freein,status='old',err=93)
      open(70,file=path(1:index(path,' ')-1)//enthpyin,status='old',err=93)

c   Reads in data from dangle.dg(d) -- free energies
      do itemp=1,8
         jtemp=0
         do while(index(rowheader(itemp,jtemp),'<--').eq.0)
            jtemp=jtemp+1
            read(60,1030)rowheader(itemp,jtemp)
         enddo

         read(60,1040)(freeen(itemp,jtemp),jtemp=1,16)

      enddo

c   Reads in data from dangle.dh(d) -- enthalpies
      do itemp=1,8
         rowskip = '                    '
         do while(index(rowskip,'3'' <-- 5''').eq.0)
            read(70,1050)rowskip
         enddo

         read(70,1040) (enthalpies(itemp,jtemp),jtemp=1,16)

      enddo

      close(60)
      close(70)

c   Makes sure value read in is number, then calculates output
      do itemp=1,8
         do jtemp=1,16
            if(index(freeen(itemp,jtemp),' . ').eq.0)then
               read(freeen(itemp,jtemp),1060,err=94)fenum
               read(enthalpies(itemp,jtemp),1060,err=95)enthnum
               outnum = (1./310.15)*(((37.0-real(temper))*enthnum) +
     .          ((273.15+real(temper))*fenum)) + 0.5*salt_corr2
               call num2char(outnum,output(itemp,jtemp))
            else
               output(itemp,jtemp) = '   . '
            endif
         enddo
      enddo

      outfile(1:7) = 'dangle.'
      read(chtemp,1020)outfile(8:10)
      open(80,file=outfile,status='unknown')

      do itemp=1,8

         jtemp=0
         do while(index(rowheader(itemp,jtemp),'<--').eq.0)
            jtemp=jtemp+1
            write(80,1030)rowheader(itemp,jtemp)
         enddo

         write(80,1040)(output(itemp,jtemp),jtemp=1,16)

      enddo
      
      close(80)
      write(6,*) '              ''''               ',outfile

      return
      
 93   write(6,*) 'ERROR -- could not open ',freein
      call exit(1)

 94   write(6,*) 'ERROR -- could not open ',enthpyin
      call exit(1)

 95   write(6,*) 'Error reading enthalpies; itemp=',itemp,' jtemp=',jtemp
      write(6,*) 'enthalpies(itemp,jtemp)=',enthalpies(itemp,jtemp)
      call exit(1)
      end
c------------------------------------------------------------------------
c   Subroutine to convert numbers to characters
      subroutine num2char(num,charac)

      character charac*6
      integer place(5),i
      real num

      write(charac,10) num
 10   format(f6.2)
c     Check for overflow
      if(num.le.-100.0.or.num.ge.1000.0) then
         write(6,*) 'ERROR in newtemp in num2char subroutine. Value = ',num
      endif   
      return
      end
c------------------------------------------------------------------------------
      subroutine loop

      character header(20)*70,gfile*8,path*80,outchar(100,3)*6,outfile*8,
     .          fechars(100,3)*6,numrow*70
      integer i,j,max
      real fenums(100,3),outnum

      character*3 chtemp,gsuffix,hsuffix
      real salt_corr1,salt_corr2
      integer temper
      common /block1/chtemp,gsuffix,hsuffix
      common /block2/salt_corr1,salt_corr2
      common /block3/temper

 2020 format(a3)
 2030 format(a70)
 2040 format(14x,a6,12x,a6,11x,a6)
 2050 format(f6.2)
 2060 format(i2,12x,a6,12x,a6,11x,a6)

      gfile(1:5) = 'loop.'
      read(gsuffix(1:3),2020)gfile(6:8)
      
c     Get path for energy data 
      call getenv('MFOLDLIB',path)
      in = index(path,' ')
      if (path.eq.'     ') then
         call getenv('MFOLD',path)
         in = index(path,' ')
         path(in:in+4) = '/dat/'
      else
         path(in:in) = '/'
      endif
      open(80,file=path(1:index(path,' ')-1)//gfile,status='old',err=96)

c   Reads in data from loop.dg(d)
      i=0
      do while(index(header(i),'-----------').eq.0)
         i=i+1
         read(80,2030)header(i)
      enddo
      
c   Sets numrow to '......' so can enter do while loop
      numrow = '.....'

      i=0
      do while(index(numrow,'.').gt.0)
         read(80,2030,end=99)numrow
         if(index(numrow,'.').gt.0)then
            i=i+1
            read(numrow,2040)(fechars(i,j),j=1,3)
            do j=1,3
               if(index(fechars(i,j),' . ').eq.0)then
                  read(fechars(i,j),2050)fenums(i,j)
               else
                  fenums(i,j)=999.9
               endif
            enddo
         endif
      enddo

 99   max=i

c   Calculates the free energy at the specified temp and stores it in outchar
      do i=1,max
         do j=1,3
            if(fenums(i,j).lt.999.9)then
               outnum=(1./310.15)*((273.15+real(temper))*fenums(i,j))
c     Salt correction for bulge and interior loops only up to size 10
c     Note that there is no mechanism for salt correction in the
c     logarithmic extrapolation used by mfold for interior loops of
c     size > 30.
               if (j.lt.3) outnum = outnum + salt_corr1 + 
     .                      0.5*float(min0(i,10))*salt_corr2
               call num2char(outnum,outchar(i,j))
            else
               outchar(i,j) = '   . '
            endif
         enddo
      enddo

c   Writes the output to loop."temper"
      outfile(1:5) = 'loop.'
      read(chtemp,2020)outfile(6:8)
      open(81,file=outfile,status='unknown')

      i=0
      do while(index(header(i),'----------').eq.0)
         i=i+1
         write(81,2030)header(i)
      enddo

      do i=1,max
         write(81,2060)i,(outchar(i,j),j=1,3)
      enddo

      close(80)
      close(81)

      write(6,*) '              ''''               ',outfile
      return
 
 96   write(6,*) 'ERROR opening ',gfile
      call exit(1)
      end
c----------------------------------------------------------------------------
      subroutine miscloop

      character gfile*12,hfile*12,text(12,0:10)*72,temp*72,path*80,outfile*12
      integer i,case,gail,tempi
      real energy(18),enthalpy(18),output(18)

      character*3 chtemp,gsuffix,hsuffix
      real salt_corr1,salt_corr2
      integer temper
      common /block1/chtemp,gsuffix,hsuffix
      common /block2/salt_corr1,salt_corr2
      common /block3/temper
     
 4020 format(a3)
 4030 format(f6.3)
 4040 format(f6.2)
 4050 format(4f6.2)
 4060 format(4x,f6.2,7x,f6.2,12x,f6.2)
 4070 format(a72)
 4080 format(i1)

c   Initialize text and temp
      do i = 1,12
         text(i,0) = '                      '
      enddo
      temp = '                      '

c   Salt correction not needed here.
c   Opens and reads in data from miscloop.dg(d)
      gfile(1:9) = 'miscloop.'
      read(gsuffix(1:3),4020)gfile(10:12)
      hfile(1:9) = 'miscloop.'
      read(hsuffix(1:3),4020)hfile(10:12)

c     Get path for energy data 
      call getenv('MFOLDLIB',path)
      in = index(path,' ')
      if (path.eq.'     ') then
         call getenv('MFOLD',path)
         in = index(path,' ')
         path(in:in+4) = '/dat/'
      else
         path(in:in) = '/'
      endif
      open(85,file=path(1:index(path,' ')-1)//gfile,status='old',err=98)
      open(87,file=path(1:index(path,' ')-1)//hfile,status='old',err=99)

      i=0
      do while(index(text(1,i),'-->').eq.0)
         i=i+1
         read(85,4070) text(1,i)
      enddo
      read(87,4070) temp
      do while(index(temp,'-->').eq.0)
         read(87,4070) temp
      enddo
      read(85,*) energy(1)
      read(87,*) enthalpy(1)

      i=0
      do while(index(text(2,i),'-->').eq.0)
         i=i+1
         read(85,4070) text(2,i)
      enddo
      read(87,4070) temp
      do while(index(temp,'-->').eq.0)
         read(87,4070) temp
      enddo
      read(85,*) energy(2)
      read(87,*) enthalpy(2)

      i=0
      do while(index(text(3,i),'-->').eq.0)
         i=i+1
         read(85,4070) text(3,i)
      enddo
      read(87,4070) temp
      do while(index(temp,'-->').eq.0)
         read(87,4070) temp
      enddo
      read(85,*) (energy(i),i=3,6)
      read(87,*) (enthalpy(i),i=3,6)

      i=0
      do while(index(text(4,i),'-->').eq.0)
         i=i+1
         read(85,4070) text(4,i)
      enddo
      read(87,4070) temp
      do while(index(temp,'-->').eq.0)
         read(87,4070) temp
      enddo
      read(85,*) (energy(i),i=7,9)
      read(87,*) (enthalpy(i),i=7,9)

      i=0
      do while(index(text(5,i),'-->').eq.0)
         i=i+1
         read(85,4070) text(5,i)
      enddo
      read(87,4070) temp
      do while(index(temp,'-->').eq.0)
         read(87,4070) temp
      enddo
      read(85,*) (energy(i),i=10,12)
      read(87,*) (enthalpy(i),i=10,12)

      do case = 6,12
         i=0
         do while(index(text(case,i),'-->').eq.0)
            i=i+1
            read(85,4070) text(case,i)
         enddo
         read(87,4070) temp
         do while(index(temp,'-->').eq.0)
            read(87,4070) temp
         enddo
         if (case.ne.12) then
            read(85,*) energy(7+case)
            read(87,*) enthalpy(7+case)
         else
            read(85,*) gail
            read(87,*) tempi
         endif
      enddo

c   Calculates the output data
      do i=1,18
         output(i) = (1./310.15)*(((37.0-real(temper))*enthalpy(i) +
     .                           (273.15+real(temper))*energy(i)))
      enddo

      outfile(1:9) = 'miscloop.'
      read(chtemp,4020)outfile(10:12)

c   Writes the output data to the output file
      open(86,file=outfile,status='unknown')
      
      i=0
      do while(index(text(1,i),'-->').eq.0)
         i=i+1
         write(86,4070)text(1,i)
      enddo
      write(86,4030)output(1)

      i=0
      do while(index(text(2,i),'-->').eq.0)
         i=i+1
         write(86,4070) text(2,i)
      enddo
      write(86,4040)output(2)

      i=0
      do while(index(text(3,i),'-->').eq.0)
         i=i+1
         write(86,4070)text(3,i)
      enddo
      write(86,4050)(output(i),i=3,6)

      i=0
      do while(index(text(4,i),'-->').eq.0)
         i=i+1
         write(86,4070)text(4,i)
      enddo
      write(86,4060)(output(i),i=7,9)

      i=0
      do while(index(text(5,i),'-->').eq.0)
         i=i+1
         write(86,4070)text(5,i)
      enddo
      write(86,4060)(output(i),i=10,12)

      do case = 6,12
         i=0
         do while(index(text(case,i),'-->').eq.0)
            i=i+1
            write(86,4070)text(case,i)
         enddo
         if (case.ne.12) then
            write(86,4040) output(7+case)
         else
            write(86,4080) gail
         endif
      enddo
      close(85)
      close(86)
      close(87)

      write(6,*) '              ''''               ',outfile
      return

 98   write(6,*) 'ERROR opening ',gfile
      call exit(1)
 99   write(6,*) 'ERROR opening ',hfile
      call exit(1)
      end

c---------------------------------------------------------------------------
      subroutine sint2

      character mainheader(20)*50,temp*50,path*80,tbleheader(6,0:30)*144,
     .          gfile*9,hfile*9,outfile*9
      integer i,j,k
      real fearray(6,4,24),harray(6,4,24),output(6,4,24)

      character*3 chtemp,gsuffix,hsuffix
      real salt_corr1,salt_corr2
      integer temper
      common /block1/chtemp,gsuffix,hsuffix
      common /block2/salt_corr1,salt_corr2
      common /block3/temper

 3010 format(a6)
 3020 format(a3)
 3030 format(a50)
 3040 format(a144)
 3050 format(24f6.2)

      gfile(1:6) = 'sint2.'
      read(gsuffix(1:3),3020)gfile(7:9)
      hfile(1:6) = 'sint2.'
      read(hsuffix(1:3),3020)hfile(7:9)

c     Get path for energy data 
      call getenv('MFOLDLIB',path)
      in = index(path,' ')
      if (path.eq.'     ') then
         call getenv('MFOLD',path)
         in = index(path,' ')
         path(in:in+4) = '/dat/'
      else
         path(in:in) = '/'
      endif
      open(82,file=path(1:index(path,' ')-1)//gfile,status='old',err=98)
      open(83,file=path(1:index(path,' ')-1)//hfile,status='old',err=99)

c   Reads in data from sint2.dg(d)
      i=0
      do while(index(mainheader(i),'(U)').eq.0.and.index(mainheader(i),'(T)').eq.0)
         i=i+1
         read(82,3030)mainheader(i)
      enddo

      read(chtemp,3020)mainheader(2)(18:20)

      do i=1,6
         j=0
         do while(index(tbleheader(i,j),'<--').eq.0)
            j=j+1
            read(82,3040)tbleheader(i,j)
         enddo
         do j=1,4
            read(82,3050)(fearray(i,j,k),k=1,24)
         enddo
      enddo

c   Reads in data from sint2.dh(d)
      do while(index(temp,'(U)').eq.0.and.index(temp,'(T)').eq.0)
         read(83,3030)temp
      enddo

      do i=1,6
         do j=0,9
            temp(1+5*j:5+5*j) = 'BLANK'
         enddo

         do while(index(temp,'<--').eq.0)
            read(83,3030)temp
         enddo
         
         do j=1,4
            read(83,3050)(harray(i,j,k),k=1,24)
         enddo
      enddo

c   Calculates free energy at specified temp
      do i=1,6
         do j=1,4
            do k=1,24
               output(i,j,k)=( (37.0 - real(temper)) * harray(i,j,k) + 
     . (273.15 + real(temper))*fearray(i,j,k))/310.15 + salt_corr1 + salt_corr2
            enddo
         enddo
      enddo

c   Writes output to loop."temperature"
      outfile(1:6) = 'sint2.'
      read(chtemp,3020)outfile(7:9)
      open(84,file=outfile,status='unknown')

      i=0
      do while(index(mainheader(i),'(U)').eq.0.and.index(mainheader(i),'(T)').eq.0)
         i=i+1
         write(84,3030)mainheader(i)
      enddo

      do i=1,6
         j=0
         do while(index(tbleheader(i,j),'<--').eq.0)
            j=j+1
            write(84,3040)tbleheader(i,j)
         enddo
         do j=1,4
            write(84,3050)(output(i,j,k),k=1,24)
         enddo
      enddo

      write(6,*) '              ''''               ',outfile
      return

 98   write(6,*) 'ERROR opening ',gfile
      call exit(1)

 99   write(6,*) 'ERROR opening ',hfile
      call exit(1)
      end
c-------------------------------------------------------------------
      subroutine sint4
      
      character mainheader(60)*96,tableheader(36,0:25)*96,temp*96,
     .path*80,outfile*9,gfile*9,hfile*9
      integer itemp,jtemp,ktemp
      real output(36,16,16),garray(36,16,16),harray(36,16,16)

      character*3 chtemp,gsuffix,hsuffix
      real salt_corr1,salt_corr2
      integer temper
      common /block1/chtemp,gsuffix,hsuffix
      common /block2/salt_corr1,salt_corr2
      common /block3/temper

 5010 format(a6)
 5020 format(a3)
 5030 format(a96)
 5040 format(a96)
 5050 format(16f6.2)

c   Opens sint4.dg(d) and sint4.dh(d) data files
      gfile(1:6) = 'sint4.'
      read(gsuffix(1:3),5020)gfile(7:9)
      hfile(1:6) = 'sint4.'
      read(hsuffix(1:3),5020)hfile(7:9)

c     Get path for energy data 
      call getenv('MFOLDLIB',path)
      in = index(path,' ')
      if (path.eq.'     ') then
         call getenv('MFOLD',path)
         in = index(path,' ')
         path(in:in+4) = '/dat/'
      else
         path(in:in) = '/'
      endif
      open(10,file=path(1:index(path,' ')-1)//gfile,status='old',err=91)
      open(20,file=path(1:index(path,' ')-1)//hfile,status='old',err=92)

c   Reads in data and headers from sint4.dg(d)

c   Stores the mainheader in mainheader character array
      itemp=0
      do while(index(mainheader(itemp),'(UU)').eq.0.and.index(mainheader(itemp),'(TT)').eq.0)
         itemp=itemp+1
         read(10,5030)mainheader(itemp)
      enddo

c   Stores the tableheaders in the tableheader array
      do itemp=1,36
         jtemp=0
         do while(index(tableheader(itemp,jtemp),'<-----').eq.0)
            jtemp=jtemp+1
            read(10,5040)tableheader(itemp,jtemp)
         enddo
      
         do jtemp=1,16
            read(10,5050)(garray(itemp,jtemp,ktemp),ktemp=1,16)
         enddo
      enddo

c   Reads in data from sint4.dh(d)
      do while(index(temp,'(UU)').eq.0.and.index(temp,'(TT)').eq.0)
         read(20,5030)temp
      enddo

      do itemp=1,36
c   Resets the temp array so that next do while loop won't be skipped
         if(itemp.gt.1)then
            do jtemp=0,11
               temp(1+6*jtemp:6+6*jtemp) = 'BLANKS'
            enddo
         endif

c   Skips down to the numerical data
         do while(index(temp,'<-----').eq.0)
            read(20,5030)temp
         enddo

c   Reads in enthalpy data and stores it in harray
         do jtemp=1,16
            read(20,5050)(harray(itemp,jtemp,ktemp),ktemp=1,16)
         enddo
      enddo

c   Calculates output data
      do itemp=1,36
         do jtemp=1,16
            do ktemp=1,16
               output(itemp,jtemp,ktemp)=(1./310.15)*(((37.0-real(temper))*
     .          harray(itemp,jtemp,ktemp)) + ((273.15+real(temper))*
     .          garray(itemp,jtemp,ktemp))) + salt_corr1 + 2.0*salt_corr2
            enddo
         enddo
      enddo

      read(chtemp(1:3),5020)mainheader(2)(18:20)

c   Opens the output file
      outfile(1:6) = 'sint4.'
      read(chtemp,5020)outfile(7:9)
      open(30,file=outfile,status='unknown')

c   Writes the mainheader to the output file
      itemp=0
      do while(index(mainheader(itemp),'(UU)').eq.0.and.index(mainheader(itemp),'(TT)').eq.0)
         itemp=itemp+1
         write(30,5030)mainheader(itemp)
      enddo

c   Writes the data tables to the output file
      do itemp=1,36
         jtemp=0
         do while(index(tableheader(itemp,jtemp),'<-----').eq.0)
            jtemp=jtemp+1
            write(30,5040)tableheader(itemp,jtemp)
         enddo

         do jtemp=1,16
            write(30,5050)(output(itemp,jtemp,ktemp),ktemp=1,16)
         enddo
      enddo

      close(10)
      close(20)
      close(30)

      write(6,*) '              ''''               ',outfile

      return

 91   write(6,*) 'ERROR -- could not open ',gfile
      call exit(1)

 92   write(6,*) 'ERROR -- could not open ',hfile
      call exit(1)
      end

c---------------------------------------------------------------------
      subroutine sint6

      character gfile*9,hfile*9,outfile*9,
     .    path*80,header(0:40)*144,tblheader(1536,0:20)*144,temp*30
      integer itemp,jtemp,ktemp
      real garray(1536,4,24),harray(1536,4,24),output(1536,4,24)

      character*3 chtemp,gsuffix,hsuffix
      real salt_corr1,salt_corr2
      integer temper
      common /block1/chtemp,gsuffix,hsuffix
      common /block2/salt_corr1,salt_corr2
      common /block3/temper

 6010 format(a6)
 6020 format(a3)
 6030 format(a120)
 6040 format(24f6.2)
 6050 format(a30)

      gfile(1:6) = 'sint6.'
      gfile(7:9) = gsuffix
      hfile(1:6) = 'sint6.'
      hfile(7:9) = hsuffix

c     Get path for energy data 
      call getenv('MFOLDLIB',path)
      in = index(path,' ')
      if (path.eq.'     ') then
         call getenv('MFOLD',path)
         in = index(path,' ')
         path(in:in+4) = '/dat/'
      else
         path(in:in) = '/'
      endif
      open(10,file=path(1:index(path,' ')-1)//gfile,status='old',err=91)
      open(20,file=path(1:index(path,' ')-1)//hfile,status='old',err=92)

c   Reads in headers and data from sint6.dg(d)
      itemp=0
      do while(index(header(itemp),'(U)').eq.0.and.index(header(itemp),'(T)').eq.0)
         itemp=itemp+1
         read(10,6030)header(itemp)
      enddo

      do itemp=1,1536
         jtemp=0
         do while(index(tblheader(itemp,jtemp),'<--').eq.0)
            jtemp=jtemp+1
            read(10,6030)tblheader(itemp,jtemp)
         enddo

         do jtemp=1,4
            read(10,6040)(garray(itemp,jtemp,ktemp),ktemp=1,24)
         enddo
      enddo

c   Reads in data from sint6.dh(d)
      do while(index(temp,'(U)').eq.0.and.index(temp,'(T)').eq.0)
         read(20,6050)temp
      enddo

      do itemp=1,1536
         temp = '                              '
         
         do while(index(temp,'<--').eq.0)
            read(20,6050)temp
         enddo

         do jtemp=1,4
            read(20,6040,err=99)(harray(itemp,jtemp,ktemp),ktemp=1,24)
         enddo
      enddo

c   Calculates the free energy at the specified temperature
      do itemp=1,1536
         do jtemp=1,4
            do ktemp=1,24
               output(itemp,jtemp,ktemp) = (1./310.15)*(((37.0-real(temper))*
     .          harray(itemp,jtemp,ktemp)) + ((273.15+real(temper))*
     .          garray(itemp,jtemp,ktemp))) + salt_corr1 + 3.0*salt_corr2
            enddo
         enddo
      enddo

      outfile(1:6) = 'sint6.'
      read(chtemp,6020)outfile(7:9)

      open(30,file=outfile,status='unknown')

      read(chtemp,6020)header(2)(18:20)

      itemp=0
      do while(index(header(itemp),'(U)').eq.0.and.index(header(itemp),'(T)').eq.0)
         itemp=itemp+1
         write(30,6030)header(itemp)
      enddo

      do itemp=1,1536
         jtemp=0

         do while(index(tblheader(itemp,jtemp),'<--').eq.0)
            jtemp=jtemp+1
            write(30,6030)tblheader(itemp,jtemp)
         enddo

         do jtemp=1,4
            write(30,6040)(output(itemp,jtemp,ktemp),ktemp=1,24)
         enddo
      enddo

      close(10)
      close(20)
      close(30)

      write(6,*) '              ''''               ',outfile

      return

 91   write(6,*) 'ERROR -- could not open ',gfile
      call exit(1)

 92   write(6,*) 'ERROR -- could not open ',hfile
      call exit(1)

 99   write(6,*) '** ERROR ** premature end of ',hfile
      write(6,*) 'itemp=',itemp,' jtemp=',jtemp,' ktemp=',ktemp
      call exit(1)
      end
c-----------------------------------------------------------------------
      subroutine stack

      character temp*20,gfile*9,hfile*9,
     .           path*80,mainheader(50)*96,tbleheader(0:4,0:30)*96,outfile*9,
     .           outchar(4,4,16)*6,gcharray(4,4,16)*6,hcharray(4,4,16)*6
      integer itemp,jtemp,ktemp
      real outnum,gnum,hnum

      character*3 chtemp,gsuffix,hsuffix
      real salt_corr1,salt_corr2
      integer temper
      common /block1/chtemp,gsuffix,hsuffix
      common /block2/salt_corr1,salt_corr2
      common /block3/temper

 7020 format(a3)
 7030 format(a96)
 7040 format(16f6.2)
 7060 format(16a6)
 7070 format(a20)
 7080 format(f6.2)

      gfile(1:6) = 'stack.'
      gfile(7:9) = gsuffix
      hfile(1:6) = 'stack.'
      hfile(7:9) = hsuffix

c     Get path for energy data 
      call getenv('MFOLDLIB',path)
      in = index(path,' ')
      if (path.eq.'     ') then
         call getenv('MFOLD',path)
         in = index(path,' ')
         path(in:in+4) = '/dat/'
      else
         path(in:in) = '/'
      endif
      open(10,file=path(1:index(path,' ')-1)//gfile,status='old',err=91)
      open(20,file=path(1:index(path,' ')-1)//hfile,status='old',err=92)

c   Read in headers and then data from stack.dg(d)
      itemp=0
      do while(index(mainheader(itemp),'STACKING').eq.0)
         itemp=itemp+1
         read(10,7030)mainheader(itemp)
      enddo

      do itemp=1,4 
         jtemp=0
         do while(index(tbleheader(itemp,jtemp),'<--').eq.0)
            jtemp=jtemp+1
            read(10,7030)tbleheader(itemp,jtemp)
         enddo

         do jtemp=1,4
            read(10,7060)(gcharray(itemp,jtemp,ktemp),ktemp=1,16)
         enddo
      enddo

c   Reads in headers and then data from stack.dh(d)
      do while(index(temp,'STACKING').eq.0)
         read(20,7070)temp
      enddo

      do itemp=1,4
         do jtemp=0,3
            temp(1+5*jtemp:5+5*jtemp) = 'BLANK'
         enddo

         do while(index(temp,'<--').eq.0)
            read(20,7070)temp
         enddo

         do jtemp=1,4
            read(20,7060)(hcharray(itemp,jtemp,ktemp),ktemp=1,16)
         enddo
      enddo

c   Calculates the free energies at specified temperature
      do itemp=1,4
         do jtemp=1,4
            do ktemp=1,16
               if(index(gcharray(itemp,jtemp,ktemp),' . ').eq.0)then
                  read(gcharray(itemp,jtemp,ktemp),7080)gnum
                  read(hcharray(itemp,jtemp,ktemp),7080)hnum

                  outnum=(1./310.15)*(((37.0-real(temper))*hnum)+
     .             ((273.15+real(temper))*gnum)) + salt_corr1
                  call num2char(outnum,outchar(itemp,jtemp,ktemp))
               else
                  outchar(itemp,jtemp,ktemp) = '   . '
               endif
            enddo
         enddo
      enddo

c   Opens the output file:  stack."temperature"
      outfile(1:6) = 'stack.'
      read(chtemp,7020)outfile(7:9)
      open(30,file=outfile,status='unknown')

c   Writes output to output file
      itemp=0
      do while(index(mainheader(itemp),'STACKING').eq.0)
         itemp=itemp+1
         write(30,7030)mainheader(itemp)
      enddo

      do itemp=1,4
         jtemp=0
         do while(index(tbleheader(itemp,jtemp),'<--').eq.0)
            jtemp=jtemp+1
            write(30,7030)tbleheader(itemp,jtemp)
         enddo

         do jtemp=1,4
            write(30,7060)(outchar(itemp,jtemp,ktemp),ktemp=1,16)
         enddo
      enddo

      close(10)
      close(20)
      close(30)

      write(6,*) '              ''''               ',outfile

      return

 91   write(6,*) 'ERROR -- could not open ',gfile
      call exit(1)

 92   write(6,*) 'ERROR -- could not open ',hfile
      call exit(1)
      end
c-----------------------------------------------------------------
      subroutine tloop

      character temp*15,gfile*9,hfile*9,
     .           path*80,outfile*9,header(10)*15,sequence(200)*8
      integer itemp,counter
      real garray(200),harray(200),output(200)

      character*3 chtemp,gsuffix,hsuffix
      real salt_corr1,salt_corr2
      integer temper
      common /block1/chtemp,gsuffix,hsuffix
      common /block2/salt_corr1,salt_corr2
      common /block3/temper

 1010 format(a6)
 1020 format(a3)
 1030 format(a15)
 1040 format(a8,f6.2)
 1050 format(8x,f6.2)

      gfile(1:6) = 'tloop.'
      gfile(7:9) = gsuffix
      hfile(1:6) = 'tloop.'
      hfile(7:9) = hsuffix

c     Get path for energy data 
      call getenv('MFOLDLIB',path)
      in = index(path,' ')
      if (path.eq.'     ') then
         call getenv('MFOLD',path)
         in = index(path,' ')
         path(in:in+4) = '/dat/'
      else
         path(in:in) = '/'
      endif
      open(10,file=path(1:index(path,' ')-1)//gfile,status='old',err=91)
      open(20,file=path(1:index(path,' ')-1)//hfile,status='old',err=92)

c   Reads in headers and data from tloop.dg(d)
      itemp=0
      do while(index(header(itemp),'-----').eq.0)
         itemp=itemp+1
         read(10,1030)header(itemp)
      enddo

      counter=0
      do while(itemp.ge.0)
         counter=counter+1
         read(10,1040,end=93)sequence(counter),garray(counter)
      enddo

c   Reads in data from tloop.dh(d)
 93   do while(index(temp,'-----').eq.0)
         read(20,1030)temp
      enddo

      itemp=0
      do while(itemp.ge.0)
         itemp=itemp+1
         read(20,1050,end=94)harray(itemp)
      enddo

c   Calculates the free energies at the specified temperature
 94   do itemp=1,counter-1
         output(itemp)=(1./310.15)*(((37.0-real(temper))*harray(itemp)) +
     .    ((273.15+real(temper))*garray(itemp)))
      enddo

      outfile(1:6) = 'tloop.'
      read(chtemp,1020)outfile(7:9)
      open(30,file=outfile,status='unknown')

      itemp=0
      do while(index(header(itemp),'-----').eq.0)
         itemp=itemp+1
         write(30,1030)header(itemp)
      enddo

      do itemp=1,counter-1
         write(30,1040)sequence(itemp),output(itemp)
      enddo

      close(10)
      close(20)
      close(30)
      
      write(6,*) '              ''''               ',outfile

      return

 91   write(6,*) 'ERROR -- could not open ',gfile
      call exit(1)

 92   write(6,*) 'ERROR -- could not open ',hfile
      call exit(1)
      end

c--------------------------------------------------------------------
      subroutine triloop

      character temp*8,gfile*11,hfile*11,outfile*11,
     .          path*80,header(10)*15,sequence(50)*8
      integer itemp,counter
      real garray(50),harray(50),output(50)

      character*3 chtemp,gsuffix,hsuffix
      real salt_corr1,salt_corr2
      integer temper
      common /block1/chtemp,gsuffix,hsuffix
      common /block2/salt_corr1,salt_corr2
      common /block3/temper

 1010 format(a8)
 1020 format(a3)
 1030 format(a15)
 1040 format(a8,f6.2)
 1050 format(8x,f6.2)

      gfile(1:8) = 'triloop.'
      gfile(9:11) = gsuffix
      hfile(1:8) = 'triloop.'
      hfile(9:11) = hsuffix

c     Get path for energy data 
      call getenv('MFOLDLIB',path)
      in = index(path,' ')
      if (path.eq.'     ') then
         call getenv('MFOLD',path)
         in = index(path,' ')
         path(in:in+4) = '/dat/'
      else
         path(in:in) = '/'
      endif
      open(10,file=path(1:index(path,' ')-1)//gfile,status='old',err=91)
      open(20,file=path(1:index(path,' ')-1)//hfile,status='old',err=92)

c   Reads in headings and data from triloop.dg(d)
      itemp=0
      do while(index(header(itemp),'-----------').eq.0)
         itemp=itemp+1
         read(10,1030)header(itemp)
      enddo
      
      counter=0
      do while(itemp.ge.0)
         counter=counter+1
         read(10,1040,end=93)sequence(counter),garray(counter)
      enddo

c   Reads in data from triloop.dh(d)
 93   do while(index(temp,'---').eq.0)
         read(20,1010)temp
      enddo

      itemp=0
      do while(itemp.ge.0)
         itemp=itemp+1
         read(20,1050,end=94)harray(itemp)
      enddo

c   Calculates free energies at specified temperature
 94   do itemp=1,counter-1
         output(itemp)=(1./310.15)*(((37.0-real(temper))*harray(itemp)) +
     .    ((273.15+real(temper))*garray(itemp)))
      enddo

c   Opens output file 
      outfile(1:8) = 'triloop.'
      read(chtemp,1020)outfile(9:11)
      open(30,file=outfile,status='unknown')

c   Writes output to output file
      itemp=0
      do while(index(header(itemp),'-----').eq.0)
         itemp=itemp+1
         write(30,1030)header(itemp)
      enddo

      do itemp=1,counter-1
         write(30,1040)sequence(itemp),output(itemp)
      enddo

      close(10)
      close(20)
      close(30)

      write(6,*) '              ''''               ',outfile

      return

 91   write(6,*) 'ERROR -- could not open ',gfile
      call exit(1)

 92   write(6,*) 'ERROR -- could not open ',hfile
      call exit(1)
      end
c-----------------------------------------------------------------
      subroutine tstack
      
      character tbleheader(0:6,0:20)*96,header(40)*96,gfile*11,hfile*11,
     .          path*80,outfile*11,outchar(4,4,16)*6,name*8,rowtemp*96,
     .          gcharray(4,4,16)*6,hcharray(4,4,16)*6
      integer itemp,jtemp,ktemp,flag
      real outnum,gnum,hnum

      character*3 chtemp,gsuffix,hsuffix
      real salt_corr1,salt_corr2
      integer temper
      common /block1/chtemp,gsuffix,hsuffix
      common /block2/salt_corr1,salt_corr2
      common /block3/temper

 1010 format(a8)
 1020 format(a3)
 1030 format(a96)
 1060 format(a96)
 1080 format(16a6)
 1090 format(f6.2)

      do counter=1,2
   
         if(counter.eq.1)then
            name(1:8) = 'tstackh.'
         else
            name(1:8) = 'tstacki.'
c            do itemp=0,6
c               do jtemp=0,20
c                  tbleheader(itemp,jtemp)(1:25)='                         '
c               enddo
c            enddo
         endif

         read(name,1010) gfile(1:8)
         gfile(9:11) = gsuffix
         read(name,1010) hfile(1:8)
         hfile(9:11) = hsuffix

c     Get path for energy data 
      call getenv('MFOLDLIB',path)
      in = index(path,' ')
      if (path.eq.'     ') then
         call getenv('MFOLD',path)
         in = index(path,' ')
         path(in:in+4) = '/dat/'
      else
         path(in:in) = '/'
      endif
         open(10,file=path(1:index(path,' ')-1)//gfile,status='old',err=91)
         open(20,file=path(1:index(path,' ')-1)//hfile,status='old',err=92)

c   Reads in data from tstacki.dg(d) or tstackh.dg(d)
         itemp=0

         do while(index(header(itemp),'STACKING').eq.0)
            itemp=itemp+1
            read(10,1060) header(itemp)
         enddo

         do itemp=1,4
            jtemp=0
            do while(index(tbleheader(itemp,jtemp),'<--').eq.0)
               jtemp=jtemp+1
               read(10,1030,end=93) tbleheader(itemp,jtemp)
            enddo
         
            do jtemp=1,4
              read(10,1080) (gcharray(itemp,jtemp,ktemp),ktemp=1,16)
            enddo
         enddo

c   Reads in data from tstackh.dh(d) or tstacki,dh(d)
 93      do while(index(rowtemp,'STACKING').eq.0)
            read(20,1030) rowtemp
         enddo

         do itemp=1,4
c   Resets rowtemp character variable
            do jtemp=0,3
               rowtemp = '               '
            enddo

            do while(index(rowtemp,'<--').eq.0)
               read(20,1030) rowtemp
            enddo

            do jtemp=1,4
               read(20,1080) (hcharray(itemp,jtemp,ktemp),ktemp=1,16)
c---------------below slows down program but useful for debugging
c               do ktemp=1,16
c                  if(index(hcharray(itemp,jtemp,ktemp),'.').lt.1)then
c                     write(6,*) '** ERROR ** reading ',hfile
c                     write(6,*) 'itemp=',itemp,' jtemp=',jtemp,' ktemp=',ktemp
c                     call exit(1)
c                  endif
c               enddo
c--------------end debug

            enddo
         enddo

c   Calculates free energies at specified temperature
         do itemp=1,4
            do jtemp=1,4
               do ktemp=1,16
                  if(index(gcharray(itemp,jtemp,ktemp),' . ').eq.0)then
                     read(gcharray(itemp,jtemp,ktemp),1090)gnum
                     read(hcharray(itemp,jtemp,ktemp),1090)hnum
                     outnum=(1./310.15)*(((37.0 - real(temper))*hnum) +
     .                ((273.15 + real(temper))*gnum))
                     call num2char(outnum,outchar(itemp,jtemp,ktemp))
                  else
                     outchar(itemp,jtemp,ktemp) = '   .  '
                  endif
               enddo
            enddo
         enddo

c   Opens output file
         read(name,1010)outfile(1:8)
         read(chtemp,1020)outfile(9:11)
         open(30,file=outfile,status='unknown')

c   Writes output to output file
         itemp=0
         do while(index(header(itemp),'STACKING').eq.0)
            itemp=itemp+1
            write(30,1060)header(itemp)
         enddo

         flag=0
         do itemp=1,4
            jtemp=0
            do while(index(tbleheader(itemp,jtemp),'<--').eq.0.and.flag.le.4)
               jtemp=jtemp+1
               write(30,1030)tbleheader(itemp,jtemp)
            enddo

            do jtemp=1,4
               write(30,1080)(outchar(itemp,jtemp,ktemp),ktemp=1,16)
            enddo
         enddo

         close(10)
         close(20)
         close(30)
         
         write(6,*) '              ''''               ',outfile
      
      enddo

      return

 91   write(6,*) 'ERROR -- could not open ',gfile
      call exit(1)

 92   write(6,*) 'ERROR -- could not open ',hfile
      call exit(1)
      end
c--------------------------------------------------------------------















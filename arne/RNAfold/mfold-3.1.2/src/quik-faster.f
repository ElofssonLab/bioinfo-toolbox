c______________________________________________________________________________
c   mfold 3.0: 
c   © Copyright 1997, Michael Zuker, Washington University in St. Louis
c------------------------------------------------------------------------------
c         - Prediction of RNA/DNA secondary structure by free energy 
c           minimization.
c         - Version 3.0
c         - Michael Zuker and Doug Turner
c
      include 'quik.inc'
      character*30 file_name
      real energy

      call testdate

c     Display information
      write(6,*) 'mfold version 3.1 by Michael Zuker & Doug Turner '

c     Initial setup for run.
5     cntrl(1) = 0
      cntrl(2) = 1
      cntrl(5) = 0
      cntrl(6) = 1
      cntrl(7) = 2
      cntrl(8) = 0

c     Read energy information 
      call enefiles
      call ergread

      mrep = 1
10    read(3,7,end=99) n,(seq(k),k=1,n)
 7    format(i2,1x,255a1)
 
c     Initialize the vst array.
      do i = 1,n
        do j = i,i+n-1
          vst((n-1)*(i-1)+j) = 0
        enddo
      enddo
      do k = 1,n
c           Non-excised bases are examined to determine their type.
c           A - type 1
c           C - type 2
c           G - type 3
c           U/T - type 4
c           anything else - type 5
c           numseq stores nucleotide type.
c           This information may be used to find invalid constraints
c           or to display the folded fragment annotated with constraints.
            numseq(k) = 5
            if (seq(k).eq.'A'.or.seq(k).eq.'a') numseq(k) = 1
            if (seq(k).eq.'C'.or.seq(k).eq.'c') numseq(k) = 2
            if (seq(k).eq.'G'.or.seq(k).eq.'g') numseq(k) = 3
            if (seq(k).eq.'U'.or.seq(k).eq.'u') numseq(k) = 4
            if (seq(k).eq.'T'.or.seq(k).eq.'t') numseq(k) = 4
      enddo
 
c    Fill the optimal energy arrays 
      call fill

c     Output
      write(7,107) mrep,float(vmin)/prec
 107  format(i9,'. dG = ',f8.1)
c
c     Multiple sequence option (CNTRL(7) = 2)
c     If sequence number (MREP) is < total number of sequences
c     (CNTRL(5)), go get another sequence.
c
      if (500*(mrep/500).eq.mrep) write(6,*) mrep
      mrep = mrep + 1
      goto 10
 99   call exit(0)
      end

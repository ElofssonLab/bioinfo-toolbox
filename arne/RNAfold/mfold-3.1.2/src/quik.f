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
      call mseq(cntrl(5))

c     Read energy information 
      call enefiles
      call ergread

      mrep = 1
10    call mseq(mrep)

      call process

      if (n*2.gt.fldmax) then
         tt = fldmax/2
c        Fragment is too long. Try again.
         write(6,*) 'STOP: Segment larger than ',tt
         call exit(1)
      endif

c        Fill the optimal energy arrays 
         call fill

c     Output
         write(6,107) mrep,float(vmin)/prec,seqlab
 107     format(/,i6,'. dG = ',f8.1,3x,a50)
c
c     Multiple sequence option (CNTRL(7) = 2)
c     If sequence number (MREP) is < total number of sequences
c     (CNTRL(5)), go get another sequence.
c
         if (mrep.eq.cntrl(5)) call exit(0)
         mrep = mrep + 1
         goto 10

      end

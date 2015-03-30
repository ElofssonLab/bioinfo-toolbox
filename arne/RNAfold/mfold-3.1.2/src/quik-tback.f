c         - Prediction of RNA/DNA secondary structure by free energy 
c           minimization.
c     This version reads in a file containing many sequences.
c     A simple format is used and multid is bypassed
c     The fill algorithm is executed for only the 1 to n region.
c     A single optimal fold is created.
c     The traceback algorithm has been simplfied a bit to accomplish this.

      include 'quik.inc'
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
 10   read(5,*,end=99) seqlab
      if (seqlab(1:1).eq.">") seqlab(1:1) = ' '
      write(6,*) 'Got here:',seqlab
      read(5,*,end=99) seqrec
      n = index(seqrec,' ') - 1
      if(n.le.0) then
         call exit(0)
      endif
      do i = 1,n
         seq(i) = seqrec(i:i)
      enddo

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
 
c     Fill the optimal energy arrays 
      call fill

c     Compute an optimal folding.

      call trace(error)

c     Output
      write(6,107) mrep,float(vmin)/prec,seqlab
 107  format(i9,'. dG = ',f8.1,' for ',a50)
      write(6,108) n,float(vmin)/prec,seqlab
 108  format(i5,1x,'dG = ',f7.1,4x,a50)
      do i = 1,n
         il = i-1
         iu = i+1
         if (i.eq.n) iu = 0
         write(6,110) i,seq(i),il,iu,basepr(i),i
 110     format(i5,1x,a1,4i6)
      enddo
      if (500*(mrep/500).eq.mrep) write(6,*) mrep
      mrep = mrep + 1
      goto 10
 99   call exit(0)
      end

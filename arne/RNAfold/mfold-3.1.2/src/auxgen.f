c     This program reads a sequence file and a constraint file prepared 
c     by the mfold server 
c     The output is a listing of the sequence annotated with the constraints.

      parameter (nmax=200000)
      implicit integer (a-z)
      integer inc(5,5),seqn(0:nmax),pk_test(100,3)
      character*1 seq(nmax),aux(nmax),record(80),seq_name(40),paren1,paren2
      character*2 choices(9),type
      character*40 file_name,in_name,aux_name,out_name
      real r1,r2,r3,r4
      data choices/'EN','SF','DF','CE','OE','SP','DP','RP','BF'/
      data inc/0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0/
      data numrp/16/

c     'choices' corresponds to the nafold menu
c     Energy modification, single force, double force, closed excision,
c     open excision, single prohibit, double prohibit, begin folding
c     and range prohibit

c     Get file name root from command line
      if (iargc().ne.1) then
         write(6,*) 'STOP: auxgen must be called with an argument'
         call exit(1)
      endif
      file_name = '                                        '
      call getarg(1,file_name)
      in_name = file_name
      aux_name = file_name
      out_name = file_name
      last = index(file_name,' ')
      in_name(last:last+3) = '.seq'
      aux_name(last:last+3) = '.con'
      out_name(last:last+3) = '.pnt'

c     Open files for reading and writing
      open(2,file=in_name,status='old',err=2)
      go to 4
 2    in_name(last:last+3) = '.sqd'
      open(2,file=in_name,status='old',err=3)
      go to 4
 3    in_name(last:last+3) = '    '
      open(2,file=in_name,status='old',err=94)
 4    open(3,file=aux_name,status='old',err=5)
 5    open(4,file=out_name,status='unknown',err=96)

c     Read sequence
      read(2,103,end=97) seq_name
 103  format(40a1)
      do while (seq_name(1).eq.';')
         read(2,105,end=97) seq_name
      enddo
c     'record' contains the sequence name at this point
      write(4,103) seq_name
      n = 0
      do while (1.eq.1)
         read(2,105,end=20) record
 105     format(80a1)
         do i=1,80
            if (record(i).ne.' ') then
               if (record(i).eq.'1') go to 20
               n = n + 1
               seq(n) = record(i)
               if (seq(n).eq.'a'.or.seq(n).eq.'A') then
                  seqn(n) = 1
               else if (seq(n).eq.'b'.or.seq(n).eq.'B') then
                  seqn(n) = 1
               else if (seq(n).eq.'c'.or.seq(n).eq.'C') then
                  seqn(n) = 2
               else if (seq(n).eq.'z'.or.seq(n).eq.'Z') then
                  seqn(n) = 2
               else if (seq(n).eq.'g'.or.seq(n).eq.'G') then
                  seqn(n) = 3
               else if (seq(n).eq.'h'.or.seq(n).eq.'H') then
                  seqn(n) = 3
               else if (seq(n).eq.'t'.or.seq(n).eq.'T') then
                  seqn(n) = 4
               else if (seq(n).eq.'u'.or.seq(n).eq.'U') then
                  seqn(n) = 4
               else if (seq(n).eq.'v'.or.seq(n).eq.'V') then
                  seqn(n) = 4
               else if (seq(n).eq.'w'.or.seq(n).eq.'W') then
                  seqn(n) = 4
               else
                  seqn(n) = 5
               endif
               aux(n) = ' '
            endif
         enddo
      enddo
      seqn(0) = 5
      seqn(n+1) = 0

c     Read constraints and fill in aux array
 20   num_aux = 0
      num_df = 0
      do while (1.eq.1)
         read(3,*,end=30,err=25) l
         num_aux = num_aux + 1
         type = choices(l)
         if (type.eq.'EN') then
            read(3,*,end=30,err=25) dummy
            read(3,*,end=30,err=25) dummy
         endif
         if (type.eq.'BF') then
            read(3,*,end=30,err=25) dummy
         endif
         if (type(1:1).eq.'S') then
            read(3,*,end=30,err=25) i,k
            if (i+k-1.gt.n) then
               write(4,205) type(2:2),i,0,k
 205           format('Invalid constraint: ',a1,3i6)
               write(6,*) 'STOP: Invalid constraint used.'
               call exit(1)
            endif   
            do m=0,k-1
               aux(i+m) = type(2:2)
            enddo
         else if (type(1:1).eq.'D') then
            read(3,*,end=30,err=25) i,j,k
            if (i+k-1.ge.j-k+1) then
               write(4,205) type(2:2),i,j,k
               write(6,*) 'STOP: Invalid constraint used.'
               call exit(1)
            endif   
            if (type(2:2).eq.'F') then
               num_df = num_df + 1
               pk_test(num_df,1) = i
               pk_test(num_df,2) = j
               pk_test(num_df,3) = k
               paren1 = '('
               paren2 = ')'
            else if (type(2:2).eq.'P') then
               paren1 = '{'
               paren2 = '}'
            endif
            do m=0,k-1
               test = inc(seqn(i+m-1),seqn(j-m+1))
     .              + inc(seqn(i+m+1),seqn(j-m-1)) 
               if (inc(seqn(i+m),seqn(j-m)).eq.0.or.test.eq.0) then
                  aux(i+m) = '!'
                  aux(j-m) = '!'
               else
                  aux(i+m) = paren1
                  aux(j-m) = paren2
               endif
            enddo
         else if (type.eq.'RP') then
            read(3,*,end=30,err=25) i,j,k,l
            numrp = numrp + 1
            if (numrp.eq.6) numrp = 7
            if (numrp.eq.16) numrp = 17
            if (numrp.gt.26) numrp = 1
            do m = i,j
               aux(m) = char(96+numrp)
            enddo
            do m = k,l
               aux(m) = char(64+numrp)
            enddo
         endif
      enddo
c     Abort reading constraint file on first error.
 25   continue

c     Test for forced pseudoknot

 30   if (num_df.gt.1) then
         do k1 = 1,num_df-1
            do k2 = k1+1,num_df
               i1 = pk_test(k1,1)
               i2 = pk_test(k1,2)
               i3 = pk_test(k2,1)
               i4 = pk_test(k2,2)
               r1 = float(i1)
               r2 = float(i2)
               r3 = float(i3)
               r4 = float(i4)
               if ((r1-r3)*(r1-r4)*(r2-r3)*(r2-r4).le.0.0) then
                  write(4,*) 'Invalid attempt to force a pseudoknot!'
                  write(4,206) 'F',i1,i2,pk_test(k1,3)
                  write(4,206) 'F',i3,i4,pk_test(k2,3)
 206              format(a1,1x,2i5,1x,i3)
                  write(6,*) 'STOP: Invalid constraint used.'
                  call exit(1)
               endif
            enddo
         enddo
      endif

c     Output
      write(4,104) n
 104  format('#BASES= ',i6)
      nrows = (n-1)/50 + 1
      do j = 1,nrows
         last = min0(n,50*j)
         lastn = 10*((last + 9)/10)
         write(4,100) (i,i=50*(j-1)+10,lastn,10)
 100     format(5(1x,i10))
         write(4,102) (seq(i),i=1+50*(j-1),last)
         if (num_aux.gt.0) write(4,102) (aux(i),i=1+50*(j-1),last)
 102     format(5(1x,10a1))
      enddo

c     Normal termination
      call exit(0)
 94   write(6,*) 'STOP: Cannot open sequence file' 
      call exit(1)
 96   write(6,*) 'STOP: Cannot open file for output'
      call exit(1)
 97   write(6,*) 'STOP: Unexpected end of sequence file'
      call exit(1)
c 98   write(6,*) 'STOP: Error in reading constraint file'
c      call exit(1)
      end

c     This is a simplified redo of John Jaeger's sscount.f
c     A series of RNA/DNA structures in .ct format are read
c     via standard input.
c     It is assumed that all the structures are on the same sequence.
c     Standard output: The first record contains the number of foldings.
c     Record i+1 contains: hstnum(i) j, where j is the number of times 
c     the ith base is single stranded in all the  foldings.

      parameter (maxsiz=20000)
      integer count(maxsiz),nbases,hstnum(maxsiz),check
      character*1 seq(maxsiz)
      character*80 ctrec

      do n = 1,maxsiz
         count(n) = 0
      enddo
c     k counts the nummber of foldings
      k = 0
      do while (1.eq.1)
         k = k + 1
         if (k.eq.1) then
            read(5,1010,end=99) nbases
         else
            read(5,1010,end=98) check
            if (check.ne.nbases) then
               write(6,*) 'STOP: Incorrect header in ct file.'
               call exit(1)
            endif
         endif
 1010    format(i5)
         do i = 1,nbases
            read(5,1015,end=97) ctrec
 1015       format(a80)
            read(ctrec,1020) seq(i)
            read(ctrec(8:80),*) itmp,itmp,j,hstnum(i)
 1020       format(6x,a1)
            if (j.eq.0) count(i) = count(i) + 1
         enddo
      enddo

 99   write(6,*) 'STOP: No structures.'
      call exit(1)
 97   write(6,*) 'STOP: Truncated .ct file.'
      call exit(1)
 98   k = k - 1
      write(6,*) k
      do i = 1,nbases
         write(6,*) hstnum(i),count(i),' ',seq(i)
      enddo
      call exit(0)
      end

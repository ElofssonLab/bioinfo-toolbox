      program scorer
** Compute the number of equal helixes between two foldings,
** using the CT file generated.  Use the definition of
** 'helix' defined in the function CONTIG.  Allow multiple
** CT data sets in each file, and output the data to either
** a file or the screen.
** This routine uses the code found in CTIN.FOR and MYLIB.FOR
** (found below).
**               - John Jaeger 13Nov87
        parameter (maxsiz=5000)
*
        integer pair1(maxsiz),pair2(maxsiz),size1,size2
        integer unit1,unit2,outfil,qlen
        character*1 base(maxsiz)
        character*50 seqlab
        character*80 query,outnam
        logical open1,open2,isold
        real eng
*
        data open1/.false./open2/.false./unit1/1/unit2/2/
        data isold/.false./
        data query/'Output'/qlen/6/
*
** Open CT files, get the first data set.
*
        do while (1.eq.1)
        open1 = .false.
        open2 = .false.
        write (6,5)
        call ctin(unit1,open1,size1,seqlab,base,pair1,eng)
        call ctin(unit2,open2,size2,seqlab,base,pair2,eng)
        if (size1.ne.size2) then
                write (6,1)
                stop
        endif
*
        call getfil (qlen,query,outnam,isold)
        outfil=12
        open (unit=outfil,file=outnam,status='unknown')
*
        write (outfil,3) seqlab
*
        do while (open2)
                call score (outfil,size1,pair1,pair2,eng)
                call ctin(unit2,open2,size2,seqlab,base,pair2,eng)
        enddo
        enddo
*
1       format (1x,'** You''re trying to compare two',
     1' dissimilar sequences.**')
3       format (/,a50)
5       format (1x,'Input Reference CT file, then calculated CT file.',/,
     1        1x,'Calculated file may have several CT files in it.')
        end
      
      subroutine score (outfil,n,bp1,bp2,eng)
*     * See if two foldings have equivalent helixes.  Print out
*     * the helix score (correct # of bp) and distance per helix.
*     * ASSUME that the correct folding is the FIRST one (BP1).
*     
      parameter (maxsiz=5000)
*     
      integer bp1(maxsiz),bp2(maxsiz),n,outfil
*     
      integer ptr,realbp,size,right,notbp,q(maxsiz/4,3)
      integer hstart,ptrval,in,out
      integer dist
      character*80 line
      logical helix,contig,first
      real eng
*     
      save first,line
*     
      common /queue/in,out,q
*     
      data first/.true./
*     
      in=0
      out=0
      helix=.false.
      size=0
      right=0
      notbp=0
*     * Over all bases...
      do ptr=1,n
         ptrval=ptr
         realbp=bp1(ptr)
*     * Check to see if we've looked at this base pair
         if (realbp.gt.ptr) then
*     * Are we looking at a helix already?
            if (helix) then
*     * Check to see if the helix continues..
               if (contig(bp1,ptr)) then
                  notbp=0
                  size=size+1
                  if (realbp.eq.bp2(ptr)) right=right+1
               else
*     * Otherwise, save SIZE and RIGHT and reset variables,
*     *  because we found another helix...
                  call push (right,size,dist(bp1,bp2,hstart,ptr-3))
                  size=1
                  right=0
                  if (realbp.eq.bp2(ptr)) right=1
                  helix=.true.
                  hstart=ptr
               endif
*     * We were not looking at a helix, but we found a new one
            else
               size=size+1
               if (realbp.eq.bp2(ptr)) right=right+1
               helix=.true.
               hstart=ptr
            endif
*     * We may be in the middle of a 1 or 2 base interuption of a helix
         elseif (realbp.eq.0.and.helix) then
            notbp=notbp+1
            if (notbp.eq.3) then
               call push(right,size,dist(bp1,bp2,hstart,ptr-3))
               notbp=0
               size=0
               right=0
               helix=.false.
            endif
         endif
      enddo
      if (helix) then
         call push (right,size,dist(bp1,bp2,hstart,ptrval-3))
      endif
      if (in.eq.0) then
         write (6,3)
      else
         if (first) then
            first=.false.
            call header(outfil,line,in)
         endif
         call dump(outfil,eng)
         write (outfil,1) line
      endif
      return
 1    format (a)
 3    format (1X,'No helixes in reference structure.')
      end
      
      logical function contig(bp,ptr)
*     * See if we have a contigous helix, according to the definition..
*     * "A helix is a continuous set of base pairs, not interupted
*     *  by more than 2 unpaired bases in any form."
*     
*     * Look BACKWARDS to see if this helix is part of the last one.
*     * We need this func. to look at BOTH strands, not just the 5'.
      parameter (maxsiz=5000)
*     
      integer bp(maxsiz),ptr
*     
      logical t1,t2,t3,temp
*     
*     * quick part for a continuous bp
*     
      if (bp(ptr).eq.bp(ptr-1)-1) then
         contig=.true.
      else
         temp=.false.
         if (ptr.gt.2) then
            t1=((bp(ptr).eq.bp(ptr-1)-2).or.
     1           (bp(ptr).eq.bp(ptr-1)-3))
            if (ptr.gt.3) then
               t2=((bp(ptr).eq.bp(ptr-2)-1).or.
     1              (bp(ptr).eq.bp(ptr-2)-2))
               if (ptr.gt.4) then
                  t3=(bp(ptr).eq.bp(ptr-3)-1)
               endif
            endif
         endif
         contig=t1.or.t2.or.t3
      endif
      return
      end
      
      subroutine push (right,size,dist)
*     * Push the parameters of a new helix on the queue.  Make
*     * sure the helix has more than two base pairs (if not, then
*     * discard it).  Assume less than 99 base pairs per helix.
*     
      parameter (maxsiz=5000)
*     
      integer size,right,q(maxsiz/4,3),in,out,dist
*     
      common/queue/in,out,q
*     
      if (size.gt.2) then
         in=in+1
         q(in,1)=right
         q(in,2)=size
         q(in,3)=dist
      endif
      return
      end
      
      subroutine dump(outfil,eng)
*     * Remove the parameters from the queue (FIFO) and print.
      parameter (maxsiz=5000)
*     
      parameter (maxwid=10)
      integer q(maxsiz/4,3),in,out,outfil
      real eng
      integer lines,left,case,ok,i
*     
      common/queue/in,out,q
*     
      ok=0
      lines=in/10
      left=mod(in,10)
      if (left.eq.0) then
         lines=lines-1
         left=10
      endif
      if (lines.ne.0) then
         case=1
         do i=1,lines
            call onelin(outfil,eng,ok,case)
            case=2
         enddo
         case=3
         call onelin(outfil,eng,ok,case)
      else
         case=4
         call onelin(outfil,eng,ok,case)
      endif
      return
      end
      
      
      subroutine onelin (outfil,eng,ok,case)
      parameter (maxsiz=5000)
*     
      integer outfil,ok,case
      real eng
*     
      integer ptr,times,i,in,out,q(maxsiz/4,3)
      character*60 fbuf,dbuf
      character*6 temp
*     
      common/queue/in,out,q
*     
      fbuf=' '
      dbuf=' '
      if (case.lt.3) then
         times=10
      else
         times=in-out
      endif
      ptr=1
*     Make one line of output (really 2 - the helix score and
*     the helix distance).
      do i=1,times
         write(temp,1) q(i+out,1),q(i+out,2)
         fbuf(ptr:ptr+5)=temp
         write(temp,2) q(i+out,3)
         dbuf(ptr:ptr+5)=temp
         ptr=ptr+6
         if (q(i+out,1)+2.ge.q(i+out,2)) ok=ok+1
      enddo
      out=out+10
*     output
      if (case.eq.1) then
         write (outfil,3) eng,fbuf
      elseif (case.eq.2) then
         write (outfil,5) fbuf
      elseif (case.eq.3) then
         write (outfil,7) fbuf,ok,in
      else
         write (outfil,9) eng,fbuf(1:times*6),ok,in
      endif
      write (outfil,11) dbuf
*     
      return
 1    format (i2,'/',i2,' ')
 2    format (i5,1x)
 3    format (f7.1,' : ',a)
 5    format (8x,': ',a)
 7    format (8x,': ',a,': ',i3,'/',i3)
 9    format (f7.1,' : ',a,': ',i2,'/',i2)
 11   format ('DISTANCE: ',A)
      end
      
      subroutine header(outfil,line,in)
      integer outfil,nhyph,in,i,pos,tempin
      character*80 line,outbuf
      character*2 anint
*     
      if (in.gt.10) then
         tempin=10
      else
         tempin=in
      endif
      nhyph=10+tempin*6+7
      line=' '
      do i=1,nhyph
         line(i:i)='-'
      enddo
*     
      outbuf=line
      pos=10+tempin*3-1
      outbuf(pos:pos+4)='Helix'
      pos=10+tempin*6+3
      outbuf(pos:pos+4)='Total'
      write (outfil,1) outbuf
*     
      outbuf=' '
      outbuf(5:6)='dG'
      outbuf(pos:pos+4)='Score'
      do i=1,tempin
         write(anint,3) i
         pos=10+6*(i-1)+3
         outbuf(pos:pos+1)=anint
      enddo
      write (outfil,1) outbuf
*     
      write (outfil,1) line
      return
 1    format (a)
 3    format (i2)
      end
      
      subroutine ctin(myunit,isopen,size,seqlab,base,pairto,eng)
*     * Read in a CT file.  Use the unit MYUNIT as the file so we can
*     * read in multiple CT data sets from one big file.  ISOPEN will
*     * return .FALSE. when we reach the end of the CT data set.
*     * As it is set up now, only read SIZE (# of bases),
*     * SEQLAB (comment), BASE (the sequence) and PAIRTO (the base pairs).
*     * This routine requires the GETFIL subroutine, found in MYLIB.FOR.
*     *               -John Jaeger 28Aug87
*     
      parameter (maxsiz=5000)
      integer pairto(maxsiz),size,myunit
      character*1 base(maxsiz)
      character*50 seqlab
      logical isopen
      real eng
*     
      integer i,tsiz
      logical old
      character*80 filnam,title,ctrec
*     
      data old/.true./title/'CT'/tsiz/2/
*******************************************************
*     
*     * CT file header
*     * CTSUM  - integer - total number of bases, lines of CT data
*     * SEQLAB - character*50 - label that goes with the sequence
*     
*     zuker comments out510   FORMAT(I5,10X,F7.1,4X,A50)
 510  FORMAT(i5,a50)
*     
*******************************************************
*     
*     * CT file data section
*     * SEQCT  - character*1 - base (AGCTU)
*     * K2     - integer - who K5 is basepaired to (0=unpaired)
*     * full format : (I5,1X,A1,3X,4I5)
*     
 520  FORMAT(6X,A1)
*     
*     * if the file isn't open, then get it
*     
      if (.not.isopen) then
         call getfil (tsiz,title,filnam,old)
         open (unit=myunit,file=filnam,status='old')
         isopen=.true.
      endif
*     
*     * Make sure we can fit it in the array
*     
*     zuker comments out       READ (MYUNIT,510,END=700) SIZE,ENG,SEQLAB
      read (myunit,510,end=700) size,seqlab
      eng = 0.0
      if (size.gt.maxsiz) then
         write (6,505)
 505     format (1X,'Make MAXSIZ bigger in the program.')
         stop
      else
*     
*     * Get the stuff
*     
         do 10 i=1,size
            read (myunit,519) ctrec
 519        format(a80)
            read (ctrec,520) base(i)
            read (ctrec(8:80),*) itmp,itmp,pairto(i)
 10      continue
      endif
      goto 800
*     
 700  close(unit=myunit,status='keep')
      isopen=.false.
*     
 800  return
      end
*     * MYLIB.FOR is the library of routines written in fortran.
*     * Name         Type            Comment
*     * GetFil       Sub             Get a file name and see if it's old or new.
*     * Yes          L. Func Query and wait for the answer Y or N.
*     * UpCase       C*1 Func        Take a single char and convert to upper case.
*     *
*     *                    - John Jaeger 28Aug87
*     
      subroutine getfil (tsiz,text,fname,isold)
*     * Get the name of a file and see if it's old or new.
*     
      character*80 fname,query,text
      integer tsiz
      logical isold,nofile,yes,itexists
*     
      nofile=.true.
 103  write (6,1) text(1:tsiz)
      read (5,2) fname
*     zuker adds :
      if (fname.eq.'   ') stop
*     
*     * Try to open the old FName to see if it already exists.
*     * This is to prevent writing over a previous data file by mistake
*     
      inquire (file=fname,exist=itexists)
      if (isold) then
         if (itexists) then
            nofile=.false.
         else
            write (6,5)
         endif
      else
         if (fname.eq.' ')       then
            fname='sysoutput'
            return
         endif
         if (itexists) then
            query='File already exists.  Continue'
            if (yes(31,query)) nofile=.false.
         else
            nofile=.false.
         endif
      endif
      if (nofile) goto 103
*     
 1    FORMAT (1X,A,' file name? ')
 2    FORMAT (A)
 5    FORMAT (1X,'File does not exist.')
      return
      end
 
        logical function yes(qsiz,query)
*
        character*80 query
        integer qsiz,ich
        character*1 inbuf,y,n,upcase
        data y/'Y'/n/'N'/
*
10      write (*,1) query(1:qsiz)
        read (*,2) inbuf
        ich=ichar(upcase(inbuf))
        if (ich.ne.ichar(y).and.ich.ne.ichar(n)) goto 10
        yes=ich.eq.ichar(y)
1       format (1x,a,' (Y or N)? ')
2       format (a)
        return
        end
 
       character*1 function upcase (ch)
** Change Ch to upper case if necessary
       character*1 ch
       character*1 lowera,lowerz,uppera
*
       data lowera/'a'/lowerz/'z'/uppera/'A'/
*
       if (lle(lowera(1:1),ch(1:1)).and.lge(lowerz(1:1),ch(1:1))) then
         upcase=char(ichar(uppera)+ichar(ch)-ichar(lowera))
       else
         upcase=ch
       endif
       return
       end
 
        integer function dist(bp1,bp2,from,to)
** Mike Zuker's formula to show the "distance" between
** two RNA foldings.  BP1 and BP2 have the pairing such
** that I base pairs with BP1(I) in one structure, BP2(I)
** in the other.
** This function needs the routine DPRIME
**                   -John Jaeger 28Aug87
*
        parameter (maxsiz=5000)
*
        integer bp1(maxsiz),bp2(maxsiz),from,to,dprime
*
        dist=max(dprime(bp1,bp2,from,to),dprime(bp2,bp1,from,to))
        return
        end
      
      integer function dprime(bp1,bp2,from,to)
*     * Code courtesy of MXZ to calculate one possible distance
*     * between two RNA foldings.
*     *                 - John Jaeger 28Aug87
*     
      parameter (maxsiz=5000)
*     
      integer bp1(maxsiz),bp2(maxsiz),from,to
*     
      integer d,i1,j1,d1,i2,j2,d2
*     
      n=maxsiz
      d=0
      do 30 i1=from,to
         j1=bp1(i1)
         if (j1.gt.i1.and.bp2(i1).ne.j1) then
            d1=to
            do 10 i2 = from,to
               j2=bp2(i2)
               if (j2.gt.i2) then
                  d2=max0(iabs(i2-i1),iabs(j2-j1))
                  d1=min0(d1,d2)
               endif
 10         continue
            d=max0(d,d1)
         endif
 30   continue
      dprime=d
      return
      end
      



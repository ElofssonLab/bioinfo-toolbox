       program batgen
*
** Batch file generator for RNA.
** Design:  As each subroutine is called, it asks for it's own
**          particular parameter.  Before exiting, it writes out
**          it's answers (data, file name, whatever) to file 20,
**          the batch file.  At the same time, these values are
**          stored in the program, so that if the user runs again,
**          the defaults are the answers from the previous run.
**          Parameters between main, high level and intermediate
**          level routines are passed in COMMON blocks and by name
**          to low and basement level routines (exception:GetRunType)
**
**          Tried to stay as close as possible to variable
**          definitions listed in MFOLD.DOC ( RNA.DOC in version 1)
*
** Synopsis of routines (sub=subroutine, fun=function)
*
** HIGH LEVEL ROUTINES
**  blank Block Data - initial default set up
**  Sub INIT         - Intro!, gets batch file name, header info, and
**                     whether Circular or Linear folding program.
**  Sub GetRunType   - gets run type (Regular, Save, Continuation)
**  Sub RegularRun   - gets run mode (nbest or multiple), and calls
**  Sub Finish       - KP.  Closes files.  Could do more.
** INTERMEDIATE LEVEL ROUTINES
**  Sub SaveRun      - gets save run file name, then the other data
**  Sub ContinuationRun - gets old save file name, then other info
**  Sub OneSeq       - Regular run, one sequence.  Get stuff.
**  Sub MultiSeq     - Regular run, multiple sequence.  Go to it.
** LOW LEVEL ROUTINES
**  Sub SeqInfo      - gets seq (from GetSeq, which is FORMID with
**                     different vars passed back), and begin/end
**                     if necessary.
**  Sub ForceInfo    - Gets force info, EXCEPT for EParam stuff.
**  Sub EFiles       - Gets energy file names.
**  Sub EParams      - EParam menu, separated from Force menu.
**  Sub NBest        - get the 3 nbest parameters: how many % within minimum
**                     energy, max # of structures and distance.
**  Sub OutputFiles  - inquires of 3 output files: Structure,CT and
**                     Details
** THE BASEMENT
**  Sub GetInt       - gets an integer within limits, returns default if
**                     user hits return.  Has default.
**  Fun Convt        - converts a Char*5 string to an integer, if possible
**  Fun YesSir       - returns .TRUE. if user types Y.  Has default.
**  Fun UpCase       - "uppercases" character sent to it.
**  Sub WriteStr     - writes a string to desired unit, leaves out
**                     trailing blanks
**  Sub GetStr       - gets a string from user.  Has default.
** NEEDS THIS OTHER ROUTINE (not in this file).
**  Sub GetSeq       - altered version of FORMID to return file name and
**                     seq number.
**
**                  -John Jaeger Dec88
**                  -modified for version 2 by Michael Zuker   Nov89
**                  -modified for Personal Iris by Michael Zuker Feb90
**                  -no header file in UNIX version
*
** Files:
** Unit   Description
** 20     Batch file - where all this stuff goes
*
** Main program: Figures out RunType (Regular, Save or Continuation)
**               and calls proper routine.
*
      integer inunit,outunit,runtype
      parameter (inunit=5,outunit=6)
*
      integer cntrl(10),eparam(9),list(500,4),listsz,seqlimits(3)
      character*80 batch,sequence,struc,ct,details,savefile
      character*80 asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,
     + sint6,stack,tloop,triloop,tstckh,tstcki
*
      common/ints/cntrl,eparam,list,listsz,seqlimits 
      common/fnams/batch,sequence,struc,ct,details,
     1       savefile
      common/enfiles/asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,
     + sint6,stack,tloop,triloop,tstckh,tstcki
*
** Open batch file, put in header info, get linear or circular run
*
      call init
        runtype=cntrl(1)
        call getruntype(runtype)
        cntrl(1)=runtype
        if (runtype.eq.0) then
          call regularrun
          write (20,2)
        elseif (runtype.eq.1) then
          call saverun
          write (20,2)
        elseif (runtype.eq.2) then
          call continuationrun
        endif
      call finish
2     format ('10')
5     format (a7/1x)
      end
 
      block data
*
** Initializes variables in COMMON blocks.  Change these to suit
** your own tastes.
*
      integer cntrl(10),eparam(9),list(500,4),listsz,seqlimits(3)
      character*80 batch,sequence,struc,ct,details,savefile
      character*80 asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,
     +sint6,stack,tloop,triloop,tstckh,tstcki
*
      common/ints/cntrl,eparam,list,listsz,seqlimits 
      common/fnams/batch,sequence,struc,ct,details,
     1       savefile
      common/enfiles/asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,
     + sint6,stack,tloop,triloop,tstckh,tstcki
*
      data cntrl/1,1,80,0,0,1,1,0,0,0/eparam/0,0,0,0,46,4,30,30,1/
      data listsz/0/seqlimits/1,1,999/
      data batch/'fold.com'/
      data asint1x2/'asint1x2.dat'/asint2x3/'asint2x3.dat'/dangle/'dangle.dat'/
      data loop/'loop.dat'/miscloop/'miscloop.dat'/sint2/'sint2.dat'/
      data sint4/'sint4.dat'/sint6/'sint6.dat'/stack/'stack.dat'/
      data tloop/'tloop.dat'/triloop/'triloop.dat'/tstckh/'tstackh.dat'/
      data tstcki/'tstacki.dat'/
      data sequence/'fold.seq'/struc/'fold.out'/ct/'fold.ct'/
      data details/'fold.det'/savefile/'fold.sav'/
*
      end
 
*    High level routines - called by main
 
      subroutine init
*
** Greet user.  Get batch file and open as unit 20.  
*
      integer inunit,outunit
      parameter (inunit=5,outunit=6)
*
      integer cntrl(10),eparam(9),list(500,4),listsz,seqlimits(3)
      character*80 batch,sequence,struc,ct,details,savefile
      character*80 asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,sint6,
     +stack,tloop,triloop,tstckh,tstcki
*
      common/ints/cntrl,eparam,list,listsz,seqlimits 
      common/fnams/batch,sequence,struc,ct,details,
     1       savefile
      common/enfiles/asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,
     + sint6,stack,tloop,triloop,tstckh,tstcki
*
** Intro.
*
      write (outunit,100)
      write (outunit,101)
      write (outunit,104)
      write (outunit,105)
*
** Get Batch file name
*
20    write (outunit,3) batch
      call getstr(batch)
      open (unit=20,status='unknown',err=20,file=batch)
c     write (20,*) '   ' Zuker removes this for version 3.0

3     format (1x,'Batch file name? Default:',a40)
100   format (3x,'This program generates batch files for the nucleic acid')
101   format (1x,'folding program naview (version 3.0).')
104   format (3x,'At each prompt, you may type in your own response,')
105   format (1x,'or hit return to use the default.',/)
      end
 
      subroutine getruntype(rt)
*
** Gets Run type (regular, save or continuation), and returns
** the right value (Reg=0,Save=1,Cont=2).
*
      integer inunit,outunit
      parameter (inunit=5,outunit=6)
      integer rt
      write (outunit,1)
      call getint(rt,0,2,rt)
      write (20,2) rt
      return
1     format (1x,'Enter run type',/,
     1        5x,'0  Regular run',/,
     2        5x,'1  Save run',/,
     3        5x,'2  Continuation run')
2     format (i1)
      end
 
      subroutine regularrun
*
** User wants a regular run, but there are two modes: single
** and multiple (can't run plots from batch).  Get mode type,
** (nbest=1, multiple=2) write to file, then call appropriate
** subroutine.
*
      integer inunit,outunit
      parameter (inunit=5,outunit=6)
*
      integer cntrl(10),eparam(9),list(500,4),listsz,seqlimits(3)
      character*80 batch,sequence,struc,ct,details,savefile
      character*80 asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,sint6,
     +stack,tloop,triloop,tstckh,tstcki
*
      common/ints/cntrl,eparam,list,listsz,seqlimits 
      common/fnams/batch,sequence,struc,ct,details,
     1       savefile
      common/enfiles/asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,
     + sint6,stack,tloop,triloop,tstckh,tstcki
*
** Get Run Mode a.k.a. Cntrl(7).  Either 1 (N best) or 2 (Multiple)
*
      write (outunit,1)
      call getint(cntrl(7),1,2,cntrl(7))
      write (20,2) cntrl(7)
*
      if (cntrl(7).eq.1) then
        call oneseq
      else if (cntrl(7).eq.2) then
        call multiseq
      endif

      return
1     format (1x,'Enter run mode',/,
     1        5x,'1   N best',/,
     2        5x,'2   Multiple molecules')
2     format (i1)
      end
 
      subroutine finish
*
** Close shop
*
      close (unit=10,status='keep')
      close (unit=20,status='keep')
      return
      end
 
* Medium level routines
      subroutine saverun
*
** User wants a save run.  Ask for save file name, write it out,
** then call appropriate subroutines to get the rest of the data.
*
      integer inunit,outunit
      parameter (inunit=5,outunit=6)
*
      integer cntrl(10),eparam(9),list(500,4),listsz,seqlimits(3)
      character*80 batch,sequence,struc,ct,details,savefile
      character*80 asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,sint6,
     +stack,tloop,triloop,tstckh,tstcki
*
      common/ints/cntrl,eparam,list,listsz,seqlimits 
      common/fnams/batch,sequence,struc,ct,details,
     1       savefile
      common/enfiles/asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,
     + sint6,stack,tloop,triloop,tstckh,tstcki
*
**  Get the Save file name
*
      write (outunit,1) savefile
      call getstr(savefile)
      call writestr(20,savefile)
*
** Ask for all the right info
*
      call seqinfo(sequence,seqlimits,cntrl)
      call efiles
      call forceinfo(seqlimits,list,listsz)
      call eparams(eparam)
      return
1     format (1x,'What is the save file name? Default:',a40)
      end
 
      subroutine continuationrun
*
** User wants a continuation run.  Get NBest parameters, the
** old save file name, put the continuation dump in the log file,
** then ask for the other output files.
*
      integer inunit,outunit
      parameter (inunit=5,outunit=6)
      logical itexists,yessir
      character*1 ch
*
      integer cntrl(10),eparam(9),list(500,4),listsz,seqlimits(3)
      character*80 batch,sequence,struc,ct,details,savefile
      character*80 asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,sint6,
     +stack,tloop,triloop,tstckh,tstcki
*
      common/ints/cntrl,eparam,list,listsz,seqlimits 
      common/fnams/batch,sequence,struc,ct,details,
     1       savefile
      common/enfiles/asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,
     + sint6,stack,tloop,triloop,tstckh,tstcki
*
** Since we're writing a batch file, assume "n best" mode
** for the continuation run
*
      write (20,3)
      call nbest(cntrl)
*
** Get the save file name, and see if it exists.
*
10    write (outunit,1) savefile
      call getstr(savefile)
      inquire(file=savefile,exist=itexists)
      if (.not.itexists) then
        write (outunit,*) 'That save file doesn''t exist.'
        write (outunit,*) 'Continue anyway?'
        ch='N'
        if (.not.yessir(ch)) goto 10
      endif
      call writestr(20,savefile)
*
      write (20,2)
      call outputfiles(struc,ct,details,cntrl)
      return
1     format (1x,'What is the save file name? Default:',a40)
2     format (/,'N')
3     format ('1')
      end
 
      subroutine oneseq
*
** User wants a regular run, one sequence.
*
      integer cntrl(10),eparam(9),list(500,4),listsz,seqlimits(3)
      character*80 batch,sequence,struc,ct,details,savefile
      character*80 asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,sint6,
     +stack,tloop,triloop,tstckh,tstcki
*
      common/ints/cntrl,eparam,list,listsz,seqlimits 
      common/fnams/batch,sequence,struc,ct,details,
     1       savefile
      common/enfiles/asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,
     + sint6,stack,tloop,triloop,tstckh,tstcki
*
      call nbest(cntrl)
      call seqinfo(sequence,seqlimits,cntrl)
      call efiles
      call outputfiles(struc,ct,details,cntrl)
      call forceinfo(seqlimits,list,listsz)
      call eparams(eparam)
      return
      end
 
      subroutine multiseq
*
** Regular run, many sequences.
*
      integer cntrl(10),eparam(9),list(500,4),listsz,seqlimits(3)
      character*80 batch,sequence,struc,ct,details,savefile
      character*80 asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,sint6,
     +stack,tloop,triloop,tstckh,tstcki
*
      common/ints/cntrl,eparam,list,listsz,seqlimits 
      common/fnams/batch,sequence,struc,ct,details,
     1       savefile
      common/enfiles/asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,
     + sint6,stack,tloop,triloop,tstckh,tstcki
*
      call seqinfo(sequence,seqlimits,cntrl)
      call nbest(cntrl)
      call efiles

      call outputfiles(struc,ct,details,cntrl)
      call eparams(eparam)
      return
      end
 
*     Low level routines
      subroutine seqinfo(seq,seqlimits,cntrl)
*
** Get sequence file name, sequence #, and folding limits (5' to 3').
** Use a slightly altered version of FORMID to do most of the dirty
** work.
*
      integer inunit,outunit,ptr
      parameter (inunit=5,outunit=6)
*
      character*80 seq,seqnumbuf
      integer seqlimits(3),cntrl(10),i,temp
*
      if (cntrl(1).eq.0.and.cntrl(7).eq.2) then
*
**   Just get seq, because it's a multi-seq run
*
        write (outunit,3) seq
        call getstr(seq)
        call writestr(20,seq)
      else
*
** Get seq, number and limits.
*
        call getseq(seq,seqlimits(3),seqlimits(1))
        write (outunit,1)
        call getint(1,1,seqlimits(3),seqlimits(2))
        write (outunit,2)
        temp=seqlimits(3)
         call getint(temp,seqlimits(2),temp,seqlimits(3))
        call writestr(20,seq)
*
** Have to encode and WriteStr the seq #, since FORMID is picky about
** this sort of thing.
*
        seqnumbuf=' '
        write (seqnumbuf,5) seqlimits(1)
        ptr=1
        do while (seqnumbuf(ptr:ptr).eq.' ')
                ptr=ptr+1
        enddo
        seqnumbuf(1:80-ptr+1)=seqnumbuf(ptr:80)
        call writestr(20,seqnumbuf)
                write (20,6) (seqlimits(i),i=2,3)
      endif
      return
1     format (1x,'5'' end?')
2     format (1x,'3'' end?')
3     format (1x,'Sequence file? Default:',a40)
5     format (i5)
6     format (i5,/,i5)
      end
 
      subroutine forceinfo(seqlimits,list,listsz)
*
** Get the information found in the FORCE menu - except option 1
** (change energy parameters), which is a separate routine.  Menu
** is slighly changed from RNA - check format 2.
*
      integer inunit,outunit
      parameter (inunit=5,outunit=6)
*
      integer seqlimits(3),list(500,4),listsz,i,j,k
      logical yessir
      character*1 ch
      ch='N'
      i=10
*
** Ask user if they have any info.  If they don't just go.
*
      write (outunit,1)
      if (yessir(ch)) then
10              i=10
        write(outunit,2)
         call getint(i,2,12,i)
*
** it's ugly but it works
*
        if (i.eq.10) then
*
** write info to batch file, exit.
*
                do i=1,listsz
                        if (list(i,1).eq.3.or.list(i,1).eq.7) then
                                write (20,12) (list(i,j),j=1,4)
                        else
                                write (20,12) (list(i,j),j=1,3)
                        endif
                        do j=1,4
                                list(i,j)=0
                        enddo
                enddo
                listsz=0
                return
        elseif (i.eq.3.or.i.eq.7) then
*
** Double force, double prohibit
*
                listsz=listsz+1
                list(listsz,1)=i
30              write (outunit,5)
                read (inunit,*,err=10) i,j,k
                if (i.lt.seqlimits(2).or.i.gt.j.or.
     1          j.lt.i.or.i+k.gt.j) then
                        write (outunit,4) seqlimits(2),seqlimits(3)
                        goto 30
                endif
                list(listsz,2)=i
                list(listsz,3)=j
                list(listsz,4)=k
        else if (i.lt.10) then
* note that 3,7 were taken care of above...
                listsz=listsz+1
                list(listsz,1)=i
20              if (i.eq.4.or.i.eq.5) then
*
** Excision info
*
                        write (outunit,7)
                        read (inunit,*,err=10) i,j
                        if (i.lt.seqlimits(2).or.i.gt.seqlimits(3)) then
                                write (outunit,4) seqlimits(2),seqlimits(3)
                                goto 20
                        endif
                elseif (i.eq.2.or.i.eq.6) then
*
** Single force, single prohibit.
*
                        write (outunit,3)
                        read (inunit,*,err=10) i,j
                        if (i.lt.seqlimits(2).or.i.gt.seqlimits(3).or.
     1             i+j.gt.seqlimits(3)) then
                                write (outunit,4) seqlimits(2),seqlimits(3)
                                goto 20
                        endif
                endif
                list(listsz,2)=i
                list(listsz,3)=j
                list(listsz,4)=0
                if (i.eq.8) then
c                  Prohibit Range
                   write (outunit,33)
                   read (inunit,*,err=10) i,j,k,l
                   list(listsz,1) = -i
                   list(listsz,2) =  j
                   list(listsz,3) =  k
                   list(listsz,4) =  l
                 else if (i.eq.9) then
c                     Maximum Distance between base pairs
                      write (outunit,44)
                      read (inunit,*,err=10) i
                      list(listsz,2) =  i
                 endif

       else if (i.eq.11) then
*
** Output info to the user
*
                do  j=1,listsz
                   if (list(j,1).eq.2) then
                      write (outunit,*) 'Single Force    ',list(j,2),list(j,3)
                   elseif (list(j,1).eq.3) then
                      write (outunit,*) 'Double Force    ',list(j,2),list(j,3),
     .                              list(j,4)
                   elseif (list(j,1).eq.4) then
                      write (outunit,*) 'Closed Excision ',list(j,2),list(j,3)
                   elseif (list(j,1).eq.5) then
                      write (outunit,*) 'Open Excision   ',list(j,2),list(j,3)
                   elseif (list(j,1).eq.6) then
                      write (outunit,*) 'Single Prohibit ',list(j,2),list(j,3)
                   elseif (list(j,1).eq.7) then
                      write (outunit,*) 'Double Prohibit ',list(j,2),list(j,3),
     .                              list(j,4)
                   elseif (list(j,1).eq.8) then
                      write (outunit,*) 'Prohibit Range  ',-list(j,1),'-',
     .                     list(j,1),' ',list(j,3),'-',list(j,4)
                   elseif (list(j,1).eq.9) then
                      write (outunit,*) 'Maximum Distance',list(j,2)
                   endif
              enddo
        else
*
** clear the arrray
*
                do 99 i=1,listsz
                        list(i,1)=0
                        list(i,2)=0
                        list(i,3)=0
99              continue
                listsz=0
        endif
        goto 10
      endif
1     format (1x,'Do you have Force or Excision information?')
2     format(/,
     .  10x,'2  Single Force        7 Double Prohibit ',/,
     .  10x,'3  Double Force        8 Prohibit Range  ',/,
     .  10x,'4  Closed Excision     9 Maximum Distance',/,
     .  10x,'5  Open Excision      10 Exit Force Info ',/,
     .  10x,'6  Single Prohibit    11 Show current    ',/,
     .  32x,'12 Clear Current'/1x,'Choice?')
 3    format (1x,'Starting location and length?')
 33   format (1x,'Bounds? i1-i2 i3-i4')
 4    format (1x,'Out of bounds.  Range:',i5,' to ',i5)
 44   format (1x,'Enter maximum distance between base pairs: ')
 5    format (1x,'5'' base, 3'' base and length?')
 6    format (1x,a15,'  :',i5,2x,i5,2x,i5)
 7    format (1x,'5'' base and 3'' base?')
 12   format (i1,/,i5,2(1x,i5))
      end
 
      
      subroutine efiles
*
** Get the energy file info.  Ask the user if they want default
** files.  If so, then write out defaults.  Otherwise get new
** energy file names (WHICH BECOME THE NEW DEFAULTS!!) and write
** them out.
*
      integer inunit,outunit
      parameter (inunit=5,outunit=6)
      character*80 asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,sint6,
     +stack,tloop,triloop,tstckh,tstcki,path
      logical yessir
      character*1 ch
     
      common/enfiles/asint1x2,asint2x3,dangle,loop,miscloop,sint2,sint4,
     + sint6,stack,tloop,triloop,tstckh,tstcki
*
      ch='Y'
      write (outunit,1)
      if (.not.yessir(ch)) then

c----------begin modifications-----------------------------
        write(outunit,3)
        write(outunit,2) 'AsInt1x2',asint1x2
        call getstr(asint1x2)
        write(outunit,2) 'AsInt2x3',asint2x3
        call getstr(asint2x3)
        write(outunit,2) 'Dangle',dangle
        call getstr(dangle)
        write(outunit,2) ' Loop',loop
        call  getstr(loop)
 50     write(outunit,2) ' Miscloop',miscloop
        call getstr(miscloop)
        write(outunit,2) 'SyInt2',sint2
        call getstr(sint2)
        write(outunit,2) 'SyInt4',sint4
        call getstr(sint4)
        write(outunit,2) 'SyInt6',sint6
        call getstr(sint6)
        write(outunit,2) 'Stack',stack
        call getstr(stack)
        write(outunit,2) ' TLoop',tloop
        call getstr(tloop)
        write(outunit,2) 'TriLoop',triloop
        call getstr(triloop)
        write(outunit,2) 'TStckH',tstckh
        call getstr(tstckh)
        write(outunit,2) 'TStckI',tstcki
        call getstr(tstcki)

      endif

      call writestr(20,asint1x2)
      call writestr(20,asint2x3)
      call writestr(20,dangle)
      call writestr(20,loop)
      call writestr(20,miscloop)
      call writestr(20,sint2)
      call writestr(20,sint4)
      call writestr(20,sint6)
      call writestr(20,stack)
      call writestr(20,tloop)
      call writestr(20,triloop)
      call writestr(20,tstckh)
      call writestr(20,tstcki)
c----------------end----------------------------
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

      open(32,file=miscloop,status='old',err=55)
      return
55    open(32,file=path(1:index(path,' ')-1)//miscloop,status='old',err=60)
      return
60    write(6,4) miscloop
      go to 50
1     format (1x,'Use default energy files?')
2     format (1x,a6,' file name? Default:',a40)
3     format (/,1x,'The files you type in here become',
     1         ' the new defaults.')
4     format (' ',a80/' does not exist. Try again.'/)
      end
 
      subroutine eparams(eparam)
*
** Get the energy parameters.  Usually this is part of
** the Force/Excision menu, but I separated them to
** speed things up a bit.
*
      integer inunit,outunit
      parameter (inunit=5,outunit=6)
      integer eparam(9),i,j,oldep(10),start,fin,convt
      logical yessir,find
      character*20 inbuf
      character*5  token
      character*1 ch
      real a,b,c
 
c     Get misc loop info - discard all but the EPARAM data.
      if(find(32,3,'-->')) stop 'Premature end of Miscloop file.'
      read (32,*) a
      if(find(32,3,'-->')) stop 'Premature end of Miscloop file.'
      read (32,*) a
      if(find(32,3,'-->')) stop 'Premature end of Miscloop file.'
      read (32,*) a,b,c,a
      if(find(32,3,'-->')) stop 'Premature end of Miscloop file.'
      read (32,*) a,b,c
      close (32)
      eparam(5)=nint(a*10)
      eparam(6)=nint(b*10)
      eparam(9)=nint(c*10)
      do i=1,10
        oldep(i)=eparam(i)
      enddo
      ch='N'
      write (outunit,1)
      if (yessir(ch)) then
        write (20,4)
10      write (outunit,2) (eparam(i),i=1,10)
        inbuf=' '
        read (inunit,7,err=10) inbuf
*
        start=1
        fin=20
        call parse (start,fin,inbuf,token)
        if (start-1.eq.fin) goto 20
        i=convt(token)
        call parse (start,fin,inbuf,token)
        j=convt(token)
*
        if (i.lt.1.or.i.gt.10) goto 10
        eparam(i)=j
        write (20,3) i,j
        goto 10
*
** TERMINATOR=CR right here...
*
20    write (20,5)
      endif
      do i=1,10
       eparam(i)=oldep(i)
      enddo
      return
1     format (1x,'Change any of the energy parameters?')
2     format(/,
     .  10x,'* All energy values are 10X integers.  For example',/,
     .  10x,'* -3.2 is -32 in this table.',//,
     .  10x,' 1 Extra stack energy                        [',i5,']',/,
     .  10x,' 2 Extra bulge energy                        [',i5,']',/,
     .  10x,' 3 Extra loop energy (interior)              [',i5,']',/,
     .  10x,' 4 Extra loop energy (hairpin)               [',i5,']',/,
     .  10x,' 5 Extra loop energy (multi)                 [',i5,']',/,
     .  10x,' 6 Multi loop energy/single-stranded base    [',i5,']',/,
     .  10x,' 7 Maximum size of interior loop             [',i5,']',/,
     .  10x,' 8 Maximum lopsidedness of an interior loop  [',i5,']',/,
     .  10x,' 9 Bonus Energy                              [',i5,']',/,
     .  10x,'10 Multi loop energy/closing base-pair       [',i5,']',//,
     .  1x,'Enter Parameter and new value (<return> to main menu)')
3     format (i2/i5)
4     format ('1')
5     format (' ')
6     format ()
7     format (a)
      end
 
      subroutine nbest(cntrl)
*
** Get the three NBest parameters - %,# structures, and distance.
** Stick them in Cntrl, where they belong, and write them out.
*
      integer inunit,outunit
      parameter (inunit=5,outunit=6)
*
      integer cntrl(10)
*
      write (outunit,1)
      call getint(cntrl(8),0,100,cntrl(8))
      write (outunit,2)
      call getint(cntrl(6),1,9999,cntrl(6))
      write (outunit,3)
      call getint(cntrl(9),0,100,cntrl(9))
*
      write (20,4) cntrl(8),cntrl(6),cntrl(9)
*
      return
1     format (1x,'Within how many percent of the minimum?')
2     format (1x,'Maximum number of structures found?')
3     format (1x,'Distance between structures must be greater than:')
4     format (i4,2(/,i4))
      end
 
      subroutine outputfiles(struc,ct,details,cntrl)
*
** Ask the user for the kind of output files they
** want.  Defaults here are imbedded in the code.
*
      integer inunit,outunit
      parameter (inunit=5,outunit=6)
      integer cntrl(10)
      character*80 struc,ct,details
      logical yessir
      character*1 ch
*
** Structure file?
*
      ch='Y'
      write (outunit,1)
      if (yessir(ch)) then
*
** If a structure file is desired, don't stick the output
** into the log file.
*
        write (20,9) 'Y'
        write (20,9) 'N'
        write (outunit,4) struc
        call getstr(struc)
        call writestr(20,struc)
      else
        write (20,9) 'N'
      endif
*
** CT file?
*
      ch='N'
      write (outunit,2)
      if (yessir(ch)) then
        write (20,9) 'Y'
        write (outunit,7) ct
        call getstr(ct)
        call writestr (20,ct)
      else
        write (20,9) 'N'
      endif
*
** details file?
*
      ch='N'
      write (outunit,3)
      if (yessir(ch)) then
        write (20,9) 'Y'
        write (outunit,8) details
        call getstr(details)
        call writestr(20,details)
      else
        write (20,9) 'N'
      endif
      return
*
1     format (1x,'Make a structure file?')
2     format (1x,'Make a CT file?')
3     format (1x,'Make a details file?')
4     format (1x,'Structure file name? Default:',a40)
6     format (i3)
7     format (1x,'CT file name? Default:',a40)
8     format (1x,'details file name? Default:',a40)
9     format (a1)
      end
*
** BASEMENT ROUTINES (should be in a tool box some where...)
*
      subroutine getint(default,minimum,maximum,i)
*
** Get an integer from the user.  Make sure it's in
** bounds (Minimum <= I <= Maximum).  If the user
** hits return, use the default.  Uses the routine
** Convt, stolen from RNA.
*
      integer inunit,outunit
      parameter (inunit=5,outunit=6)
      integer default,minimum,maximum,i,convt
      character*5 inbuf
*
      inbuf=' '
10    write (outunit,1) default
      read (inunit,2,err=10) inbuf
      if (inbuf.eq.' ') then
        i=default
      else
        i=convt(inbuf)
      endif
      if (i.lt.minimum.or.i.gt.maximum) then
        write (outunit,3) minimum,maximum
        goto 10
      endif
      return
*
1     format (1x,'Default=',i4,' ->')
2     format (a5)
3     format (1x,'Please enter an integer between ',i4,' and ',i4)
      end
 
      integer function convt(str)
*
** Convert a buffer STR to an integer.  Stolen from RNA
*
      integer i,place
      character*5 str
      logical neg
 
      neg = .false.
      place = 0
      convt = 0
 
      do 10 i = 5,1,-1
        if (str(i:i).eq.'-') then
          neg = .true.
        else
          if (str(i:i).ge.'0'.and.str(i:i).le.'9') then
             convt = convt + 10**place * (ichar(str(i:i)) - ichar('0'))
             place = place+1
          endif
        endif
 10     continue
      if (neg) convt = -1 * convt
      return
      end
 
 
      logical function yessir(default)
*
** Ask for, and get either a Y or N response.  If Y,
** this function returns .TRUE., .PHALSE. otherwise.
*
      integer inunit,outunit
      parameter (inunit=5,outunit=6)
      character*1 default,ch,upcase
      ch=' '
10    write (outunit,1) default
      read (inunit,2,err=10) ch
      if (ch.eq.' ') then
        yessir=(default.eq.'Y')
      else
        ch=upcase(ch)
        yessir=(ch.eq.'Y')
                if (ch.ne.'Y') then
                if (ch.ne.'N') goto 10
        endif
      endif
1     format (1x,'(Y,N) Default:',a1,'->')
2     format (a1)
      return
      end
 
       character*1 function upcase (ch)
*
** Stolen from MYLIB.FOR
** Change Ch to upper case if necessary
*
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
 
      subroutine writestr(unit,buffa)
*
** Write out what's ever in the buffer to UNIT.  Kill
** trailing spaces.  Special format for screen (unit=6)
** since it requires line control characters.
*
      integer outunit
      parameter (outunit=6)
      integer unit,i
      character*80 buffa
*
      i=80
10    if (buffa(i:i).ne.' ') goto 20
      i=i-1
      goto 10
*
20    if (unit.eq.outunit) then

        write (unit,1) buffa(1:i)
      else

        write (unit,2) buffa(1:i)
      endif

      return
1     format (1x,a)
2     format (a)
      end
 
      subroutine getstr(default)
*
** Get string from the user.  If it bombs, try again.
*
      integer inunit,outunit
      parameter (inunit=5,outunit=6)
      character*80 default,inbuf
*
10    inbuf=' '
      write (outunit,1)
      read (inunit,2,err=10) inbuf

      if (inbuf.ne.' ') default=inbuf
      return
1     format (1x,'->')
2     format (a)
      end
        subroutine parse(start,fin,inbuf,token)
        integer start,fin,i,tokenstart,tokenfin
        character*20 inbuf
        character*5 token
*
        i=start
        token=' '
        do while (i.le.fin.and.inbuf(i:i).eq.' ')
          i=i+1
        enddo
        tokenstart=i
        do while (i.le.fin.and.inbuf(i:i).ne.' ')
                i=i+1
        enddo
        tokenfin=i-1
      if (tokenfin-tokenstart+1.gt.5) then
        start=fin
          write (*,1) inbuf(tokenstart:tokenfin)
      else
        token(1:tokenfin-tokenstart+1)=inbuf(tokenstart:tokenfin)
        endif
        start=i
        return
1       format (1x,'Number too big to convert: ',a20)
        end
c     Used in reading the energy files.
c     Locates markers in the energy files so that data can be read
c     properly.
      function find(unit,len,str)
      implicit integer (a-z)
      logical find,flag
      character*20 str
      character*80 inrec
 
      find = .false.
      flag = .false.
      do  while(.not.flag)
         read(unit,100,end=200) inrec
         count = 1
         do i = 1,80-len+1
           if (inrec(i:i).eq.str(count:count)) then
             count = count + 1
             if (count.gt.len) flag = .true.
             if (inrec(i+1:i+1).ne.str(count:count)) count = 1
           endif
         enddo
      enddo
      return
100   format(a80)
200   find = .true.
      return
      end

*****
 
** MYLIB.FOR is the library of routines written in fortran.
** Name         Type            Comment
** GetFil       Sub             Get a file name and see if it's old or new.
** Yes          L. Func Query and wait for the answer Y or N.
** UpCase       C*1 Func        Take a single char and convert to upper case.
**
**                    - John Jaeger 28Aug87
*
       subroutine getfil (tsiz,text,fname,isold)
** Get the name of a file and see if it's old or new.
*
        character*80 fname,query,text
        integer tsiz
        logical isold,nofile,yes,itexists
*
        nofile=.true.
103       write (6,1) text(1:tsiz)
      read (5,2) fname
*
** Try to open the old FName to see if it already exists.
** This is to prevent writing over a previous data file by mistake
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
                                fname='SYSOUTPUT'
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
 1      format (1x,a,' file name? ')
 2    format (a)
 5     format (1x,'File does not exist.')
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

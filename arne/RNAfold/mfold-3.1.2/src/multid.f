c     Subroutine multid: adaptation of formid for multiple alignments
c
c     integer:  ftype     file format code:   1  -  GENBANK
c                                             2  -  EMBL
c                                             3  -  FASTA
c                                             4  -  STANFORD/IG
c                                             5  -  GCG
c               ftypes     ARRAY of sequence format types (as above)
c               gcglen     ARRAY of (GCG) sequence lengths
c               idno       counter for sequence identifiers found
c               line       counter for records read from a file
c               maxlen  *  maximum sequence length specified by user
c               pline      ARRAY of pointers to records that 
c                          immediately precede sequence code
c               pointr     index of a sequence chosen by user
c               seqlen **  actual length of the retrieved sequence
c               seqnum  *  index of requested sequence; if  0  then
c                            new sequence file is opened and seqnum
c                      **    becomes the no. of sequences in new file 
c               start      starting column for sequence code
c               stop       expected length of chosen sequence
c
c  characters:  chr        single character of sequence code
c               chrlen     string of 5 - contains (GCG) sequence length  
c               filnam     string of 40 - user supplied input file name
c               nullid     string of 50 - unidentified sequence label
c               rec        string of 80 - a single record from a file
c               seq    **  ARRAY of single characters - the sequence
c               seqid      ARRAY of strings of 50 - all sequence id's
c               ssid   **  string of 50 - single sequence identifier
c     logical:  used       if true, input file has already been scanned
c
c                      **  sent to main from subroutine
c                       *  from main;  returned unchanged
c _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
c
      subroutine multid(ssid,seq,seqlen,maxlen,seqnum)
c
      implicit integer  (a-z)
      character*1    seq(maxlen), chr
      character*5    chrlen, choice
      character*50   seqid(25000), ssid, nullid
      character*512   filnam
      character*80   rec
      logical        used
      dimension      ftypes(25000), gcglen(25000), pline(25000)
c
      data  inseq/1/, inpt/5/, iout/6/, ftype/0/
      data  used/.false./, nullid/'unidentified sequence         '/
c
c     Test user specified max sequence length for value > zero
      if (maxlen.le.0) then
        write (iout,101) ' * ERROR: max sequence length must be > zero '
        write (iout,101) '  (4th argument in call to subroutine MULTID)'
        write (iout,'(/)')
        call exit(1)
      endif
c
c     If requested sequence file already open, skip to retrieval
      if (used.and.seqnum.gt.0) goto 35
c
c     Obtain name of input file and open file or exit session
    1 write (iout,103) ' Input sequence file name  (q to quit)  '
      read (inpt,'(a)',err=1,end=99) filnam
      if (filnam(1:5).eq.'q    ') then
          write (iout,105) ' - End of session . . . good bye'
          write (iout,'(/)')
          stop
      endif
      open(unit=inseq,file=filnam,status='old',err=90)
      used = .true.
c
c     Initialize 
      idno = 0
      line = 0
      ftype = 0
      do 3 i = 1,25000
      ftypes(i) = 0
      gcglen(i) = 0
      pline(i) = 0
      seqid(i) = '                              '
    3 continue
c
c     Identify file type; locate sequence identifier(s) and record(s)
c     that are assumed immediately to precede actual sequence code
   10 read (inseq,107,err=92,end=25) rec
      line = line + 1
      if (index(rec,'..').gt.0) then
          indx = index(rec,'Length:')
          if (indx.gt.0) then                                      ! GCG
              if (idno.eq.0) idno = 1
              if (pline(idno).gt.0) then
c                 Avoid recording same sequence twice (GCG files only)
                  if ((line-pline(idno)).gt.3) idno = idno + 1
              endif
              ftype = 5
              ftypes(idno) = 5                                       
              pline(idno) = line + 1
              ssid = rec(1:(indx-2))
              do 12 i = 1,indx-2
              if (ssid(i:i).ne.' ') then
                  seqid(idno) = ssid(i:(indx - 2))
                  goto 15
              endif
   12         continue
              seqid(idno) = nullid
   15         chrlen = rec((indx+7):(indx+12))
              read (chrlen,*) gcglen(idno)
          endif
      elseif (rec(1:5).eq.'LOCUS') then                        ! GENBANK
          idno = idno + 1
          ftype = 1                                           
          ftypes(idno) = 1                                           
          ssid = rec(13:62)
          seqid(idno) = ssid(1:index(ssid,' '))
      elseif (rec(1:6).eq.'ORIGIN'.and.ftype.eq.1) then 
          if (pline(idno).gt.0) then
              idno = idno + 1
              ftypes(idno) = 1                                           
              seqid(idno) = nullid
          endif
          pline(idno) = line
      elseif (rec(1:5).eq.'ID   ') then                           ! EMBL
          idno = idno + 1
          ftype = 2                                           
          ftypes(idno) = 2
          ssid = rec(6:55)
          seqid(idno) = ssid(1:index(ssid,' '))
      elseif (rec(1:6).eq.'SQ   S'.and.ftype.eq.2) then
          if (pline(idno).gt.0) then
              idno = idno + 1
              ftypes(idno) = 2
              seqid(idno) = nullid
          endif
          pline(idno) = line
      elseif (rec(1:1).eq.'>') then                                ! FASTA
          idno = idno + 1
          ftype = 3
          ftypes(idno) = 3                                           
          pline(idno) = line
          ssid = rec(2:51)
          seqid(idno) = ssid
      elseif (rec(1:1).eq.';') then                           ! STANFORD
   20     read (inseq,107,err=92,end=25) rec
          line = line + 1
          if (rec(1:1).eq.';') goto 20
          idno = idno + 1
          ftype = 4
          ftypes(idno) = 4                                           
          pline(idno) = line
          seqid(idno) = rec(1:50)
      endif
      goto 10
c
c     EOF: if no recognizable sequence found, prompt for new file
   25 if (idno.gt.0) then
          seqnum = idno
          do 27 i = 1,idno
          if (pline(i).gt.0) return
   27     continue
      endif
      write (iout,101) ' * No recognizably formatted sequence found '
      close(inseq)
      goto 1
c
c     Check user's choice of sequence 
   35 pointr = seqnum
      if (pointr.gt.idno.or.pline(pointr).eq.0) then
          write (iout,101) ' * Chosen sequence non-existent or empty *'
          write (iout,'(/)')
          call exit(1)
      endif
c
c     Rewind file
      rewind(inseq,err=94)
c
c     Choose starting column and set limit on length of sequence
      start = 1
      ftype = ftypes(pointr)
      if (ftype.eq.1.or.ftype.eq.5) start = 11
      stop = maxlen
      if (ftype.eq.5) stop = min0(gcglen(pointr),stop)
      ssid = seqid(pointr)
c     Read through all records preceding chosen squence
      do 37 i = 1,pline(pointr)
      read (inseq,107,err=92,end=50) rec
   37 continue
      seqlen = 0
c
c     Read a record containing sequence code; check for terminator
   40 read (inseq,107,err=92,end=50) rec
      if (rec(1:1).eq.'/') return
c
c     Retrieve sequence one character at a time from current record
      do 45 i = start,80
      chr = rec(i:i)
c     Ignore blank characters
      if (chr.eq.' ') goto 45
c     Check for sequence terminators
      if (chr.eq.'>'.or.chr.eq.'1'.or.chr.eq.'2'.or.chr.eq.'*') return
c     Add acceptable character to sequence
      seqlen = seqlen + 1
      seq(seqlen) = chr
c     Check sequence length
      if (seqlen.eq.stop) return
   45 continue
      goto 40
c
   50 write (iout,101) ' * End of file encountered * '
      return
c
c     File handling errors
   90 write (iout,109) ' * Could not open file:   ',filnam
      goto 1
   92 write (iout,105) ' * ERROR in reading file * '
      write (iout,'(/)')
      call exit(1)
   94 write (iout,105) ' * ERROR in rewinding file * '
      write (iout,'(/)')
      call exit(1)
 99   write (iout,*) ' No sequences read. EOF encountered. '
      call exit(1)
c
  101 format (/ a45)
  103 format (/ a40)
  105 format (/ a34)
  107 format (a80)
  109 format (/ a26,a40 /)
      end

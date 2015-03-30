c     ct to nbp conversion
c     ct file --> new bp format - the Madison format, version 1.0
c     Standard in/Standard out
c     First historical number in first folding is used do determine
c     the offset (assumed constant) for all foldings
c
      implicit integer (a-z)
      parameter (max_n=100000,max_helix=10000)
      character*1 seq(max_n),sequ
      character*3 type
      character*60 rec
      character*80 ctrec
      integer helix_l(max_helix),helix_i(max_helix),helix_j(max_helix)

      if (iargc().eq.0) then
         type = 'RNA'
      else
         call getarg(1,type)
         if (type.ne.'DNA') type = 'RNA'
      endif

      n_struct = 0
      offset = 0
c     This program can handle many structures in a single file.
      do while (1.eq.1)
         read(5,1010,end=999) n,rec
         n_struct = n_struct + 1
         n_helix = 0
         n_mol = 0
         in_helix = 0
         do i = 1,n
            read(5,1015,end=998) ctrec
            if (n_struct.eq.1) then
c     jm1 and jp1 are not used in this version
               read(ctrec,1020,err=998) seq(i)
               read(ctrec(8:80),*,err=998) jm1,jp1,j,k
               if (i.eq.1) offset = k - 1
            else
               read(ctrec,1020,err=998) sequ
               read(ctrec(8:80),*,err=998) jm1,jp1,j
               if (seq(i).ne.sequ) then
                  write(6,*) 'STOP: Sequences do not match.'
                  call exit(1)
               endif
            endif
            if (i.lt.j) then
               if (in_helix.gt.0) then
                  if (helix_j(n_helix)-helix_l(n_helix).eq.j) then
                     helix_l(n_helix) = helix_l(n_helix) + 1
                  else   
                     in_helix = 1
                     n_helix = n_helix + 1
                     helix_l(n_helix) = 1
                     helix_i(n_helix) = i
                     helix_j(n_helix) = j
                  endif
               else
                  in_helix = 1
                  n_helix = n_helix + 1
                  helix_l(n_helix) = 1
                  helix_i(n_helix) = i
                  helix_j(n_helix) = j
               endif
            else
               in_helix = 0
            endif
         enddo
         if (n_struct.eq.1) then
            lab_start = index(rec,']') + 1
            if (lab_start.eq.1) lab_start = index(rec,'=') + 9
            if (lab_start.eq.9) lab_start = 7
            write(6,*) 'Ensemble_name: ', rec(lab_start:60), ';'
            write(6,*) 'Number_of_molecules: 1 ;'
            write(6,*) 'Molecule_name: ', rec(lab_start:60), ';'
            write(6,*) 'Molecule_type: ', type, ';'
            write(6,*) 'Sequence: ', (seq(i),i=1,n), ';'
            write(6,*) 'Sequence_length: ', n, ';'
            write(6,*) 'Sequence_starting_number: ', offset + 1, ';'
            write(6,*) 'Sequence_offset: ', offset, ';'
         endif
         write(6,*) '# Structure ', n_struct, ' ', rec(1:lab_start-1)
         write(6,*) 'Structure: {'
         do k = 1,n_helix
            write(6,1040) helix_i(k),helix_j(k),helix_l(k)
         enddo
         write(6,*) '}'
      enddo


 1010 format(i5,1x,a60)
 1015 format(a80)
 1020 format(6x,a1)
 1040 format(' H ',2i6,i4,' ;')
 998  write(6,*) 'STOP: Premature end or error in ct file.'
      call exit(1)
 999  call exit(0)
      end

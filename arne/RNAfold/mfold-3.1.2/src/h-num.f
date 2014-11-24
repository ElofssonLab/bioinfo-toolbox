	parameter (nmax=200000)
	integer level,length,istart,jstart,energy,p_num(nmax),h_num
        integer emin
	character*80 record

c	initialize

	do i = 1,nmax
		p_num(i) = 0
	enddo
	emin = 100000

c	first sweep through helix file : compute p_num

	read(21,1010) record
1010	format(a80)
	do while (1.eq.1)
		read(21,1020,end=100) level,length,istart,jstart,energy
1020		format(5i7)
		do k = 0,length-1
			p_num(istart+k) = p_num(istart+k) + 1
			p_num(jstart-k) = p_num(jstart-k) + 1
		enddo
		if (energy.lt.emin) emin = energy
	enddo

c	second pass : compute h_num

100	rewind(21)
	read(21,1010) record
	record(32:37) = 'h-num '
	write(22,1010) record
	do while (1.eq.1)
		read(21,1020,end=200) level,length,istart,jstart,energy

		h_num = 0
		do k = 0,length-1
			h_num = h_num + p_num(istart+k) + p_num(jstart-k) - 1
		enddo
		write(22,1021) level,length,istart,jstart,float(h_num)/float(length)
1021		format(4i7,f7.1)
	enddo

200	stop
	end

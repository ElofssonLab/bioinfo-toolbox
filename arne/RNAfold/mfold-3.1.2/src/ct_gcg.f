C     Generates a CT file. (Richard Feldmann)
      subroutine ct(r)
      include 'rna.inc'
      real r
 
      write(21,100) seqlab,hstnum(1),hstnum(n),n,r
      do k = 1,n
        k1 = k+1
        if (k.eq.n) k1 = 0
        write (21,200) k, seq(hstnum(k)),k-1,k1,basepr(k),hstnum(k)
      enddo
      return
100    format('FOLD of: ',a50,' Check: 0 from: ',i5,' to: ',i5,/,
     .        'Length: ',i5,' Energy: ',f7.1,/,'..')
200   format(i5,1x,a1,4i6)
      end

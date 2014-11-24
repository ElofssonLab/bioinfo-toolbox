      subroutine bpsrch(nb1,bp1,nb2,bp2,d)
c
c     (FS383)
c
      integer*4 d,d1,d2, bp1(1000),bp2(1000)
c
c     OBTAIN LARGEST MINIMAL DISTANCE WHERE DISTANCE IS THE GREATER
c     OF THE DIFFERENCES BETWEEN THE I's AND THE J's, RESPECTIVELY
c
      do i1 = 1,nb1
           d1 = nb1
           j1 = bp1(i1)
           if (j1.gt.i1.and.bp2(i1).ne.j1) then
              do i2 = 1,nb2
                   j2 = bp2(i2)
                   if (j2.gt.i2) then
                      d2 = max0(iabs(i1-i2),iabs(j1-j2))
                      d1 = min0(d1,d2)
                      if (d1.eq.d2) then
c                       SAVE INDICES OF CRITICAL BASE PAIRS
                        is1 = i1
                        js1 = j1
                        is2 = i2
                        js2 = j2
                      end if
                   end if
              end do
              d = max0(d,d1)
              if (d.eq.d1) then
c               SAVE INDICES OF CRITICAL BASE PAIRS
                ip1 = is1
                jp1 = js1
                ip2 = is2
                jp2 = js2
              end if
           end if
      end do
c
c     WRITE OUT INDICES OF CRITICAL BASE PAIRS AND RETURN DISTANCE
c     VALUE TO MAIN PROGRAM
c
!     WRITE (4,101) IP1,JP1,IP2,JP2,D
  101 format (/ 2i5,5x,2i5,5x,i5)
      return
      end

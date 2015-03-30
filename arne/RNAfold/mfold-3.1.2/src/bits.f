c     Marks a base-pair I,J.
c     Assumes that 1 <= I <= J <= N.
c     The information is stored in a single bit in the MARKS
c     array.
c     The conversion from double dimension to single is through the
c     transformation  I,J  ==>  (J-1)*J/2 + I .
      subroutine smark(i,j)
      include 'rna.inc'
      integer*2 bit
 
      posn = (((j-1)*j)/2) + i
      word = (posn+15) / 16
      bit = mod(posn,16)
      marks(word) = iibset(marks(word),bit)
      return
      end
 
c     sfce(i,j) is true if and only if there is a k such
c     that i < k < j and k is forced to be double stranded.
c     Note that the input base pairs belong to the extended
c     set containing both i.j and j.i+n for 1 <= i <= j <= n.
c     The information is stored in a single bit in the FORCE2
c     array.
c     The conversion from double dimension to single is through the
c     transformation  i,j  ==>  (j-i)*n + i

      subroutine sfce(i,j)
      include 'rna.inc'
      integer*2 bit
 
      posn = (j-i)*n + i 
      if (i.gt.n) posn = posn - n
      word = (posn+15) / 16
      bit = mod(posn,16)
 
      force2(word) = iibset(force2(word),bit)
      return
      end
 
c     Retrieves information on whether or not the base-pair I,J
c     has been marked by a traceback passing through or close to
c     this pair.
      logical function mark(i,j)
      include 'rna.inc'
      integer*2 bit,one
 
      one = 1
      posn = (((j-1)*j)/2) + i
      word = (posn+15) / 16
      bit = mod(posn,16)
 
      set = iibits(marks(word),bit,one)
      mark = .false.
      if (set.ne.0) mark = .true.
      return
      end
 
 
c     Retrieves information on whether or not the open interval (i,j)
c     contains a base that must be single stranded
      logical function fce(i,j)
      include 'rna.inc'
      integer*2 bit,one
  
      one = 1
      posn = (j-i)*n + i
      if (i.gt.n) posn = posn - n
      word = (posn+15) / 16
      bit = mod(posn,16)
 
      set = iibits(force2(word),bit,one)
      fce = .false.
      if (set.ne.0) fce = .true.
      return
      end
 

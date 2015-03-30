      subroutine testdate
      integer day,month,year,i1,i2,i3
      data day/32/,month/13/,year/2999/
      call system('date +"%d %m %Y" > date.test')
      open (unit=2,file='date.test',status='unknown')
      read(2,*) i1,i2,i3
      close(2)
      call system('\\rm -f date.test')
      if (i3.lt.year) then
         return
      elseif (i3.gt.year) then
         call notice
      else
         if (i2.lt.month) then
            return
         elseif (i2.gt.month) then
            call notice
         else
            if (i1.lt.day) then
               return
            elseif (i1.gt.day) then
               call notice
            else
               call system('echo "Your license expires at 24:00." >/dev/tty')
            endif
         endif
      endif   
      return
      end
      subroutine notice
      call system('echo "Sorry, but your license has expired." >/dev/tty')
      call exit(1)
      end

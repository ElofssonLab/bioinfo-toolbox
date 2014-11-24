      subroutine getstring(propt,userst)
      include '/usr/include/gl/fgl.h'
      include '/usr/include/gl/fdevice.h'
      character*(*) propt
      character*(*) userst
      logical keyboardwasqueued
      integer*4 filex,filey,fileyhi,textx,texty,curstrlen,maxwid,maxxval,
     .     oldmode,xori,yori,wxsize,wysize,promptlen,maxlen,dev
      integer*2 c,mask1,mask2,mask3,mask4
      external getdra,drawmo

      character CTRL_H, CTRL_J, CTRL_M, CTRL_U

      parameter (CTRL_H = char(8),
     .           CTRL_J = char(10),
     .           CTRL_M = char(13),
     .           CTRL_U = char(21) )


      promptlen = len(propt)
      maxlen = len(userst)
      filex = 5
      filey = 15
      fileyhi = (30+filey)
      textx = (filex+5)
      texty = (filey+10)
      userst(1:maxlen) = ' '
      call pushma
      oldmode = getdra()
      call getscr(mask1,mask2,mask3,mask4)
      keyboardwasqueued = isqueu(KEYBD)
      call drawmo(PUPDRW)
      call getori(xori,yori)
      call getsiz(wxsize,wysize)
      call ortho2(-0.5,float(wxsize),-.5,float(wysize))
      call gconfi
      maxxval = wxsize + xori
      maxwid = (wxsize - 11) - (filex + strwid(propt,promptlen))
      call scrmas(filex,wxsize - 6, filey, fileyhi)
      curstrlen = 0
      call color(white)
      call clear
      call color (pupblk)
      call linewi(2)
      call recti(filex+2,filey+2,wxsize-8,fileyhi-1)
      call linewi(1)
      call cmov2i(textx,texty)
      call charst(propt,promptlen)
      call qdevic(KEYBD)
      do while(.true.)
         dev = qread(c)
         if (dev.eq.KEYBD) then
            if (char(c).eq.CTRL_U) then
               curstrlen = 0
               call color(white)
               call clear
               call color (pupblk)
               call linewi(2)
               call recti(filex+2,filey+2,wxsize-8,fileyhi-1)
               call linewi(1)
               call cmov2i(textx,texty)
               call charst(propt,promptlen)
            else if ((char(c).eq.CTRL_M).or.(char(c).eq.CTRL_J)) then
               if (.not.keyboardwasqueued) call unqdev(KEYBD)
               call scrmas(mask1+0,mask2+0,mask3+0,mask4+0)
               call drawmo(PUPDRW)
               call color(PUPCLR)
               call clear
               call drawmo(oldmode)
               call popmat
               return
            else if (char(c).eq.'\b') then
               if (curstrlen.gt.0) then
                  curstrlen = curstrlen - 1
                  userst(curstrlen+1:curstrlen+1) = ' '
                  call color(white)
                  call clear
                  call color (pupblk)
                  call linewi(2)
                  call recti(filex+2,filey+2,wxsize-8,fileyhi-1)
                  call linewi(1)
                  call cmov2i(textx,texty)
                  call charst(propt,promptlen)
                  call charst(userst,curstrlen)
               endif
            else
               if (curstrlen.lt.maxlen) then
                  userst(curstrlen+1:curstrlen+1) = char(c)
                  call charst(userst(curstrlen+1:curstrlen+1),1)
                  curstrlen = curstrlen + 1
               else
                  call ringbe
               endif
            endif
         endif
      end do
      return
      end

      SUBROUTINE alert(string,time)
      INCLUDE '/usr/include/gl/fgl.h'
      INCLUDE '/usr/include/gl/fdevice.h'
      CHARACTER*(*) string
      INTEGER*4 xmin,xmax,ymin,ymax,xori,yori,wxsize,wysize,oldmode,stlen,
     .     textx,texty,time
      INTEGER*2 mask1,mask2,mask3,mask4
      stlen = len(string)
      call pushma
      oldmode = getdra()
      call getscr(mask1,mask2,mask3,mask4)
      call drawmo(PUPDRW)
      call getori(xori,yori)
      call getsiz(wxsize,wysize)
      textx = strwid(string,stlen)
      xmin = (wxsize - xori - (textx + 20))/2
      xmax = (wxsize - xori + (textx + 20))/2
      ymin = yori + 50
      ymax = ymin + 34
      call scrmas(xmin,xmax,ymin,ymax)
      call color(white)
      call clear
      call color(pupblk)
      call linewi(2)
      call recti(xmin+2,ymin+2,xmax-2,ymax-2)
      call linewi(1)
      xmin = (wxsize - xori - textx)/2
      ymin = ymin + 10
      call cmov2i(xmin,ymin)
      call charst(string,stlen)
      call sleep(time)
      call scrmas(mask1+0,mask2+0,mask3+0,mask4+0)
      call drawmo(PUPDRW)
      call color(PUPCLR)
      call clear
      call drawmo(oldmode)
      call popmat
      return
      end

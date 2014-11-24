      SUBROUTINE dotplt(iret,jret,jump)
c
c     by ROLAND GABOURY
c     
c     FEBRUARY/MARCH, 1990
c     ORIGINAL ALGORITHM BY M. ZUKER
      
      INCLUDE 'rna.inc'
      INCLUDE '/usr/include/gl/fgl.h'
      INCLUDE '/usr/include/gl/fdevice.h'
      INCLUDE 'locals.inc'
      
      INTEGER*4 menval,x,y,inv,curs,score,iret,jret
      
      INTEGER*2 val,curmap(128)
      CHARACTER*1 char,label(60)
      CHARACTER*30 plotfile,choice
      LOGICAL pnumon,coords_entered,quitsub,used,texton,pnumused,iconpnu
      SAVE filter,used,curs,curmap,quitsub,texton,score
      DATA used/.false./,pnumon/.false./,pnumused/.false./
      windwid = 1279
      windheit = 1023
      if (.not.used) then                       !have we been through once already?
         used = .true.                          !make sure we don't come through again
         texton = .false.                       !set up the appropriate variables for the first run through
         wins(1) = 0                            !the wins array is for window id's.  If a window is not assigned,
         wins(2) = 0                            !the variable MUST be set to 0.
         wins(3) = 0
         score = 0
         iret = 0
         jret = 0
         vinc = 0
         colors = 1
         call setup                             !code below...initial setup (only stuff that's done once)
         call pupsetup                          !code below...set up popup menus.
         call cmapinit                          !code below...initialize our color map.
         call setupq                            !code below...set up the event queue
         call initplot                          !code below...first plotting routine (does border and axes)
         call doplot                            !code below...do the plot.
         call dostats(score,iret,jret)          !code below...any variables to be written are done here.
      end if
      quitsub = .false.                         !flag to allow us to exit the main loop (back to the program)
      coords_entered = .false.                                                  
      do while(.not.quitsub)                    !loop until user quits or computes structure.
         do while(qtest().ne.0)                 !as long as there's something on the event queue.
            inv = qread(val)                    !read the item from the event queue.
            if (inv.eq.WINFRE) then             !if the window is iconified...
               call winset(wins(1))             !make sure we go back to main window.
            else if(inv.eq.QKEY) then                                           
               call gexit
               stop
            else if (inv.eq.redraw) then        !if REDRAW token:  set to the appropriate window and redraw it.
               if (val.eq.wins(1)) then
                  call winset(wins(1))
                  call reshap                   !in case the window has been resized
                  call color(240)
                  call clear
                  call initplot
                  call doplot
                  call dostats(score,iret,jret)
               else if(val.eq.wins(2)) then     !these apply to the p-num plot window.
                  call winset(wins(2))
                  call reshap
                  call pnumplot
                  call winset(wins(1))
               end if
            else if (inv.eq.WINSHU) then        !has the close box on p-num plot been pressed...
               if (val.eq.wins(2)) then
                  call winset(wins(2))
                  call winclo(wins(2))          !close the window and set the pnumon flag (tells whether it's on or not)
                  call winset(wins(1))
                  pnumon = .false.
                  wins(2) = 0
               end if
            else if ((inv.eq.middle).or.(inv.eq.leftmo)) then !if the middle or left mouse button has been pressed.
               call winset(wins(1))             !reset to main plot window
               if (val.eq.1) then
                  x = getval(mousex)            !do the routines to pick a point and display the historical
                  y = getval(mousey)            !i,j location ... i,j goes into dostats to be printed.
                  call getpixel(x,y,iret,jret,score)
                  if (iret.ne.-1) then          !make sure that only valid coordinates set this flag...lets
                     coords_entered = .true.    !us exit the routine with i,j
                  else
                     coords_entered = .false.
                  endif
                  call dostats(score,iret,jret)
               endif
            else if (inv.eq.rightm) then        !the menu button (right mouse) has been chosen.
               menval = dopup(pup)              !display the popup menu.
               if(menval.eq.1) then             !get the menu selection and act appropriately.
                  call winset(wins(1))          !get vinc has been chosen
                  call getinc                   !code below...reads in the maximum energy increment
                  if (pnumon) then              !if p-numplot is active, redraw it with the new increment
                     call winset(wins(2))
                     call pnumplot
                  end if
                  call qenter(redraw,wins(1))   !redraw the main screen with the new value.
               else if(menval.eq.2)then         !toggle p-numplot chosen?
                  if (.not.pnumon) then         !either turn on or off the p-numplot.
                     pnumon = .true.
                     if (wins(2).eq.0) then     !if there is no window assigned to wins(2), open a new window
                        call minsiz(200,200)
                        wins(2) = winope('P-num plot',10)
                     else                       !otherwise just move to that window.
                        call winset(wins(2))
                     end if
                     call pnumplot
                  else
                     call winclo(wins(2))       !shut off pnumplot...close the window, change the flag and reset 
                     pnumon = .false.           !the window id number.
                     call winset(wins(1))
                     wins(2) = 0
                     call alert('P-num plot is now off',3) !print a box that say's that you turned it off (in case it was
                  end if                        !behind something.
               else if (menval.eq.4) then       !item 4: turn on/off the textport.
                  if (texton) then
                     call tpoff
                     texton = .false.
                  else
                     call tpon
                     texton = .true.
                  endif
               else if (menval.eq.5) then       !if coordinates have been entered, send i,j back to main program.
                  if (coords_entered) quitsub = .true.
               else if (menval.eq.6) then                                      
                  call fileplot2                !code below...make a plot file for hard copy
                  call alert('Finished fileplot',3) !then print an alert box that say's when it's done
               else if (menval.eq.99) then
                  call gexit
                  stop
               else if ((menval.ge.21).and.(menval.le.27)) then !set the number of colors and then redraw
                  colors = menval - 20
                  call color(240)
                  call qenter(redraw,wins(1))
               else if (menval.eq.7) then
                  call qenter(redraw,wins(1))
               endif
            endif
         end do
      end do
      return
      end
      
      
      SUBROUTINE setup                          !subroutine to initialize routine and open windows, etc.
      INCLUDE 'locals.inc'
      INCLUDE '/usr/include/gl/fgl.h'
      INCLUDE '/usr/include/gl/fdevice.h'
      INTEGER*4 curs
      INTEGER*2 curmap(128)
      
      call foregr()                             !make the program run in the foreground
      call nobord                               !main window has no borders so it can't be resized/moved, etc.
      call prefpo(0,windwid,0,windheit)         !set preferred position for the window.
      wins(1) = winope('DOTPLOT',7)             !open the main plot window
      curs = 1
      call drawmo(cursdr)                       !set drawmode to change the cursor
      call cursty(ccross)                       !set a cross-hair cursor
      call defcur(curs,curmap)                  !remap the cursor (shape is automatic for cross-hair)
      call curori(curs,0,0)                     !set the origin of the cursor to the center.
      call setcur(curs,0,0)                     !enable the new cursor.
      call drawmo(NORMDR)                       !put us back into normal drawing mode.
      call textpo (600,1200,750,1000)           !set the size of the textport.
      call color(240)                           !clear the screen
      call clear
      return
      end
      
      SUBROUTINE pupsetup                       !procedure to initialize the popup menus.
      INCLUDE 'locals.inc'
      INCLUDE '/usr/include/gl/fgl.h'
      INCLUDE '/usr/include/gl/fdevice.h'
      
      cpup = newpup()                           !first we need two popup menus: one for the color-rollover
      pup = newpup()                            !and one for the main menu.
      call addtop(cpup,"COLORS %t|1 color %x21|2 colors %x22",36) !then we add the appropriate entries to the menus.
      call addtop(cpup,"3 colors %x23|4 colors %x24",27) !the %x changes the value returned by the popup menu to
      call addtop(cpup,"5 colors %x25|6 colors %x26",27) !the value that follows the %x.
      call addtop(cpup,"7 colors %x27",13)                                     
      call addtop(pup,"MENU %t|Enter new increment %x1",31) !%t is for the title of the menu.
      call addtop(pup,"Toggle pnumplot %x2",19)
      call addtop(pup,"Toggle textport on/off %x4",26)
      call addtop(pup,"Colors: %m",10,cpup)     !%m causes the pup manager to look for a rolover menu...the last
      call addtop(pup,"Compute structure for last (i,j) %x5",36) !field in the call to addtop
      call addtop(pup,"Create Plot file %x6",20)
      call addtop(pup,"Redraw screen %x7",17)
      call addtop(pup,"quit %x99",9)
c     call setpup(pup,3,pupgre)
                                                !changes the appearance and function of that pup entry
      return                                    !un-comment it if that entry is to be greyed out...it will
      end                                       !not return a value. 
      
      SUBROUTINE setupq                         !procedure to enable the appropriate devices to make 
      INCLUDE '/usr/include/gl/fdevice.h'       !entries on the event queue.
      call qdevic(QKEY)
      call qdevic(WINSHU)                       !the event queue is, essentially a buffer for input events
      call qdevic(WINFRE)                       !normally, on older systems, when you have to find out if someone
      call qdevic(ESCKEY)                       !pushed a key or a button, you had to "be there when it happened"...you
      call qdevic(LEFTMO)                       !had to be waiting for it, which meant, of course that the computer had
      call qdevic(MIDDLE)                       !to be waiting along with you (?!?).
      call qdevic(RIGHTM)                       !now, the event queue receives the "events"...pushed buttons, etc., and, using
      call qdevic(REDRAW)                       !the routines which handle the event queue, you pick the events off the
      call qdevic(FKEY)                         !queue when *you* are ready for them (as long as it's in a reasonably timely
      return                                    !fashion.  as long as you are busy doing something else, the events will be
      end                                       !entered on the event queue, waiting for the program to take them off.
                                                !queue management routines include routines which tell how many items are
                                                !on the event queue, pull events (activated devices, etc) off the queue along
                                                !with any associated values, force events onto the queue, and resetting the
                                                !entire queue.
      
      
      SUBROUTINE doborder                       !procedure to draw border aroun the main plot window...
      INCLUDE '/usr/include/gl/fgl.h'           !draws multiple rectangles of increasingly darker shades
      INCLUDE 'locals.inc'                      !of grey to give a "framed" effect.
      INTEGER*4 xmin,xmax,ymin,ymax,i
      
      xmin = 0
      xmax = windwid
      ymin = 0
      ymax = windheit
      
      call ortho2(0.0,float(windwid),0.0,float(windheit))
      call gconfi
      do i = 23,0,-1
         call color(i+32)
         call recti(xmin+i,ymin+i,xmax-i,ymax-i)
      end do
      return
      end
      
      
      SUBROUTINE initplot                       !initial plotting procedure...calls doborder, then sets up axes.
      INCLUDE 'rna.inc'
      INCLUDE 'locals.inc'
      INCLUDE '/usr/include/gl/fgl.h'
      INCLUDE '/usr/include/gl/fdevice.h'
      call doborder
      call axes                                 !code below.
      return
      end
      
      SUBROUTINE axes                           !procedure to set up and label axes.
      INCLUDE 'rna.inc'
      INCLUDE '/usr/include/gl/fgl.h'
      INCLUDE '/usr/include/gl/fdevice.h'
      INCLUDE 'locals.inc'
      INTEGER*4 numtics,i,j,ct,index
      REAL xmin,xmax,ymin,ymax, vert1(2),vert2(2),vert3(2)
      CHARACTER*5 st
      DATA numtics/20/
      SAVE numtics
      
      xmin = -0.4 * n
      xmax = 1.2 * n
      ymin = 1.1 * n
      ymax = -0.2 * n
      call ortho2(xmin,xmax,ymin,ymax)          !set an othographic coordinate system just larger than 
      call gconfi                               !our plotting limits
      numtics = 20
      call linesm(smlon)                        !set smoothline for the triangle
      call color(255)
      vert1(1) = 1.0
      vert1(2) = 1.0
      vert2(1) = float(n)
      vert2(2) = 1.0
      vert3(1) = float(n)
      vert3(2) = float(n)
      call bgnlin                               !use new drawing features (bgnlin) to draw triangle.
      call v2f(vert1)
      call v2f(vert2)
      call v2f(vert3)
      call v2f(vert1)
      call endlin
      call linesm(smloff)
      call color(255)
      ct = (n - 1)/numtics
      do j = 1,n,ct                             !draw tick marks on both axes and print tick labels.
         vert1(1) = float(j)
         vert1(2) = -0.02 * n
         vert2(1) = float(j)
         vert2(2) = 1.0
         call bgnlin
         call v2f(vert1)
         call v2f(vert2)
         call endlin
         write (st,fmt='(I5)') hstnum(j)
c   Zuker corrects small label bug for IRIX 3.2.2 - extra line added on April 3, 1991.
         vert1(1) = float(j-ct)
         call cmov2(vert1(1),vert1(2))          !move graphic character pointer to (x,y)
         call charst(st,5)                      !print a character string (string,length)
      end do
      do i = 1,n,ct
         vert1(1) = 1.02 * n
         vert1(2) = i
         vert2(1) = float(n)
         vert2(2) = i
         call bgnlin
         call v2f(vert1)
         call v2f(vert2)
         call endlin
         write (st,fmt='(I5)') hstnum(i)
         call cmov2(vert1(1),float(i))
         call charst(st,5)
      end do
      
c     Zuker comments out call to dolabe. This is not really needed and causes a problem with IRIX 4.0.5
c     call dolabe(seqlab,progtitle)             !use external C routine (dolabels) to print labels (uses font manager)
c     The following lines put in lables in plain font.
      call ortho2(0.0,float(windwid),0.0,float(windheit))
      call gconfi
      call color(blue)
      call cmov2i(640,970)
      call charst(progtitle,10)
      call cmov2i(580,930)
      call charst(seqlab,50)
      return
      end
      
      SUBROUTINE cmapinit                       !initialize the colormap for the smoothline.
      INCLUDE '/usr/include/gl/fgl.h'           !WARNING: THIS PROCEDURE CHANGES ENTRIES IN THE COLOR MAP.
      INTEGER*4 i,r,g,b,n                       !USE SAVEMAP <FILENAME> BEFORE RUNNING THIS PROGRAM IF YOU 
                                                !WISH TO RETAIN THE CURRENT COLORMAP.  LOADMAP <FILENAME> RELOADS
      do i = 240,255                            !THE COLORMAP SAVED BY SAVEMAP.  OTHERWISE, RUN MAKEMAP AFTER
         n = (255 - i) * 17                     !RUNNING THIS PROGRAM TO RESET THE COLORMAP TO THE DEFAULTS.
         call mapcol(i,n,n,n)                   !(FUTURE VERSIONS MAY FIX THIS!)
      end do
      return
      end
      
      
      SUBROUTINE doplot                         !subroutine to perform actual plot
      INCLUDE 'rna.inc'
      INCLUDE 'locals.inc'
      INCLUDE '/usr/include/gl/fgl.h'
      INTEGER*4 istrt,jstrt,j,i
      REAL xmin,xmax,ymin,ymax
      
      xmin = -0.4 * n
      xmax = 1.2 * n
      ymin = 1.1 * n
      ymax = -0.2 * n
      call ortho2(xmin,xmax,ymin,ymax)          !set up a coordinate system so that we can plot actual (i,j) points.
      call gconfi
      do jstrt = 1,n
         j = jstrt + 1
         i = 0
         call vector(i,j)                       !code below...goes along diagonal to plot "vectors"
      end do
      do istrt = 2,n
         j = n + 1
         i = istrt - 1
         call vector(i,j)
      end do
      return
      end
      
      SUBROUTINE vector (i,j)                   !subroutine to travel along diagonal to do plot (actually uses points)
      INCLUDE 'rna.inc'
      INCLUDE '/usr/include/gl/fgl.h'
      INCLUDE '/usr/include/gl/fdevice.h'
      INCLUDE 'locals.inc'
      INTEGER*4 i,j,pcount,vert1(2),vert2(2),x,y,col,vbest
      
      call color(black)
      do while (i+2.le.j)
         j = j - 1
         i = i + 1
         if (v(i,j)+v(j,n+i).gt.vmin+vinc.or.mark(i,j)) then
            if (flag.ge.cntrl(6)) then
               vert2(1) = j+1
               vert2(2) = i-1
               x = vert1(1)
               y = vert1(2)
               do pcount = 0,vert1(1)-vert2(1)
                  vbest = v(y+pcount,x-pcount) + v(x-pcount,n+y+pcount)
                  call plotcol(col,vbest)       !code below (chooses color based on vbest and colors).
                  call point(x-pcount,y+pcount,col) !draws a symbol based on the color
               end do
            endif
            flag = 0
         else
            flag = flag + 1
            if (flag.eq.1) then
               vert1(1) = j
               vert1(2) = i
            endif
         endif
      enddo
      if (flag.ge.cntrl(6)) then
         vert2(1) = j+1
         vert2(2) = i-1
         x = vert1(1)
         y = vert1(2)
         do pcount = 0,vert1(1)-vert2(1)
            vbest = v(y+pcount,x-pcount) + v(x-pcount,n+y+pcount)
            call plotcol(col,vbest)             !same as above.
            call point(x-pcount,y+pcount,col)
         end do
         call endlin
      endif
      return
      end
      
      SUBROUTINE getinc                         !subroutine to enter minimum energy increment.
      INCLUDE 'rna.inc'
      INCLUDE '/usr/include/gl/fgl.h'
      INCLUDE '/usr/include/gl/fdevice.h'
      INCLUDE 'locals.inc'
      REAL temp
      CHARACTER*10 string1
      call getstring('Enter maximum energy increment in kcals: ',string1) !use standard routine (code external) to enter a character string.
      read (string1,*) temp
      vinc = nint(temp * 10.0)
      return
      end
      
      SUBROUTINE point(x,y,col)                 !subroutine to plot a point (symbol) based on the color chosen for 
      INTEGER x,y,col                           !the point.
      REAL vert1(2),vert2(2),vert3(2),vert4(2),r
      r = 0.4
      col = 0
      if (col.eq.0) then                        !draw a filled rectangle
         call rectf(float(x)-0.4,float(y)-0.4,float(x)+0.4,float(y)+0.4)
      else if (col.eq.1) then                   !draw a "X"
         call move2(float(x)+0.3,float(y)+0.3)  !uses old style drawing routines.
         call draw2(float(x)-0.3,float(y)-0.3)
         call move2(float(x)+0.3,float(y)-0.3)
         call draw2(float(x)-0.3,float(y)+0.3)
      else if (col.eq.2) then                   !draw a filled circle.
         call circf(float(x),float(y),r)
      else if (col.eq.3) then                   !draw a triangle
         vert1(1) = float(x)                    !uses new drawing routines.
         vert1(2) = float(y) + 0.4
         vert2(1) = float(x) + 0.3464
         vert2(2) = float(y) - 0.2
         vert3(1) = float(x) - 0.3464
         vert3(2) = float(y) - 0.2
         call bgnpol
         call v2f(vert1)
         call v2f(vert2)
         call v2f(vert3)
         call endpol
      else if (col.eq.4) then                   !draw a "+"
         call move2(float(x)+0.3,float(y))
         call draw2(float(x)-0.3,float(y))
         call move2(float(x),float(y + 0.3))
         call draw2(float(x),float(y - 0.3))
      else if (col.eq.5) then                   !draw a frame (rectangle)
         call rect(float(x)-0.4,float(y)-0.4,float(x)+0.4,float(y)+0.4)
      else if (col.eq.6) then                   !draw a diamond
         vert1(1) = float(x)
         vert1(2) = float(y) + 0.4
         vert2(1) = float(x) + 0.4
         vert2(2) = float(y)
         vert3(1) = float(x)
         vert3(2) = float(y) - 0.4
         vert4(1) = float(x) - 0.4
         vert4(2) = float(y)
         call bgnpol
         call v2f(vert1)
         call v2f(vert2)
         call v2f(vert3)
         call v2f(vert4)
         call endpol
      else if (col.eq.7) then                   !draw an empty circle
         call circ(float(x),float(y),r)
      end if
      return
      end
      
      SUBROUTINE plotcol(col,score)             !choose a color for the current point
      INCLUDE 'rna.inc'                         !based on the score for the point and vinc/vmin for the plot.
      INCLUDE 'locals.inc'
      INTEGER*4 col,score
      
      if (colors.eq.1) then
         col = 0
         call color(col)
         return
      else
         if (score.eq.vmin) then
            col = 0
            call color(col)
            return
         else
            col = idint( (dflotj(score - vmin)-.0001) / vinc * (colors-1)) + 1
            call color( col)
         end if
      end if
      end
      
      SUBROUTINE getpixel(x,y,iret,jret,vbest)  !subroutine to convert point picked with getval() to plot-coordinates.
      INCLUDE 'rna.inc'                         !it also performs a local check (9x9 matrix with x,y at center) to 
      INCLUDE 'locals.inc'                      !find the best local point, then returns iret and jret.
      INTEGER*4 x,y,locx,locy,iret,jret,itest,jtest,check,vbest
      x = x + 1
      y = y - 1
      locx = (((dble(x)/(windwid+1.0)))*(dble(1.6)*n) - dble(.4)*n)
      locy = ((1.3*n)*((1.0-dble(y)/(1.0+windheit))) - dble(.2)*n)
      if ((locx.ge.1).and.(locy.ge.1).and.(locx.le.n).and.(locy.le.n)) then
         iret = locy
         jret = locx
         vbest = v(iret,jret) + v(jret,n+iret)
         do jtest = locx-1,locx+1
            do itest = locy-1,locy+1
               check = v(itest,jtest) + v(jtest,n+itest)
               if (check.lt.vbest) then
                  vbest = check
                  jret = jtest
                  iret = itest
               end if
            end do
         end do
         return
      else
         iret = -1
         jret = -1
      endif
      end
      
      
      SUBROUTINE dostats(score,iret,jret)       !procedure to print out current values of variables:
      INCLUDE 'rna.inc'                         !(I,J) location, the maximum energy increment,
      INCLUDE '/usr/include/gl/fgl.h'           !the optimal score, and the score for (I,J) 
      INCLUDE 'locals.inc'
      INTEGER*4  x1,x2,y1,y2,score,cnt,x,y,iret,jret
      CHARACTER*7 number
      CHARACTER*5 xmstr,ymstr,scorst
      CHARACTER*40 string
      
      call ortho2(0.0,float(windwid),0.0,float(windheit))
      call gconfi
      x1 = 75
      x2 = 600
      y1 = 300
      y2 = 100
      call color(white)
      call rectfi(x1,y1,x2,y2)
      call color(blue)
      write (number,fmt='(F7.1)') float(vmin)/10.0
      string = 'Optimal score = '//number//' kcal/mole'
      y1 = y1 - 16
      call cmov2i(x1,y1)
      call charst(string,33)
      if (iret.ge.1) then
         x = hstnum(jret)
         y = hstnum(iret)
         write (xmstr,fmt='(I5)') x
         write (ymstr,fmt='(I5)') y
      else
         xmstr = '     '
         ymstr = '     '
      end if
      string = '(i,j) basepair = ('//ymstr//','//xmstr//')'
      y1 = y1 - 16
      call cmov2i(x1,y1)
      call charst(string,30)
      write (number,fmt='(F7.1)') float(score)/10.0
      string = '(i,j) score = '//number//' kcal/mole'
      y1 = y1 - 16
      call cmov2i(x1,y1)
      call charst(string,31)
      write (number,fmt='(F7.1)') float(vinc)/10.0
      string = 'Energy increment = '//number//' kcal/mole'
      y1 = y1 - 32
      call cmov2i(x1,y1)
      call charst(string,36)
      y1 = y1 - 32
      call cmov2i(x1,y1)
      call charst('ENERGY DOT PLOT',15)
      return
      end
      
      SUBROUTINE pnumplot                       !subroutine to plot an x-y plot of the sum of j's plotted for any i, vs. i
      INCLUDE 'rna.inc'
      INCLUDE 'locals.inc'
      INCLUDE '/usr/include/gl/fgl.h'
      INCLUDE '/usr/include/gl/fdevice.h'
      INTEGER*4 i,j,cnt,sum,maxi,sizex,sizey
      REAL xmin,xmax,ymin,ymax,x,y
      CHARACTER*5 number
      call getsiz(sizex,sizey)
      call color(white)
      call clear
      maxi = 0
      do i = 1,n
         sum = 0
         do j = 1,n
            if (j.lt.i) then
               if ((v(j,i) + v(i,j+n)).le.(vmin+vinc)) sum = sum + 1
            else
               if ((v(i,j) + v(j,i+n)).le.(vmin+vinc)) sum = sum + 1
            endif
         enddo
         if (sum.gt.maxi) maxi = sum
      end do
      xmin = -0.2 * n
      xmax = 1.2 * n
      ymin = -0.2 * maxi
      ymax = 1.2 * maxi
      call ortho2 (xmin,xmax,ymin,ymax)
      call gconfi
      call color(black)
      call rect(0.0,0.0,float(n),float(maxi))
      do y = 0.0,float(maxi),float(maxi)/(0.714*float(sizey)/35)
         call move2(-0.05*n,y)
         call draw2(0.0,y)
         call cmov2(-0.17*n,y)
         write (number,fmt='(F5.1)') y
         call charst(number,5)
      end do
      do x = 0.0,float(n),n/(.714*float(sizex)/50)
         call move2(x,-0.05*maxi)
         call draw2(x,0.0)
         call cmov2(x-n/(.714*float(sizex)/25),-0.08*maxi)
         write (number,fmt='(F5.1)') x
         call charst(number,5)
      end do
      do i = 1,n
         call move2i(i,0)
         sum = 0
         do j = 1,n
            if (j.lt.i) then
               if ((v(j,i) + v(i,j+n)).le.(vmin+vinc)) sum = sum + 1
            else
               if ((v(i,j) + v(j,i+n)).le.(vmin+vinc)) sum = sum + 1
            endif
         enddo
         call draw2i(i,sum)
      end do
      
      end
      
      SUBROUTINE fileplot2                      !subroutine to save a plot file (old style)
      INCLUDE 'rna.inc'                         !all interaction is now done using graphics string prompt (getstring)
      INCLUDE 'locals.inc'
      INCLUDE '/usr/include/gl/fgl.h'
      INCLUDE '/usr/include/gl/fdevice.h'
      INTEGER levels
      CHARACTER*10 numstring 
      CHARACTER*30 plotfile
      REAL r1,r2

      r1 = float(vmin)/10.

      call getstring(' Enter helix file name : ',plotfile)
      open(unit=7,file=plotfile,status='unknown')

      if (vinc.gt.0) then
         call getstring('Enter number of levels : ',numstring)
         read (numstring,*) levels
      else
         levels = 1
      endif
c
      write(7,43)
43    format('   level  length istart jstart energy')
      if(levels.gt.1) icrit = vinc/(levels-1)
      do diag = 1,2*n-1
         flagid = 0
         i = (diag+1)/2
         j = (diag+2)/2
         do while (i.ge.1.and.j.le.n)
            level = 0
            check = v(i,j)+v(j,n+i)
            if(levels.gt.1) then
               k = (check + icrit - 1 - vmin)/icrit + 1
               if(k.le.levels) level = k
            else
               if(check.le.vmin+vinc) level = 1
            end if
            if (level.gt.0) then
               if (flagid.ne.level) then
                  if (flagid.gt.0) call plotout(hstnum(i+1),hstnum(j-1),
     .               istart-i,flagid,v(i+1,j-1)+v(j-1,n+i+1),0)
                  flagid = level
                  istart = i
                  jstart = j
               endif
            elseif (flagid.gt.0) then
               call plotout(hstnum(i+1),hstnum(j-1),istart-i,flagid,
     .                    v(i+1,j-1)+v(j-1,n+i+1),0)
               flagid = 0
            endif
            if (i.eq.1.or.j.eq.n) then
               if (flagid.gt.0) call plotout(hstnum(i),hstnum(j),
     .                          istart-i+1,flagid,v(i,j)+v(j,n+i),0)
            endif
            i = i - 1
            j = j + 1
         enddo
      enddo
      call plotout(0,0,0,0,0,1)
      close(7)
      return
c     error return.
 99   print *,' Read error or end of file...'
      return
      end
      subroutine plotout(istart,jstart,length,level,energy,dump)
      integer energy,dump,stack(5,5),k
      data k/0/

      if (dump.eq.0) then
         k = k + 1
         stack(1,k) = level
         stack(2,k) = length
         stack(3,k) = istart
         stack(4,k) = jstart
         stack(5,k) = energy
      endif

      if (k.eq.5.or.(k.gt.0.and.dump.eq.1)) then
         write(7,fmt='(5i7)') ((stack(i,j),i=1,5),j=1,k)
         k = 0
      endif

      return
      end

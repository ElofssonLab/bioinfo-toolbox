      implicit integer (a-z)
 
      parameter (maxn=20000,infinity=16000)
      parameter (maxtloops=100)
      parameter (maxtriloops=50)
      parameter (maxsiz=20000)

      real prelog

      integer loop,break,maxpen,newnum(maxsiz),hstnum(maxn*2),
     .        force(maxn*2),numseq(maxn*2),basepr(maxn),
     .        cntrl(10),n,numtloops,numtriloops
      integer asint3(6,6,5,5,5),asint5(6,6,5,5,5,5,5),bulge(30),
     .        dangle(5,5,5,2),eparam(16),hairpin(30),inter(30),poppen(4),
     .        sint2(6,6,5,5),sint4(6,6,5,5,5,5),sint6(6,6,25,5,5,5,5),
     .        stack(5,5,5,5),tloop(maxtloops,2),triloop(maxtriloops,2),
     .        tstkh(5,5,5,5),tstki(5,5,5,5)

      common /main/ loop,break,newnum,hstnum,force,eparam,cntrl,nsave,
     . n,numseq,prelog,maxpen,numtloops,numtriloops
      common /efiles/ asint3,asint5,bulge,dangle,hairpin,inter,poppen,
     . sint2,sint4,sint6,stack,tloop,triloop,tstkh,tstki

      character*1 seq(2*maxn)
      character*50 ctlabel,ctnam
      common /seq/ seq,ctlabel,ctnam
 
      common /struct/ basepr






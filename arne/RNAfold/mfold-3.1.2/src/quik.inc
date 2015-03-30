      include 'maxn1.inc'
      implicit integer (a-z)
      parameter (fldmax=2*maxn)
      parameter (infinity=999999)
      parameter (maxtloops=100)
      parameter (maxtriloops=50)
      parameter (maxsiz=8000)

      character*200 seqrec
      character*50 seqlab
      character*1 seq(maxsiz)

      integer vst(maxn*maxn),wst(maxn*maxn)

      integer maxpen,numseq(fldmax),work(fldmax,0:2),wmb(fldmax,0:2),
     .        w5(-1:maxn),w3(maxn+2),cntrl(10),nsave(2),
     .        n,vmin,num,numtloops,numtriloops,basepr(maxn)
      integer asint3(6,6,5,5,5),asint5(6,6,5,5,5,5,5),bulge(30),
     .        dangle(5,5,5,2),eparam(16),hairpin(30),inter(30),poppen(4),
     .        sint2(6,6,5,5),sint4(6,6,5,5,5,5),sint6(6,6,25,5,5,5,5),
     .        stack(5,5,5,5),tloop(maxtloops,2),triloop(maxtriloops,2),
     .        tstkh(5,5,5,5),tstki(5,5,5,5)

      real prelog,prec/100.0/

      common /chars/ seqrec,seqlab,seq
 
      common /bigstuff/ vst,wst

      common /main/ maxpen,numseq,work,wmb,w5,w3,cntrl,
     .        nsave,n,vmin,num,numtloops,numtriloops,basepr
      common /efiles/ asint3,asint5,bulge,dangle,eparam,hairpin,inter,poppen,
     .        sint2,sint4,sint6,stack,tloop,triloop,tstkh,tstki,prelog


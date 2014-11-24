      include 'maxn.inc'
      implicit integer (a-z)
      parameter (fldmax=2*maxn)
      parameter (infinity=16000,sortmax=90000)
      parameter (sortmaxp1=sortmax+1)
      parameter (mxbits=(maxn*(maxn+1)+31)/32)
      parameter (maxtloops=100)
      parameter (maxtriloops=50)
      parameter (maxsiz=200000)

      character*1 seq(maxsiz),aux(maxn),lorc
      character*4 usage
      character*50 seqlab
 
      integer*2 vst(maxn*maxn),wst(maxn*maxn),marks(mxbits),force2(2*mxbits)

      integer maxpen,newnum(maxsiz),hstnum(fldmax),force(fldmax),
     .        strand(fldmax),numseq(fldmax),work(fldmax,0:2),wmb(fldmax,0:2),
     .        w5(-1:maxn),w3(maxn+2),cntrl(10),nsave(2),list(3000,4),
     .        listsz,basepr(maxn),heapi(sortmaxp1),heapj(sortmaxp1),
     .        n,break,vmin,num,numtloops,numtriloops,maxbp
      integer asint3(6,6,5,5,5),asint5(6,6,5,5,5,5,5),bulge(30),
     .        dangle(5,5,5,2),eparam(16),hairpin(30),inter(30),poppen(4),
     .        sint2(6,6,5,5),sint4(6,6,5,5,5,5),sint6(6,6,25,5,5,5,5),
     .        stack(5,5,5,5),tloop(maxtloops,2),triloop(maxtriloops,2),
     .        tstkh(5,5,5,5),tstki(5,5,5,5)

      real prelog,prec/10.0/

      common /chars/ seq,aux,lorc,usage,seqlab
 
      common /bigstuff/ vst,wst,marks,force2

      common /main/ maxpen,newnum,hstnum,force,strand,numseq,work,wmb,w5,
     .     w3,cntrl,nsave,list,listsz,basepr,heapi,heapj,n,break,vmin,num,
     .     numtloops,numtriloops,maxbp
      common /efiles/ asint3,asint5,bulge,dangle,eparam,hairpin,inter,poppen,
     .        sint2,sint4,sint6,stack,tloop,triloop,tstkh,tstki,prelog

c      common /list/ list,listsz
c      common /traceback/ basepr
c      common /heap/ heapi,heapj,num
c      common /bits/marks,force2

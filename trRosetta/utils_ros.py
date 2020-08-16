import numpy as np
import random
from pyrosetta import *


def add_intrachain_rst(npz,rst,tmpdir,params,minprob=0.5,LB=1,UB=15,D=20,WD=1000,WB=50,allcontacts=False):
    ########################################################
    # Distance restraints to keep the two chains together
    ########################################################
    dist,omega,theta,phi = npz['dist'],npz['omega'],npz['theta'],npz['phi']
    prob = np.sum(dist[:,:,5:], axis=-1)
    for i in range(params["seqlen1"]):
        for j in range(params["seqlen1"]+1,params["seqlen2"]+params["seqlen1"]):
            # We should limit ourself to constrains that have a probablitu
            if (prob[i,j]>minprob):
                name=tmpdir.name+"/%d.%d-harm.txt"%(i+1,j+1)
                with open(name, "w") as f:
                    # LB is Lower Bound
                    # UB is upper bound
                    # WD is
                    # WB
                    f.write('INTRACHAIN'+'\t%.3f\t%.3f'%(UB,D)+'\n')
                    #f.write('y_axis'+'\t%.3f'%stuple(dist[a,b])+'\n')
                    #f.close()
                #rst_line = 'AtomPair %s %d %s %d FADE %.5f %.5f %.5f %.5f %.5f'%('CB',i+1,'CB',j+1,LB,UB,D,WD,WB)
                #
                #FLAT_HARMONIC x0 sd tol

                #  Zero in the range of x0 - tol to x0 + tol. Harmonic
                #  with width parameter sd outside that
                #  range. Basically, a HARMONIC potential (see above)
                #  split at x0 with a 2*tol length region of zero
                #  inserted.
                  
                rst_line = 'AtomPair %s %d %s %d FLAT_HARMONIC  %.5f %.5f %.5f'%('CB',i+1,'CB',j+1,0,D,UB)
                rst['intrachain'].append([i,j,1.0,rst_line]) # Change file?
            elif (allcontacts):
                name=tmpdir.name+"/%d.%d-harm2.txt"%(i+1,j+1)
                with open(name, "w") as f:
                    # LB is Lower Bound
                    # UB is upper bound
                    # WD is
                    # WB
                    f.write('HARM2'+'\t%.3f\t%.3f'%(UB,D)+'\n')
                    # f.write('FADE'+'\t%.3f\t%.3f\t%.3f\t%.3f'%(LB,UB,D,WB)+'\n')
                    #f.write('y_axis'+'\t%.3f'%stuple(dist[a,b])+'\n')
                    #f.close()
                #rst_line = 'AtomPair %s %d %s %d FADE %.5f %.5f %.5f %.5f %.5f'%('CB',i+1,'CB',j+1,LB,UB,D,WD,WB)
                #
                #FLAT_HARMONIC x0 sd tol

                #  Zero in the range of x0 - tol to x0 + tol. Harmonic
                #  with width parameter sd outside that
                #  range. Basically, a HARMONIC potential (see above)
                #  split at x0 with a 2*tol length region of zero
                #  inserted.
                # BOUNDED lb ub sd rswitch tag
                # rst_line = 'AtomPair %s %d %s %d BOUNDED  %.5f %.5f %.5f %.5f %s '%('CB',i+1,'CB',j+1,0,UB*2,D,0.5,"Bounded")
                rst_line = 'AtomPair %s %d %s %d FLAT_HARMONIC  %.5f %.5f %.5f'%('CB',i+1,'CB',j+1,0,WD,WB)
                rst['intrachain'].append([i,j,1.0,rst_line]) 
    print("Flat harmonic restraints:  %d"%(len(rst['intrachain'])))
        
def add_interacton_areas_rst(npz,rst,tmpdir,params,minprob=0.5,UB=30,D=30,WD=1000,WB=50,allcontacts=False):
    ########################################################
    # Distance restraints to keep the two chains together - using predicted interaction residues
    ########################################################
    dist,omega,theta,phi = npz['dist'],npz['omega'],npz['theta'],npz['phi']
    prob = np.sum(dist[:,:,5:], axis=-1)
    list_i=np.zeros((params["seqlen1"]))
    list_j=np.zeros((params["seqlen1"]+params["seqlen2"]))
    for i in range(params["seqlen1"]):
        for j in range(params["seqlen1"]+1,params["seqlen2"]+params["seqlen1"]):
            # We should limit ourself to constrains that have a probablitu
            if (prob[i,j]>minprob):
                #print (i,j,prob[i,j])
                list_i[i]+=1
                list_j[j]+=1
    for i in range(params["seqlen1"]):
        for j in range(params["seqlen1"]+1,params["seqlen2"]+params["seqlen1"]):
            if (list_i[i]>0 and list_j[j]>0):

                name=tmpdir.name+"/%d.%d-harm.txt"%(i+1,j+1)
                with open(name, "w") as f:
                    # LB is Lower Bound
                    # UB is upper bound
                    # WD is
                    # WB
                    f.write('INTRACHAIN'+'\t%.3f\t%.3f'%(UB,D)+'\n')
                  
                #     f(x) = weight * limit^2 * ( 1 - e^( -(x-x0)^2/limit^2 ) )
                #W=10
                #D=8
                #UB=12
                
                rst_line = 'AtomPair %s %d %s %d FLAT_HARMONIC  %.5f %.5f %.5f'%('CB',i+1,'CB',j+1,0,D,UB)
                #rst_line = 'AtomPair %s %d %s %d TOPOUT  %.5f %.5f %.5f'%('CB',i+1,'CB',j+1,10,8,12)
                #print (rst_line)
                rst['intrachain'].append([i,j,1.0,rst_line]) # Change file?
            elif (allcontacts):
                name=tmpdir.name+"/%d.%d-harm2.txt"%(i+1,j+1)
                with open(name, "w") as f:
                    # LB is Lower Bound
                    # UB is upper bound
                    # WD is
                    # WB
                    f.write('HARM2'+'\t%.3f\t%.3f'%(UB,D)+'\n')
                    # f.write('FADE'+'\t%.3f\t%.3f\t%.3f\t%.3f'%(LB,UB,D,WB)+'\n')
                    #f.write('y_axis'+'\t%.3f'%stuple(dist[a,b])+'\n')
                    #f.close()
                #rst_line = 'AtomPair %s %d %s %d FADE %.5f %.5f %.5f %.5f %.5f'%('CB',i+1,'CB',j+1,LB,UB,D,WD,WB)
                #
                #FLAT_HARMONIC x0 sd tol

                #  Zero in the range of x0 - tol to x0 + tol. Harmonic
                #  with width parameter sd outside that
                #  range. Basically, a HARMONIC potential (see above)
                #  split at x0 with a 2*tol length region of zero
                #  inserted.
                # BOUNDED lb ub sd rswitch tag
                # rst_line = 'AtomPair %s %d %s %d BOUNDED  %.5f %.5f %.5f %.5f %s '%('CB',i+1,'CB',j+1,0,UB*2,D,0.5,"Bounded")
                # Using default from 
                rst_line = 'AtomPair %s %d %s %d FLAT_HARMONIC  %.5f %.5f %.5f'%('CB',i+1,'CB',j+1,0,WD,WB)
                #print (rst_line)
                rst['intrachain'].append([i,j,1.0,rst_line]) 
    print("Intrachain attraction restraints:  %d"%(len(rst['intrachain'])))
        
        
def gen_rst(npz, tmpdir, params):

    dist,omega,theta,phi = npz['dist'],npz['omega'],npz['theta'],npz['phi']

    # dictionary to store Rosetta restraints
    rst = {'dist' : [], 'omega' : [], 'theta' : [], 'phi' : [], 'rep' : [], 'intrachain' : []}

    ########################################################
    # assign parameters
    ########################################################
    PCUT  = 0.05 #params['PCUT']
    PCUT1 = params['PCUT1']
    EBASE = params['EBASE']
    EREP  = params['EREP']
    DREP  = params['DREP']
    PREP  = params['PREP']
    SIGD  = params['SIGD']
    SIGM  = params['SIGM']
    MEFF  = params['MEFF']
    DCUT  = params['DCUT']
    ALPHA = params['ALPHA']

    DSTEP = params['DSTEP']
    ASTEP = np.deg2rad(params['ASTEP'])

    seq = params['seq']

    ########################################################
    # repultion restraints
    ########################################################
    #cbs = ['CA' if a=='G' else 'CB' for a in params['seq']]
    '''
    prob = np.sum(dist[:,:,5:], axis=-1)
    i,j = np.where(prob<PREP)
    prob = prob[i,j]
    for a,b,p in zip(i,j,prob):
        if b>a:
            name=tmpdir.name+"/%d.%d_rep.txt"%(a+1,b+1)
            rst_line = 'AtomPair %s %d %s %d SCALARWEIGHTEDFUNC %.2f SUMFUNC 2 CONSTANTFUNC 0.5 SIGMOID %.3f %.3f\n'%('CB',a+1,'CB',b+1,-0.5,SIGD,SIGM)
            rst['rep'].append([a,b,p,rst_line])
    print("rep restraints:   %d"%(len(rst['rep'])))
    '''


    ########################################################
    # dist: 0..20A
    ########################################################
    nres = dist.shape[0]
    bins = np.array([4.25+DSTEP*i for i in range(32)])
    prob = np.sum(dist[:,:,5:], axis=-1)
    bkgr = np.array((bins/DCUT)**ALPHA)
    attr = -np.log((dist[:,:,5:]+MEFF)/(dist[:,:,-1][:,:,None]*bkgr[None,None,:]))+EBASE
    repul = np.maximum(attr[:,:,0],np.zeros((nres,nres)))[:,:,None]+np.array(EREP)[None,None,:]
    dist = np.concatenate([repul,attr], axis=-1)
    bins = np.concatenate([DREP,bins])
    i,j = np.where(prob>PCUT)
    prob = prob[i,j]
    nbins = 35
    step = 0.5
    for a,b,p in zip(i,j,prob):
        if b>a:
            name=tmpdir.name+"/%d.%d.txt"%(a+1,b+1)
            with open(name, "w") as f:
                f.write('x_axis'+'\t%.3f'*nbins%tuple(bins)+'\n')
                f.write('y_axis'+'\t%.3f'*nbins%tuple(dist[a,b])+'\n')
                f.close()
            rst_line = 'AtomPair %s %d %s %d SPLINE TAG %s 1.0 %.3f %.5f'%('CB',a+1,'CB',b+1,name,1.0,step)
            rst['dist'].append([a,b,p,rst_line])
    print("dist restraints:  %d"%(len(rst['dist'])))

    
    ########################################################
    # omega: -pi..pi
    ########################################################
    nbins = omega.shape[2]-1+4
    bins = np.linspace(-np.pi-1.5*ASTEP, np.pi+1.5*ASTEP, nbins)
    prob = np.sum(omega[:,:,1:], axis=-1)
    i,j = np.where(prob>PCUT)
    prob = prob[i,j]
    omega = -np.log((omega+MEFF)/(omega[:,:,-1]+MEFF)[:,:,None])
    omega = np.concatenate([omega[:,:,-2:],omega[:,:,1:],omega[:,:,1:3]],axis=-1)
    for a,b,p in zip(i,j,prob):
        if b>a:
            name=tmpdir.name+"/%d.%d_omega.txt"%(a+1,b+1)
            with open(name, "w") as f:
                f.write('x_axis'+'\t%.5f'*nbins%tuple(bins)+'\n')
                f.write('y_axis'+'\t%.5f'*nbins%tuple(omega[a,b])+'\n')
                f.close()
            rst_line = 'Dihedral CA %d CB %d CB %d CA %d SPLINE TAG %s 1.0 %.3f %.5f'%(a+1,a+1,b+1,b+1,name,1.0,ASTEP)
            rst['omega'].append([a,b,p,rst_line])
    print("omega restraints: %d"%(len(rst['omega'])))


    ########################################################
    # theta: -pi..pi
    ########################################################
    prob = np.sum(theta[:,:,1:], axis=-1)
    i,j = np.where(prob>PCUT)
    prob = prob[i,j]
    theta = -np.log((theta+MEFF)/(theta[:,:,-1]+MEFF)[:,:,None])
    theta = np.concatenate([theta[:,:,-2:],theta[:,:,1:],theta[:,:,1:3]],axis=-1)
    for a,b,p in zip(i,j,prob):
        if b!=a:
            name=tmpdir.name+"/%d.%d_theta.txt"%(a+1,b+1)
            with open(name, "w") as f:
                f.write('x_axis'+'\t%.3f'*nbins%tuple(bins)+'\n')
                f.write('y_axis'+'\t%.3f'*nbins%tuple(theta[a,b])+'\n')
                f.close()
            rst_line = 'Dihedral N %d CA %d CB %d CB %d SPLINE TAG %s 1.0 %.3f %.5f'%(a+1,a+1,a+1,b+1,name,1.0,ASTEP)
            rst['theta'].append([a,b,p,rst_line])
            #if a==0 and b==9:
            #    with open(name,'r') as f:
            #        print(f.read())
    print("theta restraints: %d"%(len(rst['theta'])))


    ########################################################
    # phi: 0..pi
    ########################################################
    nbins = phi.shape[2]-1+4
    bins = np.linspace(-1.5*ASTEP, np.pi+1.5*ASTEP, nbins)
    prob = np.sum(phi[:,:,1:], axis=-1)
    i,j = np.where(prob>PCUT)
    prob = prob[i,j]
    phi = -np.log((phi+MEFF)/(phi[:,:,-1]+MEFF)[:,:,None])
    phi = np.concatenate([np.flip(phi[:,:,1:3],axis=-1),phi[:,:,1:],np.flip(phi[:,:,-2:],axis=-1)], axis=-1)
    for a,b,p in zip(i,j,prob):
        if b!=a:
            name=tmpdir.name+"/%d.%d_phi.txt"%(a+1,b+1)
            with open(name, "w") as f:
                f.write('x_axis'+'\t%.3f'*nbins%tuple(bins)+'\n')
                f.write('y_axis'+'\t%.3f'*nbins%tuple(phi[a,b])+'\n')
                f.close()
            rst_line = 'Angle CA %d CB %d CB %d SPLINE TAG %s 1.0 %.3f %.5f'%(a+1,a+1,b+1,name,1.0,ASTEP)
            rst['phi'].append([a,b,p,rst_line])
            #if a==0 and b==9:
            #    with open(name,'r') as f:
            #        print(f.read())

    print("phi restraints:   %d"%(len(rst['phi'])))

    return rst

def set_random_dihedral(pose):
    nres = pose.total_residue()
    for i in range(1, nres):
        phi,psi=random_dihedral()
        pose.set_phi(i,phi)
        pose.set_psi(i,psi)
        pose.set_omega(i,180)

    return(pose)


#pick phi/psi randomly from:
#-140  153 180 0.135 B
# -72  145 180 0.155 B
#-122  117 180 0.073 B
# -82  -14 180 0.122 A
# -61  -41 180 0.497 A
#  57   39 180 0.018 L
def random_dihedral():
    phi=0
    psi=0
    r=random.random()
    if(r<=0.135):
        phi=-140
        psi=153
    elif(r>0.135 and r<=0.29):
        phi=-72
        psi=145
    elif(r>0.29 and r<=0.363):
        phi=-122
        psi=117
    elif(r>0.363 and r<=0.485):
        phi=-82
        psi=-14
    elif(r>0.485 and r<=0.982):
        phi=-61
        psi=-41
    else:
        phi=57
        psi=39
    return(phi, psi)


def read_fasta(file):
    fasta=""
    with open(file, "r") as f:
        for line in f:
            if(line[0] == ">"):
                continue
            else:
                line=line.rstrip()
                fasta = fasta + line;
    return fasta


def remove_clash(scorefxn, mover, pose):
    for _ in range(0, 5):
        if float(scorefxn(pose)) < 10:
            break
        mover.apply(pose)


def add_rst_chain2(pose, rst, sep1, sep2, params, nogly=False):

    pcut=params['PCUT']
    seq = params['seq']
    seqlen=params['seqlen1']
    inter=params['interchain'] # THis can be "
    
    array=[]

    if nogly==True:
        array += [line for a,b,p,line in rst['intrachain'] if seq[b]!='G' and seq[a]!='G' and seq[b]!='G' and p>=pcut]
        if inter:
            array += [line for a,b,p,line in rst['dist'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut]
            if params['USE_ORIENT'] == True:
                array += [line for a,b,p,line in rst['omega'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut+0.5] #0.5
                array += [line for a,b,p,line in rst['theta'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut+0.5] #0.5
                array += [line for a,b,p,line in rst['phi'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut+0.6] #0.6
        else:
            array += [line for a,b,p,line in rst['dist'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut  and np.sign(a-seqlen)==np.sign(b-seqlen) ] #0.5
            if params['USE_ORIENT'] == True:
                array += [line for a,b,p,line in rst['omega'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut+0.5  and np.sign(a-seqlen)==np.sign(b-seqlen) ] #0.5
                array += [line for a,b,p,line in rst['theta'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut+0.5  and np.sign(a-seqlen)==np.sign(b-seqlen) ] #0.5
                array += [line for a,b,p,line in rst['phi'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut+0.6 and np.sign(a-seqlen)==np.sign(b-seqlen) ] #0.6
    else:
        array += [line for a,b,p,line in rst['intrachain'] if  p>=pcut]
        if inter:
            array += [line for a,b,p,line in rst['dist']  if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut ]
            if params['USE_ORIENT'] == True:
                array += [line for a,b,p,line in rst['omega'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5]
                array += [line for a,b,p,line in rst['theta'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5]
                array += [line for a,b,p,line in rst['phi'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.6] #0.6
        else:
            array += [line for a,b,p,line in rst['dist'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut  and np.sign(a-seqlen)==np.sign(b-seqlen)  ]
            if params['USE_ORIENT'] == True:
                array += [line for a,b,p,line in rst['omega'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5  and np.sign(a-seqlen)==np.sign(b-seqlen)  ]
                array += [line for a,b,p,line in rst['theta'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5 and np.sign(a-seqlen)==np.sign(b-seqlen) ]
                array += [line for a,b,p,line in rst['phi'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.6 and np.sign(a-seqlen)==np.sign(b-seqlen) ] #0.6



    if len(array) < 1:
        return

    random.shuffle(array)

    # save to file
    tmpname = params['TDIR']+'/minimize.cst'
    with open(tmpname,'w') as f:
        for line in array:
            f.write(line+'\n')
        f.close()

    # add to pose
    constraints = rosetta.protocols.constraint_movers.ConstraintSetMover()
    constraints.constraint_file(tmpname)
    constraints.add_constraints(True)
    #print (pose,tmpname)
    constraints.apply(pose)

    os.remove(tmpname)

def add_intra_rst(pose, rst, sep1, sep2, params, nogly=False):

    pcut=params['PCUT']
    seq = params['seq']
    seqlen=params['seqlen1']
    inter=params['interchain'] # THis can be "
    
    array=[]

    if nogly==True:
        array += [line for a,b,p,line in rst['intrachain'] if seq[b]!='G' and seq[a]!='G' and seq[b]!='G' and p>=pcut]
        array += [line for a,b,p,line in rst['dist'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut  and np.sign(a-seqlen)==np.sign(b-seqlen) ] #0.5
        if params['USE_ORIENT'] == True:
            array += [line for a,b,p,line in rst['omega'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut+0.5  and np.sign(a-seqlen)!=np.sign(b-seqlen) ] #0.5
            array += [line for a,b,p,line in rst['theta'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut+0.5  and np.sign(a-seqlen)!=np.sign(b-seqlen) ] #0.5
            array += [line for a,b,p,line in rst['phi'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut+0.6 and np.sign(a-seqlen)!=np.sign(b-seqlen) ] #0.6
    else:
        array += [line for a,b,p,line in rst['intrachain'] if  p>=pcut]
        array += [line for a,b,p,line in rst['dist'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut  and np.sign(a-seqlen)!=np.sign(b-seqlen)  ]
        if params['USE_ORIENT'] == True:
            array += [line for a,b,p,line in rst['omega'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5  and np.sign(a-seqlen)!=np.sign(b-seqlen)  ]
            array += [line for a,b,p,line in rst['theta'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5 and np.sign(a-seqlen)!=np.sign(b-seqlen) ]
            array += [line for a,b,p,line in rst['phi'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.6 and np.sign(a-seqlen)!=np.sign(b-seqlen) ] #0.6



    if len(array) < 1:
        return

    random.shuffle(array)

    # save to file
    tmpname = params['TDIR']+'/minimize.cst'
    with open(tmpname,'w') as f:
        for line in array:
            f.write(line+'\n')
        f.close()

    # add to pose
    constraints = rosetta.protocols.constraint_movers.ConstraintSetMover()
    constraints.constraint_file(tmpname)
    constraints.add_constraints(True)
    #print (pose,tmpname)
    constraints.apply(pose)

    os.remove(tmpname)

def add_rst(pose, rst, sep1, sep2, params, nogly=False):

    pcut=params['PCUT']
    seq = params['seq']

    array=[]

    if nogly==True:
        array += [line for a,b,p,line in rst['dist'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut]
        if params['USE_ORIENT'] == True:
            array += [line for a,b,p,line in rst['omega'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut+0.5] #0.5
            array += [line for a,b,p,line in rst['theta'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut+0.5] #0.5
            array += [line for a,b,p,line in rst['phi'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut+0.6] #0.6
    else:
        array += [line for a,b,p,line in rst['dist'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut]
        if params['USE_ORIENT'] == True:
            array += [line for a,b,p,line in rst['omega'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5]
            array += [line for a,b,p,line in rst['theta'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5]
            array += [line for a,b,p,line in rst['phi'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.6] #0.6


    if len(array) < 1:
        return

    random.shuffle(array)

    # save to file
    tmpname = params['TDIR']+'/minimize.cst'
    with open(tmpname,'w') as f:
        for line in array:
            f.write(line+'\n')
        f.close()

    # add to pose
    constraints = rosetta.protocols.constraint_movers.ConstraintSetMover()
    constraints.constraint_file(tmpname)
    constraints.add_constraints(True)
    constraints.apply(pose)

    os.remove(tmpname)


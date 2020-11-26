#!/usr/bin/env python3
import math
import argparse
import numpy  as np
import tensorflow as tf
import multiprocessing as mp 
from functools import partial
from Bio.PDB import *


"""
In order to transform coordinates of the ligand (receptor is always considered to be intact) use transformation-rotation matrix and apply the formulas (they are used in aace_gramm):

          xx=xlig_init(k); yy=ylig_init(k); zz=zlig_init(k);
          xlig_model(k)=t(1)+u(1,1)*xx+u(1,2)*yy+u(1,3)*zz
          ylig_model(k)=t(2)+u(2,1)*xx+u(2,2)*yy+u(2,3)*zz
          zlig_model(k)=t(3)+u(3,1)*xx+u(3,2)*yy+u(3,3)*zz

where  xlig_init(k), ylig_init(k) zlig_init(k) are Cartesian coordinates pf the k-th atom for the original ligand position (which you submit to GRAMM fro docking) and
xlig_mid(k), ylig_mod(k) zlig_mod(k) are transformed Cartesian coordinates.

In order to get translation vector t(3) and rotation matrix u(3,3) from the GRAMM output, use the following transformations:

!-----------------------------------------------------------------------------                                            
!               Transforming Euler angles into rotation matrix                                                            
!-----------------------------------------------------------------------------                                            
       alpha=pi*irot1/180.d0;  beta=pi*irot2/180.d0;  gamma=pi*irot3/180.d0;
       cosa=dcos(alpha);cosb=dcos(beta); cosg=dcos(gamma)
       sina=dsin(alpha);sinb=dsin(beta); sing=dsin(gamma)
       u(1,1)=cosg*cosa
       u(1,2)=cosg*sina
       u(1,3)=-sing
       u(2,1)=sinb*sing*cosa-cosb*sina
       u(2,2)=sinb*sing*sina+cosb*cosa
       u(2,3)=cosg*sinb
       u(3,1)=cosb*sing*cosa+sinb*sina
       u(3,2)=cosb*sing*sina-sinb*cosa
       u(3,3)=cosg*cosb

Here irot1, irot2 and irot3 are Euler rotation angles given by GRAMM.

!-----------------------------------------------------------------------------                                            
!                 Tranforming translation vector                                                                          
!-----------------------------------------------------------------------------                                            
       t(1)=tr1+xligcm-u(1,1)*xligcm-u(1,2)*yligcm-u(1,3)*zligcm
       t(2)=tr2+yligcm-u(2,1)*xligcm-u(2,2)*yligcm-u(2,3)*zligcm
       t(3)=tr3+zligcm-u(3,1)*xligcm-u(3,2)*yligcm-u(3,3)*zligcm

Here tr1, tr2 and tr3 are translations given by GRAMM and xligcm, yligcm and zligcm are Cartesian coordinates of the center of mass for the ligand in its original position.
"""
weights = {'C':12.011, 'N':14.007, 'O':15.999, 'S':32.06}

allowed_atoms = ['N','CA','C','O','OXT','CB',
                 'CZ','CZ1','CZ2','CE','CE1','CE2',
                 'CD','CD1','CD2','CG','CG1','CG2','CG3',
                 'CH2','OE','OE1','OE2','OD','OD1','OD2',
                 'OH','OG','OG1','OG2','NZ','NE','NE1','NE2',
                 'NH','NH1','NH2','ND','ND1','ND2',
                 'SG','SD','SE']

def ligand_center(atom_list, geometric=False):
    masses = []
    positions = [ [], [], [] ] # [ [X1, X2, ..] , [Y1, Y2, ...] , [Z1, Z2, ...] ]

    for atom in atom_list:
        if atom.get_id() in allowed_atoms: 
            masses.append(weights[atom.get_id()[0]])

        for i, coord in enumerate(atom.get_coord()): 
            positions[i].append(coord)

    if geometric:
        return [sum(coord_list)/len(masses) for coord_list in positions]
    else:       
        w_pos = [ [], [], [] ]
        for atom_index, atom_mass in enumerate(masses):
            w_pos[0].append(positions[0][atom_index]*atom_mass)
            w_pos[1].append(positions[1][atom_index]*atom_mass)
            w_pos[2].append(positions[2][atom_index]*atom_mass)

        return [sum(coord_list)/sum(masses) for coord_list in w_pos]

def get_rototranslation(data):
    result_batch = []
    with tf.Session(graph=get_rt, config=config) as sess:
        for pose, m1, m2 in data:
            r_mat, t_vec = sess.run([r, t], feed_dict = {gr:m1, gt:m2, lcm:cm})
            result_batch.append([pose, r_mat, t_vec])
    return result_batch

def rototranslate_coord(data, c_only=True):
    count = 0
    result_batch = []
    with tf.Session(graph=rt_comp, config=config) as sess:
        for res in data:
            mat = res[2]
            vec = res[3]
            for atom in str2[0]['B'][res[0]]:
                name = atom.get_id()
                coord = atom.get_coord()
                main_c = get_main_coord(str2[0]['B'][res[0]])
                if c_only and name != main_c: continue
                rt_coord = sess.run(rtcoord, feed_dict = {r_mat:mat, t_vec:vec, xyz:coord})
                result_batch.append([res[0], res[1], name, rt_coord])
                count += 1
                if count % (1000*res_num) == 0: print (count/res_num)
    return result_batch

def get_main_coord(res):
    atoms = [atom.get_id() for atom in res]
    if 'CB' in atoms: return 'CB'
    elif 'CA' in atoms: return 'CA'
    else: return None

def split_jobs(job_list, cores):
    batch_list = []
    bs = math.floor(len(job_list)/cores)
    for n in range(0,cores-1): batch_list.append(job_list[n*bs:(n+1)*bs])
    batch_list.append(job_list[bs*(cores-1):])
    return batch_list

    
if __name__ == "__main__":
    p = argparse.ArgumentParser(description = '- Plot PPV stats for a Cmap or a list of Cmaps')
    p.add_argument('-g', required= True, help='gramm output')
    p.add_argument('-c', required= True, help='trRosetta predictions')
    p.add_argument('-s1', required= True, help='structure file 1')
    p.add_argument('-s2', required= True, help='structure file 2')
    ns = p.parse_args()

    with np.load(ns.c) as data: cmap = data['dist'] 
    #get top contacts

    ##### parse structures #####
    p = PDBParser(QUIET=True)
    str1 = p.get_structure('', ns.s1)
    str2 = p.get_structure('', ns.s2)
    str3 = p.get_structure('', ns.s2)
    lig_res = Selection.unfold_entities(str2, 'R')
    lig_atoms = Selection.unfold_entities(str2, 'A')
    res_num = len(lig_res)
    ##### calculate lcm #####
    cm = np.array(ligand_center(lig_atoms), dtype=np.float32)

    ##### graph to elaborate GRAMM output #######################################
    with tf.Graph().as_default() as get_rt:                                     #
        with tf.name_scope('input0'):                                           #
            pi = tf.constant(math.pi/180.0)                                     #
            gr = tf.placeholder(dtype=tf.float32, shape=(3))                    #
            gt = tf.placeholder(dtype=tf.float32, shape=(3))                    #
            lcm = tf.placeholder(dtype=tf.float32, shape=(3))                   #
                                                                                #
        gr1 = pi*gr                                                             #
        sina = tf.math.sin(gr1[0])                                              #
        sinb = tf.math.sin(gr1[1])                                              #
        sing = tf.math.sin(gr1[2])                                              #
        cosa = tf.math.cos(gr1[0])                                              #
        cosb = tf.math.cos(gr1[1])                                              #
        cosg = tf.math.cos(gr1[2])                                              #
                                                                                #
        r = tf.stack([tf.stack([cosg*cosa,                                      #
                                cosg*sina,                                      #
                                -sing]),                                        #
                      tf.stack([sinb*sing*cosa-cosb*sina,                       #
                                sinb*sing*sina+cosb*cosa,                       #
                                cosg*sinb]),                                    #
                      tf.stack([cosb*sing*cosa+sinb*sina,                       #
                                cosb*sing*sina-sinb*cosa,                       #
                                cosg*cosb])])                                   #
                                                                                #
        r_lcm = tf.linalg.matvec(r, lcm)                                        #
        t = gt+(lcm-r_lcm)                                                      #

    ##### graph to roto-translate atom coordinates ##############################
    with tf.Graph().as_default() as rt_comp:                                    #
        with tf.name_scope('input1'):                                           #
            t_vec = tf.placeholder(dtype=tf.float32, shape=(3))                 #
            r_mat = tf.placeholder(dtype=tf.float32, shape=(3, 3))              #
            xyz = tf.placeholder(dtype=tf.float32, shape=(3))                   #
                                                                                #
        rtcoord = tf.math.add(t_vec, tf.linalg.matvec(r_mat, xyz))              #

    cores = mp.cpu_count()-1
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True

    ##### compute rototranslations #####
    mat_jobs = []
    poses = False
    for line in open(ns.g):
        if '[match]' in line: 
            poses = True
            continue
        if not poses: continue
        line = [float(el.strip()) for el in line.split() if el.strip() != '']
        grammr = np.array([line[2], line[3], line[4]], dtype=np.float32)
        grammt = np.array([line[5], line[6], line[7]], dtype=np.float32)
        mat_jobs.append([int(line[0]), grammr, grammt])

    pool = mp.Pool(processes = cores)
    job_list = split_jobs(mat_jobs, cores)
    results = pool.map(get_rototranslation, job_list)
    pool.close()
    pool.join()

    ##### rototranslate CBs coordinates #####
    rtlist = []
    for batch in results: rtlist.extend(batch) 
    res_jobs = [ [res.get_id()]+mat for mat in rtlist for res in lig_res ]

    pool = mp.Pool(processes = cores)
    job_list = split_jobs(res_jobs, cores)
    results = pool.map(rototranslate_coord, job_list)
    pool.close()
    pool.join()

    rtcdic = {}
    for batch in results:
        for res in batch:
            if res[1] not in rtcdic: rtcdic[res[1]]={}
            if res[0] not in rtcdic[res[1]]: rtcdic[res[1]][res[0]] = {}
            rtcdic[res[1]][res[0]][res[2]] = res[3]

    #io = PDBIO()
    #for name, coord in coords: str3[0]['B'][res][name].set_coord(coord)
    #io.set_structure(str3)
    #io.save('try.pdb')
        
        





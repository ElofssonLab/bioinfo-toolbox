#!/usr/bin/env python3
import argparse
from Bio.PDB.vectors import rotaxis, calc_angle, calc_dihedral
from Bio.PDB.Polypeptide import is_aa
from math import pi
import numpy as np
import scipy.stats as st
import pandas as pd

def virtual_cb_vector(residue):
    # get atom coordinates as vectors
    n = residue['N'].get_vector()
    c = residue['C'].get_vector()
    ca = residue['CA'].get_vector()
    # center at origin
    n = n - ca
    c = c - ca
    # find rotation matrix that rotates n -120 degrees along the ca-c vector
    rot = rotaxis(-pi*120.0/180.0, c)
    # apply rotation to ca-n vector
    cb_at_origin = n.left_multiply(rot)
    # put on top of ca atom
    cb = cb_at_origin + ca
    return cb


if __name__ == "__main__":
    arg_parser = argparse.\
        ArgumentParser(
                description="Convert a pdb/mcif to trRosetta distances/angles")

    #in_group = arg_parser.add_mutually_exclusive_group(required=True)
    #in_group.add_argument("-p", "--pdb_file", type=argparse.FileType('r'))
    #in_group.add_argument("-m", "--mmCIF_file", type=argparse.FileType('r'))

    arg_parser.add_argument("pdbA", type=str)
    arg_parser.add_argument("pdbB", type=str)
    arg_parser.add_argument("ispredA", type=str)
    arg_parser.add_argument("ispredB", type=str)
    arg_parser.add_argument("npz_name", type=str)
    #arg_parser.add_argument("-c", "--chain", type=str, default='A')
    # arg_parser.add_argument("-s", "--std", default=1, type=float,
    #                         help="Standard deviation in Ångström")
    args = arg_parser.parse_args()
    std = 1
    ispredA=pd.read_csv(args.ispredA,sep="\t",header=None)
    ispredB=pd.read_csv(args.ispredA,sep="\t",header=None)
    print (ispredA[5],ispredB[5])
    #sys.exit()
    #if args.pdb_file:
    from Bio.PDB.PDBParser import PDBParser
    bio_parser = PDBParser(PERMISSIVE=1)
    structure_fileA = args.pdbA
    structure_idA = args.pdbA[:-4]
    structure_fileB = args.pdbB
    structure_idB = args.pdbB[:-4]
    #else:
    #    from Bio.PDB.MMCIFParser import MMCIFParser
    #    bio_parser = MMCIFParser()
    #    structure_file = args.mmCIF_file
    #    structure_id = args.mmCIF_file.name[:-4]

    # Load structures
    structureA = bio_parser.get_structure(structure_idA, structure_fileA)
    structureB = bio_parser.get_structure(structure_idB, structure_fileB)

    # Get residues and length of proteins
    residuesA = []
    for chain in structureA[0]:
        for residue1 in structureA[0][chain.id]:
            if not is_aa(residue1):
                continue
            residuesA.append(residue1.get_resname())
    seqlenA = len(residuesA)
    residuesB = []
    for chain in structureB[0]:
        for residue1 in structureB[0][chain.id]:
            if not is_aa(residue1):
                continue
            residuesB.append(residue1.get_resname())
    seqlenB = len(residuesB)

    dist_rst=np.array([0.249,
                                0.,0.,0.,0.,0.,0.,
                                0.,0.,0.,0.,0.,0.03,
                                0.03,0.03,0.03,0.03,0.03,0.03,
                                0.03,0.03,0.03,0.03,0.03,0.03,
                                0.03,0.03,0.03,0.03,0.03,0.03,
                                0.03,0.03,0.03,0.03,0.03,0.03
    ],dtype=np.float32)

    # Setup bins and step for the final matrix
    DIST_STEP = 0.5
    OMEGA_STEP = 15
    THETA_STEP = 15
    PHI_STEP = 15

    z_per_bin = std/DIST_STEP
    z_step = z_per_bin/2
    angle_z_step = 1

    plen=seqlenA+seqlenB
    minvalue = 0.1
    dist_wanted_bins = (20 - 2)/DIST_STEP
    omega_wanted_bins = 360/OMEGA_STEP
    theta_wanted_bins = 360/THETA_STEP
    phi_wanted_bins = 180/PHI_STEP
    cumm_cutoff = 0.899
    # cb_lst = get_cb_coordinates(open(args.pdb_file, 'r'), "A")
    # contact_mat = get_cb_contacts(len(residues))
    dist_mat = np.full((plen, plen, 37), minvalue/36)
    omega_mat = np.full((plen, plen, 25), minvalue/24)
    theta_mat = np.full((plen, plen, 25), minvalue/24)
    phi_mat = np.full((plen, plen, 13), minvalue/12)

    dist_bins = [i for i in np.arange(2, 2 + (dist_wanted_bins)*0.5, 0.5)]
    omega_bins = [i for i in np.arange(-180, -180 +
                  (omega_wanted_bins)*OMEGA_STEP, OMEGA_STEP)]
    theta_bins = [i for i in np.arange(-180, -180 +
                  (theta_wanted_bins)*THETA_STEP, THETA_STEP)]
    phi_bins = [i for i in np.arange(0, (phi_wanted_bins)*PHI_STEP, PHI_STEP)]
    # print(dist_bins)
    # print(len(dist_bins))
    # print(omega_bins)
    # print(len(omega_bins))
    # sys.exit()
    dist_num_bins = len(dist_bins)
    omega_num_bins = len(omega_bins)
    theta_num_bins = len(theta_bins)
    phi_num_bins = len(phi_bins)

    # Iterate over all residues and calculate distances
    i = 0
    j = 0

    structure=structureA
    for chain in structureB[0]:
        chain.id = 'X'
        structure[0].add(chain)
    
    for chain in structure[0]: # We assume only one chain
        for residue1 in structure[0][chain.id]:
        # Only use real atoms, not HET or water
            if not is_aa(residue1):
                continue
            
            # If the residue lacks CB (Glycine etc), create a virtual
            if residue1.has_id('CB'):
                c1B = residue1['CB'].get_vector()
            else:
                c1B = virtual_cb_vector(residue1)
            
            j = 0
            for chain in structure[0]:
                for residue2 in structure[0][chain.id]:
                    symm = False
                    if not is_aa(residue2):
                        continue
                    # print(i,j)
                    if i == j:
                        dist_mat[i, j, 0] = 0.9
                        omega_mat[i, j, 0] = 0.9
                        theta_mat[i, j, 0] = 0.9
                        phi_mat[i, j, 0] = 0.9
                        j += 1
                        continue
                    
                    if i > j:
                        dist_mat[i, j] = dist_mat[j, i]
                        omega_mat[i, j] = omega_mat[j, i]
                        symm = True
                    # If the residue lacks CB (Glycine etc), create a virtual
                    if residue2.has_id('CB'):
                        c2B = residue2['CB'].get_vector()
                    else:
                        c2B = virtual_cb_vector(residue2)
                    ###############################################
                    dist = (c2B-c1B).norm()
                    
                    if dist > 20 or np.sign(i-seqlenA)!=np.sign(j-seqlenA):
                        dist_mat[i, j, 0] = 0.9
                        omega_mat[i, j, 0] = 0.9
                        theta_mat[i, j, 0] = 0.9
                        phi_mat[i, j, 0] = 0.9
                    else:
                        # Dist and omega are symmetrical and have already been copied
                        if not symm:
                            ix = np.digitize(dist, dist_bins)
                            cum_prob = 0
                            b_step = 0
                            while cum_prob < cumm_cutoff:
                                bin_prob = st.norm.cdf(b_step*-z_step) -\
                                    st.norm.cdf(-z_step*(1+b_step))
                                dist_mat[i, j, np.min([ix+b_step, dist_num_bins])]\
                                    += bin_prob
                                dist_mat[i, j, np.max([ix-b_step, 1])] += bin_prob
                                cum_prob += bin_prob*2
                                b_step += 1
                        ###############################################
                        # # Omega
                        c1A = residue1['CA'].get_vector()
                        c2A = residue2['CA'].get_vector()
                    
                        if not symm:
                            raw_omega = calc_dihedral(c1A, c1B, c2B, c2A)
                            omega = (raw_omega*180)/pi
                    
                            ix = np.digitize(omega, omega_bins)
                            cum_prob = 0
                            b_step = 0
                            while cum_prob < cumm_cutoff:
                                bin_prob = st.norm.cdf(b_step*-angle_z_step) -\
                                    st.norm.cdf(-angle_z_step*(1+b_step))
                                omega_mat[i, j, np.min([ix+b_step, omega_num_bins])]\
                                    += bin_prob
                                omega_mat[i, j, np.max([ix-b_step, 1])] += bin_prob
                                cum_prob += bin_prob*2
                                b_step += 1
                    
                        ###############################################
                        # # Theta
                        N1 = residue1['N'].get_vector()
                    
                        raw_theta = calc_dihedral(N1, c1A, c1B, c2B)
                        theta = (raw_theta*180)/pi
                    
                        ix = np.digitize(theta, theta_bins)
                        cum_prob = 0
                        b_step = 0
                        while cum_prob < cumm_cutoff:
                            bin_prob = st.norm.cdf(b_step*-angle_z_step) -\
                                st.norm.cdf(-angle_z_step*(1+b_step))
                            theta_mat[i, j, np.min([ix+b_step, theta_num_bins])]\
                                += bin_prob
                            theta_mat[i, j, np.max([ix-b_step, 1])] += bin_prob
                            cum_prob += bin_prob*2
                            b_step += 1
                    
                        ###############################################
                        # # Phi
                    
                        raw_phi = calc_angle(c1A, c1B, c2B)
                        phi = (raw_phi*180)/pi
                    
                        ix = np.digitize(phi, phi_bins)
                        cum_prob = 0
                        b_step = 0
                        while cum_prob < cumm_cutoff:
                            bin_prob = st.norm.cdf(b_step*-angle_z_step) -\
                                st.norm.cdf(-angle_z_step*(1+b_step))
                            phi_mat[i, j, np.min([ix+b_step, phi_num_bins])]\
                                += bin_prob
                            phi_mat[i, j, np.max([ix-b_step, 1])] += bin_prob
                            cum_prob += bin_prob*2
                            b_step += 1
                    j += 1
            i += 1


    # set all predicted surface residue pairs to be 12A apart
    for i in range(seqlenA):
        for j in range(seqlenB):
            if (ispredA[5][i]=="I"  and ispredB[5][j]=="I"   ):
                dist_mat[i, j+seqlenA] = dist_rst
                dist_mat[j+seqlenA,i] = dist_rst

    np.savez_compressed(args.npz_name,
                        dist=dist_mat,
                        omega=omega_mat,
                        theta=theta_mat,
                        phi=phi_mat)
    

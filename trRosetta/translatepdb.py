#!/usr/bin/env python3
import argparse
from Bio.PDB.vectors import rotaxis, calc_angle, calc_dihedral
from Bio.PDB.Polypeptide import is_aa
from math import pi
import math
import numpy  as np
from Bio.PDB import PDBIO

# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Module with assorted geometrical functions on
macromolecules.
"""

from Bio.PDB import Entity

def center_of_mass(entity, geometric=False):
    """
    Returns gravitic [default] or geometric center of mass of an Entity.
    Geometric assumes all masses are equal (geometric=True)
    """
    
    # Structure, Model, Chain, Residue
    if isinstance(entity, Entity.Entity):
        atom_list = entity.get_atoms()
    # List of Atoms
    elif hasattr(entity, '__iter__') and [x for x in entity if x.level == 'A']:
        atom_list = entity
    else: # Some other weirdo object
        raise ValueError("Center of Mass can only be calculated from the following objects:\n"
                            "Structure, Model, Chain, Residue, list of Atoms.")
    
    masses = []
    positions = [ [], [], [] ] # [ [X1, X2, ..] , [Y1, Y2, ...] , [Z1, Z2, ...] ]
    
    for atom in atom_list:
        
        masses.append(atom.mass)
        
        for i, coord in enumerate(atom.coord.tolist()):
            positions[i].append(coord)

    # If there is a single atom with undefined mass complain loudly.
    if 'ukn' in set(masses) and not geometric:
        raise ValueError("Some Atoms don't have an element assigned.\n"
                         "Try adding them manually or calculate the geometrical center of mass instead.")
    
    if geometric:
        return [sum(coord_list)/len(masses) for coord_list in positions]
    else:       
        w_pos = [ [], [], [] ]
        for atom_index, atom_mass in enumerate(masses):
            w_pos[0].append(positions[0][atom_index]*atom_mass)
            w_pos[1].append(positions[1][atom_index]*atom_mass)
            w_pos[2].append(positions[2][atom_index]*atom_mass)

        return [sum(coord_list)/sum(masses) for coord_list in w_pos]


    
if __name__ == "__main__":
    arg_parser = argparse.\
        ArgumentParser(                description="Translate and rotate a pdb file")
    in_group = arg_parser.add_mutually_exclusive_group(required=True)
    in_group.add_argument("-p", "--pdb_file", type=argparse.FileType('r'))
    in_group.add_argument("-m", "--mmCIF_file", type=argparse.FileType('r'))

    arg_parser.add_argument("-r", "--rotation", required=False, help="12 parameters for translation and rotation (Note you need a space before a - sign)",nargs=12,type=float)
    arg_parser.add_argument("-a", "--angles", required=False, help="6 parameters for translation and rotation (Note you need a space before a - sign)",nargs=6,type=float)

    arg_parser.add_argument("-o","--outfile", type=str)
    #arg_parser.add_argument("-c", "--chain", type=str, default='A')
    #arg_parser.add_argument("-s", "--std", default=1, type=float,
    #                         help="Standard deviation in Ångström")
    args = arg_parser.parse_args()
    std = 1

    if args.pdb_file:
        from Bio.PDB.PDBParser import PDBParser
        bio_parser = PDBParser(PERMISSIVE=1)
        structure_file = args.pdb_file
        structure_id = args.pdb_file.name[:-4]
    else:
        from Bio.PDB.MMCIFParser import MMCIFParser
        bio_parser = MMCIFParser()
        structure_file = args.mmCIF_file
        structure_id = args.mmCIF_file.name[:-4]

    # Load structure
    structure = bio_parser.get_structure(structure_id, structure_file)

    
    if (args.rotation):
        rotation_matrix = np.array((
            (args.rotation[3],args.rotation[4],args.rotation[5]),
            (args.rotation[6],args.rotation[7],args.rotation[8]),
            (args.rotation[9],args.rotation[10],args.rotation[11])
        ),'f')
        translation_matrix = np.array((args.rotation[0],args.rotation[1],args.rotation[2]),'f')
        structure.transform(rotation_matrix, translation_matrix)
    elif (args.angles):
        print ("Input",args.angles)
        alpha=args.angles[0]*pi/180
        beta=args.angles[1]*pi/180
        gamma=args.angles[2]*pi/180
        
        for chain in structure[0]:
            ca_atoms = [atom for atom in chain.get_atoms() if atom.name=="CA"]
            com=center_of_mass(ca_atoms)
            #print (com)
            break
        rotation_matrix=np.array((
            (math.cos(alpha)*math.cos(gamma),
             math.sin(alpha)*math.cos(gamma),
             -1*math.sin(gamma)),
            (math.cos(alpha)*math.sin(beta)*math.sin(gamma)-math.cos(beta)*math.sin(alpha),
             math.sin(alpha)*math.sin(beta)*math.sin(gamma)+math.cos(beta)*math.cos(alpha),
             math.cos(gamma)*math.sin(beta)),
            (math.cos(alpha)*math.cos(beta)*math.sin(gamma)+math.sin(beta)*math.sin(alpha),
             math.sin(alpha)*math.cos(beta)*math.sin(gamma)-math.sin(beta)*math.cos(alpha),
             math.cos(gamma)*math.cos(beta))
        ),'f')
        #rotation_matrix2=np.array((
        #    (math.cos(alpha)*math.cos(beta),
        #     math.cos(alpha)*math.sin(beta)*math.sin(gamma)-math.sin(alpha)*math.cos(gamma),
        #     math.cos(alpha)*math.sin(beta)*math.cos(gamma)+math.sin(alpha)*math.sin(gamma)),
        #    (math.sin(alpha)*math.cos(beta),
        #     math.sin(alpha)*math.sin(beta)*math.sin(gamma)+math.cos(alpha)*math.cos(gamma),
        #     math.sin(alpha)*math.sin(beta)*math.cos(gamma)-math.cos(alpha)*math.sin(gamma)),
        #    (-1*math.sin(beta),
        #     math.cos(beta)*math.sin(gamma),
        #     math.cos(beta)*math.cos(gamma))
        #),'f')
        print ("center of mass",com)
        translation_matrix = np.array((
            args.angles[3]+com[0]-rotation_matrix[0,0]*com[0]-rotation_matrix[0,1]*com[1]-rotation_matrix[0,2]*com[2],
            args.angles[4]+com[1]-rotation_matrix[1,0]*com[0]-rotation_matrix[1,1]*com[1]-rotation_matrix[1,2]*com[2],
            args.angles[5]+com[2]-rotation_matrix[2,0]*com[0]-rotation_matrix[2,1]*com[1]-rotation_matrix[2,2]*com[2]
        ),'f')
        pos=[0,0,0]
        
        #for atom in structure.get_atoms():
        #    cord=atom.get_vector()
        #    pos=[0,0,0]
        #    #print ("test-org",cord)
        #    for i in range(3):
        #        pos[i]=translation_matrix[i]+rotation_matrix[i,0]*cord[0]+rotation_matrix[i,1]*cord[1]+rotation_matrix[i,2]*cord[2]
        #        #print (i,cord[i],pos[i])
        #    #atom.transform(rotation_matrix2, translation_matrix)
        #    #print ("test-rot",atom.get_vector())
        #    print (atom.coord)
        #    atom.coord = np.dot(atom.coord, rotation_matrix.T) + translation_matrix
        #    print (atom,atom.coord,pos)
        #    atom.set_coord(pos)
        #for atom in structure.get_atoms():
        structure.transform(rotation_matrix.T, translation_matrix)

            
        #t(1)=tr1+xligcm-u(1,1)*xligcm-u(1,2)*yligcm-u(1,3)*zligcm
        #t(2)=tr2+yligcm-u(2,1)*xligcm-u(2,2)*yligcm-u(2,3)*zligcm
        #t(3)=tr3+zligcm-u(3,1)*xligcm-u(3,2)*yligcm-u(3,3)*zligcm
       
    else:
        sys.die("No rotation information")
        
    print (rotation_matrix)
    print (translation_matrix)


        
    #resl=0
    #   LIGTRANSFORM: do k=1,Natomlig
    #      xx=xlig_init(k); yy=ylig_init(k); zz=zlig_init(k);
    #      xlig_model(k)=t(1)+u(1,1)*xx+u(1,2)*yy+u(1,3)*zz
    #      ylig_model(k)=t(2)+u(2,1)*xx+u(2,2)*yy+u(2,3)*zz
    #      zlig_model(k)=t(3)+u(3,1)*xx+u(3,2)*yy+u(3,3)*zz
    #   end do LIGTRANSFORM


    io = PDBIO()
    io.set_structure(structure)
    io.save(args.outfile)

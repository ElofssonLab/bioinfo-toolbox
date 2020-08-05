#!/usr/bin/env python3
import argparse
from Bio.PDB.vectors import rotaxis, calc_angle, calc_dihedral
from Bio.PDB.Polypeptide import is_aa
from math import pi
import numpy 
from Bio.PDB import PDBIO


if __name__ == "__main__":
    arg_parser = argparse.\
        ArgumentParser(                description="Translate and rotate a pdb file")
    in_group = arg_parser.add_mutually_exclusive_group(required=True)
    in_group.add_argument("-p", "--pdb_file", type=argparse.FileType('r'))
    in_group.add_argument("-m", "--mmCIF_file", type=argparse.FileType('r'))

    arg_parser.add_argument("-r", "--rotation", required=True, help="12 parameters for translation and rotation (Note you need a space before a - sign)",nargs=12,type=float)

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

    rotation_matrix = numpy.array((
        (args.rotation[3],args.rotation[4],args.rotation[5]),
        (args.rotation[6],args.rotation[7],args.rotation[8]),
        (args.rotation[9],args.rotation[10],args.rotation[11])
        ),'f')
    translation_matrix = numpy.array((args.rotation[0],args.rotation[1],args.rotation[2]),'f')

    print (rotation_matrix)
    print (translation_matrix)
    for atom in structure.get_atoms():
        atom.transform(rotation_matrix, translation_matrix)



    io = PDBIO()
    io.set_structure(structure)
    io.save(args.outfile)

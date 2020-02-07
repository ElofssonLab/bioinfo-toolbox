"""\
Purpuse: With this module you can get PPV, print distances between residues
         and set a custom "cutoff". It's works on intra-chain contacts and
         with oligomer inter-chain contacts.

Preconditions:
                - Files for input
                    * Contact map
                    * Fasta file
                    * PDB file
                - Libraries
                    * Biopython
                    * Numpy

Arguments:    -Positional
                * [fasta_file]      // file
                * [contact_file]    // file
                * [pdb_file]        // file
              -Optional
                * "-d", "--cb_cutoff"    // float, default=8.0
                * "-o", "--outfile"      // string, default=""
                * "-f", "--factor"       // float, default=1.0
                * "-s", "--score"        // float, default=-1.0
                * "--chain1"             // string, default=""
                * "--chain2"             // string, default=""
                * "--noalign"            // boolean, action="store_true"
                * "--name"               // string, default=""
                * "--min_dist"           // int, default=5
                * "--print_dist"         // boolean, action="store_true"

Observations:
            - This is a first step on getting the real distances from
              each predicted contact. A lot of work need to be done
              from now.

            - PPV, TP and FP are not working yet.

            - If you don't put the "--name" flag, you won't save the result
              like this: (name, PPV, TP, FP)
              Otherwise, you will get this final function result:
              (pdb_filename, PPV, TP, FP)

            - If you put the "--outfile" flag when run ppv_with_olig.py,
              you need to add also the "--print_dist" flag, otherwise it won't
              active the printing function.


TODO: - help section on argparse
      - description on argparse
      - lot of stuff :)

"""

# import the module sys
import sys

# import the module argparse
import argparse

# import an all builtin module, math.
#from math import *

# import the module pairwise2 from package Biopython
from Bio import pairwise2

# import the package Numpy as "np"
import numpy as np

# From the folder "parsing" import three modules
from parsing import parse_contacts
from parsing import parse_fasta
from parsing import parse_pdb


def get_cb_contacts(gapped_cb_lst):
    """ Return a numpy matrix of all cb vs cb
    intra-chain distance contacts. """

    seqlen = len(gapped_cb_lst)
    dist_mat = np.zeros((seqlen, seqlen), np.float)
    dist_mat.fill(float("inf"))

    # Go simulstaneously over sequence's index and values like this:
    # for idx, val in enumerate(lst):
    for i, cb1 in enumerate(gapped_cb_lst):
        # In this case every element is an array
        if cb1[0] == "-":
            continue
        for j, cb2 in enumerate(gapped_cb_lst):
            # In this case every element is an array
            if cb2[0] == "-":
                continue
            diff_vec = cb1 - cb2
            dist_mat[i, j] = np.sqrt(np.sum(diff_vec * diff_vec))
    return dist_mat


def get_cb_contacts_PPI(gapped_cb_chain1_lst, gapped_cb_chain2_lst):
    """ Return a numpy matrix of all cb vs cb
    inter-chain distance contacts. """

    seqlen1 = len(gapped_cb_chain1_lst)
    seqlen2 = len(gapped_cb_chain2_lst)
    dist_mat_ppi = np.zeros((seqlen1, seqlen2), np.float)
    dist_mat_ppi.fill(float("inf"))

    # Go simulstaneously over sequence's index and values like this:
    # for idx, val in enumerate(lst):
    for i, cb1 in enumerate(gapped_cb_chain1_lst):
        # In this case every element is an array
        if cb1[0] == "-":
            continue
        for j, cb2 in enumerate(gapped_cb_chain2_lst):
            # In this case every element is an array
            if cb2[0] == "-":
                continue
            diff_vec = cb1 - cb2
            dist_mat_ppi[i, j] = np.sqrt(np.sum(diff_vec * diff_vec))
    return dist_mat_ppi


def print_distances(
    contacts_x, contacts_y, scores, dist_mat, atom_seq_ali_chain1, outfilename
):
    """ Return a file with all intra-chain distances
    and scores per pair of contact. """

    num_c = len(contacts_x)
    outstr = ""
    for i in range(num_c):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if atom_seq_ali_chain1:
            if atom_seq_ali_chain1[c_x] == "-":
                continue
            if atom_seq_ali_chain1[c_y] == "-":
                continue
        if outfilename:
            outstr += "%s %s %s %s\n" % (
                c_x + 1,
                c_y + 1,
                scores[i],
                dist_mat[c_x, c_y]
            )
        else:
            print(c_x + 1, c_y + 1, scores[i], dist_mat[c_x, c_y])
    if outfilename:
        with open(outfilename, "w") as outf:
            outf.write(outstr)
    return outfilename


def print_distances_PPI(
    contacts_x,
    contacts_y,
    scores,
    dist_mat_chain1_vs_chain2,
    atom_seq_ali_chain1,
    atom_seq_ali_chain2,
    outfilename
):
    """ Return a file with all inter-chain distances
    and scores per pair of contact. """

    num_c = len(contacts_x)
    outstr = ""
    for i in range(num_c):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if len(atom_seq_ali_chain1) == len(atom_seq_ali_chain2):
            if atom_seq_ali_chain1[c_x] == "-" or atom_seq_ali_chain2[c_x] == "-":
                continue
            if atom_seq_ali_chain1[c_y] == "-" or atom_seq_ali_chain2[c_y] == "-":
                continue
        if outfilename != "":
            # Append the 2 distances contact between A & B
            # that is: resA_chain1--->resB_chain2 and resA_chain2--->resB_chain1
            outstr += "%s %s %s %s\n" % (
                c_x + 1,
                c_y + 1,
                scores[i],
                dist_mat_chain1_vs_chain2[c_x, c_y],
            )
            outstr += "%s %s %s %s\n" % (
                c_y + 1,
                c_x + 1,
                scores[i],
                dist_mat_chain1_vs_chain2[c_y, c_x],
            )
        else:
            print(c_x + 1, c_y + 1, scores[i], dist_mat_chain1_vs_chain2[c_x, c_y])
            print(c_y + 1, c_x + 1, scores[i], dist_mat_chain1_vs_chain2[c_y, c_x])
    if outfilename != "":
        with open(outfilename, "w") as outf:
            outf.write(outstr)
    return outfilename


def get_separator(c_filename):
    """ Return a string separator.
    Guessing separator of constraint file. """

    with open(c_filename, "r") as contact_filename:
        for line in contact_filename:
            if len(line.split(",")) != 1:
                sep = ","
            elif len(line.split(" ")) != 1:
                sep = " "
            else:
                sep = "\t"
    return sep


def get_scores_from_contacts(c_filename, min_dist, factor_value, min_score, ref_len):
    """ Return a tupla unpacking of three lists,
    [contacts_x], [contacts_y], [scores]. """

    # get separator from c_filename
    sep = get_separator(c_filename)

    # get a list ranked predicted contacts for those
    # carbon-beta (CB) that are 5 residues separated.
    # In the function "parse_contacts.parse()", min_dist is 5 for default.
    # This returns parse_contacts.parse: [(score, resA_CB, resB_CB)]
    contacts = parse_contacts.parse(c_filename, sep, min_dist)

    # Build a list for each residue numparseber
    # e.g.: from resA_CB to resA_numN ---> contacts_x = [resA_num1, ..., resA_numN]
    # e.g.: from resB_CB to resB_numN ---> contacts_y = [resB_num1, ..., resB_numN]
    contacts_x = []
    contacts_y = []
    # Build a score list
    scores = []

    num_c = len(contacts)
    count = 0
    for i in range(num_c):
        score = contacts[i][0]
        # It use "- 1" because the calling gives the real biological position
        # of the residue and python start counting at zero.
        c_x = contacts[i][1] - 1
        c_y = contacts[i][2] - 1

        # Calculate the distance in the sequence position between resA and resB
        pos_diff = abs(c_x - c_y)
        # Boolean declaration with the distances.
        # Check if those are less than 5 residues far from each other.
        too_close = pos_diff < min_dist

        if not too_close:
            contacts_x.append(c_x)
            contacts_y.append(c_y)
            scores.append(score)
            count += 1

        # Check if the contact predicted is below than min_score
        # and the count are grater or equal than ref_len*factor_value (default, len*1.0)
        if min_score == -1.0 and count >= ref_len * factor_value:
            break
        if score < min_score:
            break
    return contacts_x, contacts_y, scores


def get_global_align_from_pdb(pdb_filename, chain, seq):
    """ Return a list with 1 tupla of 5 elements:
    [(pdb_aligned_seq,fasta_seq, float, float, int)]. """

    # Generate the atom sequence from input chain.
    # Default values in get_atom_seq(pdbfile, chain="", model=1, return_lines=False)
    atom_seq_chain = parse_pdb.get_atom_seq(pdb_filename, chain)

    # Align seq from fasta with seq from pdb
    # 2: match, -1: missmatch, -0.5: open gap, -0.1: extend gap
    # For Homo-oligomer we should use the two chains, e.g:
    # atom_seq_chain1 and atom_seq_chain2.
    # The result is a list with 1 tupla of 5 elements:
    # [(pdb_aligned_seq,fasta_seq, float, float, int)]
    align = pairwise2.align.globalms(atom_seq_chain, seq, 2, -1, -0.5, -0.1)

    return align


def get_gapped_cb_lts(pdb_filename, chain, seq, cb_lst):
    """ Return a list with gapped cb from
    a PDB_chain and fasta_seq. """

    # Get a list with 1 tupla of 5 elements:
    # [(pdb_aligned_seq,fasta_seq, float, float, int)]
    align = get_global_align_from_pdb(pdb_filename, chain, seq)
    # Take the sequence from the pdb to use it
    atom_seq_ali = align[-1][0]
    # Take the sequence from the fasta to use it
    seq_ali = align[-1][1]

    j = 0  # <------- What is it for this counter?
    # Create a gapped carbon beta list to be used in "get_cb_contacts()"
    gapped_cb_lst = []
    seqlen = len(atom_seq_ali)
    for i in range(seqlen):
        if atom_seq_ali[i] == "-":
            gapped_cb_lst.append(["-"])
        elif seq_ali[i] == "-":
            j += 1
            continue
        else:
            gapped_cb_lst.append(cb_lst[j])
            j += 1
    return gapped_cb_lst


def get_ppv_helper(contacts_x, contacts_y, ref_contact_map, atom_seq_ali):
    """ Return a tupla of tree numbers
    PPV, TP and FP. """

    num_c = len(contacts_x)
    TP = 0.0
    FP = 0.0
    PPV = 0.0
    for i in range(num_c):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if atom_seq_ali:
            # Check this error from below: TypeError: string indices must be integers
            # I think it is the "c_x" that could be a string and expect a interger.
            # Other check: see if atom_seq_ali is not a dict, otherwise should work fine!
            if atom_seq_ali[c_x] == "-":
                continue
            if atom_seq_ali[c_y] == "-":
                continue
        # ref_contact_map is a boolean matrix.
        # Should be a distance matrix instead?
        if (ref_contact_map[c_x, c_y] == True):
            TP += 1.0 / num_c
        else:
            FP += 1.0 / num_c
    # Compute the PPV calculation
    if TP > 0:
        PPV = TP / (TP + FP)

    return (PPV, TP, FP)


def get_ppv(
    fasta_filename,
    contact_filename,
    pdb_filename,
    factor_value,
    cb_cutoff,
    min_score,
    chain1,
    chain2,
    outfilename,
    name,
    noalign,
    min_dist,
    print_dist
):
    """ Return a tupla of 1 str and 3 floats,
    (pdb_filename, PPV, TP, FP). """

    # From a dictionary builded from a file, could be a simple fasta,
    # a3m, or a MSA, get first seq. The dictionary has this structure:
    # {header:[seq1,seq2,...,seqN], header2:[seq1,seq2,...,seqN]}    <--- Is it correct?
    # The KEYs are the query IDs.
    print("In get_ppv function")
    seq = list(parse_fasta.read_fasta(fasta_filename).values())[0][0]
    ref_len = len(seq)

    # Get all the scores from contacts that satisfies the arguments
    contacts_x, contacts_y, scores = get_scores_from_contacts(contact_filename, min_dist, factor_value, min_score, ref_len)

    # Create a carbon-beta list from chain1 like this:
    # [array_cb1([x1, y1, z1]), array_cb2([x2, y2, z2], ... array_cbN([xN, yN, zN])
    cb_chain1_lst = parse_pdb.get_cb_coordinates(pdb_filename, chain1)

##########################################
# Using Biopython to get the coordinates #
##########################################
# from Bio.PDB.PDBParser import PDBParser
# parser = PDBParser(PERMISSIVE=1) # See 11.7  Common problems in PDB files  --> http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc150 
# structure_id = "1bih"
# filename = "<PATH-TO-PDB-FILE>XXXX.pdb"
# structure = parser.get_structure(structure_id, filename)
# 
#    for model in structure.get_list():
#        for chain in model.get_list():
#            if chain == model["B"]:
#                for residue in chain.get_list():
#                    if residue.has_id("CB"):
#                        cb = residue["CB"]
#                        print(cb.get_coord())
#
#[46.114 29.797 48.287]


    # Here "noalign" is always setting FALSE in function definition. Why?
    if noalign:
        dist_mat = get_cb_contacts(cb_chain1_lst)
        #       cb_cutoff = 8
        # Check if those are less than "cb_cutoff" angstrom far from each other.
        ref_contact_map = dist_mat < cb_cutoff  # cb_cutoff is 8 in default mode
        # Get all the PPV, TP and FP results.
        PPV, TP, FP = get_ppv_helper(
            contacts_x, contacts_y, ref_contact_map, atom_seq_ali=[]
        )
    else:

        # Check if PPI is need it
        if (chain2 != "" and chain2 != chain1):
            print("In the PPI branch of get_ppv")
            # Create a carbon-beta list from chain2
            cb_chain2_lst = parse_pdb.get_cb_coordinates(pdb_filename, chain2)

            gapped_cb_chain1_lst = get_gapped_cb_lts(pdb_filename, chain1, seq, cb_chain1_lst)
            gapped_cb_chain2_lst = get_gapped_cb_lts(pdb_filename, chain2, seq, cb_chain2_lst)

            # Get the distance matrix from chain1 only.
            # I do not use this distance matrix but could
            # be useful if intra-chain is also need it to be printed.
            # dist_mat_chain1 = get_cb_contacts(gapped_cb_chain1_lst)

            # Get the distance matrix from chain1 vs chain2. Used in PPI
            dist_mat_chain1_vs_chain2 = get_cb_contacts_PPI(gapped_cb_chain1_lst, gapped_cb_chain2_lst)

            #            cb_cutoff = 8
            # Check if those are less than "cb_cutoff" angstrom far from each other.
            # This create a boolean matrix called: ref_contact_map
            # cb_cutoff is 8 in default mode
            ref_contact_map = dist_mat_chain1_vs_chain2 < cb_cutoff

            # Get atoms seq aligned from a PDB_chain.
            # atom_seq_ali it is a string.
            atom_seq_ali_chain1 = get_global_align_from_pdb(pdb_filename, chain1, seq)[-1][0]
            atom_seq_ali_chain2 = get_global_align_from_pdb(pdb_filename, chain2, seq)[-1][0]

            ###################################################
            ##  Which atom_seq_ali_chainX we have to use for ##
            ##  get_ppf_helper() ? Both? or that one with    ##  <-----  LOOK !
            ##  the best aligned pdb sequence to fasta_seq ? ##
            ###################################################
            print("Getting PPV, TP and FP values...")
            # Get all the PPV, TP and FP results.
            PPV, TP, FP = get_ppv_helper(
                contacts_x,
                contacts_y,
                ref_contact_map,
                atom_seq_ali_chain1
            )

            # Check if print is need it.
            if print_dist:
                print("Printing PPI's distance results...")
                print_distances_PPI(
                    contacts_x,
                    contacts_y,
                    scores,
                    dist_mat_chain1_vs_chain2,
                    atom_seq_ali_chain1,
                    atom_seq_ali_chain2,
                    outfilename
                )

        else:
            print("In the monomer branch of get_ppv")
            gapped_cb_chain1_lst = get_gapped_cb_lts(pdb_filename, chain1, seq, cb_chain1_lst)

            # Get the distance matrix from chain1 only. Could be useful if
            # intra-chain is also need it to be printed.
            dist_mat = get_cb_contacts(gapped_cb_chain1_lst)

            #            cb_cutoff = 8
            # Check if those are less than "cb_cutoff" angstrom far from each other.
            # This create a boolean matrix called: ref_contact_map
            # cb_cutoff is 8 in default mode.
            ref_contact_map = dist_mat < cb_cutoff

            # Get atoms seq aligned from a PDB_chain.
            # atom_seq_ali it is a string.
            atom_seq_ali_chain1 = get_global_align_from_pdb(pdb_filename, chain1, seq)[-1][0]

            # Get the PPV, TP and FP results.
            PPV, TP, FP = get_ppv_helper(
                contacts_x,
                contacts_y,
                ref_contact_map,
                atom_seq_ali_chain1
            )

            # Check if print is need it.
            if print_dist:
                print_distances(
                    contacts_x,
                    contacts_y,
                    scores,
                    dist_mat,
                    atom_seq_ali_chain1,
                    outfilename
                )

    # Here "name" is always empty (so False) by default.
    if name:
        print("%s\n" % ("----------------------------------"))
        print("%s %s %s %s" % (name, PPV, TP, FP))
    else:
        print("Finished")
#        print("%s %s %s %s %s" % (fasta_filename, contact_filename, PPV, TP, FP))
    return (pdb_filename, PPV, TP, FP)


if __name__ == "__main__":

    p = argparse.ArgumentParser(description="Plot protein residue contact maps.")
    # Positional arguments
    p.add_argument("fasta_filename", help="Path to Fasta file.")  # , required=True)
    p.add_argument("contact_filename", help="Path to Contact file.")  # , required=True)
    p.add_argument("pdb_filename", help="Path to PDB file.")
    # Optional arguments
    p.add_argument("-d", "--cb_cutoff", default=8.0, type=float, help="Bla bla")
    p.add_argument("-f", "--factor_value", default=1.0, type=float, help="Bla bla")
    p.add_argument("-o", "--outfilename", default="", help="Bla bla")
    p.add_argument("-s", "--min_score", default=-1.0, type=float, help="Bla bla")
    p.add_argument("--chain1", default="", help="Bla bla")
    p.add_argument("--chain2", default="", help="Bla bla")
    p.add_argument("--noalign", action="store_true", help="Bla bla")
    p.add_argument("--name", default="", help="Bla bla")
    p.add_argument("--min_dist", default=5, type=int, help="Bla bla")
    p.add_argument("--print_dist", action="store_true", help="Bla bla")

#    args = vars(p.parse_args(sys.argv[1:]))

    # From:
    # https://www.peterbe.com/plog/vars-argparse-namespace-into-a-function
    # https://stackoverflow.com/questions/35822477/unpacking-arguments-from-argparse
    args = p.parse_args()
    get_ppv(**vars(args))

#    get_ppv(
#        fasta_filename=args["fasta_file"],
#        c_filename=args["contact_file"],
#        pdb_filename=args["pdb"],
#        factor_value=args["factor"],
#        cb_cutoff=args["cb_cutoff"],
#        chain1=args["chain1"],
#        chain2=args["chain2"],
#        outfilename=args["outfile"],
#        noalign=args["noalign"],
#        min_score=args["score"],
#        name=args["name"],
#        min_dist=args["min_dist"],
#        print_dist=args["print_dist"])

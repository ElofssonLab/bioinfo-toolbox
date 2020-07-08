#!/usr/bin/env python3
import numpy as np
import os
import sys
from collections import namedtuple


cwd = os.getcwd()
if len(sys.argv) != 3:
    if len(sys.argv) != 4:
        print("Usage: {} trRosetta_contacts.npz seq.fasta".format(os.path.basename(__file__)))
        sys.exit()
    else:
        fasta2_path = os.path.join(cwd,sys.argv[3])
npz_path = os.path.join(cwd,sys.argv[1])
fasta_path = os.path.join(cwd,sys.argv[2])

Bin_values = namedtuple("Bin_values", ["bin_step", "min_bin_value", "max_bin_value", "ending"])
default_bin_values = {"dist":  Bin_values(0.5, 2, 16, ".rr"),     # In Angstrom
                      "omega": Bin_values(15, -180, 180, ".omega"),  # Dihedral angle in degrees
                      "theta": Bin_values(15, -180, 180, ".theta"),  # Dihedral angle in degrees
                      "phi":   Bin_values(15,    0, 180, ".phi")}  # Planar angle in degrees

# fasta_path = os.path.join(cwd,sys.argv[2])
# 
# try:
#     input_file = np.load(npz_path)
# except:
#     print("Input file {} not found or not a npz file".format(npz_path), file=sys.stderr)
#     sys.exit()
# try:
#     seq = ''.join(open(fasta_path).read().split('\n')[1:]).strip()
# except:
#     print("Fasta file {} not found or could not be opened", file=sys.stderr)
#     sys.exit()
# 
out_base_path = npz_path[:-4]

def npz_to_casp(input_file, info_key="dist", fasta_file=None,  fasta2_file=None, out_base_path="",
                min_sep=0, pthres=0.5):
    """
    Convert a trRosetta npz-file into casp formated restraints file.
    """

    try:
        with np.load(input_file) as npz_file:
            raw_data = npz_file[info_key]
    except:
        print("Error reading npz file: " + input_file)
        sys.exit()

    target = '.'.join(input_file.split('/')[-1].split('.')[:-1])

    header = '\n'.join(("PFRMAT RR",
                       "TARGET {}".format(target),
                       "AUTHOR pyCONFOLD",
                       "METHOD trRosetta_contact",
                       "REMARK {}".format(info_key),
                       "MODEL 1\n"))
    # Length of the protein
    nres = raw_data.shape[0]

    # If a fasta file is supplied, add in the sequence to the header
    if fasta_file:
        n=50
        try:
            with open(fasta_file) as fasta_handle:
                fasta_seq = ''.join(fasta_handle.read().split('\n')[1:]).strip()
        except:
            print("Error reading fasta file: "+fasta_file)
            sys.exit()

        if fasta2_file:
            try:
                with open(fasta2_file) as fasta2_handle:
                    fasta_seq += ''.join(fasta2_handle.read().split('\n')[1:]).strip()
            except:
                print("Error reading fasta file: "+fasta_file)
                sys.exit()


        # Assert the fasta sequence is of the same length as the resulting contacts results
        fasta_length = len(fasta_seq)
        assert nres == fasta_length, "Sequence from fasta file is not of the same length as" +\
                                       " the predictions, {} vs {}".format(fasta_length, nres)

        header += '\n'.join(fasta_seq[i:i+n] for i in range(0, len(fasta_seq), n))
        header += '\n' 

    bin_values= default_bin_values[info_key] 

    # Number of bins, minus one for the first contact binary bin
    nbins = raw_data.shape[-1] - 1

    # # How many of the bins do we want?
    wanted_bins = int((bin_values.max_bin_value - bin_values.min_bin_value)/bin_values.bin_step)
    if wanted_bins > nbins:
        print("Number of wanted bins are higher than existing bins")
        sys.exit()
    # Generate up the midpoint of each distance bin to be able to
    # calculate correct distances. Measure in Angstrom
    bins = np.array([(bin_values.min_bin_value+bin_values.bin_step/2)+bin_values.bin_step*i for i in range(wanted_bins)])

    # Pull out the predictions, first array element is 
    # probability of non-contact, only pull out as many bins as we want
    # [L, L, nbins] -> [L, L, wanted_bins]
    raw_predictions = raw_data[:,:,1:wanted_bins+1]

    # What is the total probability of contact?
    contact_probabilities = np.sum(raw_predictions, axis=-1)

    # Normalize all probabilities by dividing by the sum of the raw_predictions
    # across the last dimension, [L, L, nbins] -> [L, L, nbins]
    norm_predictions = np.divide(raw_predictions, np.sum(raw_predictions, axis=-1)[:,:,None])

    # Multiply normalized predictions with the bin (sizes) to get the correct value distribution
    values = np.multiply(norm_predictions, bins[None, None, :])

    # Make some assertions to confirm array operations har occured correctly
    # Has the normalization worked?
    assert np.array_equal(np.divide(raw_predictions[0, 5, :],np.sum(raw_predictions[0,5,:])), norm_predictions[0,5])
    # Has the multiplication with the bins worked?
    assert np.array_equal(np.multiply(np.divide(raw_predictions[0, 5, :],np.sum(raw_predictions[0, 5, :])), bins),
                          values[0,5, :])

    # Calculate mean value and errors
    mean_values = np.sum(values, axis=-1)
    variances = np.sum(
                    np.multiply(
                        np.square(
                            np.subtract(bins[None,None,:], mean_values[:,:,None])),
                            norm_predictions), axis=-1)
    sd = np.sqrt(variances)


    # Assert summation is correct
    assert np.array_equal(np.sum(values[0, 5, :]), mean_values[0, 5])
    # Asser variance is correct
    assert np.array_equal(np.sum((mean_values[0, 5]-bins)**2*norm_predictions[0,5]), variances[0,5])

    line = "{i} {j} 0 {m:.8f} {c:.8f} {sd:.8f}\n"
    data = []

    for i in range(nres-1):
        for j in range(i+1, nres):
            if info_key == "dist":
                if (abs(i - j) > min_sep) and (contact_probabilities[i][j] > pthres):  # and (fasta_seq[i] != 'G') and (fasta_seq[j] != 'G'):
                    data.append((i+1, j+1, mean_values[i,j], contact_probabilities[i,j], sd[i,j]))
            else:
                if (abs(i - j) > min_sep) and (contact_probabilities[i][j] > pthres):  # and (fasta_seq[i] != 'G') and (fasta_seq[j] != 'G'):
                    data.append((i+1, j+1, mean_values[i,j], contact_probabilities[i,j], sd[i,j]))
                if info_key == "theta":  # Theta is assymmetric, add the other side of the diagonal if applicable
                    if (abs(j - i) > min_sep) and (contact_probabilities[j][i] > pthres):  # and (fasta_seq[i] != 'G') and (fasta_seq[j] != 'G'):
                        data.append((j+1, i+1, mean_values[i,j], contact_probabilities[j,i], sd[j,i]))

    # Sort on the probability of contact
    data.sort(key=lambda x: x[3], reverse=True)

    content = [header]

    # Format each line according to CASP rr format
    for i, j, m, c, sd in data:
        content.append(line.format(i=i, j=j, m=m, c=c, sd=sd))
    content.append("END\n")

    with open(out_base_path + bin_values.ending, 'w') as contacts_handle:
        contacts_handle.write(''.join(content))
    # return ''.join(content)




if __name__ == "__main__":
    # Distances
    try:
        npz_to_casp(npz_path, "dist",  fasta2_file=fasta2_path,fasta_file=fasta_path, out_base_path=out_base_path)
        npz_to_casp(npz_path, "omega", fasta2_file=fasta2_path,fasta_file=fasta_path, out_base_path=out_base_path)
        npz_to_casp(npz_path, "theta", fasta2_file=fasta2_path,fasta_file=fasta_path, out_base_path=out_base_path)
        npz_to_casp(npz_path, "phi",   fasta2_file=fasta2_path,fasta_file=fasta_path, out_base_path=out_base_path)
    except:
        npz_to_casp(npz_path, "dist", fasta_file=fasta_path, out_base_path=out_base_path)
        npz_to_casp(npz_path, "omega", fasta_file=fasta_path, out_base_path=out_base_path)
        npz_to_casp(npz_path, "theta", fasta_file=fasta_path, out_base_path=out_base_path)
        npz_to_casp(npz_path, "phi", fasta_file=fasta_path, out_base_path=out_base_path)


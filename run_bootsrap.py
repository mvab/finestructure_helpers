#!/usr/bin/python3

"""
Tool to create alternative painting vector via bootstrapping - sampling the data with replacements.
The input is the raw chunklengths files and the donor file for a given panel.
The output is a 1001 row matrix, where row 1 is a painting vector (samples' chunklengths summed by
populations over all chromosomes), and rows 2-1001 are the bootstrapped versions of the
painting vector.

Example usage:
$ python3 run_bootstrap.py donorids.txt Europe/ 
"""
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import collections
import os.path
import random
import numpy as np


def parse_args():
    """Parse command line arguments"""
    parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('donor_ids_file', help="Column file with samples IDs and population IDs")
    parser.add_argument('sample_dir', help="Full path to directory containing chunklengths files")
    parser.add_argument('output_directory', nargs='?', default='./', help="Directory to write to")
    args = parser.parse_args()
    return args


def main():
    """ Main: create donor dictionary ; run bootstraps; save matrix to file """
    args = parse_args()

    # check inputs
    donorfile = args.donor_ids_file
    if not os.path.isfile(donorfile):
        print("Either file is missing or is not readable:", donorfile)
    outdir = args.output_directory
    sampledir = args.sample_dir

    # creating donor file  dictionary
    pop_dict = {}
    id_lines = [s.strip() for s in open(donorfile)]
    for line in id_lines:
        [sid, pid] = line.split()
        pop_dict.setdefault(pid, [])
        pop_dict[pid].append(sid)

    bs_matrix, popnames = chunklengths_to_bootstrap(sampledir, pop_dict, bootstrap_size=1000)

    # create header if want to output matrix as a tab-delimited file
    matrix_header = '\t'.join(popnames)
    np.savetxt(outdir + "bootstrap_matrix.txt", bs_matrix, fmt='%4f',
               delimiter='\t', header=matrix_header, comments='')


def chunklengths_to_bootstrap(sampledir, pop_dict, bootstrap_size=1000):
    """Create a bootstrap matrix with bootstrap_size rows,
        and columns equal to the number of donor populations"""

    chrs = list(range(1, 23))

    # creating a nested dictionary -
    #  chr: {{SID:value}, {SID:value}, ...}, chr: {{SID:value}, {SID:value}, ...}, ...
    allchr_dict = collections.defaultdict(dict)
    for i in chrs:
        chunkfile = sampledir + "/chr" + str(i) + ".chunklengths.out"
        if not os.path.isfile(chunkfile):
            print("Either file is missing or is not readable:", chunkfile)
            break
        allchr_dict[i] = chunkfile_to_dict(chunkfile)

    # iterate over every chr dict and apply make_donor_dict; save in master_donor_dict
    master_donor_dict = collections.defaultdict(dict)
    master_pop_chunk = collections.defaultdict(dict)
    for i in chrs:
        # for every chromosome chunk file create a pop:samples dictionary
        master_donor_dict[i] = make_population_chunkdata_dict(pop_dict, allchr_dict[i])
        # for every chromosome dictionary, sum chunklengths within populations
        master_pop_chunk[i] = dict((key, sum(vals)) for key, vals in master_donor_dict[i].items())

    # run bootstrap

    random.seed(1)
    bootstrap_size += 1  # number of bootstraps
    # (increase by 1 to make while loop do the exact specified number)
    # current bs; set to 1 to easier keep track
    current_bs = 1
    # matrix to hold the bs iterations
    bs_matrix = np.empty((0, len(master_pop_chunk[1].keys())), float)
    while current_bs < bootstrap_size:
        print("Bootstrap iteration " + str(current_bs) + " out of " + str(bootstrap_size - 1))
        # sample 22 chromosomes  with replacement
        chr_sample = random.choices(chrs, k=22)
        # create a single bs array
        bs_iter = bootstrap_chunkdata(master_pop_chunk, chr_sample)
        # append array to matrix
        bs_matrix = np.append(bs_matrix, bs_iter, axis=0)
        current_bs += 1
    print("\nBootstrap matrix shape after all iterations:", str(bs_matrix.shape), "\n")

    # create the combined chunklengths per population
    print("Adding the combined non-bootstrapped vector on top:")
    combined = make_painting_vector(master_pop_chunk)
    print(combined)
    # add to the bootstrap matrix as the first row
    bs_matrix = np.append(combined, bs_matrix, axis=0)
    print(bs_matrix.shape)

    # get population names (the order corresponds to columns in the matrix)
    popnames = list(master_pop_chunk[1].keys())

    return bs_matrix, popnames


def bootstrap_chunkdata(master_pop_chunk, chr_sample):
    """Create a single bootstap iteration output for population chunks"""

    chrs = list(range(1, 23))

    # get chunklength sums per chromosome
    allsums = []
    for i in chrs:
        chr_sum = sum(list(master_pop_chunk[i].values()))
        allsums.append(chr_sum)

    # bs array to fill
    bs_array = np.zeros((1, len(master_pop_chunk[1].keys())), float)

    # iterate simultaneously over lists of
    #  1) a randomly sampled chromosomes and 2) chromosomes in 1:22 order
    for rnd, chrm in list(zip(chr_sample, chrs)):
        # print (rnd, chr)

        # resample the current random chromosome
        chr_list = list(master_pop_chunk[rnd].values())
        chr_resampled = []

        for val in chr_list:
            chr_resampled.append(val * allsums[chrm - 1] / allsums[rnd - 1])
            # print("Looked at chrm"+ str(rnd) )
            # print( " ... multiplied by " + str(allsums[chrm-1]) )
            # print(" ... divided by "+ str(allsums[rnd-1]) )
            # Uncomment the print statement above to see the logic

        # add up all resampled chromosomes
        bs_array = np.sum([bs_array, chr_resampled], axis=0)
    return bs_array


def make_painting_vector(master_pop_chunk):
    """Make painting vector aka combined chunklengths file per population"""
    paintingvector = np.zeros((1, len(master_pop_chunk[1].keys())), float)
    # iterate over values in the dict, and sum them as vectors
    for i in list(range(1, 23)):
        my_chr = list(master_pop_chunk[i].values())
        paintingvector = np.sum([paintingvector, my_chr], axis=0)
    return paintingvector


def chunkfile_to_dict(chunkfile):
    """Read in a chunkfile (counts or lengths), ignore #header (if present),
       assume exactly 2 rows (recipients, samples), and create a dictionary
       with key:value as sample_id: length """

    lines = [s.strip() for s in open(chunkfile) if not s.startswith('#')]
    sample_ids = lines[0].split()[1:]
    values = [float(x) for x in lines[1].split()[1:]]

    return dict(zip(sample_ids, values))


def make_population_chunkdata_dict(popdict, chunkdict):
    """Take population dict (PID: [list of SIDs] and chunkdict (SID: length)
    and output a new dictionary - PID: [list of lengths for that PID's samples) """

    population_chunkdict = {}  # dictionary to keep chunk data per population (PIDs as keys)

    for pid in popdict.keys():  # take a population
        for sid in popdict[pid]:  # iterate over all samples in it

            if sid in chunkdict:  # take that sample in the current chunkdict
                if pid in population_chunkdict:  # check if it's already in the dict
                    # update this PID's value list
                    population_chunkdict[pid].append(float(chunkdict[sid]))
                else:
                    # create new key
                    population_chunkdict[pid] = [float(chunkdict[sid])]
            else:
                print("Sample " + sid + " not found in chunkdict")

    # remove populations whose chunklengths add up to zero
    for key in list(population_chunkdict.keys()):
        if sum(population_chunkdict[key]) == 0:
            del population_chunkdict[key]

    return population_chunkdict


if __name__ == '__main__':
    main()

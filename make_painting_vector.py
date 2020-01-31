#!/usr/bin/python3

"""
Example usage:
$ python3 make_painting_vector.py donorids.txt Europe/
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import collections
import os.path
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

    painting_vector, popnames = chunklengths_to_vector(sampledir, pop_dict)

    # create header if want to output matrix as a tab-delimited file
    matrix_header = '\t'.join(popnames)
    np.savetxt(outdir + "painting_vector.txt", painting_vector, fmt='%4f',
               delimiter='\t', header=matrix_header, comments='')


def chunklengths_to_vector(sampledir, pop_dict):
    """Prepare chunklengths to create a painting vector """

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

    combined = make_painting_vector(master_pop_chunk)
    print(combined)

    popnames = list(master_pop_chunk[1].keys())
    return combined, popnames


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

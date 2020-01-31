#!/usr/local/bin/python3

"""
Python command line chunks data parser into JSON format.
Example usage:
$ python3 convert_to_json.py donorids.txt combined.chunkcounts.out combined.chunklengths.out my_chunk_file
"""

import os.path
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import numpy as np
import json
__author__ = 'Marina Vabistsevits'


def parse_args():
    """Parse command line arguments"""
    parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('donor_ids_file', help="Column file with samples IDs (SID) and population IDs (PID)")
    parser.add_argument('chunkcounts', help="Text file of chunkcounts per sample")
    parser.add_argument('chunklengths', help="Text file of chunklengths per sample")
    parser.add_argument('out_file_prefix', nargs='?', default='chunksdata',
                        help="Output JSON file prefix (without .json)")
    parser.add_argument('output_directory', nargs='?', default='./', help="Directory to write to")
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    # check inputs
    ids_file = args.donor_ids_file
    if not os.path.isfile(ids_file):
        print("Either file is missing or is not readable:", ids_file)
    counts_file = args.chunkcounts
    if not os.path.isfile(counts_file):
        print("Either file is missing or is not readable:", counts_file)
    lenghts_file = args.chunklengths
    if not os.path.isfile(lenghts_file):
        print("Either file is missing or is not readable:", lenghts_file)
    outdir = args.output_directory
    outfile = args.out_file_prefix

    # creating PID : SIDs  dictionary
    ids_array = np.genfromtxt(ids_file, dtype='str')
    pop_dict = {}
    for sample in ids_array:
        sid = str(sample[0])
        pid = str(sample[1])
        if pid in pop_dict:  # fill in the dictionary FID : SIDs
            pop_dict[pid].append(sid)
        else:
            pop_dict[pid] = [sid]

    # create SID : count and SID : length dictionaries
    counts_dict = chunkfile_to_dict(counts_file)
    lenghts_dict = chunkfile_to_dict(lenghts_file)

    # check that the number of samples (donors (SID) is the same in all 3 supplied files)
    if len(lenghts_dict.keys()) != len(counts_dict.keys()):
        print("The number of donor samples is different between chunk files!")
    elif len(lenghts_dict.keys()) != ids_array.shape[0]:
        print("The number of donor samples in ids file does not match lenghts file!")
    elif len(counts_dict.keys()) != ids_array.shape[0]:
        print("The number of donor samples in ids file does not match counts file!")

    # create per PID counts/lenghts dictionaries
    counts_dict_values = make_population_chunkdata_dict(pop_dict, counts_dict)
    lenghts_dict_values = make_population_chunkdata_dict(pop_dict, lenghts_dict)

    # creating nested dict and saving json
    chunks_dict_output = {"chunk_counts": counts_dict_values, "chunk_lengths": lenghts_dict_values}

    with open(str(outdir + outfile + '.json'), 'w') as f:
        json.dump(chunks_dict_output, f, ensure_ascii=False)


# FUNCTIONS


def make_population_chunkdata_dict(popdict, chunkdict):
    """Take population dict (PID: [list of SIDs] and chunkdict (SID: count/length)
    and output a new dictionary - PID: [list of counts/lenghts for that PID's samples) """

    population_chunkdict = {}  # dictionary to keep chunk data per population (PIDs as keys)

    for pid in popdict.keys():  # take a population
        for sid in popdict[pid]:  # iterate over all samples in it

            if sid in chunkdict:  # take that sample in the current chunkdict
                if pid in population_chunkdict:  # check if it's already in the dict
                    population_chunkdict[pid].append(float(chunkdict[sid]))  # update this PID's value list
                else:
                    population_chunkdict[pid] = [float(chunkdict[sid])]  # create new key
    return population_chunkdict


def chunkfile_to_dict(chunkfile):
    """Read in a chunkfile (counts or lengths), ignore #header (if present),
       transpose to be 2-column array (columns: recipient, sample),
        create a dictionary with key:value as SID: count/length """

    with open(chunkfile) as f:
        lines = (line for line in f if not line.startswith('#'))
        chunk_array = np.loadtxt(lines, delimiter=' ', dtype=object)
        chunk_array = chunk_array[:, 1:]  # remove 'Sample1' and 'Recipient' column
        chunk_array_t = np.transpose(chunk_array)

    chunk_dict = dict(
        zip(chunk_array_t[:, 0], chunk_array_t[:, 1]))  # key is column 0 (SID), value is column 1 (count/length)
    return chunk_dict


if __name__ == '__main__':
    main()

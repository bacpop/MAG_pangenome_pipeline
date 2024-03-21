import argparse
import json
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def remove_ranges_from_sequence(sequence, ranges, length_list):
    # Join the sequences into one string
    full_sequence = list(''.join(sequence))

    # Remove the specified ranges
    for start, end in ranges:
        full_sequence[start:end] = ["0"] * (end - start)

    # split back up into original contigs
    list_start = 0
    result = []
    for length in length_list:
        to_append = ''.join(full_sequence[list_start: list_start + length])
        to_append = to_append.replace("0", "")
        result.append(to_append)
        list_start += length

    return result

def get_options():
    description = 'Removes COGs from presence/absence matrix based on epirical distribution of completeness.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python remove_COGs')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--input',
                    default=None,
                    required=True,
                    help='List of genome fasta files to remove sequence from. One per file. ')
    IO.add_argument('--dist',
                    default=None,
                    required=True,
                    help='Genome completeness distribution for sampling. Must be single column of floats, one per line. ')
    IO.add_argument('--prop-complete',
                    default=0,
                    type=float,
                    help='Proportion of complete genomes between 0-1. Default = 0. ')
    IO.add_argument('--avg-breaks',
                    type=int,
                    default=1,
                    help='Average number of sections to be removed from genomes. Must be >=1. Default = 1')
    IO.add_argument('--outdir',
                    default="augmented_genomes",
                    help='Output directory. Default = "augmented_genomes"')

    return parser.parse_args()

def main():
    options = get_options()
    input = options.input
    prop_complete = options.prop_complete
    dist_file = options.dist
    outdir = options.outdir
    avg_breaks = options.avg_breaks

    #testing
    # input = "data/input.txt"
    # prop_complete = 0.1
    # dist_file = "data/completeness_distribution.txt"
    # outdir = "test_output"
    # avg_breaks = 10

    # set average breaks to one less than specified to ensure that the number of blocks made is equal to that specified
    avg_breaks -= 1
    if avg_breaks < 0:
        avg_breaks = 0

    # create output directory
    try:
        os.mkdir(outdir)
    except OSError:
        pass

    # read in prop_complete, generate distribution
    dist_list = []
    with open(dist_file, "r") as f:
        for line in f:
            dist_list.append(float(line.rstrip()) / 100)

    dist_list = np.array(dist_list)

    # read in file list
    file_list = []
    with open(input, "r") as f:
        for line in f:
            file_list.append(line.rstrip())

    num_genomes = len(file_list)

    # determine complete genomes
    number_complete = round(num_genomes * prop_complete)

    # sample from distribution genomes at random
    incomplete_sample = np.random.choice(dist_list, size=(num_genomes - number_complete), replace=True)
    completeness_array = np.full(num_genomes, 1.0)

    # # set random isolates to non-complete
    incomplete_genomes = np.random.choice(num_genomes, size=incomplete_sample.size, replace=False)
    completeness_array[incomplete_genomes] = incomplete_sample

    to_remove_full = []
    output_file_list = []

    # iterate of genomes to augment
    for index, completeness in enumerate(completeness_array):
        fasta_file = file_list[index]
        basename = os.path.splitext(os.path.basename(fasta_file))[0]

        fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
        output_sequences = []

        sequence_list = []
        name_list = []
        length_list = []

        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)

            # skip if complete
            if completeness == 1.0:
                output_sequences.append(SeqRecord(Seq(sequence), id=name, description=""))
                continue

            name_list.append(name)
            sequence_list.append(sequence)
            length_list.append(len(sequence))
        
        # remove sections from genome
        to_remove = []
        if completeness < 1.0:
            genome_size = sum(length_list)
            num_breaks = np.random.poisson(lam=avg_breaks, size=1)[0]

            # testing
            #num_breaks = 1

            size_removal = round(genome_size * (1 - completeness))

            # get sizes of each break. If no breaks then just take one big block
            break_sites = []
            if num_breaks > 0:
                break_sites = np.random.choice(size_removal, size=num_breaks, replace=False).tolist()
            break_sites.extend([0, size_removal])

            # determine size of each break
            break_sizes = np.diff(np.sort(break_sites))
            np.random.shuffle(break_sizes)
            sum_breaks = sum(break_sizes)

            # search between start and furthest possible position to remove all blocks
            start_pos = 0
            end_pos = genome_size - sum_breaks
            for block in break_sizes:
                break_start = np.random.choice(range(start_pos, end_pos), size=1)[0]
                break_end = break_start + block
                
                to_remove.append((break_start, break_end))
                
                # update site values
                start_pos = break_end
                sum_breaks -= block
                end_pos = genome_size - sum_breaks

            # determine where contig breaks are
            # list_A = ['ATGC', 'AATC', 'CCTG']
            # list_B = [(2, 6), (7, 10)]
            # length_list = [4, 4, 4]
            augmented_sequence = remove_ranges_from_sequence(sequence_list, to_remove, length_list)

            for name, sequence in zip(name_list, augmented_sequence):
                if len(sequence) > 0:
                    #test = len(sequence)
                    output_sequences.append(SeqRecord(Seq(sequence), id=name, description=""))
        else:
            to_remove.append((None, None))

        to_remove_full.append(to_remove)
        SeqIO.write(output_sequences, outdir + "/" + basename + "_comp_" + str(round(completeness, 4)) + ".fasta", 'fasta')
        output_file_list.append(basename + "_comp_" + str(round(completeness, 4)) + ".fasta")
    
    # output completeness file
    with open(outdir + "/completeness.tsv", "w+") as f:
        f.write("File_name\tCompleteness\tRemoved\n")
        for output_file, completeness, removed in zip(output_file_list, completeness_array, to_remove_full):
            f.write(output_file + "\t" + str(completeness) + "\t" + ', '.join([f'({item[0]}, {item[1]})' for item in removed]) + "\n")
            

if __name__ == "__main__":
    main()

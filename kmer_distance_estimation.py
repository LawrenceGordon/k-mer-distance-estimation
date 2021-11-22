#! /usr/bin/env python3

import argparse, Bio.SeqIO, math, itertools, re
import matplotlib.pyplot as plt
from statistics import pstdev

# parses command-line arguments with fasta files, kmer size, and an output file
def parse_args():
    parser = argparse.ArgumentParser(description="K-mer Based Distance Estimation")
    parser.add_argument("fasta", nargs="+", help="fasta files containing nucleic acid sequences")
    parser.add_argument("-k", type=int, default=1, help="k-mer size")
    parser.add_argument("-m", action="store_true", 
        help="if specified, returns mahanalobis distance; defaults to euclidean distance")
    parser.add_argument("--out", "-o", help="kmer data table output")
    return parser.parse_args()

'''
# extracts organisms from fasta headers using regex terms
def get_org(fasta):
    organisms = []
    for f in fasta:
        with open(f, "r") as handle:
            for line in handle:
                if line.startswith(">"):
                    regex = "[A-Z](\.?|[a-z]+) [a-z]+"
                    organisms.append(re.search(regex, line))
                    break
    return organisms
'''

# generates kmer dictionary from fasta formatted data
def parse_fasta(fasta, k):
    kmer_dict = {}

    # parses fasta file using SeqRecord iterator
    for record in Bio.SeqIO.parse(fasta, "fasta"):
        seq = record.seq

        # Iterates over all combinations of k continuous nucleotides and adds the kmers to
        # a dictionary with a counter of the number of instances of the kmers
        for nuc_idx in range(len(seq) - (k - 1)):
            kmer = seq[nuc_idx:(nuc_idx + k)]
            kmer_dict.setdefault(kmer, 0)
            kmer_dict[kmer] += 1
            
    return kmer_dict

# compares the keys in two dictionaries, and adds the missing keys
# from each into one another
def compare_kmer(dict_1, dict_2):
    set_1 = set(dict_1.keys())
    set_2 = set(dict_2.keys())
    diff_keys = list(set_1.symmetric_difference(set_2))
    for key in diff_keys:
        dict_1.setdefault(key, 0)
        dict_2.setdefault(key, 0)

# computes euclidean distance based on two kmer dictionaries
def euclid_kmer(dict_1, dict_2):
    euclid_sum = 0
    for coord in dict_2:
        euclid_sum += (dict_2[coord] - dict_1[coord])**2
    return math.sqrt(euclid_sum)

# computes mahanalobis distance based on two kmer dictionaries
def maha_kmer(dict_1, dict_2):
    maha_sum = 0
    for coord in dict_2:
        dev = pstdev((dict_2[coord], dict_1[coord]))
        if dev == 0:
            continue
        #print("C1:", dict_2[coord], "C2:", dict_1[coord], "Dev:", dev)
        maha_sum += ((dict_2[coord]/dev) - (dict_1[coord]/dev))**2
    return math.sqrt(maha_sum)

# makes an output table comparing kmer distance for fasta samples
def make_table(samples, k, maha, output):
    result_cache = {}
    distances = []
    # distances and sub_dist lists used to create heatmap
    ### SAMPLES = [Sample1, Sample2, ... SampleN]
    # generates unique combinations of all samples
    sample_combinations = list(itertools.combinations(samples, r=2))
    table = "{0}-mer comparison\t{1}".format(k, ("\t").join(samples))
    for sample1 in samples:
        sub_dist = []
        table += "\n{0}\t".format(sample1)
        for sample2 in samples:
            if sample1 == sample2:
                table += "0\t"
                sub_dist.append(0)
                continue
            
            # runs kmer distance pipleine to generate table data
            if (sample1, sample2) in sample_combinations:
                kmer1 = parse_fasta(sample1, k)
                kmer2 = parse_fasta(sample2, k)
                compare_kmer(kmer1, kmer2)
                dist = maha_kmer(kmer1, kmer2) if maha else euclid_kmer(kmer1, kmer2)
                sub_dist.append(dist)
                # caches result from two samples in a dictionary
                result_cache[(sample1, sample2)] = dist
                table += "{0}\t".format(dist)
            # prevents duplicate running of sample comparison by accessing
            # cached result in the dictionary
            else:
                table += "{0}\t".format(result_cache[(sample2, sample1)])
                sub_dist.append(result_cache[(sample2, sample1)])
        distances.append(sub_dist)

    if output is None:
        print(table)
    else:
        with open(output, "w") as out_data:
            out_data.write(table + "\n")
    
    print(distances)
    # generates heatmap from the distances calculated above
    fig, ax = plt.subplots()

    ax.set_xticks(range(len(samples))); ax.set_yticks(range(len(samples)))
    ax.set_xticklabels(samples); ax.set_yticklabels(samples)

    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # displays heatmap wit
    heatmap = ax.imshow(distances, cmap="coolwarm")
    fig.colorbar(heatmap)

    plt.show()
    fig.savefig(output, format="jpg")
            

            

#test1 = {"A": 4, "B": 16}
#test2 = {"A": 7, "B": 12, "C": 14}
#compare_kmer(test1, test2)
#print(euclid_kmer(test1, test2))


args = parse_args()
make_table(args.fasta, args.k, args.m, args.out)

#PA = parse_fasta(args.f[0], args.k)
#KP = parse_fasta(args.f[1], args.k)
#compare_kmer(PA, KP)
#print(euclid_kmer(PA, KP))
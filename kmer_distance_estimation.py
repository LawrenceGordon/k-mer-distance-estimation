#! /usr/bin/env python3

import argparse, Bio.SeqIO, math, re
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

# generates kmer dictionary from fasta formatted data
def parse_fasta(fasta, k):
    kmer_dict = {}

    # parses fasta file using SeqRecord iterator
    for record in Bio.SeqIO.parse(fasta, "fasta"):
        seq = record.seq
        regex = "([A-Z])[a-z]+( [a-z]+)"
        match = re.search(regex, record.description)
        genus = match.group(1)
        species = match.group(2)
        organism = genus + "." + species
        #organism = re.search(regex, record.description).group(0)

        # Iterates over all combinations of k continuous nucleotides and adds the kmers to
        # a dictionary with a counter of the number of instances of the kmers
        for nuc_idx in range(len(seq) - (k - 1)):
            kmer = seq[nuc_idx:(nuc_idx + k)]
            kmer_dict.setdefault(kmer, 0)
            kmer_dict[kmer] += 1
            
    print(fasta + ":" + organism)
    return (organism, kmer_dict)    

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
        maha_sum += ((dict_2[coord]/dev) - (dict_1[coord]/dev))**2
    return math.sqrt(maha_sum)

def make_table(samples, k, maha, output):
    distance_dict = {}
    orgs = []

    for sample in samples:
        org, distance_dict[sample] = parse_fasta(sample, k)
        if org not in orgs:
                orgs.append(org)

    store_compare = {}
    distances = []
    for sample1 in samples:
        sub_dist = []
        for sample2 in samples:
            if (sample2, sample1) in store_compare:
                sub_dist.append(store_compare[sample2, sample1])
                continue
            if sample1 == sample2:
                sub_dist.append(0.0)
                continue

            kmer1 = distance_dict[sample1]
            kmer2 = distance_dict[sample2]

            compare_kmer(kmer1, kmer2)
            dist = maha_kmer(kmer1, kmer2) if maha else euclid_kmer(kmer1, kmer2)
            store_compare[(sample1, sample2)] = dist
            sub_dist.append(dist)

        distances.append(sub_dist)

    print(distances)
    fig, ax = plt.subplots()

    ax.set_xticks(range(len(orgs))); ax.set_yticks(range(len(orgs)))
    ax.set_xticklabels(orgs); ax.set_yticklabels(orgs)
    ax.set_title("{0}-mer Distance Estimation".format(k))
    #ax.minorticks_on
    #ax.grid(True, which='minor', color='w', linestyle='-', linewidth=3)

    for org1 in range(len(orgs)):
        for org2 in range(len(orgs)):
            data = ax.text(org1, org2, str(round(distances[org1][org2], 1)),
                ha="center", va="center", color="w", fontsize=16)

    #ax.invert_yaxis()
    ax.xaxis.tick_top()

    # displays heatmap wit
    heatmap = ax.imshow(distances, cmap="coolwarm")
    fig.colorbar(heatmap)

    #plt.setp(ax.get_yticklabels(), rotation=90, ha="center", 
        #va="center", rotation_mode="anchor", fontsize=10)
    plt.setp(ax.get_yticklabels(), fontsize=10)
    plt.setp(ax.get_xticklabels(), fontsize=10)
    #plt.tight_layout()
    #plt.show()
    #fig.savefig(output, format="jpg")

    fig.set_size_inches(10, 8)
    plt.savefig(output + ".jpg", bbox_inches="tight")

args = parse_args()
make_table(args.fasta, args.k, args.m, args.out)
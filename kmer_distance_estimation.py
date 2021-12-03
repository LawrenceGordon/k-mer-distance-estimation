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
    parser.add_argument("-n", action="store_true", help="if specified, normalizes distance to percent")
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
        # calculates stdev from two kmer coordinates
        dev = pstdev((dict_2[coord], dict_1[coord]))
        if dev == 0:
            continue
        maha_sum += ((dict_2[coord]/dev) - (dict_1[coord]/dev))**2
    return math.sqrt(maha_sum)

# normalizes distances in the output matrix to 1
def normalize(distances):
    all_values = []
    for sub in distances:
        for val in sub:
            all_values.append(val)

    minimum = min(all_values)
    maximum = max(all_values)

    # sets values in sub_lists to normalized values
    for sub in distances:
        for idx, val in enumerate(sub):
            norm_val = (val - minimum)/(maximum - minimum)
            sub[idx] = norm_val

    return distances

def matrix_output(samples, k, maha):
    ### {sample_N: {sample_N kmer dict}, ...} ###
    distance_dict = {}
    orgs = []

    # pairs kmer dictionaries with organism names
    for sample in samples:
        org, distance_dict[sample] = parse_fasta(sample, k)
        # adds organism to list for output header
        if org not in orgs:
                orgs.append(org)

    ### {(sample_x, sample_y): dist, ...} ###
    # stores previously computed distance values so program only
    # has to run unique pairs of orgs
    store_compare = {}

    # distance matrix output; list of lists
    distances = []
    for sample1 in samples:
        sub_dist = []
        for sample2 in samples:
            # if distance already computed, grab cached value
            if (sample2, sample1) in store_compare:
                sub_dist.append(store_compare[sample2, sample1])
                continue
            # if samples are the same, dist = 0
            if sample1 == sample2:
                sub_dist.append(0.0)
                continue

            # grabs kmer dictionaries from distance_dict
            kmer1 = distance_dict[sample1]
            kmer2 = distance_dict[sample2]

            compare_kmer(kmer1, kmer2)
            # calculates distance using eithe mahalanobis distance or euclidean
            # distance based on command-line input
            dist = maha_kmer(kmer1, kmer2) if maha else euclid_kmer(kmer1, kmer2)
            # adds dist to cache dictionary
            store_compare[(sample1, sample2)] = dist
            sub_dist.append(dist)

        distances.append(sub_dist)
    return orgs, distances

def make_table(samples, k, maha, norm, output):
    ### {sample_N: {sample_N kmer dict}, ...} ###
    distance_dict = {}
    orgs = []

    # pairs kmer dictionaries with organism names
    for sample in samples:
        org, distance_dict[sample] = parse_fasta(sample, k)
        # adds organism to list for output header
        if org not in orgs:
                orgs.append(org)

    ### {(sample_x, sample_y): dist, ...} ###
    # stores previously computed distance values so program only
    # has to run unique pairs of orgs
    store_compare = {}

    # distance matrix output; list of lists
    distances = []
    for sample1 in samples:
        sub_dist = []
        for sample2 in samples:
            # if distance already computed, grab cached value
            if (sample2, sample1) in store_compare:
                sub_dist.append(store_compare[sample2, sample1])
                continue
            # if samples are the same, dist = 0
            if sample1 == sample2:
                sub_dist.append(0.0)
                continue

            # grabs kmer dictionaries from distance_dict
            kmer1 = distance_dict[sample1]
            kmer2 = distance_dict[sample2]

            compare_kmer(kmer1, kmer2)
            # calculates distance using eithe mahalanobis distance or euclidean
            # distance based on command-line input
            dist = maha_kmer(kmer1, kmer2) if maha else euclid_kmer(kmer1, kmer2)
            # adds dist to cache dictionary
            store_compare[(sample1, sample2)] = dist
            sub_dist.append(dist)

        distances.append(sub_dist)

    print(distances)

    # normalizes the data in distances if user specifies
    if norm:
        distances = normalize(distances)

    print(distances)

    # generates a heatmap based on orgaism similarity
    fig, ax = plt.subplots()

    ax.set_xticks(range(len(orgs))); ax.set_yticks(range(len(orgs)))
    ax.set_xticklabels(orgs); ax.set_yticklabels(orgs)
    ax.set_title("{0}-mer Distance Estimation".format(k))
    #ax.minorticks_on
    #ax.grid(True, which='minor', color='w', linestyle='-', linewidth=3)

    # adds value from distance matrix to correct square
    for org1 in range(len(orgs)):
        for org2 in range(len(orgs)):
            data = ax.text(org1, org2, str(round(distances[org1][org2], 2)),
                ha="center", va="center", color="w", fontsize=16)

    #ax.invert_yaxis()
    ax.xaxis.tick_top()

    # creates heatmap from axes and matrix
    heatmap = ax.imshow(distances, cmap="coolwarm")
    fig.colorbar(heatmap)

    plt.setp(ax.get_xticklabels(), rotation=45, ha="left", 
        va="top", rotation_mode="anchor", fontsize=10)

    #for tick in ax.xaxis.get_majorticklabels():
    #tick.set_horizontalalignment("left")
    plt.setp(ax.get_yticklabels(), fontsize=10)
    #plt.setp(ax.get_xticklabels(), fontsize=10)
    #plt.tight_layout()
    #plt.show()
    #fig.savefig(output, format="jpg")

    fig.set_size_inches(10, 8)
    plt.savefig(output + ".jpg", bbox_inches="tight")

args = parse_args()
make_table(args.fasta, args.k, args.m, args.n, args.out)
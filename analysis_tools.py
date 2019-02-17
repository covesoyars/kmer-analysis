from itertools import product
import numpy as np
import regex as rx

""" This program conatins functions for the reading and analysis of fasta files. parse_fasta opens and reads in a fasta
    file. filter_fasta reads in a fasta file, filters it and rewrites it. get_kmers finds all possible kmers of a given
    alphabet and k. kmer_composition computes the kmer composition of a sequence given a list of kmers. kmer_distance
    computes the euclidiean distance between two kmer compositions. write_csv writes a dictionary to a csv file
    @author: Cove Soyars
"""


def parse_fasta(fasta_file, grab_all=False):
    """
    Parses fasta file and returns dictionary with its contents.
    :param fasta_file:  file name of file to be parsed
    :param grab_all:    indicates what part of header to grab
    :return: a dictionary in the format: {gene_name: "sequence"}
    """
    genes = {}  # create dictionary to hold gene names and their sequences
    start = True  # set boolean to handle start of file

    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            line = line.strip()
            if line.startswith(">"):

                if start:  # if this is the first header,
                    start = False  # change boolean

                else:  # if this isn't the first header,

                    genes[name] = ''.join(sequence)  # add sequence to dictionary with its name

                sequence = []  # reset sequence

                if not grab_all:
                    name = line.split(' ')[0][1:]  # grab gene name (ignoring '>')
                else:
                    name = line  # grab entire header

            else:  # if this line isn't a header,
                sequence.append(line)  # add its sequence to the entire sequence

    return genes


def get_kmers(k, alphabet="ACGT"):
    """

    :param k: length of kmers
    :param alphabet: alphabet to make kmers
    :return: list of possible kmers in lexigraphical order
    """
    return [''.join(kmer) for kmer in product(alphabet, repeat=k)]


def kmer_composition(kmers, seq):
    """
    Computes the kmer composition of a sequence with a given k
    :param k: length of kmers in composition
    :param seq: sequence
    :return: numpy array of kmer composition where array[i] = kmers[i]
    """

    # loop through sequence and compute kmer composition
    kmer_comp = [0] * len(kmers)

    for kmer in range(0, len(kmers)):
        kmer_comp[kmer] = len(rx.findall(kmers[kmer], seq, overlapped=True))

    # turn kmer_comp (list) into a numpy array and return it
    return np.array(kmer_comp)


def kmer_distance(a1, a2):
    """
    Computes the euclidean distance between two kmer composition arrays
    :param a1: first array
    :param a2: second array
    :return: euclidean distance between arrays
    """

    # return the euclidean distance of the two arrays
    distance = 0
    for i in range(0, len(a1)):
        distance += (a1[i] - a2[i]) ** 2

    return distance


def filter_fasta(prefiltered, postfiltered, bad_words):
    """
    Reads fasta and filters out unwanted genes and writes remaining genes to new fasta
    :param prefiltered: fasta file to be filtered
    :param postfiltered: file name for filtered file
    :param bad_words: words that indicate a gene should be ignored
    :return: the name of the filtered file
    """

    pre = parse_fasta(prefiltered, grab_all=True)  # create dictionary: {entire header: "sequence"}

    with open(postfiltered, 'w') as post:  # open file to write to
        for gene in pre:

            write = True
            # if any bad words are in header,
            for word in bad_words:
                if gene.lower().find(word.lower()) != -1:
                    write = False  # the gene is not written to postfiltered
            if write:
                post.write(gene)
                post.write('\n')
                post.write(pre[gene])
                post.write('\n')

    return postfiltered


def write_csv(dictionary, filename):
    """
    Writes a dictionary to a csv file
    :param dictionary: dictionary to be written
    :param filename:    name of file to be created
    :return: csv file
    """
    with open(filename, 'w') as out:
        for k, v in dictionary.items():
            out.write(str(k[0]) + ",")
            out.write(str(k[1]) + ",")
            out.write(str(dictionary[k]) + ",\n")

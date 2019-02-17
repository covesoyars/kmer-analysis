import analysis_tools as at
from heapq import nlargest, nsmallest

""" This program writes to 4 separate csv files the 10 most similar and 10 most different genes by their kmer composition 
    when k = 2 and k = 4. It uses functions from analysis_tools.py and adapts them to be called on entire dictionaries
    of data.
    @author: Cove Soyars
"""

def get_kmer_comp(genes_dict, kmers):
    """
    Returns kmer composition
    :param genes_dict: gene sequene dictionary
    :param kmers:  list of possible kmers
    :return: dictionary: {gene_name: kmer composition}
    """

    kmer_comp = {}

    for gene in genes_dict:
        kmer_comp[gene] = at.kmer_composition(kmers, genes_dict[gene])

    return kmer_comp


def get_distances(comp_dict):
    """
    Returns euclidiean distance of all possible pairs in a kmer comp dictionary
    :param comp_dict: composition dictionary
    :return: dictionary: {gene1, gene2 : distance}
    """

    distance_dict = {}
    for gene1 in comp_dict:
        for gene2 in comp_dict:

            # make sure to ignore redunant/self matching pairs
            if gene1 == gene2:
                pass
            elif (gene1, gene2) in distance_dict:
                pass
            elif (gene2, gene1) in distance_dict:
                pass
            # calculate distance
            else:
                distance_dict[(gene1, gene2)] = at.kmer_distance(comp_dict[gene1], comp_dict[gene2])

    return distance_dict


def get_top(distances, most_similar=True):
    """
    Return ten most similar or different gene pairs
    :param distances: distances of gene pairs
    :param most_similar: if True: 10 most similar, if False: 10 least similar
    :return: dictionary filtered for 10 most interesting pairs
    """

    # get ten most/least similar pairs
    if most_similar:
        top = nsmallest(10, distances, key = distances.get)
    else:
        top = nlargest(10, distances, key = distances.get)

    # wrap most interesting into dictionary
    top_dict = {}
    for key in top:
        top_dict[key] = distances[key]

    return top_dict


def main():
    genes = at.parse_fasta(at.filter_fasta("lambda.fna.txt", "filtered_lambda_genome.fna",
                                           ["tRNA", "rRNA", "ribosomal", "predicted", "cDNA", "transcript variant"]))

    # get possible kmers for k = 2, 4
    two_mers = at.get_kmers(2)
    four_mers = at.get_kmers(4)

    # get kmer compositions of genes for k = 2, 4
    two_mer_comp = get_kmer_comp(genes, two_mers)
    four_mer_comp = get_kmer_comp(genes, four_mers)

    # calcuate kmer composition differences between genes for k = 2, 4
    two_mer_dist = get_distances(two_mer_comp)
    four_mer_dist = get_distances(four_mer_comp)

    # get ten most similar, most different pairs for k = 2, 4 and write to csv:

    at.write_csv(get_top(two_mer_dist,most_similar=True),"lambda_most_similar_2mer.csv")

    at.write_csv(get_top(two_mer_dist, most_similar=False), "lambda_least_similar_2mer.csv")

    at.write_csv(get_top(four_mer_dist, most_similar=True), "lambda_most_similar_4mer.csv")

    at.write_csv(get_top(four_mer_dist, most_similar=False), "lambda_least_similar_4mer.csv")


if __name__ == "__main__":
    main()

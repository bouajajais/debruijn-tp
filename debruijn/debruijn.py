#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

import matplotlib.pyplot as plt

__author__ = "Ismail Bouajaja"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Ismail Bouajaja"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Ismail Bouajaja"
__email__ = "bouajajais@eisti.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()

def read_fastq(fastq_file):
    with open(fastq_file, 'rt') as f:
        for _ in f:
            line = next(f)[:-1]
            yield line
            next(f)
            next(f)

def cut_kmer(read, kmer_size):
    return (read[index:index+kmer_size] for index in range(len(read) - kmer_size + 1))
    # index = 0
    # while index + kmer_size <= len(read):
    #     yield read[index:index+kmer_size]
    #     index += 1

def build_kmer_dict(fastq_file, kmer_size):
    reads = read_fastq(fastq_file=fastq_file)
    kmer_dict = {}
    for read in reads:
        kmers = cut_kmer(read=read, kmer_size=kmer_size)
        for kmer in kmers:
            if kmer not in kmer_dict.keys():
                kmer_dict[kmer] = 1
            else:
                kmer_dict[kmer] += 1
    return kmer_dict

def build_graph(kmer_dict):
    kmer_graph = nx.DiGraph()
    for kmer, weight in kmer_dict.items():
        kmer_size = len(kmer)
        kmer_graph.add_edge(kmer[:kmer_size-1], kmer[1:], weight=weight)
    return kmer_graph

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    pass

def std(data):
    pass

def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    pass

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    pass

def get_sink_nodes(graph):
    pass

def get_contigs(graph, starting_nodes, ending_nodes):
    pass

def save_contigs(contigs_list, output_file):
    pass

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    # args = get_arguments()
    # lines = read_fastq(fastq_file=args.fastq_file)
    # kmers = cut_kmer(read=next(lines), kmer_size=5)
    # kmer_dict = build_kmer_dict(fastq_file=args.fastq_file, kmer_size=args.kmer_size)
    # kmer_graph = build_graph(kmer_dict=kmer_dict)

if __name__ == '__main__':
    main()
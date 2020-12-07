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

from networkx.algorithms.simple_paths import all_simple_paths
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

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)

def read_fastq(fastq_file):
    with open(fastq_file, 'rt') as f:
        for _ in f:
            yield next(f)[:-1]
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
    starting_nodes = []
    for node in graph.nodes:
        if len(list(graph.predecessors(node))) == 0:
            starting_nodes.append(node)
    return starting_nodes

def get_sink_nodes(graph):
    ending_nodes = []
    for node in graph.nodes:
        if len(list(graph.successors(node))) == 0:
            ending_nodes.append(node)
    return ending_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for starting_node in starting_nodes:
        for ending_node in ending_nodes:
            paths = list(all_simple_paths(graph, starting_node, ending_node))
            for path in paths:
                contig = path[0]
                for node in path[1:]:
                    contig += node[-1]
                contigs.append((contig, len(contig)))
    return contigs

def save_contigs(contigs_list, output_file):
    with open(output_file, 'wt') as f:
        for i, contig_tuple in enumerate(contigs_list):
            contig, length = contig_tuple
            f.write(f'>contig_{i} len={length}\n')
            f.write(fill(f'{contig}\n'))

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # lines = read_fastq(fastq_file=args.fastq_file)
    # # for line in lines: print(line)
    # kmers = cut_kmer(read=next(lines), kmer_size=args.kmer_size)
    # # for kmer in kmers: print(kmer)
    kmer_dict = build_kmer_dict(fastq_file=args.fastq_file, kmer_size=args.kmer_size)
    kmer_graph = build_graph(kmer_dict=kmer_dict)
    starting_nodes = get_starting_nodes(kmer_graph)
    ending_nodes = get_sink_nodes(kmer_graph)
    contigs_list = get_contigs(kmer_graph, starting_nodes, ending_nodes)
    filename = args.fastq_file.split('/')[-1].split('.')[0]
    # save_contigs(contigs_list=contigs_list, output_file=f'{filename}.fna')

if __name__ == '__main__':
    main()
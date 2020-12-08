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
from numpy.lib.function_base import select
random.seed(9001)
from random import randint
import statistics

import matplotlib.pyplot as plt
from networkx.algorithms.lowest_common_ancestors import lowest_common_ancestor

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
    '''returns a generator that yields lines containing a sequence'''
    with open(fastq_file, 'rt') as f:
        for _ in f:
            yield next(f)[:-1]
            next(f)
            next(f)

def cut_kmer(read, kmer_size):
    '''returns a generator that yields kmer of a sequence one by one'''
    return (read[index:index+kmer_size] for index in range(len(read) - kmer_size + 1))
    # index = 0
    # while index + kmer_size <= len(read):
    #     yield read[index:index+kmer_size]
    #     index += 1

def build_kmer_dict(fastq_file, kmer_size):
    '''returns a dictionary which keys are the kmer and values are their frequencies in all the sequences'''
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
    '''returns graph which edges are kmer weighted with their frequencies'''
    kmer_graph = nx.DiGraph()
    for kmer, weight in kmer_dict.items():
        kmer_size = len(kmer)
        kmer_graph.add_edge(kmer[:kmer_size-1], kmer[1:], weight=weight)
    return kmer_graph

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    '''
    removes all the nodes present in 'path_list' from 'graph'
    except entring nodes or ending nodes if their booleans are False
    Note : a node is an entring node if it's the first node of the path
    '''
    new_graph = nx.DiGraph(graph)
    for path in path_list:
        i_start = 0 if delete_entry_node else 1
        i_end = len(path) if delete_sink_node else -1
        new_graph.remove_nodes_from(path[i_start:i_end])
    return new_graph

def std(data):
    return statistics.stdev(data)

def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    '''
    find the best path in path_list and remove all the other ones from the graph
    '''
    # search for best path
    # keep heaviest paths
    max_weight = max(weight_avg_list)
    best_paths = [path_list[i] for i in range(len(path_list)) if weight_avg_list[i] == max_weight]
    path_length = [path_length[i] for i in range(len(path_length)) if weight_avg_list[i] == max_weight]
    # keep longest paths
    max_length = max(path_length)
    best_paths = [best_paths[i] for i in range(len(best_paths)) if path_length[i] == max_length]
    # keep random path
    rand_index = randint(0, len(best_paths) - 1)
    best_path = best_paths[rand_index]
    # remove the best path from path_list then use remove_paths
    path_list.remove(best_path)
    new_graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    new_graph.add_edges_from([(n1, n2, graph[n1][n2]) for (n1, n2) in zip(best_path[:-1], best_path[1:])])
    return new_graph

def path_average_weight(graph, path):
    '''computes the mean of all the weights of the edges of a path'''
    mean = 0
    for node1, node2 in zip(path[:-1], path[1:]):
        mean += graph[node1][node2]['weight']
    return mean / (len(path) - 1)

def solve_bubble(graph, ancestor_node, descendant_node):
    '''find the best path between two nodes and remove the remaining'''
    path_list = list(all_simple_paths(graph, ancestor_node, descendant_node))
    path_length = [len(path) for path in path_list]
    weight_avg_list = [path_average_weight(graph, path) for path in path_list]
    return select_best_path(graph, path_list, path_length, weight_avg_list)

def simplify_bubbles(graph):
    '''
    iterate over all pairs of nodes until a bubble is found
    remove the bubble
    repeat over the new graph until all pairs of nodes have already been tested
    '''
    new_graph = nx.DiGraph(graph)
    done_pairs = []
    while True:
        nodes = list(new_graph.nodes)
        n_nodes = len(nodes)
        if n_nodes <= 2: return new_graph
        i1, i2 = 0, 1
        while True:
            node1, node2 = nodes[i1], nodes[i2]
            if (node1, node2) not in done_pairs and len(list(all_simple_paths(new_graph, node1, node2))) > 1:
                new_graph = solve_bubble(new_graph, node1, node2)
                done_pairs.append((node1, node2))
                break
            if i2 == n_nodes - 1:
                if i1 == n_nodes - 2:
                    break
                i1 += 1
                i2 = i1 + 1
            else:
                i2 += 1
        if i1 == n_nodes - 2 and i2 == n_nodes - 1:
            break
    return new_graph

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    '''returns nodes that don't have predecessors'''
    starting_nodes = []
    for node in graph.nodes:
        if len(list(graph.predecessors(node))) == 0:
            starting_nodes.append(node)
    return starting_nodes

def get_sink_nodes(graph):
    '''returns nodes that don't have successors'''
    ending_nodes = []
    for node in graph.nodes:
        if len(list(graph.successors(node))) == 0:
            ending_nodes.append(node)
    return ending_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    '''returns a list of tuples containing each contig and its length'''
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
    '''saves the contigs into a file'''
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
    print(starting_nodes)
    print(ending_nodes)
    contigs_list = get_contigs(kmer_graph, starting_nodes, ending_nodes)
    filename = args.fastq_file.split('/')[-1].split('.')[0]
    # save_contigs(contigs_list=contigs_list, output_file=f'{filename}.fna')

if __name__ == '__main__':
    main()
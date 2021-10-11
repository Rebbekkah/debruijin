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

__author__ = "Goulancourt Rebecca & Jamay Théo"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Goulancourt Rebecca & Jamay Théo"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Goulancourt Rebecca & Jamay Théo"
__email__ = "rgoulancourt@yahoo.com, jamay.theo@gmail.com"
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
						default=22, help="K-mer size (default 21)")
	parser.add_argument('-o', dest='output_file', type=str,
						default=os.curdir + os.sep + "contigs.fasta",
						help="Output contigs in fasta file")
	parser.add_argument('-f', dest='graphimg_file', type=str,
						help="Save graph as image (png)")
	return parser.parse_args()


def read_fastq(fastq_file):
	with open(fastq_file, "r") as fasta_file : 
		lines = fasta_file.readlines()
		for i in range(1, len(lines), 4) :
			#print(lines)
			yield lines[i].strip()


def cut_kmer(read, kmer_size):
	'''
	for i in range(len(read) - len(kmer_size) - 1) :
		yield read[i:i + kmer_size]
	'''
	
	'''
	start = 0
	for i in read :
		print(read)
	for i in range(len(list(read))) :
	#for i in enumerate(read) :
		k_mer = sequence[start:kmer_size]
		start = start + kmer_size
		k_mer = ''.join(map(str, k_mer)) 
		print(k_mer)
		yield k_mer
	'''



	l = len(list(read))
	#l = len(read)
	print("longueur = {}".format(l))
	#print(read)
	print(type(read))
	#read_list = list(read)
	#read_list.append(read)
	#print(read_list)

	#print(type(read))
	#for r in read :
	#	print(r)

	#row_read = next(read)
	#r = list(read.next())
	#print(r)
	for i in range(0, l - kmer_size, kmer_size) : 
		#kmer = read.next()
		#kmer = next(read)
		#print(kmer)
		#r = list(read.next())
		kmer = read[i:i + kmer_size]
		#kmer = read_list[i:i + kmer_size]
		#kmer = row_read[i:i + kmer_size]
		yield kmer
  
  
def build_kmer_dict(fastq_file, kmer_size):
	kmer_dict = {}
	read = read_fastq(fastq_file)
	cut = cut_kmer(fastq_file, args.kmer_size)

	for i in range(len(read)) :
		if k not in kmer_dict :
			kmer_dict[k] = 1
		else :
			kmer_dict[k] += 1

	print(type(kmer_dict))
	yield kmer_dict





def build_graph(kmer_dict):
	cle = list(kmer_dict.keys())
	l = len(kmer_dict.keys())
	weight = list(kmer_dict.values())
	print(l)
	print(type(weight))
	print(weight)
	empty_digraph = nx.DiGraph()

	digraph = empty_digraph.add_nodes_from(cle)
	for i in range(l - 1) :
		val1 = cle[i]
		val2 = cle[i + 1]
		digraph.add_edge(val1, val2, weight[i])  
		#digraph.add_edge(val1, val2, [, weight[i]])        

	yield digraph

'''
def get_starting_nodes(graph):
	nodes = list(graph.nodes)
	print(nodes)
	edges = list(graph.edges)
	print(edges)

	start_node = graph.in_edge()

'''
'''
	for node in nodes :
		if len(graph.adj[node]) > 1 :
			start_node.append(node)

	yield start_node
'''
'''
def get_sink_nodes(graph):
	end_node = graph.out_edge()
	yield end_node
'''
'''
def get_contigs(graph, in_node, out_node):
	contig = []

	for i in in_node :
		contig[i] = in_node[i]

	for i + 1 in graph.nodes - 1 :
		contig.append(graph.adj[i])

	for i in out_node :
		contig.append(out_node[i])

	contig = tuple(contig)
	yield contig
'''

	'''
	for in_n in graph.nodes :
		contig[i] = in_node[i]



	for out_n in out_node :
	'''





#https://stackoverflow.com/questions/56918460/how-to-fix-error-object-of-type-generator-has-no-len-python/56918487

#def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
	

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
	nodes = list(graph.nodes)
	print(nodes)
	edges = list(graph.edges)
	print(edges)

	start_node = graph.in_edge()
		yield start_node


def get_sink_nodes(graph):
	end_node = graph.out_edge()
	yield end_node

def get_contigs(graph, starting_nodes, ending_nodes):
	contig = []

	for i in starting_nodes :
		contig[i] = starting_nodes[i]

	for i + 1 in graph.nodes - 1 :
		contig.append(graph.adj[i])

	for i in ending_nodes :
		contig.append(ending_nodes[i])

	contig = tuple(contig)
	yield contig



def save_contigs(contigs_list, output_file):
	with open(output_file, "w") as filout :
		for contig in contigs_list :
			filout.write(contigs_list[i].strip())

	fast = fill(output_file)
	yield fast


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


def save_graph(graph, graph_file):
	"""Save the graph with pickle
	"""
	with open(graph_file, "wt") as save:
			pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
	"""
	Main program function
	"""
	# Get arguments
	args = get_arguments()
	seq = read_fastq(args.fastq_file)
	kmer = cut_kmer(seq, args.kmer_size)
	dico = build_kmer_dict(seq, args.kmer_size)
	digraph = build_graph(dico)
	in_nodes = get_starting_nodes(digraph)
	out_nodes = get_sink_nodes(digraph)
	contig = get_contigs(digraph, in_nodes, out_nodes)
	save_contig = save_contigs(contig, args.output_file)
	fasta = fill(save_contig)

	#for i in seq : 
	#	print(i)
	#for i in kmer :
	#	print(i)
	#print(seq)
	#print(kmer)
	# Fonctions de dessin du graphe
	# A decommenter si vous souhaitez visualiser un petit 
	# graphe
	# Plot the graph
	# if args.graphimg_file:
	#     draw_graph(graph, args.graphimg_file)
	# Save the graph in file
	# if args.graph_file:
	#     save_graph(graph, args.graph_file)


if __name__ == '__main__':
	main()
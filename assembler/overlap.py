#!/usr/bin/python

import sys
import string
import random

import networkx as nx

import matplotlib

# suggestion from http://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined
matplotlib.use('Agg')	# call this before calling pyplot to avoid Xwindow errors

import matplotlib.pyplot as plt

'''
import reportlab

from reportlab.pdfgen import canvas, pathobject
from reportlab.pdfgen.canvas import *
from reportlab.pdfgen.pathobject import *

from reportlab.graphics import renderPDF
from reportlab.graphics.shapes import *

from reportlab.lib.units import cm
from reportlab.lib.colors import *
'''

# cur_version = int(sys.version_info[0])
# print (cur_version)

def readDataFromFile(fileName):
	infile = open(fileName, 'r')
	tmpDict = {}

	for line in infile:
		toks = line.strip().split(" ")

		seqNum = int(toks[0].strip())
		seq = toks[1].strip()

		tmpDict[seqNum] = seq

	print (tmpDict)

	return tmpDict


def meanLength(fileName):
	# readDataFromFile('genome_assembly.txt')
	aDict = readDataFromFile(fileName)
	
	totalLen = 0

	for aKey in aDict.keys():
		aSeq = aDict[aKey]
		totalLen += len(aSeq)

	print(round(float(totalLen / len(aDict.keys()))))
	return round(float(totalLen / len(aDict.keys())))


def getOverlap(left, right, start, kmer_ind):
	overlaps = []

	left = left.strip().upper()
	right = right.strip().upper()

	chunk = right[start:kmer_ind]
# 	print (chunk)

	if len(right) > len(left):
		return overlaps

	if left.find(chunk) >= 0:
		overlaps.append(chunk)

	return overlaps


def overlap():
	kmer_ind = 1
	start = 0

	orig_left = input("Please provide the first sequence: ")
	orig_right = input("Please provide the second sequence: ")

	left = orig_left
	right = orig_right

	# left = "CGATTCCAGGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTC"
	# right = "GGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTCGTCCAGACCCCTAGC"

	end = len(left)

	right = right[start:kmer_ind]

	ovs = []

	while start <= len(left):
		tmp = getOverlap(left, right, start, kmer_ind)
	# 	ovs.append(tmp)

		start += 1
		kmer_ind += 1

		right = right[start:kmer_ind]

	kmer_ind = 0

# 	left = "CGATTCCAGGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTC"
# 	right = "GGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTCGTCCAGACCCCTAGC"

	# reset
	left = orig_left
	right = orig_right

	end = len(left)

	left = left[kmer_ind:end]

	while end >= 0:
		tmp = getOverlap(left, right, kmer_ind, end)
		ovs.append(tmp)

		end = end-1

		right = right[kmer_ind:end]

	# print(ovs)

	longest = max(ovs)[0]
	print ("Max overlap: " + longest)


def assembly(fname):
	
	chunks = []

	kmer_size = 6

	infile = open(fname, 'r')
	seq = ""

	'''
	for line in infile:
		toks = line.strip().split(" ")

		tmp_toks = toks[1:]
		
		seq += "".join(tmp_toks)

	seq = seq.upper()
	'''

	# debug - 1 sequence
	#seq = "CGATTCCAGGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTC"
	seq = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	#print (seq)

	i = 0

	deBruinDict = {}

	# now split the genome into chunks of 100
	while i <= len(seq):
		chunk = seq[i:i+kmer_size]
		deBruinDict[i] = chunk

		if len(chunk) == kmer_size:
			i += 1
		else:
			break			

	revDict = {v: k for k, v in deBruinDict.items()}
	#print (revDict)
	
	chunks = list(deBruinDict.values())

	# starting point - an unordered list of contigs
	random.shuffle(chunks)
	#print (chunks)

	tmpDict = {}
	tmp_chunks = chunks
	
	i = 0
	
	#print (tmp_chunks)
	last = ""

	for ind in deBruinDict:
		chunk = deBruinDict[ind]
		#print ("Outer loop chunk: " + chunk)
		#left = chunk[0:len(chunk)-1]
		#print (left)
		right = chunk[1:]
		#print (right)

		for tc in revDict:
			#print (tc)

			if len(tc) == kmer_size:
				if tc[1:] == right:
					#print ("Chunk: " + chunk)
					#print ("Right: " + tc[1:])
					tmpDict[i] = chunk
					#tmpDict[i+1] = tc
					#print ("dict so far ")
					#print (tmpDict)
					chunk = tc
					#right = tc[1:]

					i += 1
					break
		
			else:
				last = tc

		#print ("HERE NOW!!")
		#print (chunk)

	#print (i)
	#print (last)

	if len(last) == kmer_size:
		tmpDict[i] = last

	gr=nx.DiGraph()
	#gr.add_nodes_from(tmpDict.items())

	#for t in gr.nodes():	# tuple
	for edge in tmpDict.keys():
		node = tmpDict[edge]
		gr.add_node(node)

		'''
		print (t)
		edge = int(t[0])
		node = t[1]
		'''

		print (edge)
		print (node)

		if edge > 0:
			#print (tmpDict[edge-1])
			gr.add_edge(tmpDict[edge-1], tmpDict[edge])
		else:
			gr.add_edge(tmpDict[edge], tmpDict[edge+1])

	'''
	print(gr.nodes())
	print(gr.number_of_edges())
	print(gr.edges())
	'''

	nx.draw_circular(gr, with_labels = True)
	plt.savefig("test.png")

	#nx.write_gpickle(gr,"test.gpickle")

# overlap()
assembly("genome_assembly.txt")

#!/usr/bin/python

import sys
import string
import random

# cur_version = int(sys.version_info[0])
# print (cur_version)

'''
https://kaspermunch.wordpress.com/2013/11/29/exercise-genome-assembly/

In genome assembly many short sequences (reads) from a sequencing machine is assembled into long sequences – ultimately chromosomes. This is done by ordering overlapping reads so they together represent genomic sequence. For example given these three reads: AGGTCGTAG, CGTAGAGCTGGGAG, GGGAGGTTGAAA ordering them based on their overlap like this:

AGGTCGTAG
    CGTAGAGCTGGGAG
             GGGAGGTTGAAA
produces the genomic sequence:

AGGTCGTAGAGCTGGGAGGTTGAAA
So very briefly the task is to read in the sequence reads, identify overlaps between them, find the right order of reads and then reconstruct the genomic sequence from the overlapping reads.

Real genome assembly is of course more sophisticated than what we do here, but the idea is the same. We make two important simplifying assumptions.

There are no sequencing errors implying that overlaps between reads are perfect matches.
No reads are nested in other reads. I.e. always overlaps of
this type:

XXXXXXXXXX
    XXXXXXXXXXX
never this type:

XXXXXXXXXXXXXX
   XXXXXXXX

1.
The first task is to read and parse the input data. Write a function

def readDataFromFile(fileName):

that takes a string fileName as argument. The function must return a dictionary where keys are read names and values are the associated read sequences. Both keys and values must be strings.

Example usage: readDataFromFile('genome_assembly.txt') should return a dictionary with the following content (maybe not with key-value pairs in that order)

{'1': 'GGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTCGTCCAGACCCCTAGC',
 '3': 'GTCTTCAGTAGAAAATTGTTTTTTTCTTCCAAGAGGTCGGAGTCGTGAACACATCAGT',
 '2': 'CTTTACCCGGAAGAGCGGGACGCTGCCCTGCGCGATTCCAGGCTCCCCACGGG',
 '5': 'CGATTCCAGGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTC',
 '4': 'TGCGAGGGAAGTGAAGTATTTGACCCTTTACCCGGAAGAGCG',
 '6': 'TGACAGTAGATCTCGTCCAGACCCCTAGCTGGTACGTCTTCAGTAGAAAATTGTTTTTTTCTTCCAAGAGGTCGGAGT'}

2.
Often there are too many reads to look at them manually. We want to know what the mean length of the reads is.

Write a function

def meanLength(fileName):

that takes a string fileName containing the name of an input file like the one above. The function must return a float beting the mean sequence length.

3.
Next thing is to figure out which reads overlap each other. To do that we need a function that takes two read sequences and computes their overlap. In the input data none of the reads are completely nested in another read. So given two reads that overlap, the 3′ end of the left read overlaps the 5′ end of the right read.

We know that there are no sequencing errors, so in the overlap the sequence match will be perfect. You need to loop over all possible overlaps honoring that one sequence is the left one and the other is the right one. In the for loop, start with the largest possible overlap (len(left)) and evaluate smaller and smaller overlaps until you find an exact match.

Hint: In each iteration of the loop you want to to compare the last len(left)-i bases of the left sequence with the first len(left)-i bases of the right sequence. You can slice out the latter using right[:len(left)-i].

Write a function

def getOverlap(left, right):

that takes two string arguments, left and right, each containing a read sequence. The function must return the overlapping sequence. If there is no overlap it should return an empty string. So with these two reads:

s1 = "CGATTCCAGGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTC"
s2 = "GGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTCGTCCAGACCCCTAGC"
getOverlap(s1, s2) should return 'GGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTC' and
getOverlap(s2, s1) should return 'C'

So in that case it seems that s1 and s2 overlap and that s1 is the left one and s2 is the right one. Treating s2 as the left one and s1 as the right one only gives an overlap of one base (we expect a few bases of overlap even for unrelates sequences).

4.
Now that we can evaluate the overlap between two reads in some orientation we can compute overlaps between all pairs of reads in both left-right and right-left orientations.

Write a function:

def getAllOverlaps(reads):

that takes a dictionary argument reads as produced by readDataFromFile. The function must return a dictionary of dictionaries containing the number of overlapping bases for a pair of reads in a specific orientation. Overlap of a read to itself is meaningless and must not be included. If the overlap between read 2 in left position and 5 in right position is 21 then d['2']['5'] must be 21 (if the dictionary was called d). Make sure you understand how this data structure represents the overlaps before you go on.

Hint: d['2'] will be a dictionary where keys (d['2'].keys()) are the names of reads that have an overlap with read ‘2’ when ‘2’ is put in the left position, and values (d[‘2’].values()) are the number of overlapping bases for those reads.

Example usage: getAllOverlaps(reads) should return a dictionary containing (maybe not with key-value pairs in that order):

{'1': {'3': 0, '2': 1, '5': 1, '4': 0, '6': 29},
 '3': {'1': 0, '2': 0, '5': 0, '4': 1, '6': 1},
 '2': {'1': 13, '3': 1, '5': 21, '4': 0, '6': 0},
 '5': {'1': 39, '3': 0, '2': 1, '4': 0, '6': 14},
 '4': {'1': 1, '3': 1, '2': 17, '5': 2, '6': 0},
 '6': {'1': 0, '3': 43, '2': 0, '5': 0, '4': 1}}
Hint: To generate all combinations of reads you need two for-loops. One looping over reads in left positions and another (inside the first one) looping over reads in right position.

5.
The dictionary returned by getAllOverlaps is a little messy to look at. We want to print it in a nice matrix-like format so we can better see which pairs overlap in what orientations.

Write a function

def prettyPrint(overlaps):

that takes a dictionary argument overlaps containing a dictionary of dictionaries produced by getAllOverlaps. The function should not return anything but must print a matrix exactly as shown in the example below with nicely aligned and right-justified columns. First column must hold names of reads in left conformation. The top row holds names of reads in right conformation. Remaining cells each hold the number of overlapping bases for a left-right read pair. The diagonal corresponds to overlaps to the read itself. You must put dashes in these cells.

Example usage: prettyPrint(overlaps) should print exactly

   1  2  3  4  5  6 
1  -  1  0  0  1 29 
2 13  -  1  0 21  0
3  0  0  -  1  0  1
4  1 17  1  -  2  0 
5 39  1  0  0  - 14 
6  0  0 43  1  0  -

Hint: To print something padded to fill three characters use string formatting like this:

print "% 3d" % 13,
The comma and the end suppresses the trailing newline that is added by default.

Notice that the overlaps are not either zero or a large number. A lot of the overlaps are 1 or 2. This is because you often find a few bases of random overlap between any two reads. So in the following we must distinguish true (significant) overlaps from random ones. We decide that true overlaps are the ones with an overlap larger than two (>2).

6.
Now that we know how the reads overlap we can chain them together pair by pair from left to right to get the order in which they represent the genomic sequence. To do this we take the first (left-most) read and identify which read has the largest overlap to its right end. Then we take that read and find the read with the largest overlap to the right end of that – and so on until we reach the rightmost (last) read.

The first thing you need to do is to identify the first (leftmost) read so we know where to start. This read is identified as the one that only has a significant (>2) overlap to its right end (it only has a good overlap when positioned to the left of other reads). In the example output from prettyPrint above the first read would be read ‘4’ because the ‘4’ column has no significant overlaps (no one larger than two).

Write a function

def findFirstRead(overlaps):
that takes a dictionary argument overlaps returned from getAllOverlaps. It must return a string containing the name of the first read.

Example usage: findFirstRead(reads) should return 4

7.
Now that we have the first read we can recursively find the correct ordering of reads. We want a list with the read names in the right order. To see how this can be solved using recursion, consider that the task can be broken into small identical tasks of finding the read that has the largest overlap to the right end of the current read (that is the recursive case). You just keep doing this until you reach the last (right-most) read that does not have any significant (>2) overlap to its right end (that is the base case).

In the recursive case you need to identify which read that is next. I.e. which read has the largest overlap to the right end of the current read. Here we use our dictionary of overlaps. If the first read is ‘4’ then overlaps[‘4’] is a dictionary of reads with overlap to the right end of read ‘4’. So to find the name of the read with the largest overlap you must write a function that finds the key associated with the largest value in a dictionary. We do that first:

Write a function

def findKeyForLargestValue(d):
that takes a dictionary argument d and returns the key associated with the largest value. Use that function as a tool in the next function.

Write a function

def findOrder(name, overlaps):

that takes a string argument name containing the name of the first read identified in problem 6, and a dictionary argument overlaps returned from getAllOverlaps. The function must return a list of read names in the order in which they represent the genomic sequence.

Hint:
Base case should return [name].
Recursive case should return [name] + findOrder(nextName) where nextName is the name of the read that has the largest overlap to the right end of the current read (name).

Example usage:
findOrder(‘4’, reads): should return: ['4', '2', '5', '1', '6', '3']
Make sure you understand why this is the right list of read names before you try to implement the function.

8.
Now that you have the number of overlapping bases between reads and the correct order of the reads you can reconstruct the genomic sequence.

Write a function

def assembleGenome(readOrder, reads, overlaps):
that takes a list argument readOrder containing the order of reads, a dictionary argument reads returned from readDataFromFile and a dictionary argument overlaps returned from getAllOverlaps. The function must return a string with the genomic sequence.

Hint: iterate over the reads in order and use the overlap information to extract and join the appropriate parts of the reads.

'''

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

# 	left = "CGATTCCAGGCTCCCoriCACGGGGTACCCATAACTTGACAGTAGATCTC"
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

	infile = open(fname, 'r')
	seq = ""

	for line in infile:
		toks = line.strip().split(" ")

		tmp_toks = toks[1:]
		
		seq += "".join(tmp_toks)

	seq = seq.upper()
	
# 	print (seq)

	i = 0

	# now split the genome into chunks of 100
	while i <= len(seq):
		chunks.append(seq[i:i+100])
		i += 100
	
	# and scramble the chunks
	random.shuffle(chunks)	

	kmer_size = 7
	
	tmp_chunks = chunks
	
	deBruinDict = {}
	nextKey = 0

	for chunk in chunks:
# 		print (chunk)
		print ("LAST")
		chunk_last = chunk[len(chunk)-kmer_size:len(chunk)]
		print (chunk_last)
				
		for tc in tmp_chunks:
			if tc != chunk:
				chunk_first = tc[0:kmer_size]

				print ("1")
				print (chunk_first)
					
				if chunk_first == chunk_last:
					print (tc)
					deBruinDict[nextKey] = tc
					deBruinDict[nextKey+1] = chunk
					
					nextKey += 2

# 	print (deBruinDict)

# 	print (chunks)
# 	print(len(chunks))

	

# overlap()
assembly("genome_assembly.txt")
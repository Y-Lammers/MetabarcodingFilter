#!/usr/bin/env python

# Contact: youri.lammers@gmail.com

# This tools process a set of metabarcoding libraries sequenced at the same time
# and removes all sequences that could be caused due to tag switching events. The
# tool requires that the data has been demulitplexed (ngsfilter) and collapsed
# (obiuniq) by the OBITools software package.

# Usage: TagSwitchCorrect.py -i [library link file]

# The link file is in the following format:
# FASTA sequence file {tab} NGSfilter file {tab} Number of raw reads in the library
# were each library and associated files are placed on a single line

# load a bunch of modules
import sys, json, argparse, os, itertools
from collections import defaultdict

# set argument parser
parser = argparse.ArgumentParser(description = """Process a set of metabarcoding
	libraries that were sequenced at the same time and remove all
	sequences that could be caused due to tag switching events""")

parser.add_argument('-i', '--input', metavar='input overview file', dest='inputLib',
			type=str, help="""TSV file containing the libraries, tag
					information and read counts""", required=True)
parser.add_argument('-s', '--swaprate', metavar='index swap rate', dest='swaprate',
			type=float, help="""The index swap rate that the tool will
					use, default = 0.01""", default=0.01)
args = parser.parse_args()


def parse_input(inputfiles):

	# parse through the input file and obtain the 
	# file paths to the libraries, NGS filter files
	# and the raw read counts. Return this info in a
	# nested list

	# create return list
	libraryInfo = [[],[],[],[]]

	# parse through the input file
	for line in open(inputfiles):

		# split the file based on tabs
		line = line.strip().split('\t')

		# add the info to the list
		libraryInfo[0].append(line[0])
		libraryInfo[1].append(line[1])
		libraryInfo[2].append(int(line[2]))

	# calculate the proportion ratio for each
	# library based on the raw read counts, this
	# is needed to work out the swap rates later.
	for i in range(0,len(libraryInfo[0])):
		libraryInfo[3].append(libraryInfo[2][i]/float(
		sum([count for count in libraryInfo[2]])))

	# return the results
	return libraryInfo


def read_library(library_paths):

	# parse through a library and extract the sequence and
	# sample information and store in a python dictionary
	# named sequences. Every library is store in a the
	# same dictionary following the following format:
	# {sequence:{sample_name:count}}

	# create the sequences variable
	sequences = defaultdict(dict)

	# loop through the library paths and open each path
	# to extract the sequence and sample information
	for path in library_paths:

		# open the sequence library
		seq_file = open(path)

		# parse through the fasta file and obtain the sequences
		seq_groups = (x[1] for x in itertools.groupby(seq_file,
			key=lambda line: line[0] == '>'))
		for header in seq_groups:

			# get the fasta header and parse out the
			# sample name and count information
			header = header.next().strip()
			descrip = header.split("merged_sample=")[1]
			descrip = descrip.split(";")[0]
			samples = json.loads(descrip.replace("\'","\""))

			# get the fasta sequence
			sequence = ''.join(seq_line.strip() for seq_line in seq_groups.next())

			# try to add the sample information for the
			# sequence to the sequence dictionary
			# if it fails, add a new entry
			for sample in samples:
				sequences[sequence]['_'.join([path,sample])] = \
					samples[sample]

		# close the sequence file
		seq_file.close()

	# return the sequence dictionary
	return sequences


def read_tags(inputdata):

	# parse through the NGSfilter tag files and extract
	# the PCR tag - sample name information per library.
	# The information is stored in a dictionary with the
	# following format:
	# {tag:{library:sample_name}}
	
	# create the tag and sample variable
	tags, samples = defaultdict(dict), defaultdict(list)

	# loop through the tag and sequence paths and extract the
	# tag - sample name combinations for each library.
	# The tag files and library files are connected based
	# on their position in the file lists (tag file 1 belongs
	# to library file 1, etc)
	for position in range(len(inputdata[1])):

		# parse through the tag file
		for line in open(inputdata[1][position]):

			# skip the comment line
			if line[0] == '#': continue

			# split the line based on tabs
			line = line.strip().split('\t')

			# add the sample names (2nd column) and tag
			# sequence (3rd column) to the tag dictionary
			# in combination with the library path.
			# create the tag item if not already present.
			tags[line[1]][inputdata[0][position]] = \
				['_'.join([inputdata[0][position],line[0]]),
				inputdata[3][position]]

			# add the sample names to a dictionary
			# containing the library names, if the
			# library is not present, add it.
			samples[inputdata[0][position]].append(
				'_'.join([inputdata[0][position],line[0]]))

	# return a list containing the tag information and sample names
	# for each library
	return [tags, samples]


def analyse_tagswitch(tags, sequences):

	# analyse each tag across all libraries and
	# compare read numbers

	# universal swap rate, currently set to 1%
	swaprate = args.swaprate

	# filtered sequence dictionary which will
	# store the final results.
	#sequence_dict = {}
	sequence_dict = defaultdict(dict)

	# go through the tag list
	for tag in tags:
	#for tag in ['GAGCTTAC:GAGCTTAC','GCTTGTGAC:GCTTGTGAC','GTCTGTTCG:GTCTGTTCG','ACAACCGA:ACAACCGA']:

		# loop through the sequences and print
		# the sequence + occurences for the samples
		for sequence in sequences:
			
			# set the temporary storage variable for
			# the sample and sample numbers
			sample_dict = {}

			# for each sample name that was used for
			# a given tag, see if it is present for the
			# sequence, if so add it to the temporary dict
			for sample in tags[tag]:
				try:
					sample_dict[tags[tag][sample][0]] = \
						[sequences[sequence][tags[tag][sample][0]]
						,tags[tag][sample][1]]
				except:
					pass

			# ignore sequences that have no reads for a
			# certain tag.
			if len(sample_dict) == 0: continue

			# create a new variable for the filtered
			# sample numbers and calculate the total
			# number of reads across al samples for a
			# given sequence and tag
			tot_sum = sum([sample_dict[sample][0] 
				for sample in sample_dict])

			# for each sample, check if the number of reads
			# is higher than the proportion of the library
			# times the swap rate. If the number of reads is
			# lower remove it, if it is higher add it to the
			# new filtered dictionary.
			for sample in sample_dict:
				prop = float(sample_dict[sample][0])/tot_sum
				if sample_dict[sample][0] >= (swaprate*tot_sum*
					sample_dict[sample][1]):
					sequence_dict[sequence][sample] = \
						sample_dict[sample]

	# return the filtered dictionary
	return sequence_dict


def write_results(library, library_dict):

	# write the filtered results to a new library file.
	# the file will replace the fasta extension with tagswap.fasta
	
	# open the outputfile
	output = open(os.path.splitext(library)[0] + '.tagswap.fasta', 'w')	

	# set the output name that will be included in the fasta
	# header along with the sequence number
	name, count = os.path.basename(library).split('.')[0], 1

	# go through the sequences in the library and write them
	# to the output file
	for seq in library_dict:

		# check if there are sequences present, if not skip
		if len(library_dict[seq]['sample']) >= 1:

			# prepare the fasta header
			header = '>{0}:{1} merged_sample={2}; count={3};\n'.format(
					name,count,library_dict[seq]['sample'],
					library_dict[seq]['count'])

			# write the header and sequence to the
			# output file
			output.write(header)
			output.write(seq + '\n')

			# increment the sequence count
			count += 1

	# close the outputfile
	output.close()


def output_results(sequence_dict, samples):

	# print the filtered sequence for each library

	# go through each library and get the samples and sequence counts
	for inputfile in samples:

		# set the library dictionary that will hold all the info
		# for each library	
		library_dict = {}		

		# parse through the list of sequences in the tagswap
		# filtered dictionary
		for sequence in sequence_dict:

			# for each sequence, three things will be kept:
			# the samples that have the sequence, the number of
			# read for the sequence per sample and the sum
			# reads for the library
			library_dict[sequence] = {'count': 0, 'sample': {}}

			# go through the samples for a sequence
			for hit in sequence_dict[sequence]:

				# if the sample is from the current library:
				if hit in samples[inputfile]:
					# increase the total count with the
					# number of reads for the sample
					library_dict[sequence]['count'] += \
						sequence_dict[sequence][hit][0]
					# add the sample + count
					library_dict[sequence]['sample'][hit.replace(
						inputfile+"_",'')] = \
						sequence_dict[sequence][hit][0]

		# when the library is complete, write it to the output file
		write_results(inputfile, library_dict)
		

inputdata = parse_input(args.inputLib)
tags, samples = read_tags(inputdata)
sequences = read_library(inputdata[0])
sequence_dict = analyse_tagswitch(tags, sequences)
output_results(sequence_dict,samples)

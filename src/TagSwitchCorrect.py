#!/usr/bin/env python

# Contact: youri.lammers@gmail.com

# This tools process a set of metabarcoding libraries sequenced at the same time
# and removes all sequences that could be caused due to tag switching events.

# Usage: TagSwitchCorrect.py -fasta [list of fasta files] -tag [list of ngsfilter files]
# note: the fasta files and the ngsfilter files need to be in the same order


# load a bunch of modules
import sys, json, argparse, os, itertools

# set argument parser
parser = argparse.ArgumentParser(description = """Process a set of metabarcoding
	libraries that were sequenced at the same time and remove all
	sequences that could be cuased due to tag switching event""")

parser.add_argument('-f', '--fasta', nargs='+', metavar='fasta files', dest='fasta',
			type=str, help='The dereplicated fasta files.', required=True)
parser.add_argument('-t', '--tag', nargs='+', metavar='tag files', dest='tag',
			type=str, help='The NGSfilter tag files.', required=True)
args = parser.parse_args()


def read_library(library_paths):

	# parse through a library and extract the sequence and
	# sample information and store in a python dictionary
	# named sequences. Every library is store in a the
	# same dictionary following the following format:
	# {sequence:{sample_name:count}}

	# create the sequences variable
	sequences = {}

	# loop through the library paths and open each path
	# to extract the sequence and sample information
	for path in library_paths:

		# open the sequence library
		seq_file = open(path)

		# parse through the fasta file and obtain the sequences
		seq_groups = (x[1] for x in itertools.groupby(seq_file, key=lambda line: line[0] == '>'))
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
				try:
					sequences[sequence][sample] = samples[sample]
				except:
					sequences[sequence] = {sample: samples[sample]}

		# close the sequence file
		seq_file.close()

	# return the sequence dictionary
	return sequences


def read_tags(tag_paths, library_paths):

	# parse through the NGSfilter tag files and extract
	# the PCR tag - sample name information per library.
	# The information is stored in a dictionary with the
	# following format:
	# {tag:{library:sample_name}}
	
	# create the tag variable
	tags = {}

	# loop through the tag and sequence paths and extract the
	# tag - sample name combinations for each library.
	# The tag files and library files are connected based
	# on their position in the file lists (tag file 1 belongs
	# to library file 1, etc)
	for position in range(len(tag_paths)):

		# parse through the tag file
		for line in open(tag_paths[position]):

			# skip the comment line
			if line[0] == '#': continue

			# split the line based on tabs
			line = line.strip().split('\t')

			# add the sample names (2nd column) and tag
			# sequence (3rd column) to the tag dictionary
			# in combination with the library path.
			# create the tag item if not already present.
			try:
				tags[line[2]][library_paths[position]] = line[1]
			except:
				tags[line[2]] = {library_paths[position]: line[1]}

	return tags


def analyse_tagswitch(tags, sequences):

	# analyse each tag across all libraries and
	# compare read numbers

	count2 = 0

	# go through the tag list
	for tag in tags:

		# set the sample dictionary
		sample_dict = {}

		count = 0

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
					sample_dict[tags[tag][sample]] = \
						sequences[sequence][tags[tag][sample]]
				except:
					pass

			if len(sample_dict) >= 2:

				# create a new variable for the filtered
				# sample numbers and calculate the total
				# number of reads across al samples for a
				# given sequence and tag
				filter_dict, tot_sum = {}, sum([sample_dict[a] 
					for a in sample_dict])

				# for each sample, check if the number of reads
				# is higher than a 1000th of the total read pool
				# minus the sample readsm. If the number of reads
				# is lower than the sum - sample, it is potentially
				# caused by tagswitching and removed, if it is
				# higher, it is stored in the new filter variable
				for item in sample_dict:
					if "not" in item:
						notused = 1
						prop = float(sample_dict[item])/tot_sum# - sample_dict[item])
						swap = 1/float(0.038/prop)
						print '{0}\t{1}\t{2}\t{3}'.format(sequence, str(prop), str(swap),str(sample_dict[item]))#str(float(sample_dict[item])/(tot_sum - sample_dict[item])))
						if sample_dict[item] >= round((tot_sum - 
							float(sample_dict[item]))/1000):
							filter_dict[item] = sample_dict[item]

tags = read_tags(args.tag, args.fasta)
sequences = read_library(args.fasta)
analyse_tagswitch(tags, sequences)

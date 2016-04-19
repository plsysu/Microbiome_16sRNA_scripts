#!/usr/bin/env python

__author__ = "Peng Li"
__copyright__ = "Copyright 2016."
__credits__ = ["Peng Li"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Peng Li"
__email__ = "Peng-Li@outlook.com"

import sys
import os
import os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

"""
Author : Peng Li
Extract reads from otu map and write to multiple fasta files which corresponding to different OTUs.

Executable format  in cmd:

python extract_seqs_from_otu_map.py <OTU map file> <seqs.fna> <output dir> <OTU>

"""
def check_file_exists(file_path):
	if os.path.isfile(file_path) and os.access(file_path, os.R_OK):
		return(True)
	else:
		return(False)

def safe_make_dir(dirname):
	try:
		os.makedirs(dirname)
	except OSError:
		if os.path.exists(dirname):
			# We are nearly safe
			pass
		else:
			# There was an error on creation, so make sure we know about it
			raise
			
			
def main():
	'''make command line interface'''
	
	parser = argparse.ArgumentParser(prog = 'extract_seqs_from_otu_map.py', description = 'Extract reads from otu map and write to multiple fasta files which corresponding to different OTUs.')
	parser.add_argument('--otu_map', action = "store",required=True, help = 'The OTU and corresponding sequence IDs file output by QIIME, one line one OTU, seperate by tab. ')
	parser.add_argument('--seq', action = "store",required=True, help = 'The fasta sequences used to analysis. Output by QIIME\'s split_library.py script. ')
	parser.add_argument('--out_dir', action = "store",required=True, help = 'The output directory. ')
	parser.add_argument('--otu', action = "store",required=True, help = 'The interesting OTUs to extract. one OTU one line ')
	args = parser.parse_args()
	# parse otu map file
	otu_map = args.otu_map	
	if(not check_file_exists(otu_map)):
		print ("Either OTU_MAP file is missing or is not readable")
		sys.exit()
	
	# parse seq file
	seqs=args.seq
	if(not check_file_exists(otu_map)):
		print ("Either seqs file is missing or is not readable")
		sys.exit()
	# parse output dir
	out_dir=args.out_dir
	safe_make_dir(out_dir)

	# parse OTUs interested
	otu=args.otu
	if(not check_file_exists(otu)):
		print ("Either otu file is missing or is not readable")
		sys.exit()
	
	otu_array = []
	print("Reading otu file... ...\n")
	with open(otu, "r") as f:
		for line in f:
			otu_array.append(line.rstrip())
	print("Read "+ str(len(otu_array)) +" OTUs.\n")
	
	print(otu_array)
	# read seqs fasta file
	seq_record_dict = SeqIO.index(seqs, "fasta")
	
	# read otu_map file
	otu_map_dict={}
	with open(otu_map, "r") as in_otu_map:
		for line in in_otu_map:
			tmp=line.rstrip().split("\t")
			otu_map_dict[tmp[0]]=tmp[1:]
	#print( otu_map_dict.items() )
	# loop find otu's corresponding sequences and write to file named by otu name.
	for otu_tmp in otu_array:
		if(otu_tmp in otu_map_dict):
			seq_list=otu_map_dict[otu_tmp]
			output_handle = open(out_dir+"/"+otu_tmp+".fasta", "w")
			for seq in seq_list:
				if(seq in seq_record_dict):
					out_seq_record_obj=seq_record_dict[seq]
					SeqIO.write(out_seq_record_obj, output_handle, "fasta")
			output_handle.close()

	seq_record_dict.close()


if __name__ == '__main__':
    main()


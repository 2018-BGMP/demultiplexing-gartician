#!/usr/env/bin python3

# Title: Demultiplex
# Description: To demultiplex Illumina reads based on barcodes to quantitatively infer rate of index hopping. 

import argparse
import gzip  
from collections import defaultdict
import numpy as np

##########################################################
# Table of Contents                                      #
##########################################################
# f1 == R1 Reads from 1294_S1_L008_R1_001.fastq.gz       #
##########################################################
# f2 == R2 Barcodes from 1294_S1_L008_R2_001.fastq.gz    #
##########################################################
# f3 == R3 Barcodes from 1294_S1_L008_R3_001.fastq.gz    #
##########################################################
# f4 == R4 Reads from 1294_S1_L008_R4_001.fastq.gz       #
##########################################################

def get_arguments():
	parser = argparse.ArgumentParser(description="Finds average quality score per bp of a FASTQ file.")
	parser.add_argument("-r1", "--file1", type=str, help="Select a standard FASTQ file (gzipped) to demultiplex (raw reads #1)")
	parser.add_argument("-r2", "--file2", type=str, help="Select a standard FASTQ file (gzipped) to demultiplex (barcodes #1)")
	parser.add_argument("-r3", "--file3", type=str, help="Select a standard FASTQC file (gzipped) to demultiplex (barcodes #2)")
	parser.add_argument("-r4", "--file4", type=str, help="Select a standard FASTQC file (gzipped) to demultiplex (raw reads #2)")
	parser.add_argument("-c", "--cutoff", type=int, help="Cutoff score using Phred33 scale")
	return parser.parse_args()

args = get_arguments()
file1 = args.file1
file2 = args.file2
file3 = args.file3
file4 = args.file4
cutoff = args.cutoff

complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
def reverse_complement(seq):
	comp = "".join([complement[x] for x in seq])
	comp = comp[::-1]
	return(comp)
	
def convert_phred(x):
	'''Converts character to quality score using 33 scale for Illumina sequencing.'''
	phred_score = ord(x) - 33 
	return phred_score

YesIndexHop = {'GTAGCGTA': 0, 'CGATCGAT': 0, 'GATCAAGG': 0, 'AACAGCGA': 0, 'TAGCCATG': 0, 'CGGTAATC': 0, 'CTCTGGAT': 0, 'TACCGGAT': 0, 'CTAGCTCA': 0, 'CACTTCAC': 0, 'GCTACTCT': 0, 'ACGATCAG': 0, 'TATGGCAC': 0, 'TGTTCCGT': 0, 'GTCCTAAG': 0, 'TCGACAAG': 0, 'TCTTCGAC': 0, 'ATCATGCG': 0, 'ATCGTGGT': 0, 'TCGAGAGT': 0, 'TCGGATTC': 0, 'GATCTTGC': 0, 'AGAGTCCA': 0, 'AGGATAGC': 0}
NoIndexHop = {'GTAGCGTA': 0, 'CGATCGAT': 0, 'GATCAAGG': 0, 'AACAGCGA': 0, 'TAGCCATG': 0, 'CGGTAATC': 0, 'CTCTGGAT': 0, 'TACCGGAT': 0, 'CTAGCTCA': 0, 'CACTTCAC': 0, 'GCTACTCT': 0, 'ACGATCAG': 0, 'TATGGCAC': 0, 'TGTTCCGT': 0, 'GTCCTAAG': 0, 'TCGACAAG': 0, 'TCTTCGAC': 0, 'ATCATGCG': 0, 'ATCGTGGT': 0, 'TCGAGAGT': 0, 'TCGGATTC': 0, 'GATCTTGC': 0, 'AGAGTCCA': 0, 'AGGATAGC': 0}
summary_dict = {'ContainsAmbiguousNucleotides' : 0, 'ContainsBadQualityScore' : 0, 'NotReverseComplements' : 0}

LN = 0 # line counter 

print("Commence DEMULTIPLEX!")

with gzip.open(file1,"rt") as f1, gzip.open(file2,"rt") as f2, gzip.open(file3,"rt") as f3, gzip.open(file4,"rt") as f4:
	while LN < 1452986940: # integer refers to line number of one file. 
		
		# WORKING STATUS
		LN += 4 # line counter
		if LN % 100000 == 0: # per 100K line...
			print("Processing read number: ", LN/4)

		# READ LINES ACCORDING TO FILE ONE-BY-ONE USING READLINE.STRIP METHOD:
		header1_1 = f1.readline().strip()
		sequence2_1 = f1.readline().strip()
		plus3_1 = f1.readline().strip()
		quality_score4_1 = f1.readline().strip()
		
		header1_2 = f2.readline().strip()
		sequence2_2 = f2.readline().strip()
		plus3_2 = f2.readline().strip()
		quality_score4_2 = f2.readline().strip()
		
		header1_3 = f3.readline().strip()
		sequence2_3 = f3.readline().strip()
		plus3_3 = f3.readline().strip()
		quality_score4_3 = f3.readline().strip()
		
		header1_4 = f4.readline().strip()
		sequence2_4 = f4.readline().strip()
		plus3_4 = f4.readline().strip()
		quality_score4_4 = f4.readline().strip()
		
		# print(sequence2_1, sequence2_2, sequence2_3, sequence2_4, sep=" | ")
		
		# FILTER BAD BARCODES WITH AMBIGUOUS NUCLEOTIDE 'N' FROM INDICES 
		
		if ('N' in sequence2_2) or ('N' in sequence2_3):
			summary_dict['ContainsAmbiguousNucleotides'] += 1
			with open("Read1_ambiguous.fq", "a") as r1_ambiguous, open("Read2_ambiguous.fq", "a") as r2_ambiguous:
				print(header1_1 + "_" + sequence2_2 + "_" + sequence2_3 + "\n" + sequence2_1 + "\n" + "+" + "\n" + quality_score4_1, file=r1_ambiguous) # print to read1_ambiguous 
				print(header1_4 + "_" + sequence2_2 + "_" + sequence2_3 + "\n" + sequence2_4 + "\n" + "+" + "\n" + quality_score4_4, file=r2_ambiguous) # print to read2_ambiguous
				continue # omission of continue will print information to r1_ambiguous and r2_ambiguous n times where n is amount of 'N' in a read.
		
		# FILTER BAD QUALITY SCORES FROM INDICES
		
			# Create a temporary numpy array (index1_QS, index2_QS) filled with quality score of barcode2 and barcode3. This is important because logic flow is controlled by a single 'if quality score < 30' conditional. 
		
		index1_QS = np.zeros(8)
		index2_QS = np.zeros(8)
		QS_counter1 = 0
		QS_counter2 = 0
		
		for x in quality_score4_2: # append barcode2 into index1_QS
			index1_QS[QS_counter1] += int(convert_phred(x))
			QS_counter1 += 1
		for y in quality_score4_3:
			index2_QS[QS_counter2] += int(convert_phred(y))
			QS_counter2 += 1
			
		# print(sequence2_2, index1_QS, index2_QS, sequence2_3, LN) 
		
		if np.amin(index1_QS) < cutoff: # if the minimum score of a read is less than cutoff, print information to r1_ambiguous/r2_ambiguous.
			summary_dict['ContainsBadQualityScore'] += 1
			with open("Read1_ambiguous.fq", "a") as r1_ambiguous, open("Read2_ambiguous.fq", "a") as r2_ambiguous:
				print(header1_1 + "_" + sequence2_2 + "_" + sequence2_3 + "\n" + sequence2_1 + "\n" + "+" + "\n" + quality_score4_1, file=r1_ambiguous) # print to read1_ambiguous 
				print(header1_4 + "_" + sequence2_2 + "_" + sequence2_3 + "\n" + sequence2_4 + "\n" + "+" + "\n" + quality_score4_4, file=r2_ambiguous) # print to read2_ambiguous
				continue
		if np.amin(index2_QS) < cutoff:
			summary_dict['ContainsBadQualityScore'] += 1
			with open("Read1_ambiguous.fq", "a") as r1_ambiguous, open("Read2_ambiguous.fq", "a") as r2_ambiguous:
				print(header1_1 + "_" + sequence2_2 + "_" + sequence2_3 + "\n" + sequence2_1 + "\n" + "+" + "\n" + quality_score4_1, file=r1_ambiguous) # print to read1_ambiguous 
				print(header1_4 + "_" + sequence2_2 + "_" + sequence2_3 + "\n" + sequence2_4 + "\n" + "+" + "\n" + quality_score4_4, file=r2_ambiguous) # print to read2_ambiguous
				continue

		# AFTER FILTERING 'N's AND BAD QUALITY SCORES, IF BARCODE 2 AND 3 ARE REVERSE COMPLEMENTARY SEQUENCES:

		temp_read1 = "Read1_{}_{}.fq".format(sequence2_2, sequence2_3) # Assigning a temporary variable so that output settings look cleaner. 
		temp_read2 = "Read2_{}_{}.fq".format(sequence2_2, sequence2_3) # These variable refer to output files. 
		
		# print(temp_read1, temp_read2, sequence2_2, reverse_complement(sequence2_3))
		
		if (sequence2_2 == reverse_complement(sequence2_3)): # test if barcodes are reverse complement. No need to check if the complement argument is true because this situation is only true or false. 
			if sequence2_2 in list(NoIndexHop.keys()): # assess whether a barcode is the one we provided. 
				NoIndexHop[sequence2_2] += 1 # update that specific barcode by one. 
				with open(temp_read1, "a") as temp_read1, open(temp_read2, "a") as temp_read2: # then output all R1/R4 information to appropriate files with respect to their barcodes.
					print(header1_1 + "_" + sequence2_2 + "_" + sequence2_3, "\n" + sequence2_1 + "\n" + "+" + "\n" + quality_score4_1, file=temp_read1)
					print(header1_4 + "_" + sequence2_2 + "_" + sequence2_3, "\n" + sequence2_4 + "\n" + "+" + "\n" + quality_score4_4, file=temp_read2)
		else: # if barcodes are not reverse complements...
			if sequence2_2 in list(YesIndexHop.keys()): # assess whether a barcode is the one we provided. 
				YesIndexHop[sequence2_2] += 1 # update that specific barcode by one. 
				with open("Read1_ambiguous.fq", "a") as r1_ambiguous, open("Read2_ambiguous.fq", "a") as r2_ambiguous:
					print(header1_1 + "_" + sequence2_2 + "_" + sequence2_3, "\n" + sequence2_1 + "\n" + "+" + "\n" + quality_score4_1, file=r1_ambiguous) # print to read1_ambiguous
					print(header1_4 + "_" + sequence2_2 + "_" + sequence2_3, "\n" + sequence2_4 + "\n" + "+" + "\n" + quality_score4_4, file=r2_ambiguous) # print to read2_ambiguous

# SUMMARY STATISTICS 
import operator
SortedYesIndexHop = sorted(YesIndexHop.items(), key=operator.itemgetter(1), reverse=True)
TotalIndexHop = sum(YesIndexHop.values())
NoIndexHopTotal = sum(NoIndexHop.values())
IndexHopProportion = sum(YesIndexHop.values())/(LN/4) # LN/4 refers to total amount of reads. 

print("\n" + "SUMMARY STATISTICS" + "\n")
print("Bad Reads: ", summary_dict['ContainsAmbiguousNucleotides'])
print("Bad Quality Score: ", summary_dict['ContainsBadQualityScore'])
print("Not Reverse Complements: ", summary_dict['NotReverseComplements'])

print("Total Index Hopping: ", sum(YesIndexHop.values()))
print("Proportion of Index Hopping to Total Reads: ", IndexHopProportion)
print("Percent of Index Hopping to Total Reads: ", IndexHopProportion * 100, "%", "\n")
print("SUMMARY STATISTICS", "\n")

print("All barcodes that have index hopped are in file 'BarcodesWithIndexHop'")
with open("BarcodesWithIndexHop.txt", "w") as BarcodesWithIndexHop:
	print(SortedYesIndexHop, file=BarcodesWithIndexHop)
# SUMMARY STATISTICS
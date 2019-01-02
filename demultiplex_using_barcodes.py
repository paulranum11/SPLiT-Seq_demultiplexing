#!/usr/bin/python
#
# Description
# 1. Pre-process barcode data files based on number of errors allowed, creating dictionaries with the allowed variant barcodes as the keys, and the actual barcodes as the values.
# 2. Walk the fastq data file in segments matching the expected barcode length, performing O(1) dictionary lookup on each segment to get a list of matching barcodes from each barcode file.
# 3. If the read has matching barcodes from at least 1 barcode in barcode file 1, 2, and 3 - save the read in the corresponding fastq results files.

import sys
import getopt
import os
import logging

from datetime import datetime
from itertools import combinations

#####################
# GLOBAL VIARIABLES #
#####################
#
# Track how many reads are added to a given result file. Files not satisfying a minimum criteria will be filtered.
#
statistics = {}
#
# Keys of dictionaries will be values that constitute a match
# Values of dictionaries will be actual barcodes that the matching potentially variant barcodes correspond to.
#
round1barcodeDictionary = {}
round2barcodeDictionary = {}
round3barcodeDictionary = {}
#
# chars that can appear in the data line of the read, used for insertion and substitution operations
#
possibleChars = ['A', 'C', 'G', 'T', 'N']

###########
# CLASSES #
###########
		
#############################

#
# FastQRead definition
#
# round1BarcodeMatches - Barcodes from the 1st round barcode file that match the read
# round2BarcodeMatches - Barcodes from the 2nd round barcode file that match the read
# round3BarcodeMatches - Barcodes from the 3rd round barcode file that match the read
# data - lines of data corresponding to the read (4 line read)
#
class FastQRead(object):
	def __init__(self):
		self.round1barcodeMatches = set()
		self.round2barcodeMatches = set()
		self.round3barcodeMatches = set()
		self.data = []

##################
# HELPER METHODS #
##################

#
# Populate dictionaries with all the possible barcode variations that constitute matches
# 
# barcodeDictionary - global dictionary hash map to populate
# barcodeFile - file containing barcodes to add to the hash set
#
def barcodesToDictionary(barcodeDictionary, barcodeFile):

	file = open(barcodeFile, "r")
	
	for barcode in file:
		barcode = barcode.strip()
		if "LENGTH" not in barcodeDictionary:
			barcodeDictionary["LENGTH"] = []
			for i in range(-errors, errors + 1):
				barcodeDictionary["LENGTH"].append(len(barcode) - i)
			
		addToDictionaryListValue(barcodeDictionary, barcode, barcode)
		if errors != 0:
			barcodeVariantsToDictionary(barcodeDictionary, barcode, barcode, 0, 0)
	file.close()
	
#
# Pre-process barcodes by calculating all of the variant barcodes allowed based on the errors parameter.
# Errors could include insertions, deletions, or substitutions.
# After each error, recurse to handle additional errors.
#
# barcodeDictionary - global dictionary hash map to populate
# barcode - actual barcode present in the barcode file
# variant - current altered barcode (after an insertion, deletion, or substitution operation)
# index - current character of the barcode we are evaluating
# depth - current number of operations that has been performed on the barcode
#
def barcodeVariantsToDictionary(barcodeDictionary, barcode, variant, index, depth):

	#Termination Clause
	if depth >= errors:
		return

	# Perform all possible deletion operations
	for i in range(index, len(variant)):
		newVariant = variant[:i] + '' + variant[i+1:]
		addToDictionaryListValue(barcodeDictionary, newVariant, barcode)
		barcodeVariantsToDictionary(barcodeDictionary, barcode, newVariant, i, depth + 1)
	
	# Perform all possible substitution operations
	for i in range(index, len(variant)):
		currentChar = variant[i]
		for possibleChar in possibleChars:
			if currentChar != possibleChar:
				newVariant = variant[:i] + possibleChar + variant[i+1:]
				addToDictionaryListValue(barcodeDictionary, newVariant, barcode)
				barcodeVariantsToDictionary(barcodeDictionary, barcode, newVariant, i + 1, depth + 1)
				
	# Perform all possible insertion operations
	for i in range(index, len(variant)):
		for possibleChar in possibleChars:
			newVariant = variant[:i] + possibleChar + variant[i:]
			addToDictionaryListValue(barcodeDictionary, newVariant, barcode)
			barcodeVariantsToDictionary(barcodeDictionary, barcode, newVariant, i + 1, depth + 1)
	
	
#
# Either append new value to an existing dictionary value list or initialize a new value list
#
# dictionary - dictionary to populate
# key - key in the dictionary, create if it doesn't already exist.
# value - value to add to the value list at the given key, create a new value list if key did not exist previously
#
#
def addToDictionaryListValue(dictionary, key, value):

	if key in dictionary:
		dictionary[key].add(value)
	else:
		dictionary[key] = set()
		dictionary[key].add(value)

#
# Clear all files from a given directory
#
# directory - directory to clear files from
#
def clearFiles(directory):
	for (dirpath, dirnames, filenames) in os.walk(directory):
		for filename in filenames:
			os.remove(directory + "/" + filename)
	
#
# Clean up any files that didn't meet filter criteria, having less reads than the minimum reads
#
# directory - directory to clean up files from
#
def filterFiles(directory):
	
	for (dirpath, dirnames, filenames) in os.walk(directory):
		for filename in filenames:
			if filename.endswith("fastq"):
				if statistics[filename] < minreads:
					os.remove(directory + "/" + filename)
			
#
# Count the number of fastq files in a specific directory
#
# directory - directory to count files in
# prefix - string to prepend to the console message
#		
def countFiles(directory, prefix):

	count = 0
	for (dirpath, dirnames, filenames) in os.walk(directory):
		for filename in filenames:
			if filename.endswith("fastq"):
				count = count + 1
	print(prefix + str(count) + " files found!")

#
# Fast Q Search Algorithm
#
# Assumes specific file format for matching in blocks of four lines.
# Line 1. Non-relevant data for search, Relevant data for results
# Line 2. Relevant data for search, Relevant data for results
# Line 3. Non-relevant data for search, Relevant data for results
# Line 4. Non-relevant data for search, Relevant data for results
#
def crawlFastQ(fastqr):
	
	startTime = datetime.now()
	counter = 1
	
	with open(fastqr) as infile:
		read = FastQRead()
		for line in infile:
			if len(read.data) == 0:
				read.data.append(line)
			elif len(read.data) == 1:
			
				read.round1barcodeMatches = crawlSegments(line, round1barcodeDictionary, round1barcodeDictionary["LENGTH"])
				
				if read.round1barcodeMatches:
					read.round2barcodeMatches = crawlSegments(line, round2barcodeDictionary, round2barcodeDictionary["LENGTH"])
				
					if read.round2barcodeMatches:
						read.round3barcodeMatches = crawlSegments(line, round3barcodeDictionary, round3barcodeDictionary["LENGTH"])
				
				read.data.append(line)
				
			elif len(read.data) == 2:
				read.data.append(line)
			elif len(read.data) == 3:
			
				read.data.append(line)
				
				saveRead(read)
				del read
				read = FastQRead()
				
				if (counter % 100000) == 0:
					print("Analyzed [" + "{:,}".format(counter) + "] reads in [" + str(datetime.now() - startTime) + "]")
				counter = counter + 1

#
# Helper method to compare line of text against barcode dictionaries
#
# line - line of text to evaluate
# barcodeDictionary - hash map of allowed variant barcodes and actual barcodes
#
def crawlSegments(line, barcodeDictionary, barcodeLengths):

	results = set()

	for barcodeLength in barcodeLengths:
	
		index = 0
	
		while index < len(line) - barcodeLength:
		
			segment = line[index:index+barcodeLength]
			
			if segment in barcodeDictionary:
				for barcode in barcodeDictionary[segment]:
					results.add(barcode)
				
			index = index + 1
	
	return results
	
#
# Save a block of data to permutations of the results files
#
# read - FastQRead object to save to disk
#
def saveRead(read):

	if read.round1barcodeMatches and read.round2barcodeMatches and read.round3barcodeMatches:
		for round1barcodeMatch in read.round1barcodeMatches:
			for round2barcodeMatch in read.round2barcodeMatches:
				for round3barcodeMatch in read.round3barcodeMatches:
				
					filename = round1barcodeMatch + "-" + round2barcodeMatch + "-" + round3barcodeMatch + ".fastq"
					file = open(outputdir + "/" + filename, "a+")
					
					# Track total number of blocks added to each file
					if filename in statistics:
						statistics[filename] += 1
					else:
						statistics[filename] = 1	
					
					for line in read.data:
						file.write(line)
						
					file.close()			
		
def main(argv):

	# Load runtime and default parameters
	global minreads
	minreads = 10
	
	global errors
	errors = 1
	
	global outputdir
	outputdir = "results"

	try:
		opts, args = getopt.getopt(argv, "hm:1:2:3:r:e:o:", ["minreads=", "round1barcodes=", "round2barcodes=", "round3barcodes=", "fastqr=", "errors=", "outputdir="])
	except getopt.GetoptError:
		print("python barcode_trie_crawler.py -b <comma-separated list of barcode files> -r <fastq_R file>")
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print("python barcode_trie_crawler.py -b <comma-separated list of barcode files> -r <fastq_R file>")
			sys.exit()
		elif opt in ("-m", "--minreads"):
			minreads = int(arg)
		elif opt in ("-1", "--round1barcodes"):
			round1barcodes = arg
		elif opt in ("-2", "--round2barcodes"):
			round2barcodes = arg
		elif opt in ("-3", "--round3barcodes"):
			round3barcodes = arg
		elif opt in ("-r", "--fastqr"):
			fastqr = arg
		elif opt in ("-e", "--errors"):
			errors = int(arg)
		elif opt in ("-o", "--outputdir"):
			outputdir = arg
	
	# Print execution parameters
	print("minreads [" + str(minreads) + "]")
	print("round1barcodes [" + round1barcodes + "]")
	print("round2barcodes [" + round2barcodes + "]")
	print("round3barcodes [" + round3barcodes + "]")
	print("fastqr [" + fastqr + "]")
	print("errors [" + str(errors) + "]")
	print("outputdir [" + outputdir + "]")
	
	# Create outputdir directory
	try:  
    		os.mkdir(outputdir)
	except OSError:  
    		print ("Directory already present %s" % outputdir)
	else:  
  			print ("Successfully created the directory %s" % outputdir)
	
	preprocessingStartTime = datetime.now()
	barcodesToDictionary(round1barcodeDictionary, round1barcodes)
	barcodesToDictionary(round2barcodeDictionary, round2barcodes)
	barcodesToDictionary(round3barcodeDictionary, round3barcodes)
	print("Constructured dictionary hash data structure in [" + str(datetime.now() - preprocessingStartTime) + "]")
	
	# Crawl the fastq files referencing dictionary hash for O(1) lookup
	clearFiles(outputdir)
	crawlFastQ(fastqr)
	countFiles(outputdir, "Before Filter: ")
	filterFiles(outputdir)
	countFiles(outputdir, "After Filter: ")
	
if __name__ == "__main__":
	main(sys.argv[1:])
	
		
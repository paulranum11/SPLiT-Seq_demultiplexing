#!/usr/bin/python
#
# Description
# 1. Pre-process barcode data files based on number of errors allowed, creating dictionaries with the allowed variant barcodes as the keys, and the actual barcodes as the values.
# 2. Walk the fastq data file in segments matching the expected barcode length, performing O(1) dictionary lookup on each segment to get a list of matching barcodes from each barcode file.
# 3. If the read has matching barcodes from at least 1 barcode in barcode file 1, 2, and 3 - save the read in the corresponding fastq results files.

import argparse
import sys

import splitseq_utilities

from datetime import datetime
from itertools import islice

##########################
# ADVANCED CONFIGURATION # 
##########################
# The following may move to command line arguments or regular configuration in the future
#
# Explicitly define constraints for where in the read to search for each barcode, by character index
# Will drastically speed up processing time, but could affect correctness if barcodes are
# offset by more than an expected amount.
#
# Supply one number to specify a starting location to begin searching for the barcode
# Supply two numbers separated by a hyphen to specify an inclusive range of locations to search for the 
# barcode in, i.e. 10-11 will expect the barcode to begin at index 10 or index 11
# Indices start at 0 from the left
# 
# i.e. in TAACTCCAATACAGATTCGTGGCCGCTGTTTCGCATCGGCGTACGACTCGAACTTAATCCACGTGCTTGAGAGGCCAGAGCACTGTCTCTTATA
# 'AGCACTGTCTCTTATA' is at index 78
#
# If barcodes are always in the exact same location in the read, a configuration could be, for example
# barcode1constraints = "78-78"
# barcode2constraints = "48-48"
# barcode3constraints = "10-10"
#
barcode1constraints = "10-"
barcode2constraints = "10-"
barcode3constraints = "10-"
#
# Error types to evaluate, possible types are deletion, substitution, and insertion
#
DELETION = 'deletion'
SUBSTITUTION = 'substitution'
INSERTION = 'insertion'
errorTypes = [ SUBSTITUTION ]

#####################
# GLOBAL VIARIABLES #
#####################
#
# Track how many reads are added to a given result file. Files not satisfying a minimum criteria will be filtered.
#
statistics = {}
buffers = {}
#
# Keys of dictionaries will be values that constitute a match
# Values of dictionaries will be actual barcodes that the matching potentially variant barcodes correspond to.
#
round1barcodeDictionary = {}
round2barcodeDictionary = {}
round3barcodeDictionary = {}
#
# Chars that can appear in the data line of the read, used for insertion and substitution operations
#
possibleChars = ['A', 'C', 'G', 'T', 'N']
#
# Constants
#
LENGTH = "LENGTH"
HYPHEN = "-"
FASTQ_EXTENSION = ".fastq"

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

#
# BarcodeMatch definition
#
# barcode - Barcode of the barcode match
# startIndex - Location in the read the barcode starts at
# endIndex - Location in the read the barcode ends at
#        
class BarcodeMatch(object):
    def __init__(self, barcode, startIndex, endIndex):
        self.barcode = barcode
        self.startIndex = startIndex
        self.endIndex = endIndex
    def __eq__(self, other):
        if isinstance(other, BarcodeMatch):
            return self.barcode == other.barcode
        return False
    def __hash__(self):
        return hash(self.barcode)

#
# ReadIndex definition
#
# index - Index we are currently tracking the read.  Glorified pointer.
#        
class ReadIndex(object):
    def __init__(self, index):
        self.index = index

#
# BarcodeIndexConstraint definition
#
# lo - Barcode cannot start before the specified index
# hi - Barcode cannot start after the specified index
#       
class BarcodeIndexConstraint(object):
    def __init__(self, lo, hi):
        if lo != "":
            self.lo = int(lo)
        else:
            self.lo = None
        if hi != "":
            self.hi = int(hi)
        else:
            self.hi = None
    def getLoOrDefault(self, default):
        if self.lo != None:
            return self.lo
        else:
            return default
    def getHiOrDefault(self, default):
        if self.hi != None:
            return self.hi
        else:
            return default

#
# Populate dictionaries with all the possible barcode variations that constitute matches
# 
# barcodeDictionary - global dictionary hash map to populate
# barcodeFile - file containing barcodes to add to the hash set
# errors - number of insertion, deletion, or substitution errors to allow in a match
#
def barcodesToDictionary(barcodeDictionary, barcodeFile, errors):

    file = open(barcodeFile, "r")
    
    for barcode in file:
        barcode = barcode.strip()
        if LENGTH not in barcodeDictionary:
            barcodeDictionary[LENGTH] = []
            if DELETION in errorTypes:
                for i in range(-errors, 1):
                    barcodeDictionary[LENGTH].append(len(barcode) - i)
            barcodeDictionary[LENGTH].append(len(barcode))
            if INSERTION in errorTypes:
                for i in range(0, errors + 1):
                    barcodeDictionary[LENGTH].append(len(barcode) - i)
            
        splitseq_utilities.addToDictionarySet(barcodeDictionary, barcode, barcode)
        if errors != 0:
            barcodeVariantsToDictionary(barcodeDictionary, barcode, barcode, 0, 0, errors, set())
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
# errors - number of insertion, deletion, or substitution errors to allow in a match
# memo - use dynamic programming to increase performance for high error count
#
def barcodeVariantsToDictionary(barcodeDictionary, barcode, variant, index, depth, errors, memo):

    #Termination Clause
    if depth >= errors:
        return
        
    #DP Termination Clause
    memoItem = variant + "-" + str(depth)
    if memoItem in memo:
        return
    
    memo.add(memoItem)
    
    # Support for deletions removed after we discovered it was a primary contributor to erroneous assignment of reads.
    # Can be enabled in the advanced configuration section
    # Perform all possible deletion operations and recurse
    if DELETION in errorTypes:
        for i in range(index, len(variant)):
            newVariant = variant[:i] + '' + variant[i+1:]
            splitseq_utilities.addToDictionarySet(barcodeDictionary, newVariant, barcode)
            barcodeVariantsToDictionary(barcodeDictionary, barcode, newVariant, i, depth + 1, errors, memo)
    
    # Perform all possible substitution operations and recurse
    if SUBSTITUTION in errorTypes:
        for i in range(index, len(variant)):
            currentChar = variant[i]
            for possibleChar in possibleChars:
                if currentChar != possibleChar:
                    newVariant = variant[:i] + possibleChar + variant[i+1:]
                    splitseq_utilities.addToDictionarySet(barcodeDictionary, newVariant, barcode)
                    barcodeVariantsToDictionary(barcodeDictionary, barcode, newVariant, i + 1, depth + 1, errors, memo)
    
    # Support for insertions removed after we discovered it was a primary contributor to erroneous assignment of reads.
    # Can be enabled in the advanced configuration section
    # Perform all possible insertion operations and recurse
    if INSERTION in errorTypes:
        for i in range(index, len(variant)):
            for possibleChar in possibleChars:
                newVariant = variant[:i] + possibleChar + variant[i:]
                splitseq_utilities.addToDictionarySet(barcodeDictionary, newVariant, barcode)
                barcodeVariantsToDictionary(barcodeDictionary, barcode, newVariant, i + 1, depth + 1, errors, memo)
  
#
# addBarcodeIndexConstraints
#
def addBarcodeIndexConstraints(barcodeIndexConstraints, key, constraints):
    
    if constraints != None:
        constraintsArr = constraints.split(HYPHEN, 2)
        if (len(constraintsArr) == 1):
            barcodeIndexConstraints[key] = BarcodeIndexConstraint(constraintsArr[0], None)
        else:
            barcodeIndexConstraints[key] = BarcodeIndexConstraint(constraintsArr[0], constraintsArr[1])
    
#
# FastQ Barcode Search
#
# fastqr - Relative path to the fastqR file
# barcodeIndexConstraints - Constraints for where to barcodes are allowed to appear
# outputdir - Output directory of results files
# targetMemory - Target memory of the script, in bytes
# granularity - Number of reads to evaluate before pausing to analyze memory usage and log progress
#
def crawlFastQ(fastqr, barcodeIndexConstraints, outputdir, targetMemory, granularity):

    startTime = datetime.now()
    counter = 0
    
    with open(fastqr) as f:
        
        while True:
            read = FastQRead()
            read.data = list(islice(f, 4))
            if not read.data:
                break
            
            # Assume the second line of data in a read is what we are searching for in the fastqR file
            line = read.data[1]
            
            readIndex = ReadIndex(len(line) - 1)
            
            loConstraint = barcodeIndexConstraints.get(1).getLoOrDefault(0)
            hiConstraint = barcodeIndexConstraints.get(1).getHiOrDefault(len(line)-1)
            read.round1barcodeMatches = crawlSegments(line, readIndex, loConstraint, hiConstraint, round1barcodeDictionary, round1barcodeDictionary[LENGTH])
            
            if read.round1barcodeMatches:
                
                loConstraint = barcodeIndexConstraints.get(2).getLoOrDefault(0)
                hiConstraint = barcodeIndexConstraints.get(2).getHiOrDefault(len(line)-1)
                read.round2barcodeMatches = crawlSegments(line, readIndex, loConstraint, hiConstraint, round2barcodeDictionary, round2barcodeDictionary[LENGTH])
                
                if read.round2barcodeMatches:
                    loConstraint = barcodeIndexConstraints.get(3).getLoOrDefault(0)
                    hiConstraint = barcodeIndexConstraints.get(3).getHiOrDefault(len(line)-1)
                    read.round3barcodeMatches = crawlSegments(line, readIndex, loConstraint, hiConstraint, round3barcodeDictionary, round3barcodeDictionary[LENGTH])
                    
            addToBuffers(read)

            counter += 1
            if (counter % granularity) == 0:
            
                print("Analyzed [" + "{:,}".format(counter) + "] reads in [{}]".format(datetime.now() - startTime))
                memoryUsage = splitseq_utilities.analyzeMemoryUsage()
                if memoryUsage > targetMemory:
                    print("\tCurrent memory [{}] exceeded target memory [{}]. Flushing buffers...".format(splitseq_utilities.bytesToDisplay(memoryUsage), splitseq_utilities.bytesToDisplay(targetMemory)))
                    splitseq_utilities.flushBuffers(outputdir, buffers)
                    print("Flushed buffers to disk in [" + str(datetime.now() - startTime) + "]")
  
    print("Analyzed [" + "{:,}".format(counter) + "] reads in [" + str(datetime.now() - startTime) + "]. Flushing buffers...")
    splitseq_utilities.flushBuffers(outputdir, buffers)
    print("Flushed buffers to disk in [" + str(datetime.now() - startTime) + "]")


#
# Helper method to compare line of text against barcode dictionaries
#
# line - line of text to evaluate
# readIndex - current location in the read
# loConstraint - barcode should not start before this index
# hiConstraint - barcode should not start after this index
# barcodeDictionary - hash map of allowed variant barcodes and actual barcodes
# barcodeLengths - possible lengths the barcodes could appear as
#
def crawlSegments(line, readIndex, loConstraint, hiConstraint, barcodeDictionary, barcodeLengths):

    results = set()

    # Segments of the file must be analyzed in chunks of all possible barcode lengths.
    # Lengths could vary due to the length of the original barcode as well as deletion and insertion errors.
    
    nextIndex = 0
    for barcodeLength in barcodeLengths:
    
        index = readIndex.index - barcodeLength
        # Ignore everything outside of the barcode index constraints
        if index > hiConstraint:
            index = hiConstraint
        
        while index >= loConstraint:
        
            segment = line[index:index+barcodeLength]
            if segment in barcodeDictionary:
                for barcode in barcodeDictionary[segment]:
                    results.add(BarcodeMatch(barcode, index, index + barcodeLength))
                    
                    if nextIndex < index:
                        nextIndex = index
                
            index = index - 1
    
    readIndex.index = nextIndex
    return results
                    
#
# addToBuffers - Add reads that satisfied the barcode input to the in memory buffer
#
# read - FastQRead object containing data corresponding to a read and matching barcodes for the read
#
def addToBuffers(read):

    if read.round1barcodeMatches and read.round2barcodeMatches and read.round3barcodeMatches:
        
        for round1barcodeMatch in read.round1barcodeMatches:
            for round2barcodeMatch in read.round2barcodeMatches:
                for round3barcodeMatch in read.round3barcodeMatches:
                
                    # Filter out false positives that didn't find barcodes in the proper order
                    if round2barcodeMatch.endIndex > round1barcodeMatch.startIndex:
                        continue
                    if round3barcodeMatch.endIndex > round2barcodeMatch.startIndex:
                        continue
                
                    bufferKey = round1barcodeMatch.barcode + HYPHEN + round2barcodeMatch.barcode + HYPHEN + round3barcodeMatch.barcode + FASTQ_EXTENSION
                    for line in read.data:
                        splitseq_utilities.addToDictionaryList(buffers, bufferKey, line)
                        
                    #Track total number of blocks added to each file
                    if bufferKey in statistics:
                        statistics[bufferKey] += 1
                    else:
                        statistics[bufferKey] = 1
        
def main(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--minreads', required=False, help='Minimum number of reads to keep a results file. Default is 10.', default=10, type=int)
    parser.add_argument('-1', '--round1barcodes', required=True, help='Relative path to the file containing round 1 barcodes.')
    parser.add_argument('-2', '--round2barcodes', required=True, help='Relative path to the file containing round 2 barcodes.')
    parser.add_argument('-3', '--round3barcodes', required=True, help='Relative path to the file containing round 3 barcodes.')
    parser.add_argument('-r', '--fastqr', required=True, help='Relative path the the fastqR file.',)
    parser.add_argument('-e', '--errors', required=False, help='Number of insertion, deletion, or substitution errors allowed in each barcode match. Default is 1.', default=1, type=int)
    parser.add_argument('-o', '--outputdir', required=False, help='Output directory for the results files. Default is results.', default='results')
    parser.add_argument('-t', '--targetMemory', required=False, help='Target memory of the application. If RSS memory exceeds this limit, in memory buffers will be written to disk.', default=256, type=splitseq_utilities.mbToBytes)
    parser.add_argument('-g', '--granularity', required=False, help='Number of reads to evaluate before pausing to evaluate memory usage and log progress.', default=100000, type=int)
    args = parser.parse_args()
    
    # Generate reference data structures
    preprocessingStartTime = datetime.now()
    barcodesToDictionary(round1barcodeDictionary, args.round1barcodes, args.errors)
    barcodesToDictionary(round2barcodeDictionary, args.round2barcodes, args.errors)
    barcodesToDictionary(round3barcodeDictionary, args.round3barcodes, args.errors)
    print("Pre-processing data structures completed in [{}]".format(datetime.now() - preprocessingStartTime))
    
    barcodeIndexConstraints = {}
    addBarcodeIndexConstraints(barcodeIndexConstraints, 1, barcode1constraints)
    addBarcodeIndexConstraints(barcodeIndexConstraints, 2, barcode2constraints)
    addBarcodeIndexConstraints(barcodeIndexConstraints, 3, barcode3constraints)
    
    print("\tCurrent memory [{}]...".format(splitseq_utilities.bytesToDisplay(splitseq_utilities.analyzeMemoryUsage())))
    
    # Analyze fastq file and generate results files
    splitseq_utilities.createDirectory(args.outputdir)
    splitseq_utilities.clearFilesMatchingFilter(args.outputdir, lambda f: True)
    crawlFastQ(args.fastqr, barcodeIndexConstraints, args.outputdir, args.targetMemory, args.granularity)
    
    # Analyze results files post-execution
    print("# of results files [{}]".format(splitseq_utilities.countFilesMatchingFilter(args.outputdir, lambda f: f.endswith("fastq"))))
    splitseq_utilities.clearFilesMatchingFilter(args.outputdir, lambda f: f.endswith("fastq") and statistics[f] < args.minreads)
    print("# of results files [{}]".format(splitseq_utilities.countFilesMatchingFilter(args.outputdir, lambda f: f.endswith("fastq"))))

    
if __name__ == "__main__":
    main(sys.argv[1:])
    
        

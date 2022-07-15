#!/usr/bin/env python

#####
#Description reframe the barcode extraction script to do the following.
#1. Improve speed by storing barcodes by position instead of string searching for them.
#2. Support sequence error handling by filtering and correcting barcodes on the fly using hamming distance.  
#   Filtering this way will improve speed by reducing reads carried forward to downstream steps.
#3. Store matepair information.
#4. Reduce the number of output files created. Depricate the single .fastq file per cell output format.
#####

import sys
import os
import argparse
import itertools
import json
#from Bio.Seq import Seq

# Get arguments
#parser = argparse.ArgumentParser()
#parser.add_argument('-f', '--inputFastqF', required=False, help='Input a forward .fastq format file')
#parser.add_argument('-r', '--inputFastqR', required=False, help='Input a reverse .fastq format file')
#parser.add_argument('-o', '--outputDir', required=False, help='Name of the output directory used to store the output.fastq')
#parser.add_argument('-b', '--bin', required=False, help='Number of reads to process before saving to disc. Binning helps accomodate large input files')
#parser.add_argument('-e', '--errorThreshold', required=False, help='Enter "0" or "1" if to indicate per barcode error threshold')
#parser.add_argument('-p', '--performanceMetrics', required=False, action='store_true', help='Provide -p flag to turn on performance metrics reporting', default=False)
#parser.add_argument('-t', '--readsPerCellThreshold', required=False, help='Provide a minimum reads per cell threshold for retaining a cell', default=1)
#parser.add_argument('-v', '--verbose', required=False, action='store_true', help='Provide -v flag to turn on verbose progress reporting for each bin', default=False)
#parser.add_argument('-k', '--positionDetection', required=True, type=bool, help='Use automated barcode position detection using flanking sequences')
#= parser.parse_)

def demultiplex_fun(inputFastqF, inputFastqR, outputDir, numReadsBin, errorThreshold, performanceMetrics, readsPerCellThreshold, positionDetection, verbose = False):
    #####
    # Set consistent parameters here
    if positionDetection == True:
        Round1_barcode_staticSeq = "CATTCG"
        Round2_barcode_staticSeq = "ATCCAC"
        Round3_barcode_staticSeq = "GTGGCC"
    #####

    #####
    # Define "print" function to print to stderr
    #def eprint(* **kw:
    #    eprint(* file=sys.stderr, **kw
    #####
    
    ######
    # Step1: Gather barcode information from provided barcodes.txt files
    ######
    # Read in text files containing possible barcodes store them in a list
    Round1_barcodes = []
    Round2_barcodes = []
    Round3_barcodes = []
    Eight_BP_barcode = []

    with open('Round1_barcodes_new5.txt', "r") as infile:
        for line in infile:
            Round1Barcode = line.rstrip()
            Round1_barcodes.append(line.rstrip())
            Eight_BP_barcode.append(Round1Barcode[8:16]) #Note this is the same for each round 

    with open('Round2_barcodes_new4.txt', "r") as infile:
        for line in infile:
            Round2_barcodes.append(line.rstrip())

    with open('Round3_barcodes_new4.txt', "r") as infile:
        for line in infile:
            Round3_barcodes.append(line.rstrip())

    # Create a functon to compare barcodes
    # Use hamming distance function from a stack overfow question I created "https://stackoverflow.com/questions/65258822/fast-python-short-dna-string-comparison/65259404#65259404"
    def hamming(s1, s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))


    ######
    # Create bin parameters
    ######
    # Create a value used to bin the dataset
    binIterator = int(numReadsBin * 4)
    print("binIterator is set to " + str(binIterator))
    # Define workingBin
    counter = 0
    readCounterFinal = 0
    workingBin = counter + binIterator

    # Get the number of lines in the input file
    with open(inputFastqF, "r") as infile:
        linesInInputFastq = sum(1 for line in infile)
        print("The linesInInputFastq value is set to " + str(linesInInputFastq))


    ######
    # Build a class to store information from each read.
    ######
    class FastQRead():
        def __init__(self, name, read, quality, lineNumber):
            self.name = name
            self.read = read
            self.quality = quality
            self.lineNumber = lineNumber

        def display_read(self):
            print("name = " + str(self.name) + "\n" \
                "read = " + str(self.read) + "\n" \
                "quality = " + str(self.quality) + "\n" \
                "lineNumber = " + str(self.lineNumber) + "\n", end='')

        def return_fastq(self):
            print(str(self.name) + "\n" \
                + str(self.read) + "\n" \
                + "+" + "\n" \
                + str(self.quality) + "\n", end='')

    class barcodeRead(FastQRead):
        def __init__(self, name, read, quality, lineNumber, barcode1, barcode2, barcode3, umi):
            self.name = name
            self.read = read
            self.quality = quality
            self.lineNumber = lineNumber
            self.barcode1 = barcode1
            self.barcode2 = barcode2
            self.barcode3 = barcode3
            self.umi = umi

        def display_read(self):
            print("name = " + str(self.name) + "\n" \
                "read = " + str(self.read) + "\n" \
                "quality = " + str(self.quality) + "\n" \
                "lineNumber = " + str(self.lineNumber) + "\n" \
                "Barcodes = " + str(barcode1 + "_" + barcode2 + "_" + barcode3) + "\n"\
                "UMI = " + str(umi), end='')

        def return_fastq(self):
            readID = self.name
            readID_split = readID.split("/")
            noSpaceReadID = readID_split[0].replace(" ", "")
            return(str(noSpaceReadID.strip() + "_" + self.barcode1 + self.barcode2 + self.barcode3 + "_" + self.umi) + "\n" \
                + str(self.read) + "\n" \
                + "+" + "\n" \
                + str(self.quality) + "\n")

    ######
    # Create some lists that need to be outside of the loop in order to aggregate performance metrics
    ######
    #Total_barcodes_detected = []
    #Total_barcodes_passing_minReadThreshold = []

    ######
    # Learn barcode positions from input fastqR
    ######
    print("Using default barcode positions...")
    print("Barcode position setting is " + str(positionDetection))
    #Set Default Positions
    umi_start=0
    umi_end=10
    barcode3_start=10
    barcode3_end=18
    barcode2_start=48
    barcode2_end=int(48+8)
    barcode1_start=86
    barcode1_end=int(86+8)

    if bool(positionDetection) == True:
        print("Learning barcode positions to replace defaults...")
        # Code for automated barcode position extractor based on static sequences
        line_ct_Learner = 0
        learner_bc1_list = []
        learner_bc2_list = []
        learner_bc3_list = []
        with open("position_learner_fastqr.fastq", "r") as infile:
            for line in infile:
                if (line_ct_Learner % 4 == 1):
                    if line.find(Round1_barcode_staticSeq) >= 17:
                        learner_bc1_list.append(line.find(Round1_barcode_staticSeq))
                    if line.find(Round2_barcode_staticSeq) >= 17: 
                        learner_bc2_list.append(line.find(Round2_barcode_staticSeq))
                    if line.find(Round3_barcode_staticSeq) >= 17:
                        learner_bc3_list.append(line.find(Round3_barcode_staticSeq))
                line_ct_Learner += 1
            foundPosition_Round1_barcode=max(set(learner_bc1_list), key=learner_bc1_list.count)
            foundPosition_Round2_barcode=max(set(learner_bc2_list), key=learner_bc2_list.count)
            foundPosition_Round3_barcode=max(set(learner_bc3_list), key=learner_bc3_list.count)
            print("Extracted position1 = " + str(foundPosition_Round1_barcode))
            print("Extracted position2 = " + str(foundPosition_Round2_barcode))
            print("Extracted position3 = " + str(foundPosition_Round3_barcode))
            # Use extracted static sequence positions to infer barcode positions
            umi_start=int(foundPosition_Round3_barcode - 18)
            umi_end=int(foundPosition_Round3_barcode - 8)
            print("UMI position has been extracted as " + str(umi_start) + ":" + str(umi_end))
            barcode3_start=int(foundPosition_Round3_barcode - 8) 
            barcode3_end=int(foundPosition_Round3_barcode)
            print("Barcode3 position has been extracted as " + str(barcode3_start) + ":" + str(barcode3_end))
            barcode2_start=int(foundPosition_Round2_barcode - 8)
            barcode2_end=int(foundPosition_Round2_barcode)
            print("Barcode2 position has been extracted as " + str(barcode2_start) + ":" + str(barcode2_end))
            barcode1_start=int(foundPosition_Round1_barcode + 6)
            barcode1_end=int(foundPosition_Round1_barcode + 14)
            print("Barcode1 position has been extracted as " + str(barcode1_start) + ":" + str(barcode1_end))
           


    ######
    # Step2: Iterate through input fastqs in bins.
    ######
    bin_counter = 0
    for i in range(0,int(linesInInputFastq),int(binIterator)):
        # Create dictionaries used to store the parsed read information from the fastq files
        readsF = {}
        readsR = {}
        readsF_BC_UMI_dict = {}
        
        # Create empty lists
        filteredBarcode1 = []
        filteredBarcode2 = []
        filteredBarcode3 = []
        
        #startingline = int(bin_counter * binIterator)
        if (verbose == True):
            print("Processing range " + str(i) + " - " + str(int(i + binIterator)))

        # Iterate through the forward reads
        with open(inputFastqF, "r") as infile:
            start_ct = 0
            line_ct1 = 0
            read_counter = 0 # To get the read counter to match we need to add 1. Each read = 4 lines.
            completeReadCounter = 0
            starterF = 0
            for line in itertools.islice(infile, i, int(i + binIterator)):
                if (starterF == 0 and line.startswith('@') == False):
                    continue
                if (starterF == 0 and line.startswith('@') == True):
                    starterF += 1
                    lineName=str(line[0:].rstrip())
                    lineNameSplit = lineName.split()
                    readIDF = lineNameSplit[0]
                    completeReadCounter += 1
                    line_ct1 += 1
                    continue
                if (line_ct1 % 4 == 0 and line.startswith('@')):
                    lineName=str(line[0:].rstrip())
                    lineNameSplit = lineName.split()
                    readIDF = lineNameSplit[0]
                    completeReadCounter += 1
                if (line_ct1 % 4 == 1):
                    lineRead=str(line[0:].rstrip())
                    completeReadCounter += 1
                if (line_ct1 % 4 == 3):
                    lineQuality=str(line[0:].rstrip())
                    completeReadCounter += 1
                if (completeReadCounter == 3): 
                    processedRead = FastQRead(name = lineName, \
                        read = lineRead, \
                        quality = lineQuality, \
                        lineNumber = read_counter)
                    #print("Adding Fwd Read Key " + str(str(bin_counter) + "_" + str(read_counter)))
                    readsF[readIDF]=processedRead
                    del(readIDF)
                    completeReadCounter = 0
                    read_counter += 1
                if (starterF == 1):
                    line_ct1 += 1

               # if (line.startswith('@') == False and start_ct <=3):
               #     start_ct += 1
               #     line_ct = 0
               #     continue

        # Iterate through the reverse reads
        with open(inputFastqR, "r") as infile:
            start_ct = 0
            line_ct1 = 0
            read_counter = 0
            completeReadCounter = 0
            starterR = 0
            for line in itertools.islice(infile, i, int(i + binIterator)):
                if (starterR == 0 and line.startswith('@') == False):
                    continue
                if (starterR == 0 and line.startswith('@') == True):
                    starterR += 1
                    lineName=str(line[0:].rstrip())
                    lineNameSplit = lineName.split()
                    readIDR = lineNameSplit[0]
                    completeReadCounter += 1
                    line_ct1 += 1
                    continue
                if (line_ct1 % 4 == 0 and line.startswith('@')):
                    lineName=str(line[0:].rstrip())
                    lineNameSplit = lineName.split()
                    readIDR = lineNameSplit[0]
                    completeReadCounter += 1
                if (line_ct1 % 4 == 1):
                    lineRead=str(line[0:].rstrip())
                    #lineReadUMI = lineRead[0:10]
                    lineReadUMI = lineRead[umi_start:umi_end]
                    #lineReadBarcode3 = lineRead[10:18]
                    lineReadBarcode3 = lineRead[barcode3_start:barcode3_end]
                    #lineReadBarcode2 = lineRead[48:int(48+8)]
                    lineReadBarcode2 = lineRead[barcode2_start:barcode2_end]
                    #lineReadBarcode1 = lineRead[86:int(86+8)]
                    lineReadBarcode1 = lineRead[barcode1_start:barcode1_end]
                    filteredBarcode1 = [s for s in Eight_BP_barcode if hamming(s, lineReadBarcode1) <= int(errorThreshold)]  # Match each extracted barcode to a greenlist of possible barcodes.  If a match within hamming distance of 1 is found move forward with that match (not the extracted sequence).
                    filteredBarcode2 = [s for s in Eight_BP_barcode if hamming(s, lineReadBarcode2) <= int(errorThreshold)]
                    filteredBarcode3 = [s for s in Eight_BP_barcode if hamming(s, lineReadBarcode3) <= int(errorThreshold)]
                    completeReadCounter += 1
                if (line_ct1 % 4 == 3):
                    lineQuality=str(line[0:].rstrip())
                    completeReadCounter += 1
                if (completeReadCounter == 3):     
                    if len(filteredBarcode1) == 0:  # The following if statments break the loop if a barcode does not pass the HD <=1 filter
                        line_ct1 += 1  # If the barcode observed does not match a barcode in our greenlist we escape the loop, count the line and reset the complete read counter.
                        completeReadCounter = 0
                        continue
                    elif len(filteredBarcode2) == 0:
                        line_ct1 += 1
                        completeReadCounter = 0
                        continue
                    elif len(filteredBarcode3) == 0:
                        line_ct1 += 1
                        completeReadCounter = 0
                        continue
                    else:
                        processedRead = barcodeRead(name = lineName, \
                        read = lineRead, \
                        quality = lineQuality, \
                        lineNumber = int(read_counter + int(bin_counter * binIterator)), \
                        barcode1 = filteredBarcode1[0], \
                        barcode2 = filteredBarcode2[0], \
                        barcode3 = filteredBarcode3[0], \
                        umi = lineReadUMI)
                        #print("Adding Rev Read Key " + str(str(bin_counter) + "_" + str(read_counter)))
                        readsR[readIDR]=processedRead
                        del(readIDR)
                        completeReadCounter = 0  # Reset the complete read counter to 0 because a read has been completely extracted and stored in the dictionary.      
                    read_counter += 1 # The read counter should progress even if a read does not satisfy the following criteria for being retained.
                if (starterR == 1):
                    line_ct1 += 1  # There are 4 lines for each read in the fastq file. The line counter needs to progress even if no "if" statements are satisfied.


    ######
    # Step3: Transfer barcode from reverse read to forward read. Store new combined reads in dictionary as instances of class barcodedRead.
    ######
        #Pass Barcode information from reverse read to forward read. F and R parts of each read share the same key.
        # This needs to be done because in SPLiT-Seq the reverse read only contains the barcode
        # while the forward read contains the gene expression information.
        # So we need to create a read that contains the barcode information appended to the read name 
        # and the forward read information stored as the actual read.
        #print("beginning transfer of barcodes")
        for key in readsR.keys():
            readF_BC_UMI = barcodeRead(name = readsF[key].name, \
                read = readsF[key].read, \
                quality = readsF[key].quality, \
                lineNumber = readsF[key].lineNumber, \
                barcode1 = readsR[key].barcode1, \
                barcode2 = readsR[key].barcode2, \
                barcode3 = readsR[key].barcode3, \
                umi = readsR[key].umi)
            readsF_BC_UMI_dict[key]=readF_BC_UMI


    ######
    # Step5: Write readF_BC_UMI reads to a .fastq file in outputDir directory.
    ######
#        if not os.path.exists(outputDir):
#            os.makedirs(outputDir)
        # Write the stored reads to disc by appending to the file "MergedCells_1.fastq"
        file1 = open(str(outputDir + "/MergedCells_1.fastq"), "a")
        for key in set(readsF_BC_UMI_dict.keys()):
            file1.write(readsF_BC_UMI_dict[key].return_fastq())
        file1.close() 
        
        bin_counter += 1

        # Here we can return the stored reads to standard out to confirm our read storage program is working as expected.
        #for key in set(readsF_BC_UMI_dict.keys()):
        #    #readsF[key].return_fastq()
        #    readsF_BC_UMI_dict[key].return_fastq()
        #bin_counter += 1
        
        # Flush stdout buffers
        sys.stdout.flush()


def calc_demux_results(outputDir, performanceMetrics, readsPerCellThreshold ):
    ######
    # Step6: Report final total barcodes observed and barcodes passing min read threshold.
    ######
    # Read in
    line_ct2 = 0
    final_counting_dict = {}
    if (performanceMetrics == True):
        with open(str(outputDir + "/MergedCells_1.fastq"), "r") as infile:
            for line in infile:
                if (line_ct2 % 4 == 0):
                    lineName=line.rstrip()
                    lineNameComponentsList=lineName.split("_")
                    cellBarcode=lineNameComponentsList[1]
                    cellUMI=lineNameComponentsList[2]
                    if (str(cellBarcode) not in final_counting_dict.keys()):
                        final_counting_dict[str(cellBarcode)]=[]
                        final_counting_dict[str(cellBarcode)].append(str(cellUMI))
                    else:
                        final_counting_dict[str(cellBarcode)].append(str(cellUMI))
                line_ct2 += 1

        totalFinalCellsDetected = len(final_counting_dict.keys())

        #Calculate number of cells that meet the min reads-per-cell threshold
        filtered_counting_dict = {}
        for key in final_counting_dict.keys():
            if (len(set(final_counting_dict[key])) >= int(readsPerCellThreshold)):
                filtered_counting_dict[key]=final_counting_dict[key]
        
        #Use the filtered_counting_dict keys to remove reads from cells containing < the minimum reads-per-cell threshold
        outfile = open(str(outputDir + "/MergedCells_passing.fastq"), "a")
        delete_counter = 0
        line_ct3 = 0
        with open(str(outputDir + "/MergedCells_1.fastq"), "r") as infile:
            for line in infile:
                if delete_counter == 3:
                    delete_counter = 0
                    line_ct3 += 1
                    continue
                if delete_counter == 0:
                    if (line_ct3 % 4 == 0):
                        lineName=line.rstrip()
                        lineNameComponentsList=lineName.split("_")
                        cellBarcode=lineNameComponentsList[1]
                        cellUMI=lineNameComponentsList[2]
                        if (cellBarcode not in filtered_counting_dict.keys()):
                            delete_counter += 1
                            line_ct3 += 1
                            continue
                        else:
                            outfile.write(str(line.strip() + "\n"))
                            line_ct3 += 1
                            continue
                    if (line_ct3 % 4 == 1):
                        outfile.write(str(line.strip() + "\n"))
                        line_ct3 += 1
                        continue
                    if (line_ct3 % 4 == 2):
                        outfile.write(str(line.strip() + "\n"))
                        line_ct3 += 1
                        continue
                    if (line_ct3 % 4 == 3):
                        outfile.write(str(line.strip() + "\n"))
                        line_ct3 += 1
                        continue
                else:
                    delete_counter += 1
                    line_ct3 += 1
                        
        outfile.close()
                        

        filteredFinalCellsDetected = len(filtered_counting_dict.keys())
        print("The total number of unique barcodes detected was " + str(totalFinalCellsDetected))
        print("The number of barcodes passing the minimum UMI threshold of " + str(readsPerCellThreshold) + " was " + str(filteredFinalCellsDetected))

        # Flush stdout buffers
        sys.stdout.flush()
  
if __name__ == '__main__':
    demultiplex_fun()
    calc_demux_results()

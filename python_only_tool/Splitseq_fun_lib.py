#!/usr/bin/env python

# This function splits a fastq file based on a provided "split_num". 
# It will not break apart reads.

def split_fastqF_fun (split_num, fastq_file, lengthFastq):
    # Get input fastq file dimensions
    print("Getting input fastqF file dimensions")
    length_fastq = int(lengthFastq)
    print("Fastq length is " + str(length_fastq))
    split_size = length_fastq / int(split_num)
    print("The split size before tuning is " + str(split_size))
    for i in range(10000):
        if (split_size % 4 != 0):
           print("Trying split size " + str(split_size))
           split_size += 1
    if (split_size % 4 != 0):
        print("WARNING!!! Unable to split input fastq without read loss, please try again with a different number of cores.")
    print("The split size after tuning is " + str(split_size))

    # Iterate through input fastq file writing lines to outfile in bins.
    print("Begin splitting the input fastq file into bins for parallel processing")
    counter = 0
    split_counter = 0
    split_fastq_list = []
    bin_counter = 0
    with open(fastq_file, "r") as infile:
        for line in infile:
            #print(counter)
            #if counter == 0 and line[0] != "@":
            #    continue
            if counter == 0:
                filename = str("./split_fastq_F_" + str(split_counter))
                split_fastq_list.append(filename)
                outfile = open(filename, "a")
                outfile.write(str(line.strip() + "\n"))
                counter += 1
            elif counter < split_size:          
                outfile.write(str(line.strip() + "\n"))
                counter += 1
            else:
                counter = 0
                split_counter += 1
                outfile.close()
                filename = str("./split_fastq_F_" + str(split_counter))
                split_fastq_list.append(filename)
                outfile = open(filename, "a")
                outfile.write(str(line.strip() + "\n"))
                counter += 1
        outfile.close() 


def split_fastqR_fun (split_num, fastq_file, lengthFastq):
    # Get input fastq file dimensions
    print("Getting input fastqF file dimensions")
    length_fastq = int(lengthFastq)
    print("Fastq length is " + str(length_fastq))
    split_size = length_fastq / int(split_num)
    print("The split size before tuning is " + str(split_size))
    for i in range(10000):
        if (split_size % 4 != 0):
           print("Trying split size " + str(split_size))
           split_size += 1
    if (split_size % 4 != 0):
        print("WARNING!!! Unable to split input fastq without read loss, please try again with a different number of cores.")
    print("The split size after tuning is " + str(split_size))


    # Iterate through input fastq file writing lines to outfile in bins.
    print("Begin splitting the input fastq file into bins for parallel processing")
    counter = 0
    split_counter = 0
    split_fastq_list = []
    with open(fastq_file, "r") as infile:
        for line in infile:
            #print(counter)
            #if counter == 0 and line[0] != "@":
            #    continue
            if counter == 0:
                filename = str("./split_fastq_R_" + str(split_counter))
                split_fastq_list.append(filename)
                outfile = open(filename, "a")
                outfile.write(str(line.strip() + "\n"))
                counter += 1
            elif counter < split_size:          
                outfile.write(str(line.strip() + "\n"))
                counter += 1
            else:
                counter = 0
                split_counter += 1
                outfile.close()
                filename = str("./split_fastq_R_" + str(split_counter))
                split_fastq_list.append(filename)
                outfile = open(filename, "a")
                outfile.write(str(line.strip() + "\n"))
                counter += 1
        outfile.close()

    # Get the first 1000 reads and write to file to create "position_learner_fastqr.fastq"
    filename = "position_learner_fastqr.fastq"
    outfile = open(filename, "w")
    counter = 0
    with open(fastq_file, "r") as infile:
        for line in infile:
            if (counter <= 1000):
                outfile.write(str(line.strip() + "\n"))
                counter += 1
            else:
                break
    outfile.close()


def remove_file_fun (filename):
    import os
    os.remove(filename)

def run_star_alignment_fun (numCores, starGenome, resultsDir):
    import os
    os.system(str("pushd " + resultsDir))
    command_str = str("STAR --runThreadN " + numCores + 
            " --readFilesIn MergedCells_passing.fastq " + 
            "--outFilterMismatchNoverLmax 0.05 " + 
            "--genomeDir " + starGenome +
            " --alignIntronMax 20000 " +
            "--outSAMtype BAM SortedByCoordinate")
    os.system(command_str)
    os.system("popd")

def run_featureCounts_SAF_fun (numCores, countsFile, resultsDir):
    import os
    os.system(str("pushd " + resultsDir))
    command_str = str("featureCounts -F SAF " +
            "-a " + countsFile +
            " -o gene_assigned " +
            "-R BAM Aligned.sortedByCoord.out.bam " +
            "-T " + numCores +
            " -M")
    os.system(command_str)
    os.system("popd")

def run_samtools_fun (resultsDir):
    import os
    os.system(str("pushd " + resultsDir))
    command_str = str("samtools sort Aligned.sortedByCoord.out.bam.featureCounts.bam -o assigned_sorted.bam")
    command_str2 = str("samtools index assigned_sorted.bam")
    os.system(command_str)
    os.system(command_str2)
    os.system("popd")

def run_umi_tools_fun (resultsDir):
    import os
    os.system(str("pushd " + resultsDir))
    command_str = str("umi_tools count --wide-format-cell-counts --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I assigned_sorted.bam -S counts.tsv.gz")
    os.system(command_str)
    os.system("popd")

if __name__ == '__main__':
    split_fastqF_fun()
    split_fastqR_fun()
    remove_file_fun()

#!/usr/bin/env python

# This function splits a fastq file based on a provided "split_num". 
# It will not break apart reads.

def split_fastqF_fun (split_num, fastq_file, lengthFastq):
    # Get input fastq file dimensions
    print("Getting input fastqF file dimensions")
    length_fastq = lengthFastq
    print("Fastq length is" + str(length_fastq))
    split_size = length_fastq / split_num
    while (split_size  % 4 != 0):
        split_size += 1
    print(split_size)

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
    length_fastq = lengthFastq
    print("Fastq length is" + str(length_fastq))
    split_size = length_fastq / split_num
    while (split_size  % 4 != 0):
        split_size += 1
    print(split_size)

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

if __name__ == '__main__':
    split_fastqF_fun()
    split_fastqR_fun()

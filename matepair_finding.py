import sys
from datetime import datetime
import argparse
import os
from os import fsencode
import re
from itertools import islice
import splitseq_utilities

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='Directory containing input files')
    parser.add_argument('-f', '--fastqf', required=True, help='Location of fastq_f file')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-t', '--targetMemory', required=False, help='Target memory of the application. If RSS memory exceeds this limit, in memory buffers will be written to disk.', default=256, type=splitseq_utilities.mbToBytes)
    parser.add_argument('-g', '--granularity', required=False, help='Number of reads to evaluate before pausing to evaluate memory usage and log progress.', default=100000, type=int)
    args = parser.parse_args()

    directory = fsencode(args.input)
    read_id_pattern = r'(^@[^ ]+)'
    d = {}
    i = 0

    print('> building map')

    for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
        # Skip directories that may be lingering from other steps/versions
        if os.path.isdir(filepath):
            continue
        with open(filepath, 'r') as f:
            for line in f:
                if re.match(read_id_pattern, line):
                    read_id = re.search(read_id_pattern, line).groups()[0]
                    d.setdefault(read_id, []).append(filename)

    print('> map complete')
    startTime = datetime.now()
    
    buffers = {}

    with open(args.fastqf, 'r') as f:
    	
        while True:
            batch = list(islice(f, 4))
            if not batch:
                break

            read_id = re.search(read_id_pattern, batch[0]).groups()[0]

            if read_id in d:
                for filename in d[read_id]:
                    out_path = os.path.join(args.output,
                                            '%s-MATEPAIR' % filename.decode('utf-8'))
                    for l in batch:
                        splitseq_utilities.addToDictionaryList(buffers, out_path, l)

            i += 1
            if i % args.granularity == 0:
                print("Analyzed [" + "{:,}".format(i) + "] reads in [" + str(datetime.now() - startTime) + "]")
                memoryUsage = splitseq_utilities.analyzeMemoryUsage()
                if memoryUsage > args.targetMemory:
                    print("\tCurrent memory [{}] exceeded target memory [{}]. Flushing buffers...".format(splitseq_utilities.bytesToDisplay(memoryUsage), splitseq_utilities.bytesToDisplay(args.targetMemory)))
                    splitseq_utilities.flushBuffers('', buffers)

        print("Analyzed [" + "{:,}".format(i) + "] reads in [" + str(datetime.now() - startTime) + "]")		
        splitseq_utilities.flushBuffers('', buffers)
                
if __name__ == '__main__':
    main()

sys.exit()
                


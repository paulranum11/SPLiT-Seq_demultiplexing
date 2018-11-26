import argparse
import os
import re
from itertools import islice

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='Directory containing input files')
    parser.add_argument('-f', '--fastqf', required=True, help='Location of fastq_f file')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    args = parser.parse_args()

    directory = os.fsencode(args.input)
    read_id_pattern = r'(^@[^ ]+)'
    d = {}
    i = 0

    print('> building map')

    for filename in os.listdir(directory):
        with open(os.path.join(directory, filename), 'r') as f:
            for line in f:
                if re.match(read_id_pattern, line):
                    read_id = re.search(read_id_pattern, line).groups()[0]
                    d.setdefault(read_id, []).append(filename)

    print('> map complete')

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
                    with open(out_path, 'a') as out_f:
                        for l in batch:
                            out_f.write(l)

            i += 1
            if i % 100000 == 0:
                print('> %d' % i)

if __name__ == '__main__':
    main()

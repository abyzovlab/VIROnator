#/usr/local/biotools/python/3.10.7/bin/python

import sys

def main(reads_ids_file):
    with open(reads_ids_file) as f:
        reads_ids = set(line.strip() for line in f)

    # Read from stdin (the output of samtools view)
    for line in sys.stdin:
        if line.startswith('@'):
            print(line, end='')
            continue
        # Extract the read ID from the SAM/BAM line (assuming it's the first column)
        read_id = line.split('\t')[0]
        if read_id in reads_ids:
            print(line, end='')

if __name__ == "__main__":
    main(sys.argv[1])



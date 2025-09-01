### Remove reads aligned to impLbCas12a/PCR cassette from the demux/trimmed fastq.gz file ###
### Author: Dan Lou ###
### 2024.07.19 ###

import gzip
import argparse


def remove_cas12a_cassette(input_fastq, output_fastq, exclude_list):
    # Read names to exclude
    with open(exclude_list) as f:
        exclude_reads = set(line.strip() for line in f)

    # Open input and output files
    with gzip.open(input_fastq, 'rt') as infile, gzip.open(output_fastq, 'wt') as outfile:
        while True:
            # Read four lines (one FASTQ entry)
            header = infile.readline()
            if not header:
                break  # End of file
            sequence = infile.readline()
            plus = infile.readline()
            quality = infile.readline()

            # Extract the read name from the header
            read_name = header.split()[0][1:]

            # Write the entry to the output file if the read name is not in the exclude list
            if read_name not in exclude_reads:
                outfile.write(header)
                outfile.write(sequence)
                outfile.write(plus)
                outfile.write(quality)

def main():
    parser = argparse.ArgumentParser(description='Filter a FASTQ file by excluding reads with specific names.')
    parser.add_argument('input_fastq', help='Path to the input FASTQ.gz file')
    parser.add_argument('output_fastq', help='Path to the output filtered FASTQ.gz file')
    parser.add_argument('exclude_list', help='Path to the file containing read names to exclude')

    args = parser.parse_args()

    remove_cas12a_cassette(args.input_fastq, args.output_fastq, args.exclude_list)

if __name__ == '__main__':
    main()


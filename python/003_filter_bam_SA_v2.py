### Remove HDR reads with supplementary alignment to M1M1, m1M1 of the same gene ###
### Author: Dan Lou, Shengdi Li ###
### 2024.08.22 ###
######## check for SA:Z instead of the 15th column, since the 15th column is not necessarily always the "SA:Z"#########
######## get exact the same output with Shengdi's "filter_bam_SA.pl"  ############

import argparse
import re

def parse_arguments():
    """Set up the argument parser and parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Filter BAM file')
    parser.add_argument('-in', dest='input_sam', required=True, help='Input SAM file')
    parser.add_argument('-out', dest='output_sam', required=True, help='Output SAM file path')
    parser.add_argument('-read', dest='rm_read', required=True, help='Output read file path')
    parser.add_argument('-samid', dest='sample_id', required=True, help='Sample ID')
    parser.add_argument('-label', default='M1M1', 
                        help='Comma-separated labels. Reads with these labels in their supplementary alignment '
                             'will be discarded when found in chromosomes containing the label. '
                             'Default: "m1M1,M1M1"')
    return parser.parse_args()

def process_sam_file(args):
    """Process the input SAM file according to the provided arguments."""
    labels = args.label.split(',')

    with open(args.input_sam, 'r') as infile, \
         open(args.output_sam, 'w') as outfile1, \
         open(args.rm_read, 'w') as outfile2:
        
        for read_line in infile:
            # Write header lines directly to output SAM file
            if read_line.startswith('@'):
                outfile1.write(read_line)
                continue
            
            # Process read lines
            read_columns = read_line.split('\t')
            primary_chr = read_columns[2]
            supplementary_info = next((col for col in read_columns if col.startswith('SA:Z')), None)
            
            if primary_chr == "*":  # Handle unmapped reads
                outfile1.write(read_line)
                continue
            
            # Process HDR-mapped reads
            if "_HDR" in primary_chr:
                primary_chr = primary_chr.replace("_HDR", "")
                if supplementary_info:
                    sa_chr = supplementary_info.split(':')[2]
                    if any(re.search(f"{primary_chr}_{label}", sa_chr) for label in labels):
                        outfile2.write(f"{args.sample_id}\t{read_columns[0]}\t{primary_chr}\t{sa_chr}\n")
                    else:
                        outfile1.write(read_line)
                else:
                    outfile1.write(read_line)
            else:
                outfile1.write(read_line)

def main():
    args = parse_arguments()
    process_sam_file(args)

if __name__ == "__main__":
    main()



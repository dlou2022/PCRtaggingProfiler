### consider HDR reads that are not completely mapped to the M1(90bp) region as NHEJ deletion, write to a new sam file ###
### Author: Dan Lou, Shengdi Li ###
### 2024.07.18 ###


import sys
import argparse
import pysam

def main():
    parser = argparse.ArgumentParser(description="Identify reads that are not completely mapped to M1 region or with large INDELs.")
    parser.add_argument("-in", dest="input_bam", required=True, help="Input BAM file.")
    parser.add_argument("-hdr", dest="hdr_bam", required=True, help="Output filtered HDR file.")
    parser.add_argument("-indel", dest="indel_bam", required=True, help="Output BAM file for indel reads.")
    #parser.add_argument("-indel_length", dest="indel_threshold", type=int, default=10, help="INDEL length threshold.")
    #parser.add_argument("-softclip_length", dest="soft_clip_threshold", type=int, default=10, help="Left soft-clipping length threshold.")
    args = parser.parse_args()

    # Open the input BAM file and output BAM file using pysam
    with pysam.AlignmentFile(args.input_bam, 'rb') as input_file, \
         pysam.AlignmentFile(args.hdr_bam, 'wb', template=input_file) as output_hdr, \
         pysam.AlignmentFile(args.indel_bam, 'wb', template=input_file) as output_indel:
        
        # Process each read in the input SAM file
        for read in input_file:
            # Check if the read maps to a chromosome containing "_HDR" and starts at a position > 5
            if "_HDR" in read.reference_name:
                has_large_indel = False
                for operation, length in read.cigartuples:
                    if (operation == 1 or operation == 2) and length > 10:  # 1 is insertion, 2 is deletion
                        has_large_indel = True
                        break

                if has_large_indel:
                    output_indel.write(read)
                elif read.reference_start > 10:
                    output_indel.write(read)
                elif read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > 10:  # 4 is soft clipping
                    output_indel.write(read)
                else:
                    output_hdr.write(read)
            else:
                output_hdr.write(read)

if __name__ == "__main__":
    main()



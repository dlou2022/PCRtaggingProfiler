### extract >10bp soft-clipping sequences from low MAPQ reads ###
### Author: Dan Lou###
### 2024.05.27 ###


import re
import sys
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def main(sam_file, output_fastq):
    soft_clips = []

    # Open SAM file
    with open(sam_file, "r") as sam, open(output_fastq, "w") as out_fastq:
        # Iterate over each line in the SAM file
        for line in sam:
            # Skip header lines
            if not line.startswith('@'):
                fields = line.strip().split('\t')
                read_id = fields[0]
                PA_chr = fields[2]
                PA_start = int(fields[3])
                cigar = fields[5]
                seq = fields[9]
                qual = fields[10]
                flag = int(fields[1])

                # conditions: Primary aligned to certain gene, 0 < align_start_position < 90, has secondary alignments, and flag < 2048
                if PA_chr.endswith("_HDR") and 0 < PA_start < 90 and any('XA:Z' in field for field in fields) and flag < 2048:
                    PA_chr = PA_chr.replace("_HDR", "")
                    read_id = f"{read_id}:{PA_chr}"  # add the tag name to read_id description

                    # Check for soft-clipped reads
                    if 'S' in cigar:
                        digits = re.findall('[0-9]+', cigar)
                        digits = [int(digit) for digit in digits]
                        chars = re.findall('[A-Z]+', cigar)
                        cigar_ops = list(zip(digits, chars))
                        soft_clip_start = 0
                        for length, op in cigar_ops:
                            if op == 'S':
                                soft_clip_end = soft_clip_start + length
                                soft_clip = seq[soft_clip_start:soft_clip_end]
                                if len(soft_clip) > 10:  # Check if soft-clipped sequence is longer than 10bp
                                    qual_scores = qual[soft_clip_start:soft_clip_end]
                                    record = SeqRecord(Seq(soft_clip), id=read_id, description='', letter_annotations={"phred_quality": [ord(c) - 33 for c in qual_scores]})
                                    SeqIO.write(record, out_fastq, "fastq")
                            soft_clip_start += length



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert MAPQ=0 reads into fastq.')
    parser.add_argument('-in', dest='input_sam', type=str, required=True, help='Input sam file.')
    parser.add_argument('-out', dest='output_fastq', type=str, required=True, help='Output fastq file.')
    args = parser.parse_args()

    main(args.input_sam, args.output_fastq)

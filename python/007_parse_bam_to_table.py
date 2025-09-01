import argparse
import re

def parse_cigar(cigar):
    # Parses the CIGAR string to calculate the length of the alignment
    import re
    matches = re.findall(r'(\d+)([MIDSH])', cigar)
    length = 0
    for length_str, type_char in matches:
        if type_char in 'MD':  # Consider only M, D, =, and X for the alignment length
            length += int(length_str)
    return length




def is_reverse_strand_mapped(flag):
    # Bitwise check for the reverse strand (0x10 bit in the flag)
    return (flag & 0x10) != 0

def main(input_sam, output_table, sample_id):
    with open(input_sam, 'r') as infile, open(output_table, 'w') as outfile:
        for read_line in infile:
            read_line = read_line.strip()
            if not read_line.startswith('@'):
                read_column = read_line.split('\t')
                id, description = read_column[0].split(":")
                PA_chr = read_column[2]
                flag = int(read_column[1])
                cigar = read_column[5]
                coord = int(read_column[3])
                seq = read_column[9]
                if PA_chr != "*":
                    if is_reverse_strand_mapped(flag) and flag < 2048:
                        coord = coord + parse_cigar(cigar) - 1
                        strand = '-'
                    else:
                    	strand = '+'
                    outfile.write(f"{sample_id}\t{id}\t{description}\t{PA_chr}\t{strand}\t{coord}\t{parse_cigar(cigar)}\t{len(seq)}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter SAM file and write to a table.')
    parser.add_argument('-in', dest='input_sam', type=str, required=True, help='Input SAM file.')
    parser.add_argument('-out', dest='output_table', type=str, required=True, help='Output table.')
    parser.add_argument('-samid', dest='sample_id', type=str, required=True, help='Sample ID.')
    args = parser.parse_args()

    main(args.input_sam, args.output_table, args.sample_id)

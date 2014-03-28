#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

# Extract samples from a VCF without requiring that the VCF be properly formatted
# Example usage:
# cat my.vcf | subset-vcf-samples.py --sn sample1 sample2 > samples1and2.vcf

import sys
import argparse

# read from file or default to stdin; output to file or default to stdout
parser = argparse.ArgumentParser(description='Extract samples from a VCF.')
parser.add_argument('--vcf', nargs='?', type=argparse.FileType('r'),
                   default=sys.stdin)
parser.add_argument('--out', nargs='?', type=argparse.FileType('wb'),
                   default=sys.stdout)
parser.add_argument('--sn', dest='sn', nargs='+', type=str,
                   help='Sample names (space separated)' )
args = parser.parse_args()

for line in args.vcf.readlines():
    # keep ## header lines as is
    if line[:2] == '##':
        args.out.write(line)
        continue
    # all other lines must be split on tab
    cols = line.strip().split("\t")
    # extract #CHROM line for futher processing
    if line[0] == '#':
        keep_indices = [i for i in range(len(cols)) if cols[i] in args.sn]
    # now for all rows including the #CHROM row, keep desired cols
    new_cols = cols[:9] + [cols[i] for i in keep_indices]
    new_row = "\t".join(new_cols)+"\n"
    args.out.write(new_row)



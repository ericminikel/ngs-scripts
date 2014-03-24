#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

# cptoigv.py
# Given a list of BAMs and loci, make for IGV_plotter
# a samples.txt and sample_locus.txt file corresponding
# to the cartesian product of all BAMs at all loci.

import sys
import re
import argparse
import os

parser = argparse.ArgumentParser(description='Produce sample and sample locus files.')
parser.add_argument('--bams', dest='bams', nargs='+', 
    type=str, help='BAMs (full paths, space separated)' )
parser.add_argument('--loci', dest='loci', nargs='+',
    type=str, help='Loci (chr:pos+buffer, space separated)' )

args = parser.parse_args()

sfile = open('samples.txt',mode='wb')
slfile = open('sample_locus.txt',mode='wb')

for bam in args.bams:
    samplename = os.path.basename(bam)[:-4] # remove directory and .bam extension
    sfile.write(samplename+"\t"+bam+"\n")
    for locus in args.loci:
        samplelocusid = samplename+locus
        slfile.write(samplelocusid+"\t"+locus+"\t"+samplename+"\n")

sfile.close()
slfile.close()
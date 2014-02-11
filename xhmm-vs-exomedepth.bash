#!/bin/bash

cd /humgen/atgu1/fs03/eminikel/048muscle/data
# grab muscle disease sample names
cat MD1.existent.bam.list | awk -F"[/\.]" '{print $7}' > MD1.sample.names
# turn them all into a greppable string. otherwise no way to grep -f with ^ and the XHMM file is to large to grep without ^
grepstring=`cat MD1.sample.names | awk '{printf $0"|"}' | awk '{print "^("$0"MADEUPSAMPLENAME)"}'`
# MADEUPSAMPELNAME is some very sloppy coding - the first awk prints a trailing | which breaks the grep string,
# and it seemed quicker to just add something that will never match than to fix the underlying problem.
# test it out
cat MD1.sample.names | egrep $grepstring | wc -l
# subset those samples from the xhmm calls
cat /broad/hptmp/ebanks/xhmm/xhmm.calls.PCA_normalized.filtered.sample_zscores.RD.txt | egrep $grepstring > muscle.xhmm.calls.txt

# group all the ExomeDepth CNV calls into one big file
echo -n "samplename,rowno" > muscle.ed.calls.csv # create header row
head -1 X9C_DH_1.bam.csv >> muscle.ed.calls.csv # create header row cont'd
for fname in *.bam.csv
do
    # get sample name from file name
    sname=${fname%.bam.csv} # strip .bam.csv from end
    sname=${sname:1} # strip X from front
    cat $fname | tail -n +2 | awk -v sname=$sname '{print sname","$0}' >> muscle.ed.calls.csv
done



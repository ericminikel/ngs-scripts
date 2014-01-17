#!/bin/bash

# to run:
# bsub -o /humgen/atgu1/fs03/eminikel/sandbox/runXHMMjoboutput.out "bash runXHMM.bash -v -o /humgen/atgu1/fs03/eminikel/sandbox/2 -b /humgen/atgu1/fs03/eminikel/sandbox/bam.list"

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
outdir=""
verbose=0
bamlist=""

while getopts "vo:b:" opt; do
    case "$opt" in
    v)  verbose=1
        ;;
    o)  outdir=$OPTARG
        ;;
    b)  bamlist=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

if test $verbose
  then
    echo "verbose=$verbose, output_directory='$outdir', unused arguments: $@"
fi

refdir=

# cat /humgen/gsa-hpprojects/GATK/bundle/2.8/b37/Broad.human.exome.b37.interval_list | grep ^@ > b37.chr22.interval_list
# cat /humgen/gsa-hpprojects/GATK/bundle/2.8/b37/Broad.human.exome.b37.interval_list | grep ^22 >> b37.chr22.interval_list

# for now only I'm handling *one* bam list, no parallelization
java -Xmx3072m \
    -jar /seq/software/picard/1.625/3rd_party/gatk/GenomeAnalysisTK-2.5.jar \
    -T DepthOfCoverage \
    -I $bamlist \
    -L /humgen/atgu1/fs03/eminikel/sandbox/b37.chr22.interval_list \
    -R /humgen/gsa-hpprojects/GATK/bundle/2.8/b37/human_g1k_v37.fasta \
    -dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable \
    --minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 \
    --includeRefNSites \
    --countType COUNT_FRAGMENTS \
    -o $outdir/bamDepthOfCov.DATA
# above: ~1 hr

    
    
# "combine" GATK output into one file
xhmm --mergeGATKdepths -o $outdir/DATA.RD.txt \
    --GATKdepths bamDepthOfCov.DATA.sample_interval_summary
# 2 seconds
# rows are samples, cols are genomic regions

xhmm --matrix -r $outdir/DATA.RD.txt --centerData --centerType target \
    -o $outdir/DATA.filtered_centered.RD.txt \
    --outputExcludedTargets $outdir/DATA.filtered_centered.RD.txt.filtered_targets.txt \
    --outputExcludedSamples $outdir/DATA.filtered_centered.RD.txt.filtered_samples.txt \
    --minTargetSize 10 --maxTargetSize 10000 \
    --minMeanTargetRD 10 --maxMeanTargetRD 500 \
    --minMeanSampleRD 25 --maxMeanSampleRD 200 \
    --maxSdSampleRD 150

# note excluded the GC filtering step
# --excludeTargets $outdir/extreme_gc_targets.txt --excludeTargets $outdir/low_complexity_targets.txt \

xhmm --PCA -r $outdir/DATA.filtered_centered.RD.txt --PCAfiles $outdir/DATA.RD_PCA

xhmm --normalize -r $outdir/DATA.filtered_centered.RD.txt --PCAfiles $outdir/DATA.RD_PCA \
    --normalizeOutput $outdir/DATA.PCA_normalized.txt \
    --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7
    
xhmm --matrix -r $outdir/DATA.PCA_normalized.txt --centerData --centerType sample --zScoreData \
    -o $outdir/DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
    --outputExcludedTargets $outdir/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
    --outputExcludedSamples $outdir/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
    --maxSdTargetRD 30
    
xhmm --matrix -r $outdir/DATA.RD.txt \
    --excludeTargets $outdir/DATA.filtered_centered.RD.txt.filtered_targets.txt \
    --excludeTargets $outdir/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
    --excludeSamples $outdir/DATA.filtered_centered.RD.txt.filtered_samples.txt \
    --excludeSamples $outdir/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
    -o $outdir/DATA.same_filtered.RD.txt

echo "1e-08   6       70      -3      1       0       1       3       1" > $outdir/params.txt

xhmm --discover -p $outdir/params.txt \
    -r $outdir/DATA.PCA_normalized.filtered.sample_zscores.RD.txt -R $outdir/DATA.same_filtered.RD.txt \
    -c $outdir/DATA.xcnv -a $outdir/DATA.aux_xcnv -s $outdir/DATA

xhmm --genotype -p $outdir/params.txt \
    -r $outdir/DATA.PCA_normalized.filtered.sample_zscores.RD.txt -R $outdir/DATA.same_filtered.RD.txt \
    -g $outdir/DATA.xcnv -F /humgen/gsa-hpprojects/GATK/bundle/2.8/b37/human_g1k_v37.fasta \
    -v $outdir/DATA.vcf

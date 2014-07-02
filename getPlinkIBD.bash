#!/bin/bash

# to run:
# bsub -o /humgen/atgu1/fs03/eminikel/sandbox/getPlinkIBD.out -q bweek -P $RANDOM "bash getPlinkIBD.bash -v -o /humgen/atgu1/fs03/eminikel/sandbox/5 -r /seq/dax/macarthur_muscle_disease_ALL/v1/macarthur_muscle_disease_ALL.vcf"

# The below code for parsing command line arguments was lightly modified from: 
# http://stackoverflow.com/a/14203146

# I don't fully understand what this line does.
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
outdir=""
verbose=1
pedmapprefix="" # full path and base file name to which to append .ped or .map prefix
vcf="" # -r for "raw vcf" because -v was taken by "verbose"

# Constants
# list of 5k snps for IBD, in 2-column format: chromosome and position
ibdsnplist=/humgen/atgu1/fs03/DM-Lab/ref/Purcell5k.2col.txt

while getopts "vo:p:r:" opt; do
    case "$opt" in
    v)  verbose=1
        ;;
    o)  outdir=$OPTARG
        ;;
    p)  pedmapprefix=$OPTARG
        ;;
    r)  vcf=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

if test $verbose
  then
    echo "verbose=$verbose, output_directory='$outdir', pedmapprefix='$pedmapprefix', unused arguments: $@"
fi

if test $vcf!=""
  then
    if test $verbose
      then echo "Starting from raw vcf file: "$vcf
    fi
    if test ${vcf: -3} == ".gz"
      then
        vcfflag=--gzvcf
      else
        vcfflag=--vcf
    fi
    # call vcftools to subset the 5k snps and recode to ped/map files
    vcftools $vcfflag $vcf \
             --positions $ibdsnplist \
             --plink --recode --out $outdir/$(basename $vcf)
    # now use the newly created ped/map file
    pedmapprefix=$outdir/$(basename $vcf)
    if test $verbose
      then echo "Subsetted IBD SNPs and wrote to PLINK file: "$pedmapprefix
    fi
  else
    if test $verbose
      then echo "Using user-specified PLINK file: "$pedmapprefix
    fi
fi

# Calculate allele frequencies
# this literally takes 1 second - that doesn't mean it failed.
plink --file $pedmapprefix --freq --out $pedmapprefix

# calculate number of samples. 
# see http://stackoverflow.com/questions/9458752/variable-for-number-of-lines-in-a-file
nlines=`cat $pedmapprefix.ped | wc -l`
nlines=$(($nlines + 1))
nlines=$(($nlines - 1)) # dumb way of casting to integer in bash

# if fewer than 1000 samples, you can do the genome calculation in serial
# in a few minutes and avoid burdening the LSF system.
# this block either calls PLINK once, or submits (nsamples/100 choose 2) jobs
if test $nlines -lt 1000
  then
    if test $verbose
      then
        echo "Fewer than 1000 samples, calling PLINK to calculate IBD in serial."
    fi
    plink --file $pedmapprefix --read-freq $pedmapprefix.frq --all \
          --genome --out $pedmapprefix
    if test $verbose
      then
        echo "PLINK IBD calcuations complete: " $pedmapprefix.genome
    fi
  else
    if test $verbose
      then
        echo "More than 999 samples, submitting tons of jobs to LSF"
        echo "to parallelize PLINK IBD calculations..."
    fi
    
    # split ped file into lists of 100 individuals for parallelization
    cut -d ' ' -f 1,2 $pedmapprefix.ped | split -d -a 3 -l 100 - $outdir/tmp.list
    
    # jobs will be like 1 second each. hour queue is fine.
    let i=0 a=$((nlines/100-1)) # a is number of files created by the above step, minus 1
    let j=0
    while [ $i -le $a ]
    do
        while [ $j -le $a ]
        do
            bsub -q hour -o /dev/null -e /dev/null -J plinkibd -P $RANDOM \
                plink --file $pedmapprefix \
                --read-freq $pedmapprefix.frq\
                --all \
                --genome \
                --genome-lists $outdir/tmp.list`printf "%03i\n" $i` \
                               $outdir/tmp.list`printf "%03i\n" $j` \
                --out $outdir/data.sub.$i.$j
                let j=$j+1
        done
        let i=$i+1
        let j=$i
    done
    
    # now wait for the above jobs to all finish, updating the user every 60 seconds
    while test `bjobs | grep plinkibd | wc -l` -gt 0
    do
      if test $verbose
        then
          echo "Waiting for" `bjobs | grep plinkibd | wc -l` "jobs to finish..."
      fi
      sleep 60
    done
    
    if test $verbose
      then
        echo "PLINK IBD jobs done, now re-combining into one big file..."
    fi
    
    # recombine them
    head -1 $outdir/data.sub.0.0.genome > $outdir/header # use the header only once
    cat $outdir/data.sub.*.genome | grep -v FID1 | cat $outdir/header - > $outdir/$(basename $pedmapprefix).genome
    
    if test $verbose
      then
        echo "PLINK IBD calcuations complete: " $pedmapprefix.genome
    fi
fi



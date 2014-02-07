#!/bin/bash

# runs the entire ExomeDepth CNV-calling pipeline end-to-end
# inputs: list of BAMs (-b), and a VCF (-v) or PLINK IBD .genome file (-g)
# outputs: ExomeDepth CNV calls for all BAMs

# to test run:
# cd /humgen/atgu1/fs03/eminikel/048muscle/data/
# bsub -q priority -W 00:15 -P$RANDOM -J pipetest \
#      -o fulltest/pipeline.out -e fulltest/pipeline.err \
"exomeDepthPipeline.bash \
-o /humgen/atgu1/fs03/eminikel/048muscle/data/fulltest \
-b /humgen/atgu1/fs03/eminikel/048muscle/data/muscle.20.test.list \
-e /humgen/atgu1/fs03/eminikel/048muscle/data/exclude.list \
-g /humgen/atgu1/fs03/eminikel/048muscle/data/macarthur_muscle_disease_ALL.vcf.genome"

# to really run:
# cd /humgen/atgu1/fs03/eminikel/048muscle/data/
bsub -q bweek -W 167:00 -P$RANDOM -J edpipe \
     -o fulltest/pipeline.out -e fulltest/pipeline.err \
 "exomeDepthPipeline.bash \
 -o /humgen/atgu1/fs03/eminikel/048muscle/data/fulltest \
 -b /humgen/atgu1/fs03/eminikel/048muscle/data/MD1.existent.bam.list \
 -e /humgen/atgu1/fs03/eminikel/048muscle/data/exclude.list \
 -r /humgen/atgu1/fs03/eminikel/048muscle/data/macarthur_muscle_disease_ALL.vcf"

# choose one of these two options:
# -r /humgen/atgu1/fs03/eminikel/048muscle/data/macarthur_muscle_disease_ALL.vcf
# -g /humgen/atgu1/fs03/eminikel/048muscle/data/macarthur_muscle_disease_ALL.vcf.genome


# The below code for parsing command line arguments was lightly modified from: 
# http://stackoverflow.com/a/14203146

# I don't fully understand what this line does.
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
outdir=""
verbose=1
vcf="" # -r for "raw vcf" because -v was taken by "verbose"
bamlist="" 
genome=""
excludelist=""

# Constants
# list of 5k snps for IBD, in 2-column format: chromosome and position
ibdsnplist=/humgen/atgu1/fs03/DM-Lab/ref/Purcell5k.2col.txt

while getopts "vo:b:r:g:e:" opt; do
    case "$opt" in
    v)  verbose=1
        ;;
    o)  outdir=$OPTARG
        ;;
    b)  bamlist=$OPTARG
        ;;
    r)  vcf=$OPTARG
        ;;
    g)  genome=$OPTARG
        ;;
    e)  excludelist=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

if test $verbose
  then
    echo "verbose=$verbose, 
          output_directory='$outdir', 
          vcf='$vcf', 
          genome='$genome', 
          bamlist='$bamlist', 
          excludelist='$excludelist', 
          unused arguments: $@"
fi

# temporary settings for testing
outdir=/humgen/atgu1/fs03/eminikel/048muscle/data/fulltest
bamlist=/humgen/atgu1/fs03/eminikel/048muscle/data/muscle.bam.existent.list
excludelist=/humgen/atgu1/fs03/eminikel/048muscle/data/exclude.list
vcf=/humgen/atgu1/fs03/eminikel/048muscle/data/macarthur_muscle_disease_ALL.vcf
genome=/humgen/atgu1/fs03/eminikel/048muscle/data/macarthur_muscle_disease_ALL.vcf.genome


# Step 1: counts
# Parallelize the calculation of ExomeDepth counts, then merge back together
if test bamlist == ""
	then
      echo "You didn't specify a list of BAMs using -b. Exiting w/ code 1"
      exit 1
    else
      if test $verbose
        then
          echo "Using BAM list: $bamlist"
      fi
      # parallelize by submitting a jobarray with one job per single bam
      nbam=`wc -l $bamlist | cut -d ' ' -f1`
      # make sure $nbam is an integer
      nbam=$(($nbam + 1))
      nbam=$(($nbam - 1)) # yes, this is stupid, but no, $nbam + 0 *doesn't* cast it to int.

      mkdir $outdir/clutter
      bsub -o $outdir/clutter/ExomeDepthCountJobOutput_%I.out \
           -e $outdir/clutter/ExomeDepthCountJobOutput_%I.err \
           -q bweek -P $RANDOM -J exdep[1-$nbam] -W 10:00 \
           "getExomeDepthCounts.r -s \`sed -n \${LSB_JOBINDEX}p $bamlist\` -o $outdir"
fi

# bash one-liner to test the above method
# bsub -o %I.out -e %I.err -q bhour -P $RANDOM -J jtest[1-5] -W 00:01 "echo \$LSB_JOBINDEX > \$LSB_JOBINDEX.internal.out"



# Step 2: relatedness
# Check if we have a .genome file. If not, create one.
if test genome == ""
	then
	  if test $verbose
        then
          echo "Creating a .genome file..."
      fi
      # since the getPlinkIBD.bash step is already as parallel as it's going to get,
      # just run it directly from this script rather than submitting a job.
	  bash getPlinkIBD.bash -v -o $outdir -r $vcf
	  # now grab the genome file created by that script
      genome=$outdir/$(basename $vcf)
      if test $verbose
        then
          echo "New .genome file is: $genome"
      fi
    else
      if test $verbose
        then
          echo "Using existing .genome file: $genome"
      fi
fi


# Interlude: loop and wait to make sure the jobs from Step 1 are done.
while test $((`bjobs | grep exdep | wc -l`)) -gt 0
do
    echo -n "Waiting for" $((`bjobs | grep exdep | wc -l`)) 
    echo " ExomeDepth count jobs to finish"
	sleep 60
done

# Step 3:
# Fuse the counts back together
cd $outdir
# grab the first 5 columns which give data about exons
cut -f1-5 `ls -1 *.counts.txt | head -1` > temp.exonmeta
# now grab the 6th column of every counts file
for countsfile in *.counts.txt
do
	cut -f6 $countsfile > temp.$countsfile
done
paste temp.exonmeta `ls -1 | \grep temp | \grep counts.txt` > merged.counts

# clean up a little bit of the clutter
rm temp.exonmeta
rm temp.*.counts.txt
cd -

# now set the countsfile variable to what we used above.
countsfile=$outdir/merged.counts

# create lists of samples to call CNVs on in next step
# for now do it so each file contains 1 name.
for i in $(seq 1 $nbam)
do
    basename `sed -n ${i}p $bamlist`  | awk '{print "echo \""$1"\" > "$1".blist"}' | bash
done

cat *.blist > all.samples.blist

# Step 4: Call CNVs using the counts and .genome file
# parallelize this by each individual sample
bsub -o $outdir/ExomeDepthCNVJobOutput_%I.out \
     -e $outdir/ExomeDepthCNVJobOutput_%I.err \
     -q bhour -P $RANDOM -J exdep[1-$nbam] -W 00:45 \
     "getExomeDepthCNVs.r -c $countsfile -s \`sed -n \${LSB_JOBINDEX}p all.samples.blist\`.blist -g $genome -t .15 -e $excludelist -o $outdir"

# Step 5. Now create a read ratio matrix and combined matrix from counts
# this could go before Step 4 but since Step 4 is slow, might as well get those jobs
# submitted and then do step 5 while waiting
bsub -o $outdir/exomedepth/clutter/2/RatioMatrixJobOutput.out \
     -e $outdir/exomedepth/clutter/2/RatioMatrixJobOutput.err \
     -q bweek -P $RANDOM -J ratmat -W 05:00 \
     "ratioMatrix.r -c $countsfile -r -o $outdir"

# Outro: loop and wait to make sure the jobs from Step 4 are done.
while test $((`bjobs | grep exdep | wc -l`)) -gt 0
do
    echo -n "Waiting for" $((`bjobs | grep exdep | wc -l`)) 
    echo " ExomeDepth CNV jobs to finish"
  sleep 60
done

# Outro: loop and wait to make sure the jobs from Step 5 are also done.
while test $((`bjobs | grep exdep | wc -l`)) -gt 0
do
    echo -n "Waiting for" $((`bjobs | grep ratmat | wc -l`)) 
    echo " ExomeDepth read ratio matrix jobs to finish"
  sleep 60
done



echo "Wow, if you got here, then amazingly this whole pipeline actually worked."
echo "Now exiting with code 0"

exit 0


# below is stuff that I'm using to run some steps manually during testing.

countsfile=/humgen/atgu1/fs03/eminikel/048muscle/data/merged.counts 
genome=/humgen/atgu1/fs03/eminikel/048muscle/data/macarthur_muscle_disease_ALL.vcf.genome
excludelist=/humgen/atgu1/fs03/eminikel/048muscle/data/exclude.list
outdir=/humgen/atgu1/fs03/eminikel/048muscle/data
getExomeDepthCNVs.r -c $countsfile \
                    -g $genome \
                    -t .15 \
                    -o $outdir \
                    -e $excludelist \
                    -s 66T_NG_1.bam.list

# submitted all at 3:32p on 2/5
bsub -o $outdir/clutter/ExomeDepthCNVJobOutput_%I.out \
     -e $outdir/clutter/ExomeDepthCNVJobOutput_%I.err \
     -q bhour -P $RANDOM -J exdep[1-$nbam] -W 00:45 \
     "getExomeDepthCNVs.r -c $countsfile -s \`sed -n \${LSB_JOBINDEX}p all.samples.blist\`.blist -g $genome -t .15 -e $excludelist -o $outdir"

# submitted all at 1:57p on Thurs. 1-2 ran successfully in 2-3 mins on priority queue
bsub -o $outdir/exomedepth/clutter/2/ExomeDepthCNVJobOutput_%I.out \
     -e $outdir/exomedepth/clutter/2/ExomeDepthCNVJobOutput_%I.err \
     -q bhour -P $RANDOM -J exdep[1-$nbam] -W 00:45 \
     "getExomeDepthCNVs.r -c $countsfile -s \`sed -n \${LSB_JOBINDEX}p all.samples.blist\`.blist -g $genome -t .15 -e $excludelist -o $outdir"
# success.


# now run the ratiomatrix script
bsub -o $outdir/exomedepth/clutter/2/RatioMatrixJobOutput.out \
     -e $outdir/exomedepth/clutter/2/RatioMatrixJobOutput.err \
     -q bhour -P $RANDOM -J ratmat -W 00:50 \
     "ratioMatrix.r -c $countsfile -r $outdir -o $outdir"

# temporary script to convert .txt into .csv for Brett
cd humgen/atgu1/fs03/eminikel/048muscle/data/exomedepth
mkdir csv
for txtfile in *.bam.txt
do
  cat $txtfile | sed 's/\t/,/g' | sed 's/^start.p/,start.p/' > csv/${txtfile}.csv
done

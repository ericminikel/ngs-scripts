#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_3.0.2/bin/Rscript

# Eric Minikel
# script to run ExomeDepth to get CNVs (Step 1)
# how to run with no special options:
# getExomeDepthCNVs.r -c /humgen/atgu1/fs03/eminikel/sandbox/counts2.rda -o /humgen/atgu1/fs03/eminikel/sandbox/3/ 

# A note on design:
# We discussed having one list of samples to call CNVs on and another list of 
# samples to use as optional references. For now I am assuming only one input 
# counts.rda file because ExomeDepth does not have any native capability to merge
# count files, though we could certainly hack one if needed.  For now I figure 
# we'll generate counts on all samples in a dataset, and then only specify the
# list of samples we want CNVs for (with -s) and then the PED (-p) and/or
# PLINK IBD (-m) files will be used to exclude relateds (if provided), but other
# than that any and all samples in the counts.rda file will be fair game for
# using as references.

start_time = Sys.time()

suppressPackageStartupMessages(require(optparse)) # http://cran.r-project.org/web/packages/optparse/optparse.pdf
options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings


options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
  make_option(c("-c", "--countsfile"), action="store", default='', 
              type='character', help="Full path to counts file"),
  make_option(c("-o", "--outdir"), action="store", default='./',
              type='character', help="Output directory for CNV calls [default %default]"),
  make_option(c("-p", "--pedfile"), action="store", default='',
              type='character', help="PED file [default none]"),
  make_option(c("-m", "--relmat"), action="store", default='',
              type='character', help="Relatedness matrix in PLINK IBD format [default none]"),
  make_option(c("-s", "--callsamples"), action="store", default='',
              type='character', help="Path to list of samples to call CNVs on [default all]"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print verbose output [default %default]"),
  make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
              help="Do not print verbose output (this is the default)")
)
opt = parse_args(OptionParser(option_list=option_list))

if (opt$verbose) {
    cat(paste("Execution started at ",start_time,"\n",sep=""),file=stdout())
}

if (opt$verbose) {
    cat("Running getExomeDepthCNVs.r with the following parameters:\n",file=stdout())
    print(opt,file=stdout())
}

if (opt$verbose) {
    cat("Loading ExomeDepth package...\n",file=stdout())
}

suppressPackageStartupMessages(require(ExomeDepth)) # http://cran.r-project.org/web/packages/ExomeDepth/ExomeDepth.pdf
data(exons.hg19) # from ExomeDepth package

if (file.exists(opt$countsfile)) {
    if(opt$verbose) {
        cat(paste("Loading count data from ",opt$countsfile," ...\n",sep=""),file=stdout())
    }
    load(opt$countsfile) # this contains an S4 object called "counts"
    if(opt$verbose) {
        cat(paste("Loaded.\n",sep=""),file=stdout())
    }
    stopifnot(exists("counts")) # assert that a counts object was loaded & now exists
    if(opt$verbose) {
        cat("Successfully loaded count data.\n",file=stdout())
    }
} else {
    cat("**Count file not found. You specified this path:\n",file=stderr())
    cat(paste(opt$countsfile,"\n",sep=""),file=stderr())
}

if (file.exists(opt$outdir)) {
    setwd(opt$outdir)
    if(opt$verbose) {
        cat(paste("Will write CNV calls to: ",opt$outdir,"\n",sep=""),file=stdout())
    }
} else {
    cat(paste("You specified this directory: ",opt$outdir,"\n",sep=""),file=stdout())
    cat(paste("And it doesn't exist.\n",sep=""),file=stdout())
    stop()
}
# counts is an S4 object.
# you need to cast it to a data frame (for bin length)
# AND to a matrix (for reference.count)
# and for some reason you can't cast S4 directly to matrix, only via df
if(opt$verbose) {
    cat("Casting count data to a data frame...\n",file=stdout())
}
countdf = as.data.frame(counts)
if(opt$verbose) {
    cat("Column names in the count data frame:\n",file=stdout())
    print(colnames(countdf),file=stdout())
}


if(opt$verbose) {
    cat("Casting count data to a matrix...\n",file=stdout())
}
countmat = as.matrix(countdf[,6:dim(countdf)[2]]) # remove cols 1-5 metadata

if(opt$callsamples=='') {
    # default mode: all samples (colnames of count matrix)
    # will have CNVs called on them
    call_indices = 1:dim(countmat)[2] 
    
    if(opt$verbose) {
        cat("Default mode - calling CNVs on all samples:\n",file=stdout())
        print(colnames(countmat)[call_indices],file=stdout())
    }
    
} else {
    # user has specified a file listing the samples on which to call CNVs
    if(file.exists(opt$callsamples)) {
        cat("Reading list of samples to call CNVs on from:\n",file=stdout())
        cat(paste(opt$callsamples,"\n",sep=""),file=stdout())
        call_sample_names = read.table(opt$callsamples,header=FALSE)$V1
        cat("Sample list read successfully.\n",file=stdout())
        # figure out if they specified samples that aren't in the count data
        missing_samples = call_sample_names[!call_sample_names %in% colnames(countmat)]
        if(length(missing_samples) > 0) {
            cat("You specified some samples that aren't in the counts file:\n",file=stdout())
            print(missing_samples,file=stdout())
        }
        # get indices of ones they specified that are actually in the counts file
        call_indices = match(call_sample_names,call_colnames(countmat))
        call_indices = call_indices[!is.na(call_indices)] # remove NAs
        if(length(call_indices) > 0) {
            cat("Will be calling CNVs on the following samples:\n",file=stdout())
            print(colnames(countmat)[call_indices],file=stdout())
        } else {
            cat("**Didn't find any specified sample names in the counts file.\n",file=stdout())
        }
    }
}

# At this point we can assume that call_indices is a vector of length > 0
# which indicates the columns of countmat on which to call CNVs
# let's make some assertions about that:
stopifnot(exists("call_indices"))
stopifnot(length(call_indices) > 0)

# beta version: assume you want CNVs on all samples
for (i in call_indices) {
    sample_name = colnames(countmat)[i]
    if(opt$verbose) {
        cat(paste("Now about to call CNVs on ",sample_name,"...",sep=""),file=stdout())
    }
    reference_list = select.reference.set(test.counts = countmat[,i], 
        reference.count = countmat[,-i],
        bin.length=(countdf$end-countdf$start)/1000,
        n.bins.reduced = 10000)
    reference_set = apply(
        X = as.matrix(countdf[, reference_list$reference.choice]), 
        MAR=1, FUN=sum)
    all_exons = new('ExomeDepth', test=countmat[,i], 
        reference=reference_set,
        formula = 'cbind(test,reference) ~ 1')
    all_exons = CallCNVs(x = all_exons, transition.probability=10^-4,
        chromosome=countdf$space, start=countdf$start,
        end=countdf$end, name=countdf$names)
    write.table(all_exons@CNV.calls, file=paste(sample_name,".txt",sep=''), 
        sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
    if (opt$verbose) {
        cat(paste("..done.\n",sep=''),file=stdout())
    }
}


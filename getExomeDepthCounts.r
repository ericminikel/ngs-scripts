#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_3.0.2/bin/Rscript

# Eric Minikel
# script to run ExomeDepth to get counts (Step 1)
# how to run:
# getExomeDepthCounts.r -b bamlist.txt -o /output/path/and/filename.rda -v > runExomeDepthOutput.txt

start_time = Sys.time()

suppressPackageStartupMessages(require(optparse)) # http://cran.r-project.org/web/packages/optparse/optparse.pdf
options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
  make_option(c("-b", "--bamlist"), action="store", default='', 
              type='character', help="Path to list of BAMs"),
  make_option(c("-o", "--outpath"), action="store", default='./',
              type='character', help="Output path to store counts [default %default]"),
  make_option(c("-r", "--refexons"), action="store", default='hg19',
              type='character', help="Reference exome to use. Must be hg19. [default %default]"),
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
    cat("Running getExomeDepthCounts.r with the following parameters:\n",file=stdout())
    print(opt,file=stdout())
}


if (opt$verbose) {
    cat("Loading ExomeDepth package...\n",file=stdout())
}

suppressPackageStartupMessages(require(ExomeDepth)) # http://cran.r-project.org/web/packages/ExomeDepth/ExomeDepth.pdf
data(exons.hg19) # from ExomeDepth package

if (opt$refexons != 'hg19') {
    cat("**Stopping execution because we don't support any reference except hg19.\n",file=stderr())
    stop()
}

if (opt$verbose) {
    cat("Reading list of BAMs...\n",file=stdout())
}


# read list of BAMs
# to avoid writing a tryCatch I use file.exists, and set default to '', which is
# a file that never exists.
if (file.exists(opt$bamlist)) {
    # read bam list directly into a vector (note use of $V1)
    bams = read.table(opt$bamlist,header=FALSE)$V1
} else {
    cat("**You need to specify a valid BAM list using -b.\n",file=stderr())
    cat(paste("The filename you specified was '",opt$bamlist,"'.",sep=''),file=stderr())
    stop()
}

if (opt$verbose) {
    cat(paste("Successfully read the BAM list from ",opt$bamlist,"\n",sep=''),file=stdout())
    cat("Using the following BAMs: \n",file=stdout())
    write.table(bams,row.names=FALSE,quote=FALSE,col.names=FALSE,file=stdout())
}


# read output directory
# note right now if not existent, I stop execution. 
# an alternative is to create the dir.
# see http://stackoverflow.com/questions/4216753/check-existence-of-directory-and-create-if-doesnt-exist
outdir = dirname(opt$outpath)
filename = basename(opt$outpath)
if (file.exists(outdir)) {
    setwd(outdir)
} else {
    cat("**You need to specify a valid output directory using -o.\n",file=stderr())
    cat(paste("The directory you specified was '",opt$outdir,"'.",sep=''),file=stderr())
    stop()
}

if (opt$verbose) {
    cat("Calling ExomeDepth to compute counts for the BAMs...",file=stdout())
}

counts = getBamCounts(bed.frame = exons.hg19, bam.files = bams)

if (opt$verbose) {
    cat("Successfully calculated counts.\n",file=stdout())
    cat("Saving the counts to: \n",file=stdout())
    cat(paste(opt$outpath,"\n",sep=""),file=stdout())
}

save(counts,file=filename)

if (opt$verbose) {
    cat(paste("Successfully wrote counts to: ",opt$outpath,"\n",sep=""),file=stdout())
}

duration = format(Sys.time() - start_time)

if(opt$verbose) {
    cat(paste("Completed execution in ",duration,".\n",sep=''),file=stdout())
}


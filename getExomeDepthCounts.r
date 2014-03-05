#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_3.0.2/bin/Rscript


# Eric Minikel
# script to run ExomeDepth to get counts (Step 1)
# how to run:
# getExomeDepthCounts.r -b bamlist.txt -o /output/path/without/final/slash -v > runExomeDepthOutput.txt

# A note on design:
# ExomeDepth does not provide native capability to merge counts from multiple files
# We could hack a solution if need be, but for now I'm just doing all counts
# for a project in serial in this script.
# This script takes 5.77 minutes to do chr22 on 20 BAMs, so I predict it will take
# 5.77 * (1/.02) * (200/20) / (60*24) = ~2 days to do whole exome on 200 BAMs

start_time = Sys.time()

suppressPackageStartupMessages(require(optparse)) # http://cran.r-project.org/web/packages/optparse/optparse.pdf
options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
  make_option(c("-b", "--bamlist"), action="store", default='', 
              type='character', help="Path to list of BAMs"),
  make_option(c("-s", "--singlebam"), action="store", default='',
              type='character', help="Path to a single BAM"),
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

# try to load the ExomeDepth package in a safe and informative way
tryCatch ({
    # try
    suppressPackageStartupMessages(require(ExomeDepth)) # http://cran.r-project.org/web/packages/ExomeDepth/ExomeDepth.pdf
    data(exons.hg19) # from ExomeDepth package
    if (opt$verbose) {
        cat("Loaded successfully.\n",file=stdout())
    }
}, error = function(cond) {
    cat("Experienced an error trying to load the ExomeDepth package.\n",file=stderr())
    print(cond, file=stderr())
    stop()
})
##### can't include this warning condition or else the package load fails on warnings:
# , warning = function(cond) {
#     cat("Got a warning which might or might not mean the ExomeDepth package isn't installed.\n\n",file=stderr())
#     cat("If so, try install.packages() or check whether you're using\n",file=stderr())
#     cat("the version of R for which you already installed it.\n\n",file=stderr())
#     print(cond, file=stderr())
# })


if (opt$refexons != 'hg19') {
    cat("**Stopping execution because we don't support any reference except hg19.\n",file=stderr())
    stop()
}

# read list of BAMs
# to avoid writing a tryCatch I use file.exists, and set default to '', which is
# a file that never exists.
if (file.exists(opt$bamlist)) {
    if (opt$verbose) {
        cat(paste("Reading list of BAMs from ",opt$bamlist,"\n"),file=stdout())
    }
    # read bam list directly into a vector (note use of $V1)
    bams = read.table(opt$bamlist,header=FALSE)$V1
    if (opt$verbose) {
        cat(paste("Successfully read the BAM list from ",opt$bamlist,"\n",sep=''),file=stdout())
        cat("Using the following BAMs: \n",file=stdout())
        write.table(bams,row.names=FALSE,quote=FALSE,col.names=FALSE,file=stdout())
    }

} else if (file.exists(opt$singlebam)) {
    bams = c(opt$singlebam)
    if (opt$verbose) {
        cat(paste("Will read single BAM from: ","\n",sep=''),file=stdout())
        write.table(bams,row.names=FALSE,quote=FALSE,col.names=FALSE,file=stdout())
    }
} else {
    cat("**You need to specify a valid BAM list using -b or single BAM using -s.\n",file=stderr())
    cat(paste("The BAM list you specified was '",opt$bamlist,"'.",sep=''),file=stderr())
    cat(paste("The single BAM you specified was '",opt$singlebam,"'.",sep=''),file=stderr())
    stop()
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

# ExomeDepth sometimes spectacularly fails to find the BAM index, for no clear reason
# so help it out by testing a few likely filenames
bamindexes = bams
for (i in 1:length(bamindexes)) {
  stardotbai = gsub(".bam",".bai",bams[i])
  stardotbamdotbai = gsub(".bam",".bam.bai",bams[i])
  if (file.exists(stardotbai)) {
    bamindexes[i] = stardotbai
  } else if (file.exists(stardotbamdotbai)) {
    bamindexes[i] = stardotbamdotbai
  }
  else {
    cat(paste("Cannot find a .bai index for BAM: ",bams[i],"\n",sep=""),file=stderr())
    cat("stopping execution....",file=stderr())
    stop()
  }
}


counts = getBamCounts(bed.frame = exons.hg19, bam.files = bams)

if (opt$verbose) {
    cat("Successfully calculated counts.\n",file=stdout())
    cat("Saving the counts to: \n",file=stdout())
    cat(paste(opt$outpath,"\n",sep=""),file=stdout())
}

# old method: save counts as IRanges object in an RData file
# save(counts,file=filename)

# new method: save counts as a text file
# If operating in multi-BAM mode, save all to counts.txt
if (file.exists(opt$bamlist)) {
    countspath = paste(opt$outpath,"/counts.txt",sep="")
} else  {
    countspath = paste(opt$outpath,"/",gsub(".bam",".counts.txt",basename(opt$singlebam)),sep="")
}
write.table(counts,countspath,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

if (opt$verbose) {
    cat(paste("Successfully wrote counts to: ",countspath,"\n",sep=""),file=stdout())
}

duration = format(Sys.time() - start_time)

if(opt$verbose) {
    cat(paste("Completed execution in ",duration,".\n",sep=''),file=stdout())
}


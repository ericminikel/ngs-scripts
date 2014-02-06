#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_3.0.2/bin/Rscript

# converts an ExomeDepth counts matrix to a ratio matrix


start_time = Sys.time()

suppressPackageStartupMessages(require(optparse)) # http://cran.r-project.org/web/packages/optparse/optparse.pdf
options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings


option_list = list(
  make_option(c("-c", "--countsfile"), action="store", default='', 
              type='character', help="Full path to counts matrix"),
  make_option(c("-r", "--reflistdir"), action="store", default='./',
              type='character', help="Directory full of .reflist.txt files [default %default]"),
  make_option(c("-o", "--outdir"), action="store", default='./',
              type='character', help="Output directory [default %default]"),
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
    cat("Running ratioMatrix.r with the following parameters:\n",file=stdout())
    print(opt,file=stdout())
}

# open the counts table or stop if it can't be found
if(file.exists(opt$countsfile)) {
    countdf = read.table(opt$countsfile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
} else {
    cat("Couldn't find the counts file you specified with -c: \n",file=stderr())
    cat(paste(opt$countsfile,"\n"),file=stderr())
    stop()
}

# check if the reference directory exists
if(file.exists(opt$reflistdir)) {
    if(opt$verbose) {
        cat("Will be pulling reference lists from this directory: \n",file=stdout())
        cat(paste(opt$reflistdir,"\n"),file=stdout())
    }
} else {
    cat("Couldn't find the reference list directory you specified with -r: \n",file=stdout())
    cat(paste(opt$reflistdir,"\n"),file=stdout())
    stop()
}

# check output directory exists
if (file.exists(opt$outdir)) {
    setwd(opt$outdir)
    if(opt$verbose) {
        cat(paste("Will write data to: ",opt$outdir,"\n",sep=""),file=stdout())
    }
} else {
    cat(paste("You specified this directory: ",opt$outdir,"\n",sep=""),file=stdout())
    cat(paste("And it doesn't exist.\n",sep=""),file=stdout())
    stop()
}


# cast counts to matrix
if(opt$verbose) {
    cat("Casting count data to a matrix...\n",file=stdout())
}
countmat = as.matrix(countdf[,6:dim(countdf)[2]]) # remove cols 1-5 metadata

# create a ratio matrix of the same dimensions & names as the count matrix
ratiomat = matrix(nrow=dim(countmat)[1],ncol=dim(countmat)[2])
colnames(ratiomat) = colnames(countmat)
rownames(ratiomat) = rownames(countmat)

# create matrices for observed and expected as well
observedmat = countmat
expectedmat = ratiomat # just to initialize with same names and dimensions

# now loop through and populate the ratio matrix
for (i in 1:dim(countmat)[2]) {
    sample_name = colnames(countmat)[i]
    refpath = paste(opt$reflistdir,"/",sample_name,".reflist.txt",sep="")
    reflist = read.table(refpath,skip=2,sep="\t",header=TRUE)
    
    # check that the samples are in the same order
    stopifnot(all(colnames(reflist) == colnames(countmat)))

    # grab a boolean vector for indexing columns - "was each sample used as reference?"
    ref_used = reflist[1,] == 1

    expectedmat[,i] = rowMeans(countmat[,ref_used]) # calculate expected based on ref samples
    ratiomat[,i] = countmat[,i] / expectedmat[,i] # ratio obs/exp
}

# put back into a dataframe and write out to disk
ratiodf = cbind(countdf[,1:5],as.data.frame(ratiomat))
write.table(ratiodf,"readratio.txt",row.names=FALSE,col.names=TRUE,sep="\t")

# also try a semicolon-delimited VCF-like version
# because OBS/EXP can be more useful than NaN when expectation is low
combinedmat = matrix(paste(observedmat,expectedmat,ratiomat,sep=";"),
    nrow=dim(ratiomat)[1],ncol=dim(ratiomat)[2])
combodf = cbind(countdf[,1:5],as.data.frame(combinedmat))
write.table(combodf,"obs-exp-ratio.txt",row.names=FALSE,col.names=TRUE,sep="\t")



duration = format(Sys.time() - start_time)

if(opt$verbose) {
    cat(paste("Completed execution in ",duration,".\n",sep=''),file=stdout())
}



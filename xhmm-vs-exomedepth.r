#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_3.0.2/bin/Rscript

# compares the output of an XHMM file and an ExomeDepth file

start_time = Sys.time()

suppressPackageStartupMessages(require(optparse)) # http://cran.r-project.org/web/packages/optparse/optparse.pdf
options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
  make_option(c("-x", "--xhmmfile"), action="store", default='', 
              type='character', help="Full path to XHMM output"),
  make_option(c("-e", "--edfile"), action="store", default='',
              type='character', help="Full path to ExomeDepth output"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print verbose output [default %default]"),
  make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
              help="Do not print verbose output")
)
opt = parse_args(OptionParser(option_list=option_list))

if (opt$verbose) {
    cat(paste("Execution started at ",start_time,"\n",sep=""),file=stdout())
}

if (opt$verbose) {
    cat("Loading Bioconductor packages...\n",file=stdout())
}
suppressPackageStartupMessages(source("http://bioconductor.org/biocLite.R"))
suppressPackageStartupMessages(biocLite("IRanges"))
suppressPackageStartupMessages(library(IRanges))


# manually set opts for testing
opt=data.frame(xhmmfile="muscle.xhmm.calls.txt",edfile="muscle.ed.calls.csv",verbose=TRUE)

xhmm = read.table(opt$xhmmfile,header=TRUE,sep="\t")
ed = read.table(opt$edfile,header=TRUE,sep=',')

xhmm$samplename = xhmm$SAMPLE

# extract chr, start, stop from xhmm
extracted_strings = strsplit(xhmm$INTERVAL,split="[:-]")
xhmm$chromosome = sapply(extracted_strings,"[",1)
xhmm$start      = as.integer(sapply(extracted_strings,"[",2))
xhmm$end        = as.integer(sapply(extracted_strings,"[",3))

# make unique ids to track many:manyness
xhmm$xuid = paste(xhmm$samplename,"_",xhmm$INTERVAL,sep="")
ed$euid = paste(ed$samplename,"_",ed$id,sep="")

intersection = sqldf("
select   *
from     xhmm, ed
where    xhmm.samplename = ed.samplename
and      xhmm.chromosome = ed.chromosome
and      (xhmm.start <= ed.end and xhmm.end >= ed.start)
;")

xhmmonly = sqldf("
select   *
from     xhmm x
where    not exists
    (select   null
     from     ed e
     where    e.samplename = x.samplename
     and      e.chromosome = x.chromosome
     and      (x.start <= e.end and x.end >= e.start))
;")

edonly = sqldf("
select   *
from     ed e
where    not exists
    (select   null
     from     xhmm x
     where    e.samplename = x.samplename
     and      e.chromosome = x.chromosome
     and      (x.start <= e.end and x.end >= e.start))
;")

# write out dimensions
cat(paste("ed: ",dim(ed)[1],"\n",sep=""),file=stdout())
cat(paste("xhmm: ",dim(xhmm)[1],"\n",sep=""),file=stdout())
cat(paste("intersection: ",dim(intersection)[1],"\n",sep=""),file=stdout())
cat(paste("    xhmm unique: ",length(unique(intersection$xuid)),"\n",sep=""),file=stdout())
cat(paste("    ed unique: ",length(unique(intersection$euid)),"\n",sep=""),file=stdout())
cat(paste("xhmm only: ",dim(xhmmonly)[1],"\n",sep=""),file=stdout())
cat(paste("ed only: ",dim(edonly)[1],"\n",sep=""),file=stdout())

# are the quality of CNVs higher in the intersection than when unique to either side?

summary(ed$BF)
summary(intersection$BF)
summary(edonly$BF)
ks.test(edonly$BF, intersection$BF)

summary(xhmm$Q_SOME)
summary(intersection$Q_SOME)
summary(xhmmonly$Q_SOME)
ks.test(xhmmonly$Q_SOME, intersection$Q_SOME)

# above analysis suggests yes, but there are some really high quality variants
# that are unique to either call set.
# can we explore those?




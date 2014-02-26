#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_3.0.2/bin/Rscript

# Eric Minikel
# script to run ExomeDepth to get CNVs (Step 1)
# how to run with no special options:
# getExomeDepthCNVs.r -c /humgen/atgu1/fs03/eminikel/sandbox/counts2.rda -o /humgen/atgu1/fs03/eminikel/sandbox/3/ -e /humgen/atgu1/fs03/eminikel/sandbox/exclude.list -g /humgen/atgu1/fs03/eminikel/sandbox/sandbox.genome

# A note on design:
# We discussed having one list of samples to call CNVs on and another list of 
# samples to use as optional references. For now I am assuming only one input 
# counts.rda file because ExomeDepth does not have any native capability to merge
# count files, though we could certainly hack one if needed.  For now I figure 
# we'll generate counts on all samples in a dataset, and then only specify the
# list of samples we want CNVs for (with -s) and then PLINK .genome IBD files
# (with -g) will be used to exclude relateds (if provided), but other
# than that any and all samples in the counts.rda file will be fair game for
# using as references.

start_time = Sys.time()

suppressPackageStartupMessages(require(optparse)) # http://cran.r-project.org/web/packages/optparse/optparse.pdf
options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings


option_list = list(
  make_option(c("-c", "--countsfile"), action="store", default='', 
              type='character', help="Full path to counts file"),
  make_option(c("-o", "--outdir"), action="store", default='./',
              type='character', help="Output directory for CNV calls [default %default]"),
  make_option(c("-g", "--genome"), action="store", default='', # recommend using long name --genome; otherwise R thinks -g is specifying the gui.
              type='character', help="PLINK .genome IBD file"),
  make_option(c("-e", "--exclude"), action="store", default='',
              type='character', help="txt file listing samples to exclude from reference"),
  make_option(c("-t", "--threshold"), action="store", default=.03125,
              type='numeric', help="Maximum IBD to still use as reference [default %default]"),
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


# if (file.exists(opt$countsfile)) {
#     if(opt$verbose) {
#         cat(paste("Loading count data from ",opt$countsfile," ...\n",sep=""),file=stdout())
#     }
#     load(opt$countsfile) # this contains an S4 object called "counts"
#     if(opt$verbose) {
#         cat(paste("Loaded.\n",sep=""),file=stdout())
#     }
#     stopifnot(exists("counts")) # assert that a counts object was loaded & now exists
#     if(opt$verbose) {
#         cat("Successfully loaded count data.\n",file=stdout())
#     }
# } else {
#     cat("**Count file not found. You specified this path:\n",file=stderr())
#     cat(paste(opt$countsfile,"\n",sep=""),file=stderr())
# }

# counts is an S4 object.
# you need to cast it to a data frame (for bin length)
# AND to a matrix (for reference.count)
# and for some reason you can't cast S4 directly to matrix, only via df
if(opt$verbose) {
    cat("Casting count data to a data frame...\n",file=stdout())
}
countdf = read.table(opt$countsfile,sep="\t",header=TRUE) # as.data.frame(counts)
if(opt$verbose) {
    cat("Column names in the count data frame:\n",file=stdout())
    print(colnames(countdf),file=stdout())
}


if(opt$verbose) {
    cat("Casting count data to a matrix...\n",file=stdout())
}
countmat = as.matrix(countdf[,6:dim(countdf)[2]]) # remove cols 1-5 metadata

if(opt$verbose) {
    cat("OK, here's the head of the matrix...\n",file=stdout())
    #print(head(countmat),file=stdout())
    # it was too big and made the output hard to read
}

# prependX = function(bamnames) {
    # modified_names = bamnames
    # starts_with_number = substr(bamnames,1,1) %in% 0:9
    # modified_names[starts_with_number] = paste("X",bamnames[starts_with_number],sep="")
    # return ( modified_names )
# }


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
        # user can specify full path of BAMs OR just BAM filenames.
        # either way we just use the BAM filenames, because that's what the
        # columns are named in the counts matrix.
        call_sample_names = basename(read.table(opt$callsamples,header=FALSE)$V1)
        # note that when you read the BAMs originally in getExomeDepthCounts.r,
        # R will have modified the BAM names to prepend X before numbers and
        # remove characters that are illegal in R column names. Luckily R
        # exposes this same logic in the make.names function, so if we call this
        # on the user's input as to which BAMs to use, then it should exactly
        # match what was done in the getExomeDepthCounts.r script.
        call_sample_names = make.names(call_sample_names)
        cat("Sample list read successfully.\n",file=stdout())
        # figure out if they specified samples that aren't in the count data
        missing_samples = call_sample_names[!call_sample_names %in% colnames(countmat)]
        if(length(missing_samples) > 0) {
            cat("You specified some samples that aren't in the counts file:\n",file=stdout())
            print(missing_samples,file=stdout())
        }
        # get indices of ones they specified that are actually in the counts file
        call_indices = match(call_sample_names,colnames(countmat))
        call_indices = call_indices[!is.na(call_indices)] # remove NAs
        if(length(call_indices) > 0) {
            cat("Will be calling CNVs on the following samples:\n",file=stdout())
            print(colnames(countmat)[call_indices],file=stdout())
        } else {
            cat("**Didn't find any specified sample names in the counts file.\n",file=stdout())
        }
    } else {
        cat("Couldn't find the callsamples list you specified with -s: \n",file=stderr())
        cat(paste(opt$callsamples),"\n",file=stderr())
        stop()
    }
}

# At this point we can assume that call_indices is a vector of length > 0
# which indicates the columns of countmat on which to call CNVs
# let's make some assertions about that:
stopifnot(exists("call_indices"))
stopifnot(length(call_indices) > 0)

call_names = colnames(countmat)[call_indices]

# Now try to read .genome file
if (file.exists(opt$genome)) {
    if (opt$verbose) {
        cat(paste("Reading .genome IBD file from: ",opt$genome,
        "\nin order to exclude relateds from reference sets.\n\n",
        sep=""),file=stdout())
    }
    genome = read.table(opt$genome,header=TRUE)
    exclude = genome[genome$PI_HAT > opt$threshold,c("IID1","IID2")]
    exclude_proportion = as.numeric(dim(exclude)[1]/dim(genome)[1])
    if (opt$verbose) {
        cat(paste("Based on your threshold of ",opt$threshold*100,"% IBD, ",
            formatC(exclude_proportion*100,digits=2),"% of possible\n",
            "combinations of people will be excluded from the reference\n\n",
            sep=""),file=stdout())
    }
    # make the sample names from PLINK match the R-ified BAM names
    exclude$IID1 = paste(make.names(exclude$IID1),".bam",sep="")
    exclude$IID2 = paste(make.names(exclude$IID2),".bam",sep="")
    
    # match to the names of samples to call on
    genome_iids = unique(c(genome$IID1,genome$IID2))
    missing_names = call_names[!(call_names %in% genome_iids)]
    if (length(missing_names) > 0) {
        cat("Warning: The following samples that you want to call\n",file=stderr())
        cat("CNVs on are not in your .genome IBD file:\n",file=stderr())
        print(missing_names,file=stderr())
    }
    
    # create a named list where the name is the sample name
    # and the value is a vector of other samples to exclude
    exclude_list = list()
    for (individual in call_names) {
        exclude_list[[individual]] = unique(
            c(exclude$IID1[exclude$IID2==individual],
              exclude$IID2[exclude$IID1==individual]))
    }
    
    
} else {
    if (opt$genome=='') {
        if (opt$verbose) {
            cat("You didn't specify a .genome IBD file, so we're just going\n",file=stdout())
            cat("to assume that no individuals are related.\n",file=stdout())
        }
        exclude_list = list()
        for (individual in call_names) {
            exclude_list[[individual]] = character(0)
        }
    } else {
        cat("**Couldn't find the .genome file you specified:\n",file=stderr())
        cat(opt$genome,file=stderr())
        cat("\n\n",file=stderr())
        stop()
    }
}


if (file.exists(opt$exclude)) {
    if(opt$verbose) {
        cat(paste("Loading list to exclude from ",opt$exclude," ...\n",sep=""),file=stdout())
    }
    exclude_all = read.table(opt$exclude)$V1
    # get the exclude_all list to match the R-ified BAM names used in the counts data
    exclude_all = paste(make.names(exclude_all),".bam",sep="")

    if(opt$verbose) {
        cat(paste("Loaded.\n",sep=""),file=stdout())
    }
    if(opt$verbose) {
        cat("The following samples will be excluded from the reference for all CNV calls:\n",file=stdout())
        print(exclude_all,file=stdout())
    }
    
    for (individual in call_names) {
            exclude_list[[individual]] = c(exclude_list[[individual]], exclude_all)
    }
    
} else if (opt$exclude == '') {
    if(opt$verbose) {
        cat(paste("No exclude all list (-e) detected, continuing...\n",sep=""),file=stdout())
    }
    exclude_all = character(0) # just make it an empty vector
 } else {
    cat("**Exclude all list file not found. You specified this path:\n",file=stderr())
    cat(paste(opt$exclude,"\n",sep=""),file=stderr())
}


# now we are done reading input, so switch to output directory.
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


for (i in call_indices) {
    sample_name = colnames(countmat)[i]
    if(opt$verbose) {
        cat(paste("Now about to call CNVs on ",sample_name,"...",sep=""),file=stdout())
        cat("Excluding the following samples from reference:\n",file=stdout())
        print(exclude_list[[sample_name]],file=stdout())
    }
    
    exclude_indices = match(colnames(countmat),exclude_list[[sample_name]]) # individual excludes
    exclude_indices = c(exclude_indices,match(colnames(countmat),exclude_all)) # excludes from all samples
    exclude_indices = exclude_indices[!is.na(exclude_indices)] # remove NA
    
    # call ExomeDepth to select best reference set, excluding
    # the sample itself as well as all IBD-excluded samples from the
    # set of potential references.
    reference_list = select.reference.set(test.counts = countmat[,i], 
        reference.count = as.matrix(countmat[,c(-i,-exclude_indices)]), # as.matrix in case only 1 col
        bin.length=(countdf$end-countdf$start)/1000,
        n.bins.reduced = 10000)
    
    if(opt$verbose) {
        cat("Using the following samples as reference:\n",file=stdout())
    }
    # whether verbose or not, write to disk what the reference set was, for later use
    refpath = paste(sample_name,".reflist.txt",sep="")
    cat(paste("# SAMPLE: ",sample_name,"\n",sep=""),file=refpath)
    cat(paste("# REFERENCE: ",paste(reference_list$reference.choice,collapse=" "),"\n",sep=""),file=refpath,append=TRUE)
    cat(paste(paste(colnames(countmat),collapse="\t"),"\n",sep=""),file=refpath,append=TRUE)
    cat(paste(paste(1*(colnames(countmat) %in% reference_list$reference.choice),collapse="\t"),"\n",sep=""),file=refpath,append=TRUE)

    reference_set = apply(
        X = as.matrix(countdf[, reference_list$reference.choice]), 
        MAR=1, FUN=sum)
    
    all_exons = new('ExomeDepth', test=countmat[,i], 
        reference=reference_set,
        formula = 'cbind(test,reference) ~ 1')
    
    all_exons = CallCNVs(x = all_exons, transition.probability=10^-4,
        chromosome=countdf$space, start=countdf$start,
        end=countdf$end, name=countdf$names)
    
    # write out cnv calls
    cnvpath = paste(sample_name,".csv",sep='')
    cat(",",file=cnvpath) # to *exactly* match Fengmei's file format, need an initial comma.
    write.table(all_exons@CNV.calls, file=cnvpath, append=TRUE,
        sep=',', row.names=TRUE, col.names=TRUE, quote=FALSE)
    
    if (opt$verbose) {
        cat(paste("..done.\n",sep=''),file=stdout())
    }
}

duration = format(Sys.time() - start_time)

if(opt$verbose) {
    cat(paste("Completed execution in ",duration,".\n",sep=''),file=stdout())
}



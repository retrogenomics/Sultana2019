#!/usr/bin/env Rscript

# hotspot discovery using Poisson distribution and FDR

# import arguments from bash script
args <- commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# set base directory
setwd(args[1])

# import files into data tables or values
nb_ins_per_bin <- read.delim(args[2], header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", check.names = FALSE)
perbase_ins_rate <- as.numeric(scan(args[3]))
output_file <- as.character(args[4])
cover <- read.table(args[5],header=FALSE)

# Give names to the cover file, so that we can use it
# in the same order as nb_ins_per_bin later on
rownames(nb_ins_per_bin)<-paste(nb_ins_per_bin[,1],nb_ins_per_bin[,2],nb_ins_per_bin[,3])
rownames(cover)<-paste(cover[,1],cover[,2],cover[,3])

# Compute the effective cover corrective factor, normalized so that the absolute value of coverage in the HelaS3 run has no incidence on results.
cover$val<-cover$V6*sum(as.double(cover$V5))/sum(as.double(cover$V6))

# calculate poisson p-values and FRD-corrected p-values
delta<-0.1
poisson_p <- ppois(nb_ins_per_bin[,7]-delta, lambda = perbase_ins_rate * cover[rownames(nb_ins_per_bin),"val"],lower.tail=FALSE)
poisson_fdr <- p.adjust(poisson_p, method = "fdr")
mydata <- cbind(nb_ins_per_bin, poisson_p, poisson_fdr)

# export result table
write.table(mydata, file = output_file, sep = "\t", row.names=FALSE, quote=FALSE)


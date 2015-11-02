#! /usr/bin/env Rscript --vanilla

library(optparse)
library(infotheo)
library(minet)

# Argument handling
option_list <- list(
  make_option(c("-i", "--input"), help = "Input filename. Mandatory."),
  make_option(c("-o", "--out"), help = "Output filename. Mandatory."),
  make_option(c("-m", "--method"), default="mrnet",
      help = 'Mutial Information estimation method. One of "clr", "mrnet," "mrnetb", "aracne".'),
  make_option(c("-n", "--num"), action="store", default=100000, type="integer",
      help="maximum number of coexpressions to return")
  
)

opt.parser <- OptionParser(option_list = option_list, 
  description='Builds mutual information matrix from a table of gene expression data using "build.mim" from "minet" package in R.')

opts <- parse_args(opt.parser)

mim <- read.table(opts$input, header=TRUE, row.names=1)

miAlgorithm <- function(mim, method) {
  gene_names <- rownames(mim)
  if (method == "clr") 
    net <- clr(mim)
  else if (method == "mrnet") 
    net <- mrnet(mim)
  else if (method == "mrnetb") 
    net <- mrnetb(mim)
  else if (method == "aracne") 
    net <- aracne(mim)
#  else if (method == "c3") 
#    net <- c3(mim)
  
#  net <- net/max(net)
  pos <- as.data.frame(which(net>0, arr.ind=T))
  ans <- data.frame(A=gene_names[pos$row], B=gene_names[pos$col],
		    MI=net[net>0])
  return(ans)
}

net <- miAlgorithm(mim, opts$method)
o <- order(net$MI, decreasing=TRUE)
N <- min(opts$n, length(o))
net <- net[o[1:N],]
write.table(net, file=opts$out, quote=F, sep="\t")


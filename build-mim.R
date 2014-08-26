#! /usr/bin/env Rscript --vanilla

library(optparse)
library(infotheo)
library(minet)

# Argument handling
option_list <- list(
  make_option(c("-i", "--input"), help = "Input filename. Mandatory."),
  make_option(c("-o", "--out"), help = "Output filename. Mandatory."),
  make_option(c("-e", "--estimator"), default="pearson",
      help = 'Mutial Information estimation method. One of "pearson", "spearman", "kendall", "mi.empirical", "mi.mm", "mi.shrink", "mi.sg".'),
  make_option(c("-d", "--discretization"), default="equalfreq",
      help = 'Discretization method. One of "equalfreq", "equalwidth", "globalequalwidth".')
)

opt.parser <- OptionParser(option_list = option_list, 
  description='Builds mutual information matrix from a table of gene expression data using "build.mim" from "minet" package in R.')

opts <- parse_args(opt.parser)

M <- read.table(opts$input, header=T, sep="\t", quote="", row.names=1)
gene_names <- rownames(M)
expdata <- t(M)
rm(M)

mim <- build.mim(expdata, estimator=opts$estimator, disc=opts$discretization)

save(mim, opts, file=opts$out)

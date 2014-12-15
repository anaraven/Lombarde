#! /usr/bin/env Rscript --vanilla

library(methods)
library(igraph)
library(optparse)

# Argument handling
option_list <- list(
  make_option(c("-c", "--coexp"), help = "Filename of coexpressions."),
  make_option(c("-n", "--count"), action="store", default=1000, help="Number of arcs to swap."),
  make_option(c("-o", "--outfile"), help = "Output filename.", default="/dev/stdout")
)

opt.parser <- OptionParser(option_list = option_list, 
                           description="Usage: randomizeInfluence.R [options]")
opts <- parse_args(opt.parser)

g <- read.graph(opts$coexp, format="ncol")

write.graph(rewire(g, niter=opts$count), opts$outfile, format="ncol")


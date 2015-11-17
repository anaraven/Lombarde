#!/usr/bin/Rscript --vanilla

library(optparse)
library(methods)
library(igraph)

# Argument handling
option_list <- list(
  make_option(c("--header"), action="store_true", default=FALSE,
              help="Show field names as header"),
  make_option(c("-m","--mode"), action="store", default="all",
              help="One of all, in, out, total"),
  make_option(c("--delim","-d"), action="store", default="\t",
              help="Field delimiter in output"),
  make_option(c("-c", "--cols"), action="store", default=10, type="integer",
              help="Number of columns on the result.")
)

opt.parser <- OptionParser(option_list = option_list,
  description="Alternative usage: Rscript --vanilla degree-dist.R [options] test-graph.ncol ...")

opts <- parse_args(opt.parser, positional_arguments = TRUE)

argv <- opts$args
if(length(argv) < 1) {
  print_help(opt.parser)
  quit(save = "no", status = 1)
}

N <- opts$options$cols

for(net.file in argv)  {
    g <- read.graph(net.file, format="ncol", names=T, direct=T, weights="yes")
    d <- degree.distribution(g, mode=opts$options$mode)
    dd <- round(100*d[2:N],1)
    dd[N] <- round(100*sum(d[-(1:N)]),1)

    cat(net.file, paste(dd, sep=opts$options$delim), "\n", sep=opts$options$delim)  
} 


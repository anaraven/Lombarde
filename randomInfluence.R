#! /usr/bin/env Rscript --vanilla

library(igraph)
library(optparse)

# Argument handling
option_list <- list(
  make_option(c("-N", "--nodes"),  help = "Filename of node names."),
  make_option(c("-n", "--count"), action="store", default=10000,
              help="Number of arcs to generate."),
  make_option(c("-c", "--column"), action="store", default=2, type="integer",
              help="Column id of node names Default: 2.")
)

opt.parser <- OptionParser(option_list = option_list, 
                           description="Usage: reachability.R [options] net.file [...]")

# arg <-strsplit("--nodes Input/RDB8.1/operon_names.txt --count 10000"," ")[[1]]
# opts <- parse_args(opt.parser, args=arg)
opts <- parse_args(opt.parser)
ops <- read.table(opts$nodes, as.is=T)
ops <- ops[, opts$column]

size <- opts$count
g <- erdos.renyi.game(length(ops), size, type="gnm")
edges <- get.edges(g,E(g))
ans <- data.frame(X=ops[edges[,1]], Y=ops[edges[,2]], Z=runif(size))

write.table(ans, file="",
            sep="\t", quote=F, row.names=F, col.names=F)

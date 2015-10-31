#!/usr/bin/env Rscript --vanilla

library(methods)
library(optparse)
library(igraph)
library(parallel)

# Argument handling
option_list <- list(
  make_option(c("-o", "--out"), help = "Filename of output graph (\"ncol\" format)"),
  make_option(c("-a", "--asp-out"), help = "Filename of solutions descriptions (to process in ASP)",
	      dest="asp.out"),
  make_option(c("-w", "--wgt"), action="store_true", default=FALSE,
              help="Use weights as specified in the file. Don't recalculate weights."),
  make_option(c("-b", "--base"), action="store", default=10, type="double",
              help="Base to use in weight conversion"),
  make_option(c("-c", "--cores"), action="store", default=detectCores(), type="integer",
              help="Number of parallel process to run. Default: all cores.")
)

opt.parser <- OptionParser(option_list = option_list, 
  description="Usage: lombarde.R [options] coexp.file net.file")

opts <- parse_args(opt.parser, positional_arguments = TRUE)

argv <- opts$args
if(length(argv) < 2) {
  print_help(opt.parser)
  quit(save = "no", status = 1)
}
coexp.file <- argv[1]
net.file   <- argv[2] 
out.file   <- opts$options$out
asp.file   <- opts$options$asp.out

options(mc.cores=opts$options$cores, digits=10, scipen=3)

# Step one: read input graph
g  <- read.graph(net.file, format="ncol", names=TRUE, direct=TRUE, weights="yes")
if(! opts$options$wgt) {
  message("changing weights ",opts$options$base)
  E(g)$orig   <- E(g)$weight
  E(g)$weight <- opts$options$base^(E(g)$weight)
}
g.name <- V(g)$name

# Step two: read coexpressions
coexps <- read.table(coexp.file, as.is=TRUE, col.names=c("from", "to", "weight"))

# keep only coexpressions involving valid vertices (those in `g` graph)
is.valid.obs <- coexps$from %in% g.name & coexps$to %in% g.name # & coexps$from < coexps$to
coexps <- coexps[is.valid.obs,]
N <- nrow(coexps)

# `W` is the cost of the shortst path betwen each pair of vertices
W <- shortest.paths(g, mode="out")

# `shared_pred` is a list that has, for each pair in `coexps`, the list of
# their minimal cost common predecessors 
shared_pred <- mclapply(1:N, function(i) {
			v <- rowSums(W[,unlist(coexps[i,1:2])]);
			names(which(v==min(v)))
})

for(i in 1:N) {
    cat(unlist(coexps[i,1:2]), shared_pred[[i]], "\n", file=asp.file, append=TRUE)
}

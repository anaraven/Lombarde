#! /usr/bin/env Rscript --vanilla

library(optparse)
library(igraph)
library(parallel)

# Argument handling
option_list <- list(
  make_option(c("-N", "--net"),  help = "Filename of transcriptional regulation network. Mandatory."),
  make_option(c("-n", "--count"), action="store_true", default=FALSE,
              help="Only shows the number of reachable arcs."),
  make_option(c("-c", "--cores"), action="store", default=detectCores(), type="integer",
              help="Number of parallel process to run. Default: all cores.")
)

opt.parser <- OptionParser(option_list = option_list, 
  description="Usage: reachability.R [options] net.file [...]")

# arg <-strsplit("-v --obs Input/RDB8.1/M3D/pearson/aracne/10K.ncol Input/RDB8.1/RegulonDB/MEME/g0-3.ncol"," ")[[1]]
# opts <- parse_args(opt.parser, positional_arguments = TRUE, args=arg)
opts <- parse_args(opt.parser, positional_arguments = TRUE)
net.file <- opts$options$net  # coexpressions to explain

if(length(opts$args) <1 | is.null(net.file)) {
  print_help(opt.parser)
  quit(save = "no", status = 1)
}

options(mc.cores=opts$options$cores)

g  <- read.graph(net.file, format="ncol", names=T, direct=T, weights="yes")
nom <- V(g)$name
# `W` shows the cost of the cheapest path betwen vertices
W <- shortest.paths(g, mode="out")
  
getCoverage <- function(coexp.file, nom, W) {
  obs <- read.table(coexp.file, as.is=T, col.names=c("from", "to", "weight"))
  N0 <- nrow(obs)

  is.valid.obs <- obs$from %in% nom & obs$to %in% nom
  v.obs <- obs[is.valid.obs,]
  N1 <- nrow(v.obs)

  # `TF` is the list of common predecessors for each pair in `v.obs`
  TF <- sapply(1:N1, function(i) {v <- rowSums(W[,unlist(v.obs[i,1:2])]); is.finite(min(v))}) 
  N2 <- sum(TF)
  return(list(coexp.file=coexp.file, N0=N0, N1=N1, N2=N2))
}

all.cvg <- mclapply(opts$args, getCoverage, nom, W)
for(cvg in all.cvg) {
  with(cvg, cat(coexp.file, net.file, N0, N1, N2, 100*N1/N0, 100*N2/N0, sep="\t"))
  cat("\n")
}


#! /usr/bin/env Rscript --vanilla

library(optparse)
library(igraph)
library(parallel)

# Argument handling
option_list <- list(
  make_option(c("-O", "--obs"),  help = "Filename of Observed Coexpressed pairs. Mandatory"),
  make_option(c("-g", "--gold"), help = "Filename of Gold Standard graph. Mandatory"),
  make_option(c("-i", "--intersect"), action="store_true", default=FALSE,
              help="Only considers reachability of arcs in the intersection."),
  make_option(c("-n", "--count"), action="store_true", default=FALSE,
              help="Only shows the number of reachable arcs."),
  make_option(c("-t", "--table"), action="store_true", default=FALSE,
              help="Shows the table of degree (?)."),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Be more explicit on what is being done."),
  make_option(c("-c", "--cores"), action="store", default=detectCores(), type="integer",
              help="Number of parallel process to run. Default: all cores.")
)

opt.parser <- OptionParser(option_list = option_list, 
  description="Usage: reachability.R [options] net.file [...]")

arg <-strsplit("--intersect --obs Input/RDB8.1/M3D/pearson/mrnet/100K.ncol --gold Input/RDB8.1/gold-std.ncol Input/RDB8.1/Prodoric/MEME/g0-3.ncol "," ")[[1]]
opts <- parse_args(opt.parser, positional_arguments = TRUE, args=arg)
# opts <- parse_args(opt.parser, positional_arguments = TRUE)
coexp.file <- opts$options$obs  # coexpressions to explain
gs.file    <- opts$options$gold # gold standard filename

if(length(opts$args) <1 | is.null(coexp.file) | is.null(gs.file)) {
  print_help(opt.parser)
  quit(save = "no", status = 1)
}

if(opts$options$verbose) {
  cat("coexp.file",coexp.file,"\n")
#  cat("net.file",net.file,"\n")
  cat("gs.file",gs.file,"\n")
#  quit(save = "no", status = 0)
}

options(mc.cores=opts$options$cores)

gs <- read.graph(gs.file, format="ncol", names=T, direct=T, weights="no")
# we store the gold-standard vertices names in a vector for easy access 
gs.v.names <- V(gs)$name
gs.edges <- get.edges(gs, E(gs))
# ...trasnformed in a two-level tree of node names
X <- gs.v.names[gs.edges[,1]]
Y <- gs.v.names[gs.edges[,2]]

obs <- read.table(coexp.file, as.is=T, col.names=c("from", "to", "weight"))

for(net.file in opts$args) {
  g  <- read.graph(net.file, format="ncol", names=T, direct=T, weights="yes")
  nom <- V(g)$name
  
  is.valid.obs <- obs$from %in% nom & obs$to %in% nom # & v.obs$from < v.obs$to
  v.obs <- obs[is.valid.obs,]
  N <- nrow(v.obs)
  # cat("n.obs",N,"\n")
  if(opts$options$verbose) cat("is.valid.obs",table(is.valid.obs),"\n")

  # `W` shows the cost of the cheapest path betwen vertices
  W <- shortest.paths(g, mode="out")
  # `M` is a boolean matrix that indicates if there is a path from $i$ to $j$
  M <- is.finite(W)
  # `TF` is the list of common predecessors for each pair in `v.obs`
  TF <- mclapply(1:N, function(i) {v <- rowSums(W[,unlist(v.obs[i,1:2])]); which(v==min(v))})
  # TODO: tengo la duda si considerar sólo el que está a mínima distancia o todos
  if(opts$options$table) {
    print(table(sapply(TF,length)))
  }
  
  if(opts$options$intersect) {
    g.edges <- get.edges(g, E(g))
    # ...trasnformed in a two-level tree of node names
    included <- sapply(1:length(X), function(i) any(nom[g.edges[,1]] %in% X[i] & nom[g.edges[,2]] %in% Y[i]))
  } else {
    included <- rep(TRUE,length(X))
  }

  X.i <- X[included]
  Y.i <- Y[included]
  X.nom <- intersect(X.i, nom) # node names valid for both `gs` and `net` graph
  Y.nom <- intersect(Y.i, nom)
  
  # tests if $i$-th TF reaches each `X.nom`
  reach_TF_to_x <- function(i) apply(M[TF[[i]], X.nom, drop=F], 2, any)
  # tests if `Y.nom` reaches any node of $i$-th coexpression pair
  reach_y_to_AB <- function(i) apply(M[Y.nom, unlist(v.obs[i,1:2]), drop=F], 1, any)
  
  A <- mcmapply(reach_TF_to_x, 1:N) # rows:gs nodes x, cols:observations, cell: TF for observation j reaches arc x[i]
  B <- mcmapply(reach_y_to_AB, 1:N) # rows:gs nodes y, cols:observations, cell: arc y[i] reaches a node in onservation j 
  C <- (A %*% t(B))>0
  
  valid.X <- X.i %in% nom
  valid.Y <- Y.i %in% nom
  connected <- sapply(1:length(X.i),
      function(i) ifelse(valid.X[i] & valid.Y[i], C[X.i[i],Y.i[i]], FALSE))
  cat(coexp.file, net.file, sum(connected), sum(valid.X & valid.Y), sep="\t")
  cat("\n")
}



#for(i in 1:length(X)) {
#  connected <- ifelse(valid.X[i] & valid.Y[i], C[X[i],Y[i]], FALSE)
#  cat(X[i],Y[i], valid.X[i], valid.Y[i],connected, "\n")
#}

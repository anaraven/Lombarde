#!/usr/bin/env Rscript --vanilla

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
g  <- read.graph(net.file, format="ncol", names=T, direct=T, weights="yes")
if(! opts$options$wgt) {
  message("changing weights")
  E(g)$orig   <- E(g)$weight
  E(g)$weight <- opts$options$base^(E(g)$weight)
}
g.name <- V(g)$name

# Step two: read coexpressions
coexps <- read.table(coexp.file, as.is=T, col.names=c("from", "to", "weight"))

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
			which(v==min(v))
})

path.extremes <- function(i, shared_pred, coexps) {
# this function takes the i-th co-expressed pair of vertices and returns a list
# with all the pairs of vertices that define paths connecting each common
# predecesor to both co-expressed vertices.

  unlist(lapply(names(shared_pred[[i]]),
		function(r,v) list(c(r, v$from), c(r, v$to)), coexps[i,]),
	 recursive=F, use.names=F)
}

# determine the non-redundant set of pairs of vertices that determine the
# relevant paths to evaluate
expl.path <- unique(unlist(mclapply(1:N, path.extremes, shared_pred, coexps),
			   recursive=F, use.names=F))

# change the representation from a list of (from,to) pairs to a two level tree.
# Recycle the varaible to save memory
expl.path <- split(sapply(expl.path,`[`,2), sapply(expl.path,`[`,1))

expl.path <- mcmapply(function(src, targets) mapply(function(tgt) {
  get.all.shortest.paths(g, src, tgt, mode="out")$res
  }, targets, SIMPLIFY=FALSE), names(expl.path), expl.path, SIMPLIFY=FALSE)
  
expl.path <- mclapply(expl.path, lapply, lapply, function(l) {
  if(length(l)>1) as.numeric(E(g, path=l)) else numeric()} )

if(!is.null(out.file)) {
  write.graph(subgraph.edges(g, unique(unlist(expl.path))), out.file, format="ncol")
}

if(!is.null(asp.file)) {
  cat("n.obs",N,"\n", file=asp.file)
  vid <- 1
  for(i in 1:N){
    for(r in names(shared_pred[[i]])){
      cat("explanation",r,unlist(coexps[i,1:2]),"\n", file=asp.file, append=T)
      vv1 <- expl.path[[r]][[ coexps[[i,1]] ]]
      vv2 <- expl.path[[r]][[ coexps[[i,2]] ]]
      for(p1 in vv1) {
        for(p2 in vv2) {
          cat("vshape", vid, unlist(coexps[i, 1:2]), i, "\n", file=asp.file, append=T)
          edgelist <- c(p1,p2)
          weights <- get.edge.attribute(g, "weight", edgelist)
          sides <- get.edges(g, edgelist)
          for(j in 1:length(edgelist)) {
            cat("arcInVshape", vid, g.name[sides[j,1]], g.name[sides[j,2]], weights[j], "\n",
		file=asp.file, append=T)        
          }
          vid <- vid + 1
        }
      }
    }
  }
}

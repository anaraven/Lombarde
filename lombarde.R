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
cat("N:",N,"\n")

# `W` is the cost of the shortst path betwen each pair of vertices
W <- shortest.paths(g, mode="out")

# `shared_pred` is a list that has, for each pair in `coexps`, the list of
# their minimal cost common predecessors 
shared_pred <- mclapply(1:N, function(i) {
			v <- rowSums(W[,unlist(coexps[i,1:2])]);
			names(which(v==min(v)))
})
cat("shared_pred:",length(shared_pred),"\n")

# this function takes the i-th co-expressed pair of vertices and returns a list
# with all the pairs of vertices that define paths connecting each common
# predecesor to both co-expressed vertices.
path.extremes <- function(i, shared_pred, coexps) {
  unlist(lapply(shared_pred[[i]],
		function(r,v) list(c(r, v$from), c(r, v$to)), coexps[i,]),
	 recursive=FALSE, use.names=FALSE)
}

# determine the non-redundant set of pairs of vertices that determine the
# relevant paths to evaluate
expl.path <- unique(unlist(mclapply(1:N, path.extremes, shared_pred, coexps),
			   recursive=FALSE, use.names=FALSE))
names(expl.path) <- sapply(expl.path, paste, collapse=" ")
cat("number of extremes:",length(expl.path),"\n")
# count the number of non-trivial paths
non.trivial <- sapply(expl.path, function(e) e[1]!=e[2])
print(sum(non.trivial))

# replace each pair (a,b) for a list of all short-path-a-b
expl.path <- mcmapply(function(e) {
  get.all.shortest.paths(g, e[1], e[2], mode="out")$res}, expl.path[non.trivial])
# now each element of the list is a list of all paths between (a,b) vertices
# and "a b" is the name of the element.
# each path is a list of the vertices in the corresponding order
# which can be transformed into edges with `as.numeric(E(g, path=l))`

npath <- table(sapply(expl.path, length))
cat("number of paths:",sum(as.numeric(names(npath))*npath),"\n")
cat("complexity:", round(sum(log10(as.numeric(names(npath)))*npath)),"\n")

if(!is.null(out.file)) {
  cat("Writing output to", out.file,"\n")
  all.edges <- lapply(unlist(expl.path[non.trivial], recursive = F),
                      function(l) E(g, path=l))
  write.graph(subgraph.edges(g, unique(unlist(all.edges))), out.file, format="ncol")
}

print.arcs <- function(edgelist) {
  weights <- get.edge.attribute(g, "weight", edgelist)
  sides <- get.edges(g, edgelist)
  for(j in 1:length(edgelist)) {
    cat("arcInVshape", vid, g.name[sides[j,1]], g.name[sides[j,2]], weights[j], "\n",
        file=asp.file, append=TRUE)        
  }
}

if(!is.null(asp.file)) {
  cat("Writing structure to", asp.file,"\n")
  cat("n.obs",N,"\n", file=asp.file)
  vid <- 1
  for(i in 1:N){
# <<<<<<< HEAD
    a <- coexps[i,1]
    b <- coexps[i,2]
    for(r in shared_pred[[i]]){
      cat("explanation",r,a,b,i,"\n", sep="\t",file=asp.file, append=TRUE)
      vv1 <- expl.path[[paste(r,a)]]
      vv2 <- expl.path[[paste(r,b)]]
      if(is.null(vv1)) {
        for(p2 in vv2) {
          cat("vshape", vid, a, b, i, "\n", file=asp.file, append=TRUE)
          edgelist <- as.numeric(E(g, path=p2))
          print.arcs(edgelist)
          vid <- vid + 1
        }
      } else if(is.null(vv2)) {
        for(p1 in vv1) {
            cat("vshape", vid, a, b, i, "\n", file=asp.file, append=TRUE)
            edgelist <- as.numeric(E(g, path=p1))
            print.arcs(edgelist)
            vid <- vid + 1
        }
      } else
# =======
#     for(r in shared_pred[[i]]){
#       cat("explanation",r,unlist(coexps[i,1:2]),"\n", file=asp.file, append=TRUE)
#       vv1 <- expl.path[[r]][[ coexps[[i,1]] ]]
#       vv2 <- expl.path[[r]][[ coexps[[i,2]] ]]
# >>>>>>> feature/no_discret
      for(p1 in vv1) {
        for(p2 in vv2) {
          cat("vshape", vid, a, b, i, "\n", file=asp.file, append=TRUE)
          edgelist <- as.numeric(c(E(g, path=p1),E(g, path=p2)))
          print.arcs(edgelist)
          vid <- vid + 1
        }
      }
    }
  }
}

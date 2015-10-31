library(methods)
library(igraph)
library(parallel)

coexp.file <- "~/RLombarde/Input/RDB8.1/M3D/pearson/mrnet/10K.ncol"
net.file   <- "~/RLombarde/Input/RDB8.1/Prodoric/MEME/g0-9.ncol"
out.file   <- "/dev/stdout"
asp.file   <- "/dev/stdout"
change_weights <- TRUE
base <- 10

options(mc.cores=4, digits=10, scipen=3)

# Step one: read input graph
g  <- read.graph(net.file, format="ncol", names=TRUE, direct=TRUE, weights="yes")
if(change_weights) {
  message("changing weights ",base)
  E(g)$orig   <- E(g)$weight
  E(g)$weight <- base^(E(g)$weight)
}
g.name <- V(g)$name

# Step two: read coexpressions
coexps <- read.table(coexp.file, as.is=TRUE, col.names=c("from", "to", "weight"))

# keep only coexpressions involving valid vertices (those in `g` graph)
is.valid.obs <- coexps$from %in% g.name & coexps$to %in% g.name # & coexps$from < coexps$to
coexps <- coexps[is.valid.obs,]
N <- nrow(coexps)
print(paste("N:",N))

# `W` is the cost of the shortst path betwen each pair of vertices
W <- shortest.paths(g, mode="out")

# `shared_pred` is a list that has, for each pair in `coexps`, the list of
# their minimal cost common predecessors 
shared_pred <- mclapply(1:N, function(i) {
			v <- rowSums(W[,unlist(coexps[i,1:2])]);
			names(which(v==min(v)))
})
print(paste("shared_pred",length(shared_pred)))

path.extremes <- function(i, shared_pred, coexps) {
# this function takes the i-th co-expressed pair of vertices and returns a list
# with all the pairs of vertices that define paths connecting each common
# predecesor to both co-expressed vertices.

  unlist(lapply(shared_pred[[i]],
		function(r,v) list(c(r, v$from), c(r, v$to)), coexps[i,]),
	 recursive=FALSE, use.names=FALSE)
}

# determine the non-redundant set of pairs of vertices that determine the
# relevant paths to evaluate
expl.path <- unique(unlist(mclapply(1:N, path.extremes, shared_pred, coexps),
			   recursive=FALSE, use.names=FALSE))
names(expl.path) <- sapply(expl.path, paste, collapse=" ")
print(paste("number of extremes:",length(expl.path)))
# count the number of non-trivial paths
non.trivial <- sapply(expl.path, function(e) e[1]!=e[2])
print(table(non.trivial))

# replace each pair (a,b) for a list of all short-path-a-b
expl.path <- mcmapply(function(e) {
  get.all.shortest.paths(g, e[1], e[2], mode="out")$res}, expl.path[non.trivial])
# now each element of the list is a list of all paths between (a,b) vertices
# and "a b" is the name of the element.
# each path is a list of the vertices in the corresponding order
# which can be transformed into edges with `as.numeric(E(g, path=l))`

npath <- table(sapply(expl.path, length))
print(paste("number of paths:",sum(as.numeric(names(npath))*npath)))
print(paste("complexity:", round(sum(log10(as.numeric(names(npath)))*npath))))

if(!is.null(out.file)) {
  all.edges <- lapply(unlist(expl.path[non.trivial], recursive = F),
                      function(l) E(g, path=l))
  write.graph(subgraph.edges(g, unique(unlist(all.edges))), out.file, format="ncol")
}

if(!is.null(asp.file)) {
  cat("n.obs",N,"\n", file=asp.file)
  vid <- 1
  for(extremes in names(expl.path)) {
    for(l in 1:length(expl.path[[extremes]])) {
      edgelist <- as.numeric(E(g, path=expl.path[[extremes]][[l]]))
      weights <- get.edge.attribute(g, "weight", edgelist)
      for(j in 1:length(edgelist)) {
        cat("arcInPath\tp", vid, "\te", edgelist[j],"\t", 
            weights[j], "\n", file=asp.file, append=TRUE, sep="")
      }
      expl.path[[extremes]][[l]] <- vid
      vid <- vid+1
    }
  }
  for(i in 1:N){
    a <- coexps[i,1]
    b <- coexps[i,2]
    for(r in shared_pred[[i]]){
      cat("shared",r,a,b,i,"\n", sep="\t",file=asp.file, append=TRUE)
      for(j in expl.path[[paste(r,a)]]) {
        cat("path\t",r,"\t",a,"\tp",j,"\n", sep="",file=asp.file, append=TRUE)
      }
      for(j in expl.path[[paste(r,b)]]) {
        cat("path\t",r,"\t",b,"\tp",j,"\n", sep="",file=asp.file, append=TRUE)
      }
    }
  }
}

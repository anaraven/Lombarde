#! /usr/bin/Rscript --vanilla

### Functions ###
usage <- function() {
  message("Usage: Rscript --vanilla graph-diff.R g1.ncol g2.ncol")
}

# Argument handling
argv <- commandArgs(T)
if(length(argv) != 2) {
  usage()
  quit(save = "no", status = 1)
}

library(igraph)

g1.file  <- argv[1] 
g2.file  <- argv[2]
g1 <- read.graph(g1.file, format="ncol", names=T, direct=T, weights="yes")
g2 <- read.graph(g2.file, format="ncol", names=T, direct=T, weights="yes")

g1.v.names <- V(g1)$name
g1.edges <- get.edges(g1, E(g1))
g1.wgt   <- E(g1)$weight
names(g1.wgt) <- g1.v.names[g1.edges[,2]]
g1.wgt   <- split(g1.wgt, g1.v.names[g1.edges[,1]])
g1.edges <- split(g1.v.names[g1.edges[,2]], g1.v.names[g1.edges[,1]])


g2.v.names <- V(g2)$name
g2.edges <- get.edges(g2, E(g2))
g2.wgt   <- E(g2)$weight
names(g2.wgt) <- g2.v.names[g2.edges[,2]]
g2.wgt   <- split(g2.wgt, g2.v.names[g2.edges[,1]])
g2.edges <- split(g2.v.names[g2.edges[,2]], g2.v.names[g2.edges[,1]])

in.both <- sapply(names(g1.edges), function(x) intersect(g1.edges[[x]], g2.edges[[x]]))

in.g1   <- sapply(names(g1.edges), function(x) setdiff(g1.edges[[x]], g2.edges[[x]]))
in.g2   <- sapply(names(g2.edges), function(x) setdiff(g2.edges[[x]], g1.edges[[x]]))

for(x in names(in.both))
  for(y in in.both[[x]])
    if(g1.wgt[[x]][[y]]!=g2.wgt[[x]][[y]])
      cat("I", x, y, g1.wgt[[x]][[y]], g2.wgt[[x]][[y]],"\n")

for(x in names(in.g1))
  for(y in in.g1[[x]])
    cat("A", x, y, g1.file,"\n")

for(x in names(in.g2))
  for(y in in.g2[[x]])
    cat("B", x, y, g2.file,"\n")


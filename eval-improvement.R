#!/usr/bin/Rscript --vanilla

library(optparse)
library(methods)
library(igraph)

# Argument handling
option_list <- list(
  make_option(c("-u","--undirected"), action="store_true", default=FALSE,
              help="Input graphs are undirected"),
  make_option(c("-p","--pretty"), action="store_true", default=FALSE,
              help="print pretty names"),
  make_option(c("--ref","-r"), action="store", 
              help="Reference network of validated arcs"),
  make_option(c("-i", "--ini"), action="store", 
              help="Initial network of putative arcs")
		    
)

opt.parser <- OptionParser(option_list = option_list,
  description="Alternative usage: Rscript --vanilla eval-improvement.R [options] gold-standard.ncol test-graph.ncol ...")

opts <- parse_args(opt.parser, positional_arguments = TRUE)

argv <- opts$args
if(length(argv) < 1) {
  print_help(opt.parser)
  quit(save = "no", status = 1)
}

# function to evaluate graph intersections
analize.graph <- function(g, gs.edges, rv.edges) {
  g.v.names <- V(g)$name
  g.edges <- get.edges(g, E(g))
  g.edges <- split(g.v.names[g.edges[,2]], g.v.names[g.edges[,1]])
  
  in.both <- sum(sapply(names(g.edges),
  	function(x) length(intersect(gs.edges[[x]], g.edges[[x]]))) )
  if(!is.null(rv.edges)) {
    in.both <- in.both + sum(sapply(names(g.edges),
      function(x) length(intersect(rv.edges[[x]], g.edges[[x]]))) )
    
  }
  
  list(vertices=vcount(g), edges=ecount(g), in.both=in.both)
}

# `gs` is the gold standard graph, `g` is the graph to evaluate
gs.file   <- opts$options$ref # "Input/RDB8.1/gold-std.ncol"
gs <- read.graph(gs.file,  format="ncol", names=T, direct=T, weights="no")
# we store the gold-standard vertices names in a vector for easy access 
gs.v.names <- V(gs)$name
# and the edges as a two column matrix...
gs.edges <- get.edges(gs, E(gs))
# ...trasnformed in a two-level tree of node names
rv.edges <- if(opts$options$undirected) split(gs.v.names[gs.edges[,1]], gs.v.names[gs.edges[,2]]) else NULL
gs.edges <- split(gs.v.names[gs.edges[,2]], gs.v.names[gs.edges[,1]])


library(knitr)
if(opts$options$pretty) {
    nom <- "Validated" 
} else {
    nom <- sub(".ncol","",gs.file)
}
n.vert  <-  vcount(gs)
n.edges <-  ecount(gs)
in.gold <-  ecount(gs)
prec    <-  0
c.edges <-  0
c.in.gold <- 0

ini.file   <- opts$options$ini 
ini <- read.graph(ini.file,  format="ncol", names=T, direct=T, weights="no")
l <- analize.graph(ini, gs.edges, rv.edges)
ini.edges   <- l$edges
ini.in.gold <- l$in.both
if(opts$options$pretty) {
    nom     <- c(nom, "Initial")
} else {
    nom     <- c(nom, sub(".ncol", "", ini.file))
}
n.vert  <- c(n.vert, l$vertices)
n.edges <- c(n.edges, l$edges)
in.gold <- c(in.gold, l$in.both)
prec    <- c(prec, l$in.both/l$edges*100)
c.edges <- c(c.edges, l$edges/ini.edges*100)
c.in.gold <- c(c.in.gold, l$in.both/ini.in.gold*100)

for(net.file in argv)  {
    l <- analize.graph(read.graph(net.file, format="ncol", names=T, direct=T, weights="yes"),
             gs.edges, rv.edges)
      nom     <- c(nom, sub(".ncol", "", net.file))
      n.vert  <- c(n.vert, l$vertices)
      n.edges <- c(n.edges, l$edges)
      in.gold <- c(in.gold, l$in.both)
      prec    <- c(prec, l$in.both/l$edges*100)
      c.edges <- c(c.edges, l$edges/ini.edges*100)
      c.in.gold <- c(c.in.gold, l$in.both/ini.in.gold*100)
} 
ans <- data.frame(Graph=nom, Vertices=n.vert, Edges=n.edges, `In gold`=in.gold,
      	    Precision=prec, `Cons Edges`=c.edges, `Cons Valid`=c.in.gold)  
kable(ans, digits=2, format.args=list(zero.print="-.--"))

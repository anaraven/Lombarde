#!/usr/bin/Rscript --vanilla

library(optparse)
library(methods)
library(igraph)

# Argument handling
option_list <- list(
  make_option(c("--induce","-i"), action="store_true", default=FALSE,
              help="Use only the subgraph induced by gold-standard vertices"),
  make_option(c("--header"), action="store_true", default=FALSE,
              help="Show field names as header"),
  make_option(c("-m","--markdown"), action="store_true", default=FALSE,
              help="output table in Markdown format"),
  make_option(c("-u","--undirected"), action="store_true", default=FALSE,
              help="Input graphs are undirected"),
  make_option(c("--delim","-d"), action="store", default="\t",
              help="Field delimiter in output")
)

opt.parser <- OptionParser(option_list = option_list,
  description="Alternative usage: Rscript --vanilla precision-recall.R [options] gold-standard.ncol test-graph.ncol ...")

opts <- parse_args(opt.parser, positional_arguments = TRUE)

argv <- opts$args
if(length(argv) < 2) {
  print_help(opt.parser)
  quit(save = "no", status = 1)
}

# function to evaluate graph intersections
analize.graph <- function(g, gs.edges, rv.edges, induce=FALSE) {
  g.v.names <- V(g)$name
  if(induce){
    g <- delete.vertices(g, g.v.names[!g.v.names %in% gs.v.names])
    g.v.names <- V(g)$name  
  }
  
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

gs.file   <- argv[1] # "Input/RDB8.1/gold-std-operons.ncol"
                     # "Out1/RDB8.1/mi-pearson-aracne-133/graph-operons-3.ncol"

# `gs` is the gold standard graph, `g` is the graph to evaluate
gs <- read.graph(gs.file,  format="ncol", names=T, direct=T, weights="no")
# we store the gold-standard vertices names in a vector for easy access 
gs.v.names <- V(gs)$name
# and the edges as a two column matrix...
gs.edges <- get.edges(gs, E(gs))
# ...trasnformed in a two-level tree of node names
rv.edges <- if(opts$options$undirected) split(gs.v.names[gs.edges[,1]], gs.v.names[gs.edges[,2]]) else NULL
gs.edges <- split(gs.v.names[gs.edges[,2]], gs.v.names[gs.edges[,1]])

if(!opts$options$markdown) {
  options(digits=3)
  if(opts$options$header) {
    if(opts$options$induce) {
      cat("Graph (Induced)", "Vertices", "Edges", "In gold", "Precision", "Recall",
          "F-measure", "\n", sep=opts$options$delim)  
      cat("---------------", "--------", "-----", "-------", "---------", "------",
          "---------", "\n", sep=opts$options$delim)  
    } else {
      cat("Graph", "Vertices", "Edges", "In gold", "Precision", "Recall",
          "F-measure", "\n", sep=opts$options$delim)  
      cat("-----", "--------", "-----", "-------", "---------", "------",
          "---------", "\n", sep=opts$options$delim)  
    }
  }
  cat(gs.file, vcount(gs), ecount(gs), ecount(gs), 100, 100, 100, "\n", sep=opts$options$delim)

  for(net.file in argv[2:length(argv)])  {
      l <- analize.graph(read.graph(net.file, format="ncol", names=T, direct=T, weights="yes"),
               gs.edges, rv.edges, opts$options$induce)

      cat(net.file, l$vertices, l$edges, l$in.both,
          l$in.both/l$edges*100, l$in.both/ecount(gs)*100,
          2*l$in.both/(l$edges+ecount(gs))*100, "\n", sep=opts$options$delim)  
  } 
} else { # markdown output
  library(knitr)
  nom <- gs.file
  n.vert  <-  vcount(gs)
  n.edges <-  ecount(gs)
  in.gold <-  ecount(gs)
  prec    <-  100
  recall  <-  100
  f.score <-  100

  for(net.file in argv[2:length(argv)])  {
      l <- analize.graph(read.graph(net.file, format="ncol", names=T, direct=T, weights="yes"),
               gs.edges, rv.edges, opts$options$induce)
	nom     <- c(nom, net.file)
	n.vert  <- c(n.vert, l$vertices)
	n.edges <- c(n.edges, l$edges)
	in.gold <- c(in.gold, l$in.both)
	prec    <- c(prec, l$in.both/l$edges*100)
	recall  <- c(recall, l$in.both/ecount(gs)*100)
	f.score <- c(f.score, 2*l$in.both/(l$edges+ecount(gs))*100)
  } 
  ans <- data.frame(Graph=nom, Vertices=n.vert, Edges=n.edges, `In gold`=in.gold,
		    Precision=prec, Recall=recall, `F-measure`=f.score)  
  kable(ans, digits=2)
}

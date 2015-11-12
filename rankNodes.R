#! /usr/local/bin/Rscript --vanilla
# ranks each gene for a series of NCOL files

library(optparse)
library(igraph)

### Functions ###
valenteForeman <- function(g, what="radiality") {
  sp <- shortest.paths(g, mode="out", weight=NA)
  rd <- max(sp[is.finite(sp)])+1-sp
  rd[is.infinite(rd)] <- 0
  diag(rd) <- 0
  what <- match.arg(what, c("integration","radiality"))
  if(what=="radiality") {
    return( rowSums(rd)/(ncol(rd)-1) )
  } else {
    return( colSums(rd)/(ncol(rd)-1) )
  }
}

getIndex <- function(g, what="radiality") {
  what <- match.arg(what, c("integration","radiality","outdegree", "indegree","evcenter"))
  switch(what,
         integration = valenteForeman(g, what),
         radiality   = valenteForeman(g, what),
         outdegree   = degree(g, mode="out"),
         indegree    = degree(g, mode="in"),
         evcenter    = evcent(g, directed=T, scale=T)
  )
}

analize.graph <- function(net, what) {
  score <- getIndex(net, what=what)
  if(opts$options$rank) {
    score <- 1 + length(score) - rank(score, ties.method="max")
    score[score==max(score)] <- NA  
  }
  missing <- setdiff(nodes, V(net)$name)
  score[missing] <- NA
  score
}

### Main ###
# Argument handling
option_list <- list(
  make_option(c("-a","--asp"), action="store_true", default=FALSE,
              help="Input files are ASP output"),
  make_option(c("-r", "--rank"), action="store_true", default=FALSE,
              help="Print node ranking instead of index value."),
  make_option(c("-s", "--summary"), action="store_true", default=FALSE,
              help="Only prints summary, not gene-by-gene data."),
  make_option(c("-f", "--funct"), action="store", default="row.summary", 
	      help="Function to use for summary."),
  make_option(c("-i", "--index"), action="store", default="radiality", 
              help="Index used for ranking. Options are 'radiality', 'integration', 'indegree', 'outdegree', 'evcenter'")
)

opt.parser <- OptionParser(option_list = option_list, 
  description="Usage: Rscript --vanilla rankNodes [options] graph.ncol [...]")

opts <- parse_args(opt.parser, positional_arguments = TRUE)
# arg <-strsplit("-i outdegree --asp Prodoric/g0-5-operons-RDB81.ncol ../Out2/Prodoric/mi-pearson-mrnet-133/graph-5-operons-RDB81.optim"," ")[[1]]
# opts <- parse_args(opt.parser, positional_arguments = TRUE, args=arg)

argv <- opts$args
if(length(argv) < 1) {
  print_help(opt.parser)
  quit(save = "no", status = 1)
}

if(opts$options$asp) {
  # first file should be in NCOL format
  # TODO: document, relax
  net   <- read.graph(argv[1], format="ncol", names=T, direct=T, weights="no")
  nodes <- sort(V(net)$name)
  score <- analize.graph(net, opts$options$index)
  ans   <- as.data.frame(score[nodes])
  net.names <- argv[1]
  m <- 2
  for(net.file in argv[2:length(argv)])  {
    gg <- readLines(net.file)
    gg <- gg[substr(gg, 1, 8)=="used_arc"]
    n <- 1
    for(gl in gg) {
      arcs <- lapply(strsplit(gl," "), strsplit, '[,"()]')[[1]]
      g <- graph.data.frame(data.frame(t(sapply(arcs, function(x) c(x[3], x[6], x[8])))))
      score <- analize.graph(g, opts$options$index)
      ans <- cbind(ans,score[nodes])
      net.names[m] <- paste(net.file, n, sep="-")
      n <- n+1
      m <- m+1
    }
  }
  colnames(ans) <- net.names
} else {
  net <- read.graph(argv[1], format="ncol", names=T, direct=T, weights="no")
  nodes <- sort(V(net)$name)
  ans <- matrix(NA, nrow=length(nodes), ncol=length(argv), dimnames=list(nodes,argv))

  # TODO: agregar caso de salida ASP
  for(graph.file in argv) {
    net <- read.graph(graph.file, format="ncol", names=T, direct=T, weights="no")
    score <- analize.graph(net, opts$options$index) 
    ans[, graph.file] <- score[nodes]
  }
}

column.mean <- function(mat) {
	ans <- colMeans(mat, na.rm=TRUE)
	return(t(ans))
}

row.mean <- function(mat) {
	ans <- colMeans(mat, na.rm=TRUE)
	return(ans)
}

row.summary <- function(mat) {
	# t(as.matrix(summary(mat)))
	data.frame(min=apply(mat, 2, min, na.rm=TRUE),
		   median=apply(mat, 2, median, na.rm=TRUE),
		   mean=colMeans(mat, na.rm=TRUE),
		   sd=apply(mat, 2, sd, na.rm=TRUE),
		   max=apply(mat, 2, max, na.rm=TRUE))
}

if(opts$options$summary) {
	func <-  match.fun(opts$options$funct)
	write.table(func(ans), file="", sep="\t", quote=F)
} else {
	write.table(ans, file="", sep="\t", quote=F)
}

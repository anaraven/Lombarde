#! /usr/bin/env Rscript --vanilla

library(optparse)
library(igraph)
library(parallel)

# Argument handling
option_list <- list(
  make_option(c("-o", "--out"), help = "Filename of output"),
  make_option(c("-n", "--dry"), action="store_true", default=FALSE,
              dest="dry", help="Don't calculate, just show input arguments"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Be more explicit on what is being done"),
  make_option(c("-w", "--wgt"), action="store_true", default=FALSE,
              help="Use weights as specified in the file. Don't recalculate weights."),
  make_option(c("-b", "--base"), action="store", default=10, type="double",
              help="Base to use in weight conversion"),
  make_option(c("-c", "--cores"), action="store", default=detectCores(), type="integer",
              help="Number of parallel process to run. Default: all cores.")
)

opt.parser <- OptionParser(option_list = option_list, 
  description="Usage: list-vshapes.R --out out.file coexp.file net.file extra.net")

opts <- parse_args(opt.parser, positional_arguments = TRUE)

argv <- opts$args
if(length(argv) < 2) { # TODO: better use a flag to determine when the extra.file is included
  print_help(opt.parser)
  quit(save = "no", status = 1)
}
coexp.file <- argv[1]
net.file   <- argv[2] 
extra.file <- if(length(argv)>2) argv[3] else argv[2]
out.file   <- opts$options$out

if(opts$options$verbose) {
  cat("coexp.file",coexp.file,"\n")
  cat("net.file",net.file,"\n")
  cat("extra.file",extra.file,"\n")
  cat("out.file",out.file,"\n")
#  quit(save = "no", status = 0)
}

path.extremes <- function(i, exps, v.obs) {
  unlist(lapply(names(exps[[i]]), function(r,v) list(c(r, v$from), c(r, v$to)), v.obs[i,]), recursive=F, use.names=F)
}


options(mc.cores=opts$options$cores)

g  <- read.graph(net.file, format="ncol", names=T, direct=T, weights="yes")
ge <- read.graph(extra.file, format="ncol", names=T, direct=T, weights="yes")
if(! opts$options$wgt) {
  message("changing weights")
  E(g)$orig   <- E(g)$weight
  E(g)$weight <- 10^(E(g)$weight)
  E(g)$out.w  <- 10^(E(ge)$weight)
} else {
  message("keeping weights")
  E(g)$out.w  <- E(ge)$weight
}
nom <- V(g)$name

v.obs <- read.table(coexp.file, as.is=T, col.names=c("from", "to", "weight"))
is.valid.obs <- v.obs$from %in% nom & v.obs$to %in% nom # & v.obs$from < v.obs$to
v.obs <- v.obs[is.valid.obs,]
N <- nrow(v.obs)
cat("n.obs",N,"\n")
if(opts$options$verbose) cat("is.valid.obs",table(is.valid.obs),"\n")

# `W` shows the cost of the cheapest path betwen vertices
W <- shortest.paths(g, mode="out")
# `exps` is the list of common predecessors for each pair in `v.obs`
exps <- mclapply(1:N, function(i) {v <- rowSums(W[,unlist(v.obs[i,1:2])]); which(v==min(v))})

expl.path <- unique(unlist(mclapply(1:N, path.extremes, exps, v.obs), recursive=F, use.names=F))

expl.path <- split(sapply(expl.path,`[`,2), sapply(expl.path,`[`,1))

expl.path <- mcmapply(function(src, targets) mapply(function(tgt) {
  get.all.shortest.paths(g, src, tgt, mode="out")$res
  }, targets, SIMPLIFY=FALSE), names(expl.path), expl.path, SIMPLIFY=FALSE)
  
expl.path <- mclapply(expl.path, lapply, lapply, function(l) {
  if(length(l)>1) as.numeric(E(g, path=l)) else numeric()} )

write.graph(subgraph.edges(g, unique(unlist(expl.path))), out.file, format="ncol")


vid <- 1
for(i in 1:N){
  for(r in names(exps[[i]])){
    cat("explanation",r,unlist(v.obs[i,1:2]),"\n")
    vv1 <- expl.path[[r]][[ v.obs[[i,1]] ]]
    vv2 <- expl.path[[r]][[ v.obs[[i,2]] ]]
    for(p1 in vv1) {
      for(p2 in vv2) {
        cat("vshape", vid, unlist(v.obs[i, 1:2]), i, "\n")
        edgelist <- c(p1,p2)
        weights <- get.edge.attribute(g, "out.w", edgelist)
        sides <- get.edges(g, edgelist)
        for(j in 1:length(edgelist)) {
          cat("arcInVshape", vid, nom[sides[j,1]], nom[sides[j,2]], weights[j], "\n")        
        }
        vid <- vid + 1
      }
    }
  }
}

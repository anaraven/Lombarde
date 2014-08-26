#! /usr/local/bin/Rscript --vanilla

### Functions ###
usage <- function() {
  message("Usage: edges.ncol output.file")
}

### Main ###

# Argument handling
args <- commandArgs(T)

if(length(args) != 2) {
  usage()
  quit(save = "no", status = 1)
}

graph.file      <- args[1]
node.names.file <- args[2] # "Input/operon_names.txt"
output.file     <- args[3]

valenteForeman <- function(g, what="both") {
  sp <- shortest.paths(g, mode="out", weight=NA)
  rd <- max(sp[is.finite(sp)])+1-sp
  rd[is.infinite(rd)] <- 0
  diag(rd) <- 0
  what <- match.arg(what, c("integration","radiality","both"))
  if(what=="both"){
    return( list(integration=colSums(rd)/(nrow(rd)-1), radiality=rowSums(rd)/(ncol(rd)-1)) )
  } else if(what=="radiality") {
    return( rowSums(rd)/(ncol(rd)-1) )
  } else {
    return( colSums(rd)/(ncol(rd)-1) )
  }
}

library(igraph)
net <- read.graph(graph.file, format="ncol", names=T, direct=T, weights="no")
score <- as.data.frame(valenteForeman(net))

tmp <- read.table(node.names.file, as.is=T)
m <- as.character(tmp$V2)
n <- as.character(tmp$V1)
names(n) <- tmp$V2
names(m) <- tmp$V1
rm(tmp)

missing <- m[! m %in% rownames(score)]
score <- rbind(score, data.frame(integration=rep(0,length(missing)), radiality=0, row.names=missing))
score$gene <- n[rownames(score)]
score$rank <- 1 + length(score$radiality) - rank(score$radiality, ties.method="max")
score[score$rank==max(score$rank), "rank"] <- NA

write.table(score, file=output.file, sep="\t", quote=F)
#mc <- read.delim("~/Desktop/IRISA/branches/PWM_MI/Ecoli/martinez-collado.txt", as.is=T)
#write.table(score[m[mc$gene],], file=output.file, sep="\t", quote=F)
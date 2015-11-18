#! /usr/local/bin/Rscript --vanilla

### Functions ###
usage <- function() {
  message("Usage: Rscript --vanilla evaluation-script.R association.file upstream.file can.bind.file n.levels out.file")
}

# Argument handling
args <- commandArgs(T)
if(length(args) != 5) {
  usage()
  quit(save = "no", status = 1)
}
filename      <- args[1] # "Expression/mi-pearson-aracne-133.txt"
upstream.file <- args[2] # "RegulonDB/Meme/upstreamLog.txt"
can.bind.file <- args[3] # "RegulonDB/Meme/canBind.txt"
n.levels      <- args[4]
out.file      <- args[5]

source("funcExplain.R")
source("load-gold-standard.R")
ops <-load.operons()
gold.oo <- gold.o2o(gold, ops)
v.obs <- read.obs(filename, names(ops), 10000)
v.obs <- contract.data.frame(v.obs, ops, max)

n <- as.numeric(n.levels)
adjust.w <- function(v, n) {
  w <- 1 + max(v) - v
  discretize.weight(w, n)  
}

# choose initial graph
pg <- read.delim(upstream.file, header=F, as.is=T, col.names=c("from", "to", "p", "q"))
gp <- read.delim(can.bind.file, header=F, as.is=T, col.names=c("from", "to", "e"))
pg$p <- adjust.w(pg$p, n)
pg$q <- adjust.w(pg$q, n)
gp$e <- adjust.w(gp$e, n)

links <- gene.to.gene.list(gp, pg)
oo  <- o2o(links, ops)

#oo$weight <- adjust.w(oo$w1,3)
oo$weight <- oo$w1
all.paths <- find.all.paths(oo, v.obs)
found <- count.found(all.paths, rownames(gold.oo), oo)
found.restr <- count.found.restricted(all.paths, rownames(gold.oo))
save.image(file=out.file)

pr <- function (table) {
  precision <- table$tp/(table$tp + table$fp)
  recall <- table$tp/(table$tp + table$fn)
  precision[is.nan(precision)] <- 1
  res <- data.frame(precision, recall)
  names(res) <- c("p", "r")
  res
}
print(pr(tail(table.found(found),1)))
print(pr(tail(table.found(found.restr),1)))

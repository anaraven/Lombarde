# setwd("~/Desktop/IRISA/branches/PWM_MI/Ecoli/")

source("funcExplain.R")
source("load-gold-standard.R")
ops <-load.operons()
gold.oo <- gold.o2o(gold, ops)
filename <- "Expression/mi-pearson-aracne-133.txt"
filename <- "Expression/mi-pearson-mrnet-133.txt"
v.obs <- read.obs(filename, names(ops), 10000)
v.obs <- contract.data.frame(v.obs, ops, max)

adjust.w <- function(v, n) {
  w <- 1 + max(v) - v
  discretize.weight(w, n)  
}

# choose initial graph
pg <- read.delim("RegulonDB/Meme/upstreamLog.txt", header=F, as.is=T, col.names=c("from", "to", "p", "q"))
gp <- read.delim("RegulonDB/Meme/canBind.txt", header=F, as.is=T, col.names=c("from", "to", "e"))
pg$p <- adjust.w(pg$p, 3)
pg$q <- adjust.w(pg$q, 3)
gp$e <- adjust.w(gp$e, 3)

links <- gene.to.gene.list(gp, pg)
oo  <- o2o(links, ops)

#oo$weight <- adjust.w(oo$w1,3)
oo$weight <- oo$w1
all.paths <- find.all.paths(oo, v.obs)
found <- count.found(all.paths, rownames(gold.oo), oo)
found.restr <- count.found.restricted(all.paths, rownames(gold.oo))

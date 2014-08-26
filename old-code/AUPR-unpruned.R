setwd("~/Desktop/IRISA/branches/PWM_MI/Ecoli/")

source("funcExplain.R")
source("load-gold-standard.R")
ops <-load.operons()
gold.oo <- gold.o2o(gold, ops)

adjust.w <- function(v) {
  1 + max(v) - v
}

# choose initial graph
pg <- read.delim("RegulonDB/Meme/upstreamLog.txt", header=F, as.is=T, col.names=c("from", "to", "p", "q"))
gp <- read.delim("RegulonDB/Meme/canBind.txt", header=F, as.is=T, col.names=c("from", "to", "e"))
pg$p <- adjust.w(pg$p)
pg$q <- adjust.w(pg$q)
gp$e <- adjust.w(gp$e)

links <- gene.to.gene.list(gp, pg)
or <- order(sapply(links,'[[',"w1"))
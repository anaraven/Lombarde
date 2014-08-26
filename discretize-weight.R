#! /usr/bin/Rscript --vanilla

### Functions ###
usage <- function() {
  message("Usage: Rscript --vanilla discretize-weight.R edge-list.net n1.levels n2.levels out.file")
}

# Argument handling
argv <- commandArgs(T)
if(length(argv) != 4) {
  usage()
  quit(save = "no", status = 1)
}

net.file   <- argv[1]
n1.levels  <- as.numeric(argv[2])
n2.levels  <- as.numeric(argv[3])
out.file   <- argv[4]

discretize.weight <- function(w, n) {
  thr <- quantile(w, (0:n)/n, na.rm=T)
  return(cut(w, breaks=thr, include.lowest=T, labels=F)-1)
}

net <- read.table(net.file, as.is=T)
w1 <- discretize.weight(net$V3, n1.levels)
w2 <- discretize.weight(net$V4, n2.levels)
out <- data.frame(net[,1:2], w1+w2)
write.table(out, out.file, quote=F, sep ="\t", row.names=F, col.names=F)
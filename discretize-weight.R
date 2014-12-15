#! /usr/bin/Rscript --vanilla

library(optparse)

# Argument handling
option_list <- list(
  make_option(c("-i", "--in"), help = "Filename of input graph (\"ncol\" format)", dest="input"),
  make_option(c("-o", "--out"), help = "Filename of output graph (\"ncol\" format)"),
  make_option(c("-n", "--n1"), action="store", default=3, type="integer",
      help="number of discrete levels in the first input weight (BLAST E-value)"),
  make_option(c("-N", "--n2"), action="store", default=3, type="integer",
      help="number of discrete levels in the second input weight (MEME/FIMO p-value)"),
  make_option(c("-b", "--base"), action="store", default=10, type="double",
      help="Base to use in weight exponentiation. Use a negative value to keep the original discrete values")
)

opt.parser <- OptionParser(option_list = option_list, 
  description="Usage: Rscript --vanilla discretize-weight.R [options]")

opts <- parse_args(opt.parser)

discretize.weight <- function(w, n) {
  thr <- quantile(w, (0:n)/n, na.rm=T)
  return(cut(w, breaks=thr, include.lowest=T, labels=F)-1)
}

net <- read.table(opts$input, as.is=T)
w1 <- discretize.weight(net$V3, opts$n1)
w2 <- discretize.weight(net$V4, opts$n2)
w3 <- w1+w2
if(opts$base>0) {
  w3 <- opts$base ^ w3
}
out <- data.frame(net[,1:2], w3)
write.table(out, opts$out, quote=F, sep ="\t", row.names=F, col.names=F)

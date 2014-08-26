#! /usr/bin/Rscript --vanilla


library(optparse)

# Argument handling
option_list <- list(
  make_option(c("-o", "--out"), help = "Filename of output"),
  make_option(c("-p", "--plevels"), action="store", default=3, type="integer", help="Levels on P-value discretization"),
  make_option(c("-e", "--elevels"), action="store", default=3, type="integer", help="Levels on E-value discretization"),
  make_option(c("-m", "--mlevels"), action="store", default=3, type="integer", help="Levels on Mutual Information discretization")
)

opt.parser <- OptionParser(option_list = option_list, 
                           description="Discretize graph weights from E-valeus, p-values and Mutual Information")

opts <- parse_args(opt.parser, positional_arguments = TRUE)

argv <- opts$args
if(length(argv) < 2) {
  print_help(opt.parser)
  quit(save = "no", status = 1)
}

net.file   <- argv[1]
mim.file   <- argv[2]
n1.levels  <- opts$options$elevels
n2.levels  <- opts$options$plevels
n3.levels  <- opts$options$mlevels
out.file   <- opts$options$out

if(is.null(out.file)) {
  print_help(opt.parser)
  quit(save = "no", status = 1)
}

discretize.weight <- function(w, n) {
  thr <- quantile(w, (0:n)/n, na.rm=T)
  return(cut(w, breaks=thr, include.lowest=T, labels=F)-1)
}

net <- read.table(net.file, as.is=T)
load(mim.file)

colnames(mim) <- sapply(strsplit(colnames(mim),"_"), `[` , 2)
rownames(mim) <- sapply(strsplit(rownames(mim),"_"), `[` , 2)

X <- net[,1]
Y <- net[,2]
valid.X <- X %in% rownames(mim)
valid.Y <- Y %in% rownames(mim)

c <- sapply(1:length(X),
  function(i) ifelse(valid.X[i] & valid.Y[i], -mim[X[i],Y[i]], 0))

w1 <- discretize.weight(net$V3, n1.levels)
w2 <- discretize.weight(net$V4, n2.levels)
w3 <- discretize.weight(c, n3.levels)
out <- data.frame(net[,1:2], w1+w2+w3)
write.table(out, out.file, quote=F, sep ="\t", row.names=F, col.names=F)
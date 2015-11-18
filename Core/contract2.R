#!/usr/bin/Rscript --vanilla

library(optparse)

# Argument handling
option_list <- list()
opt.parser <- OptionParser(option_list = option_list,
  description="Usage: contract2.R operon_list.txt graph.ncol output.ncol")
opts <- parse_args(opt.parser, positional_arguments = TRUE)

argv <- opts$args
if(length(argv) < 3) {
  print_help(opt.parser)
  quit(save = "no", status = 1)
}

net.file <- argv[1]
ops.file <- argv[2]
out.file <- argv[3]

op <- read.table(ops.file, row.names=1, as.is=T)
o <- op[,1]
names(o) <- rownames(op)

net <- read.table(net.file, as.is=T)
net$V1 <- o[net$V1]
net$V2 <- o[net$V2]
net <- net[!is.na(net$V1) & !is.na(net$V2) & net$V1!=net$V2,]
idx <- paste(net$V1, net$V2)
z <- lapply(1:ncol(net), function(i) tapply(net[,i], idx, min))
net2 <- do.call(cbind.data.frame,z)
write.table(net2, out.file, quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)

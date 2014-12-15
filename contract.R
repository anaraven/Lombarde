#!/usr/bin/Rscript --vanilla

library(optparse)
library(igraph)

# Argument handling
option_list <- list()
opt.parser <- OptionParser(option_list = option_list,
  description="Usage: operon_list.txt graph.ncol")
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

net <- read.graph(net.file,  format="ncol", names=T, direct=T, weights="yes")
# we store the gold-standard vertices names in a vector for easy access 
net.v.names <- V(net)$name
#cat(vcount(net), "vertex", ecount(net), "edges\n")
print(net)

net <- delete.vertices(net, net.v.names[!net.v.names %in% rownames(op)])
net.v.names <- V(net)$name
net.v.id <- 1:length(net.v.names)
names(net.v.id) <- net.v.names
#cat(vcount(net), "vertex", ecount(net), "edges\n")
print(net)

onn <- o[net.v.names]
u.name <- unique(onn)
u.id <- 1:length(u.name)
names(u.id) <- u.name

net2 <- contract.vertices(net, u.id[onn], vertex.attr.comb="first")
V(net2)$name <- u.name
#cat(vcount(net2), "vertex", ecount(net2), "edges\n")
print(net2)

net2 <- simplify(net2, edge.attr.comb="min")
#cat(vcount(net2), "vertex", ecount(net2), "edges\n")
print(net2)

write.graph(net2, out.file, format="ncol")

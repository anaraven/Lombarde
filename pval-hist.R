op <- read.table("Input/RDB8.1/operon_list.clean.txt", row.names=1, as.is=T)
o <- op[,1]
names(o) <- rownames(op)
rm(op)

library(igraph)
q <- read.table("Input/Genes/Prodoric/MEME/g0.2val", col.names=c("from","to","blast","meme"))
net <- graph.data.frame(q, directed=T)
net.v.names <- V(net)$name
net <- delete.vertices(net, net.v.names[!net.v.names %in% names(o)])
net.v.id <- 1:length(net.v.names)
names(net.v.id) <- net.v.names
onn <- o[net.v.names]
onn <- onn[!is.na(onn)]
u.name <- unique(onn)
u.id <- 1:length(u.name)
names(u.id) <- u.name

net2 <- contract.vertices(net, u.id[onn], vertex.attr.comb="first")
V(net2)$name <- u.name
net2 <- simplify(net2, edge.attr.comb="min")
rm(net.v.id, net.v.names, onn, u.id, u.name, o, net, q)

all <- get.data.frame(net2)
rownames(all) <- paste(all$from,all$to)

gs <-read.table("Input/RDB8.1/gold-std.ncol", as.is=T, col.names=c("from","to","zero"))
rownames(gs) <- paste(gs$from, gs$to)

all$in.gs <- rownames(all) %in% rownames(gs)

library(ggplot2)

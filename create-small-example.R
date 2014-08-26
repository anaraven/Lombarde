library(igraph)

ref.graph <- read.table("DrawFig3/all2.graph", as.is=T)
ref.graph <- graph.data.frame(ref.graph[, 1:2])
ref.nodes <- V(ref.graph)$name

g0 <- read.graph("Input/RDB8.1/Prodoric/MEME/g0-9.ncol", format="ncol", names=T, direct=T, weights="yes")
g1 <- induced.subgraph(g0, V(g0)[ref.nodes])
# g2 <- graph.difference(g1, ref.graph)
write.graph(g1, 'Input/RDB8.1/Prodoric/MEME/small.ncol', format="ncol")

g0 <- read.graph("Input/RDB8.1/M3D/pearson/mrnet/100K.ncol", format="ncol", names=T, direct=F, weights="yes")
g1 <- induced.subgraph(g0, V(g0)[ref.nodes])
write.graph(g1, 'Input/RDB8.1/M3D/pearson/mrnet/small.ncol', format="ncol")


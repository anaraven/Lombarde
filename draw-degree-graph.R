library(igraph)
gs <- read.graph("Input/RDB8.1/gold-std.ncol", format="ncol", names=T, direct=T, weights="no")
g0 <- read.graph("Input/RDB8.1/Prodoric/MEME/g0-3.ncol", format="ncol", names=T, direct=T, weights="yes")
gl <- read.graph("Out1/RDB8.1/M3D/pearson/mrnet/100K/Prodoric/MEME/gl-3.ncol", format="ncol", names=T, direct=T, weights="yes")


d.gs <- degree.distribution(gs)[2:16]
d.g0 <- degree.distribution(g0)[2:16]
d.gl <- degree.distribution(gl)[2:16]

d.gs[16] <- sum(degree.distribution(gs)[-(1:16)])
d.g0[16] <- sum(degree.distribution(g0)[-(1:16)])
d.gl[16] <- sum(degree.distribution(gl)[-(1:16)])

plot(1:15, d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(1,15,7))
points(d.g0, type="o", pch=2, lty=3)
points(d.gl, type="o", pch=16)

library(ggplot2)

l.gl <- "LOMBARDE output network"
l.g0 <- "Putative TRN"
l.gs <- "Network of Validated Arcs"

a <- data.frame(degree=c(d.gl, d.g0, d.gs), graph=factor(c(rep(l.gl,15), rep(l.g0,15), rep(l.gs,15)), levels=c(l.gl,l.g0,l.gs)))
pdf(width=6, height=4)
ggplot(data=a, aes(x=rep(1:15,3), y=degree, colour=graph, shape=graph)) + 
  theme_bw() + 
  scale_x_continuous(breaks=1:15) +
  xlab("Degree") + 
  ylab("Proportion") + 
  geom_line(size=1.2) + 
  geom_point(size=3, fill="white") + 
  theme(legend.position=c(1,1), legend.justification=c(1,1)) +
  scale_shape_manual(values=c(21,22,23), name="Network") +
  scale_color_hue(name="Network")
dev.off()

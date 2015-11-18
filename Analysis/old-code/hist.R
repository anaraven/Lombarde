}
TF.id <- TF.id[order(dg2, decreasing=T)]
dg2 <- dg2[order(dg2, decreasing=T)]
r <- pmax(dg2, opts$min)
r <- r/sum(r)
theta <- rep(0,N1)
theta[2:N1] <- cumsum(r[1:(N1-1)]+r[2:N1])*pi
big <- list(x=cos(theta), y=sin(theta))
XY <- matrix(NA, nrow=vcount(g), ncol=2, dimnames=list(V(g)$name, c("x","y")))
XY[TF.id, "x"] <- big$x
XY[TF.id, "y"] <- big$y
color <- rep("SkyBlue2",vcount(g))
names(color) <- V(g)$name
color[TF.id] <- "red"
for(i in 1:N1) {
neigh <- names(dg3)[dg3==TF.id[i]]
N2 <- length(neigh)
if(N2>0) {
small <- circle(N2, r[i]*3,big$x[i], big$y[i])
XY[neigh, "x"] <- small$x
XY[neigh, "y"] <- small$y
}
g1<-delete.edges(g, edge.color=="red")
get.edges(g,edge.color=="red")
get.edges(g)
E(g)[edge.color=="red"]
g1<-delete.edges(g, E(g)[edge.color=="red"])
g2<-delete.edges(g, E(g)[edge.color!="red"])
plot(g1, layout=XY, vertex.size=1.5, vertex.label=NA, edge.arrow.mode=0, vertex.color=color,
edge.curved=opts$curved, edge.color="darkgrey")
plot(g2, layout=XY, vertex.size=1.5, vertex.label=NA, edge.arrow.mode=0, vertex.color=color,
edge.curved=opts$curved, edge.color="red")
plot(g1, layout=XY, vertex.size=1.5, vertex.label=NA, edge.arrow.mode=0, vertex.color=color,
edge.curved=opts$curved, edge.color="darkgrey")
plot(g2, layout=XY, vertex.size=1.5, vertex.label=NA, edge.arrow.mode=0, vertex.color=color,
edge.curved=opts$curved, edge.color="red",add=T)
library(optparse)
library(igraph)
?plot.igraph
source('~/.active-rstudio-document')
dev.cur()
source('~/.active-rstudio-document')
library(igraph)
ops <- read.table("Input/RDB8.1/operon_names.txt")
View(ops)
X <- sample(ops, size, replace = TRUE)
ops <- ops$V2
size <- 100000
X <- sample(ops, size, replace = TRUE)
Y <- sample(ops, size, replace = TRUE)
table(X == Y)
table(X < Y)
ops <- read.table("Input/RDB8.1/operon_names.txt", as.is=T)
ops <- ops$V2
size <- 100000
X <- sample(ops, size, replace = TRUE)
Y <- sample(ops, size, replace = TRUE)
table(X < Y)
erdos.renyi.game
g <- erdos.renyi.game(lenght(ops), size, type="gnm")
g <- erdos.renyi.game(length(ops), size, type="gnm")
get.edgelist(g) [1:10]
get.edges(g) [1:10]
get.edges(g,E(g)[1:10])
edges <- get.edges(g,E(g)[1:10])
ops[edges[,1]]
ops[edges[,2]]
ans <- data.frame(X=ops[edges[,1]], Y=ops[edges[,2]])
ans
runif(10)
write.table(ans, file="Input/Random/erdos-RDB8.1-100K.ncol",
sep="\t", quote=F, row.names=F, col.names=F)
size <- 100000
g <- erdos.renyi.game(length(ops), size, type="gnm")
edges <- get.edges(g,E(g))
ans <- data.frame(X=ops[edges[,1]], Y=ops[edges[,2]], Z=runif(size))
write.table(ans, file="Input/Random/erdos-RDB8.1-100K.ncol",
sep="\t", quote=F, row.names=F, col.names=F)
source('~/Desktop/RLombarde/randomInfluence.R')
a <- read.table("Input/Genes/Prodoric/MEME/g0.2val",as.is=T)
load("Input/Genes/M3D/pearson/matrix.Rdata")
colnames(mim) <- sapply(strsplit(colnames(mim),"_"), `[` , 2)
rownames(mim) <- sapply(strsplit(rownames(mim),"_"), `[` , 2)
nom <- rownames(mim)
X <- a[,1]
Y <- a[,2]
valid.X <- X %in% nom
valid.Y <- Y %in% nom
c <- sapply(1:length(X),
function(i) ifelse(valid.X[i] & valid.Y[i], mim[X[i],Y[i]], NA))
a$c <- c
names(c) <- paste(X,Y)
g <- read.table("Input/Genes/gold-std.ncol",as.is=T)
X <- g[,1]
Y <- g[,2]
c[paste(X,Y)]
plot(c[paste(X,Y)])
net <- read.table("Input/Genes/Prodoric/MEME/g0.2val",as.is=T)
load("Input/Genes/M3D/pearson/matrix.Rdata")
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
n1.levels  <- opts$options$elevels
n2.levels  <- opts$options$plevels
n3.levels  <- opts$options$mlevels
discretize.weight <- function(w, n) {
thr <- quantile(w, (0:n)/n, na.rm=T)
return(cut(w, breaks=thr, include.lowest=T, labels=F)-1)
}
net <- read.table("Input/Genes/Prodoric/MEME/g0.2val",as.is=T)
load("Input/Genes/M3D/pearson/matrix.Rdata")
colnames(mim) <- sapply(strsplit(colnames(mim),"_"), `[` , 2)
rownames(mim) <- sapply(strsplit(rownames(mim),"_"), `[` , 2)
X <- net[,1]
Y <- net[,2]
valid.X <- X %in% rownames(mim)
valid.Y <- Y %in% rownames(mim)
c <- sapply(1:length(X),
function(i) ifelse(valid.X[i] & valid.Y[i], -mim[X[i],Y[i]], NA))
w1 <- discretize.weight(net$V3, n1.levels)
w2 <- discretize.weight(net$V4, n2.levels)
w3 <- discretize.weight(c, n3.levels)
table(w1,w2)
table(w1,w3)
table(w2,w3)
names(w1) <- paste(X,Y)
names(w2) <- paste(X,Y)
names(w3) <- paste(X,Y)
g <- read.table("Input/Genes/gold-std.ncol",as.is=T)
paste(g[,1],g[,2])
w1[paste(g[,1],g[,2])]
gg <- paste(g[,1],g[,2])
w1[gg] + w2[gg] + w3[gg]
table(w1[gg] + w2[gg] + w3[gg])
table(w1 + w2 + w3)
plot(table(w1 + w2 + w3),table(w1[gg] + w2[gg] + w3[gg]))
barplot(table(w1 + w2 + w3),table(w1[gg] + w2[gg] + w3[gg]))
barplot(table(w1 + w2 + w3))
barplot(table(w1[gg] + w2[gg] + w3[gg]),add=T)
table(w1 + w2)
table(w1[gg] + w2[gg])
gs <- read.graph("Input/RDB8.1/gold-std.ncol", format="ncol", names=T, direct=T, weights="no")
library(igraph)
gs <- read.graph("Input/RDB8.1/gold-std.ncol", format="ncol", names=T, direct=T, weights="no")
?cluster.distribution
transitivity(gs)
transitivity(gs, type="barrat")
g0 <- read.graph("Input/RDB8.1/Prodoric/MEME/g0-3.ncol", format="ncol", names=T, direct=T, weights="yes")
transitivity(g0)
transitivity(gs)
gl <- read.graph("Out1/RDB8.1/M3D/pearson/mrnet/100K/Prodoric/MEME/gl-3.ncol", format="ncol", names=T, direct=T, weights="yes")
transitivity(gl)
triad.census(gs)
triad.census(g0)
triad.census(gl)
triad.census(gs)[2:16]
?sum
triad.census(gl)/sum(triad.census(gl), na.rm=T)
plot(triad.census(gl)/sum(triad.census(gl), na.rm=T))
lines(triad.census(g0)/sum(triad.census(g0), na.rm=T))
lines(triad.census(gs)/sum(triad.census(gs), na.rm=T))
lines(triad.census(gs)[3:16]/sum(triad.census(gs), na.rm=T))
plot(triad.census(gs)[3:16]/sum(triad.census(gs), na.rm=T))
lines(triad.census(g0)[3:16]/sum(triad.census(g0), na.rm=T))
lines(triad.census(gl)[3:16]/sum(triad.census(gl), na.rm=T))
degree.distribution(g0)
plot(degree.distribution(g0))
lines(degree.distribution(gl))
plot(degree.distribution(gl)[1:50])
lines(degree.distribution(g0)[1:50])
lines(degree.distribution(gs)[1:50])
plot(degree.distribution(gs)[1:15])
plot(degree.distribution(g0)[1:15],pch=16)
plot(degree.distribution(gs)[1:15])
points(degree.distribution(g0)[1:15],pch=16)
points(degree.distribution(gl)[1:15],pch=19)
plot(degree.distribution(g0)[1:15],pch=17)
plot(degree.distribution(g0)[1:15],pch=1)
plot(degree.distribution(g0)[1:15],pch=2)
plot(degree.distribution(g0)[1:15],pch=3)
plot(degree.distribution(gs)[1:15])
points(degree.distribution(g0)[1:15],pch=2)
points(degree.distribution(gl)[1:15],pch=3)
plot(degree.distribution(gs)[1:15], type="o")
points(degree.distribution(g0)[1:15], type="o", pch=2)
points(degree.distribution(gl)[1:15], type="o", pch=3)
plot(degree.distribution(gs)[1:15], type="o", xlab="Degree", ylab="Density")
points(degree.distribution(g0)[1:15], type="o", pch=2)
points(degree.distribution(gl)[1:15], type="o", pch=3)
plot(degree.distribution(gs)[1:15], type="o", xlab="Degree", ylab="Density", lty=2)
plot(degree.distribution(gs)[1:15], type="o", xlab="Degree", ylab="Density", lty=3)
source('~/.active-rstudio-document')
?plot.default
plot(d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(2,15,2))
plot(d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(2,15,5))
plot(d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(2,15,14))
plot(d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(2,15,15))
plot(d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(2,15,13))
source('~/.active-rstudio-document')
plot(2:15,d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(2,15,13))
plot(2:16,d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(2,15,13))
source('~/.active-rstudio-document')
plot(2:15,d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(2,15,13))
d.gs
plot(2:16,d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(2,15,13))
points(d.g0, type="o", pch=2, lty=3)
plot(2:16,d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(2,15,13))
points(2:16,d.g0, type="o", pch=2, lty=3)
points(2:16,d.gl, type="o", pch=16)
d.gs <- degree.distribution(gs)[2:15]
d.g0 <- degree.distribution(g0)[2:15]
d.gl <- degree.distribution(gl)[2:15]
d.gs[15] <- sum(degree.distribution(gs)[-(1:15)])
d.g0[15] <- sum(degree.distribution(g0)[-(1:15)])
d.gl[15] <- sum(degree.distribution(gl)[-(1:15)])
plot(1:15, d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(2,15,13))
plot(1:15, d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(1,15,14))
plot(1:15, d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(1,14,13))
source('~/.active-rstudio-document')
plot(1:16, d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(1,16,5))
plot(1:16, d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(1,16,4))
plot(1:16, d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(1,15,3))
plot(1:16, d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(1,16,3))
plot(1:16, d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(1,15,8))
plot(1:16, d.gs, type="o", xlab="Degree", ylab="Density", lty=2, xaxp=c(1,15,7))
points(d.g0, type="o", pch=2, lty=3)
source('~/Desktop/RLombarde/draw-degree-graph.R')
a <- read.table("noinduced-summary.txt", header=T)
a <- read.table("noinduced-summary.txt", header=T, sep="\t")
View(a)
?with
summary(a)
subset <- a$estimator=="pearson" & a$BS=="MEME"
table(subset)
a[subset,]
plot(edges, data=a[subset,])
plot(Edges, data=a[subset,])
with(a[subset,], plot(Edges))
with(a[subset,], plot(Edges~MI+Levels))
with(a[subset,], plot(Edges~MI))
a[a$estimator=="Initial",]
subset <- (a$estimator=="pearson" | a$estimator=="Initial") & a$BS=="MEME"
with(a[subset,], plot(Edges~MI))
meme.pearson <- a$estimator=="pearson" & a$BS=="MEME"
gs0.pearson <- a$estimator=="pearson" & a$BS=="MEME-GS-0"
meme.pearson.initial <- a[(a$estimator=="pearson" | a$estimator=="Initial") & a$BS=="MEME", ]
with(meme.pearson.initial,plot(Eedges))
with(meme.pearson.initial,plot(Edges))
with(meme.pearson.initial,plot(Edges~Levels))
with(meme.pearson.initial,plot(Edges~MI))
with(meme.pearson,plot(Edges~MI))
meme.pearson <- a[a$estimator=="pearson" & a$BS=="MEME", ]
gs0.pearson  <- a[a$estimator=="pearson" & a$BS=="MEME-GS-0", ]
with(meme.pearson,plot(Edges~MI))
with(gs0.pearson,plot(Edges~MI))
attach(a)
?plot
plot(Edges~Levels)
plot(Edges~Levels, subset=(BS=="MEME"))
plot(Edges~Levels, type="o", subset=(BS=="MEME"))
plot(Edges~Levels, type="o", subset=(BS=="MEME" & estimator="pearson"))
plot(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson"))
plot(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"))
plot(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"), xaxp=c(3,9,4))
plot(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"), xaxp=c(3,9,5))
plot(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"), xaxp=c(3,9,3))
plot(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"), xaxp=c(3,9,3), ylim=c(0,6000))
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="aracne"),lty=2,pch=2)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="c3"),lty=3,pch=3)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="clr"),lty=4,pch=4)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"), col="blue")
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="aracne"), col="blue", lty=2, pch=2)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="c3"), col="blue", lty=3, pch=3)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="clr"), col="blue", lty=4, pch=4)
plot(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"), xaxp=c(3,9,3), ylim=c(0,6500), main="Ab initio")
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="aracne"),lty=2,pch=2)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="c3"),lty=3,pch=3)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="clr"),lty=4,pch=4)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"), col="blue")
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="aracne"), col="blue", lty=2, pch=2)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="c3"), col="blue", lty=3, pch=3)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="clr"), col="blue", lty=4, pch=4)
plot(Edges~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="mrnet"), xaxp=c(3,9,3), ylim=c(0,6500), main="Ab initio")
points(Edges~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="aracne"),lty=2,pch=2)
points(Edges~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="c3"),lty=3,pch=3)
points(Edges~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="clr"),lty=4,pch=4)
points(In.gold~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="mrnet"), col="blue")
points(In.gold~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="aracne"), col="blue", lty=2, pch=2)
points(In.gold~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="c3"), col="blue", lty=3, pch=3)
points(In.gold~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="clr"), col="blue", lty=4, pch=4)
plot(Edges~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="mrnet"), xaxp=c(3,9,3), ylim=c(0,6500), main="Extended")
points(Edges~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="aracne"),lty=2,pch=2)
points(Edges~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="c3"),lty=3,pch=3)
points(Edges~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="clr"),lty=4,pch=4)
points(In.gold~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="mrnet"), col="blue")
points(In.gold~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="aracne"), col="blue", lty=2, pch=2)
points(In.gold~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="c3"), col="blue", lty=3, pch=3)
points(In.gold~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="clr"), col="blue", lty=4, pch=4)
plot(Edges~Levels, type="o", subset=(BS=="MEME" & (estimator=="pearson" | estimator="Initial") & MI=="mrnet"), xaxp=c(3,9,3), ylim=c(0,6500), main="Ab initio")
plot(Edges~Levels, type="o", subset=(BS=="MEME" & (estimator=="pearson" | estimator=="Initial") & MI=="mrnet"), xaxp=c(3,9,3), ylim=c(0,6500), main="Ab initio")
a[(BS=="MEME" & (estimator=="pearson" | estimator=="Initial") & MI=="mrnet"),]
a[(BS=="MEME" & (estimator=="pearson" | estimator=="Initial") & (MI=="mrnet" | MI=="Initial"),]
a[(BS=="MEME" & (estimator=="pearson" | estimator=="Initial") & (MI=="mrnet" | MI=="Initial")),]
plot(Edges~Levels, type="o", subset=(BS=="MEME" & (estimator=="pearson" | estimator=="Initial") & MI=="Initial"), xaxp=c(3,9,3), main="Ab initio")
plot(Edges~Levels, type="o", subset=(BS=="MEME" & (estimator=="pearson" | estimator=="Initial") & MI=="Initial"), xaxp=c(3,9,3), ylim=c(0,26500), main="Ab initio")
plot(Edges~Levels, type="o", subset=(BS=="MEME" & (estimator=="pearson" | estimator=="Initial") & MI=="Initial"), xaxp=c(3,9,3), ylim=c(0,26500), main="Ab initio", log="y")
plot(Edges~Levels, type="o", subset=(BS=="MEME" & (estimator=="pearson" | estimator=="Initial") & MI=="Initial"), xaxp=c(3,9,3), ylim=c(1,26500), main="Ab initio", log="y")
plot(Edges~Levels, type="o", subset=(BS=="MEME" & (estimator=="pearson" | estimator=="Initial") & MI=="Initial"), xaxp=c(3,9,3), ylim=c(10,26500), main="Ab initio", log="y")
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"),lty=1,pch=1)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="aracne"),lty=2,pch=2)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="c3"),lty=3,pch=3)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="clr"),lty=4,pch=4)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"), col="blue")
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="aracne"), col="blue", lty=2, pch=2)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="c3"), col="blue", lty=3, pch=3)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="clr"), col="blue", lty=4, pch=4)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="Initial" & MI=="Initial"), col="blue")
plot(Edges~Levels, type="o", subset=(BS=="MEME" &  estimator=="Initial" & MI=="Initial"), xaxp=c(3,9,3), ylim=c(10,26500), main="Ab initio", log="y", pch=16)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"),lty=1,pch=1)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="aracne"),lty=2,pch=2)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="c3"),lty=3,pch=3)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="clr"),lty=4,pch=4)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="Initial" & MI=="Initial"), col="blue", pch=16)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"), col="blue")
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="aracne"), col="blue", lty=2, pch=2)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="c3"), col="blue", lty=3, pch=3)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="clr"), col="blue", lty=4, pch=4)
plot(Edges~Levels, type="o", subset=(BS=="MEME" &  estimator=="Initial" & MI=="Initial"), xaxp=c(3,9,3), ylim=c(100,26500), main="Ab initio", log="y", pch=16)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"),lty=1,pch=1)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="aracne"),lty=2,pch=2)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="c3"),lty=3,pch=3)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="clr"),lty=4,pch=4)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="Initial" & MI=="Initial"), col="blue", pch=16)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"), col="blue")
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="aracne"), col="blue", lty=2, pch=2)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="c3"), col="blue", lty=3, pch=3)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="clr"), col="blue", lty=4, pch=4)
plot(Edges~MI, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"), xaxp=c(3,9,3), ylim=c(0,6500), main="Ab initio")
plot(Edges~MI, type="o", subset=(BS=="MEME" & estimator=="pearson"), xaxp=c(3,9,3), ylim=c(0,6500), main="Ab initio")
plot(Edges~MI, type="o", subset=(BS=="MEME" & estimator=="pearson"), ylim=c(2000,6500), main="Ab initio", xlab="Co-expr. detec. meth.")
plot(Edges~factor(MI), type="o", subset=(BS=="MEME" & estimator=="pearson"), ylim=c(2000,6500), main="Ab initio", xlab="Co-expr. detec. meth.")
plot(Edges~factor(MI, levels=c("aracne","c3","clr")), type="o", subset=(BS=="MEME" & estimator=="pearson"), ylim=c(2000,6500), main="Ab initio", xlab="Co-expr. detec. meth.")
plot(Edges~factor(MI, levels=c("aracne","c3","clr","mrnet")), type="o", subset=(BS=="MEME" & estimator=="pearson"), ylim=c(2000,6500), main="Ab initio", xlab="Co-expr. detec. meth.")
plot(Edges~factor(MI, levels=c("aracne","c3","clr","mrnet")), type="o", subset=(BS=="MEME" & estimator=="pearson"), ylim=c(2000,6500), main="Ab initio",
xlab="Co-expr. detect. meth.", ylab="Num. Arcs")
library(ggplot2)
install.packages("ggplot2", lib="/Library/Frameworks/R.framework/Versions/2.15/Resources/library")
install.packages("extrafont", lib="/Library/Frameworks/R.framework/Versions/2.15/Resources/library")
library(ggplot2)
library(extrafont)
loadfonts()
theme_xkcd <- theme(
panel.background = element_rect(fill="white"),
axis.ticks = element_line(colour=NA),
panel.grid = element_line(colour="white"),
axis.text.y = element_text(colour=NA),
axis.text.x = element_text(colour="black"),
text = element_text(size=16, family="Humor Sans")
)
a <- read.table("noinduced-summary.txt", header=T, sep="\t")
p <- ggplot(data=a,aes(x=Levels, y=Edges))
p
p <- ggplot(data=a,aes(x=Levels, y=Edges)) + geom_smooth(aes(y=Nodes), colour="gold", size=1, position="jitter", fill=NA)
p
ggplot(data=a,aes(x=Levels, y=Edges)) + geom_smooth(aes(y=Vertices), colour="gold", size=1, position="jitter", fill=NA)
ggplot(data=a,aes(x=Levels, y=Edges)) + geom_smooth(aes(y=Vertices), colour="gold", size=1, position="jitter", fill=NA)+geom_smooth(colour="white", size=3, position="jitter", fill=NA)
ggplot(data=a,aes(x=Levels, y=Edges)) + geom_smooth(aes(y=Vertices), colour="gold", size=1, position="jitter", fill=NA)+geom_smooth(colour="white", size=3, position="jitter", fill=NA)+geom_smooth(colour="dark blue", size=1, position="jitter", fill=NA)
ggplot(data=a,aes(x=Levels, y=Edges)) + geom_smooth(aes(y=Vertices), colour="gold", size=1, position="jitter", fill=NA)+geom_smooth(colour="white", size=3, position="jitter", fill=NA)+geom_smooth(colour="dark blue", size=1, position="jitter", fill=NA)+ geom_text(data=a[10, ], family="Humor Sans", aes(x=Date), colour="gold", y=20, label="Searches for clegg")
ggplot(data=a,aes(x=Levels, y=Edges)) + geom_smooth(aes(y=Vertices), colour="gold", size=1, position="jitter", fill=NA)+geom_smooth(colour="white", size=3, position="jitter", fill=NA)+geom_smooth(colour="dark blue", size=1, position="jitter", fill=NA)+ geom_text(data=a[10, ], family="Humor Sans", aes(x=Levels), colour="gold", y=20, label="Searches for clegg")
ggplot(data=a,aes(x=Levels, y=Edges)) + geom_smooth(aes(y=Vertices), colour="gold", size=1, position="jitter", fill=NA, method="loess")+geom_smooth(colour="white", size=3, position="jitter", fill=NA, method="loess")+geom_smooth(colour="dark blue", size=1, position="jitter", fill=NA, method="loess")+ geom_text(data=a[1, ], family="Humor Sans", aes(x=Levels), colour="gold", y=20, label="Searches for clegg")
warnings()
ggplot(data=a,aes(x=Levels, y=Edges)) + geom_smooth(aes(y=Vertices), colour="gold", size=1, position="jitter", fill=NA, method="loess")+geom_smooth(colour="white", size=3, position="jitter", fill=NA, method="loess")+geom_smooth(colour="dark blue", size=1, position="jitter", fill=NA, method="loess")+ geom_text(data=a[1, ], family="Comic Sans", aes(x=Levels), colour="gold", y=20, label="Searches for clegg")
warnings()
?loadfonts
fints()
fonts()
ggplot(data=a,aes(x=Levels, y=Edges)) + geom_smooth(aes(y=Vertices), colour="gold", size=1, position="jitter", fill=NA, method="loess")+geom_smooth(colour="white", size=3, position="jitter", fill=NA, method="loess")+geom_smooth(colour="dark blue", size=1, position="jitter", fill=NA, method="loess")+ geom_text(data=a[1, ], aes(x=Levels), colour="gold", y=20, label="Searches for clegg")
warnings()
ggplot(data=a,aes(x=Levels, y=Edges)) + geom_smooth(aes(y=Vertices), colour="gold", size=1, position="jitter", fill=NA, method="loess")+geom_smooth(colour="white", size=3, position="jitter", fill=NA, method="loess")+geom_smooth(colour="dark blue", size=1, position="jitter", fill=NA, method="loess")+ geom_text(data=a[2, ], aes(x=Levels), colour="gold", y=20, label="Searches for clegg")+geom_line(aes(y=xaxis), position = position_jitter(h = 0.1), colour="black")
ggplot(data=a,aes(x=Levels, y=Edges)) + geom_smooth(aes(y=Vertices), colour="gold", size=1, position="jitter", fill=NA, method="loess")+geom_smooth(colour="white", size=3, position="jitter", fill=NA, method="loess")+geom_smooth(colour="dark blue", size=1, position="jitter", fill=NA, method="loess")+ geom_text(data=a[2, ], aes(x=Levels), colour="gold", y=20, label="Searches for clegg")+geom_line(aes(y=0), position = position_jitter(h = 0.1), colour="black")
ggplot(data=a,aes(x=Levels, y=Edges)) + geom_smooth(aes(y=Vertices), colour="gold", size=1, position="jitter", fill=NA, method="loess")+geom_smooth(colour="white", size=3, position="jitter", fill=NA, method="loess")+geom_smooth(colour="dark blue", size=1, position="jitter", fill=NA, method="loess")+ geom_text(data=a[2, ], aes(x=Levels), colour="gold", y=20, label="Searches for clegg")+geom_line(aes(y=0), position = position_jitter(h = 0.1), colour="black") + theme_xkcd
library("knitr", lib.loc="/Users/anaraven/Library/R/2.15/library")
knit2html
library(markdown)
markdownToHTML
render_html
renderMarkdown
registeredRenderers
registeredRenderers()
markdownToHTML
install.packages("devtools", lib="/Library/Frameworks/R.framework/Versions/2.15/Resources/library")
library("devtools", lib.loc="/Library/Frameworks/R.framework/Versions/2.15/Resources/library")
install_github("slidify", "ramnathv")
Sys.getenv("TAR")
Sys.setenv(TAR = '/usr/bin/tar')
install_github("slidify", "ramnathv")
install_github("slidifyLibraries", "ramnathv")
markdownToHTML
renderMarkdown
?renderMarkdown
library(minet)
?minet
?infotheo
condinformation
entropy
op <- read.table("Input/RDB8.1/operon_list.clean.txt", row.names=1, as.is=T)
o <- op[,1]
names(o) <- rownames(op)
rm(op)
head(o)
net <- delete.vertices(net, net.v.names[!net.v.names %in% rownames(op)])
q<- read.table("Input/Genes/Prodoric/MEME/g0.2val", col.names=c("from","to","blast","meme"))
net <- graph.data.frame(q, directed=T)
library(igraph)
net <- graph.data.frame(q, directed=T)
net.v.names <- V(net)$name
net <- delete.vertices(net, net.v.names[!net.v.names %in% rownames(op)])
net <- delete.vertices(net, net.v.names[!net.v.names %in% names(o)])
net
net.v.id <- 1:length(net.v.names)
names(net.v.id) <- net.v.names
onn <- o[net.v.names]
table(is.na(onn))
onn <- onn[!is.na(onn)]
u.name <- unique(onn)
u.id <- 1:length(u.name)
names(u.id) <- u.name
net2 <- contract.vertices(net, u.id[onn], vertex.attr.comb="first")
net2
head(V(net2))
V(net2)$name <- u.name
head(V(net2))
net2 <- simplify(net2, edge.attr.comb="min")
net2
all <- get.data.frame(net2)
rownames(all) <- paste(all$from,all$to)
gs <-read.table("Input/RDB8.1/gold-std.ncol", as.is=T)
rownames(gs) <- paste(gs$from, gs$to)
rm(net.v.id, net.v.names, onn, u.id, u.name)
gs <-read.table("Input/RDB8.1/gold-std.ncol", as.is=T, col.names=c("from","to","zero"))
rownames(gs) <- paste(gs$from, gs$to)
all$in.gs <- rownames(all) %in% rownames(gs)
rm(net.v.id, net.v.names, onn, u.id, u.name, o, net)
library(ggplot2)
qplot(meme, data=all)
summary(all$meme)
qplot(meme, data=all, fill = I('#099D09'))
qplot(meme, data=all, fill = I('#099DD9'))
qplot(meme, data=all, fill = I('#099DD9'), color=I('black'))
qplot(meme, data=all, fill = I('#099DD9'), color=I('black'), binwidth=1e-6)
qplot(meme, data=all, fill = I('#099DD9'), color=I('black'), binwidth=1e-6, xlim=c(0,1e-4))
qplot(meme, data=all, fill = I('#099DD9'), color=I('black'), binwidth=1e-6, xlim=c(0,1e-4), geom="freqpoly")
qplot(meme, data=all, fill = I('#099DD9'), color=I('black'), binwidth=1e-6, xlim=c(0,1e-4), geom="freqpoly", color=in.gs)
qplot(meme, data=all, fill = I('#099DD9'), color=I('black'), binwidth=1e-6, xlim=c(0,1e-4), geom="freqpoly", color=1+in.gs)
qplot(meme, data=all, fill = I('#099DD9'), color=I('black'), binwidth=1e-6, xlim=c(0,1e-4), geom="freqpoly", colour=1+in.gs)
qplot(meme, ..density.., data=all, fill = I('#099DD9'), color=I('black'), binwidth=1e-6, xlim=c(0,1e-4), geom="freqpoly", colour=1+in.gs)
qplot(meme, data=all, fill = I('#099DD9'), color=I('black'), binwidth=1e-6, xlim=c(0,1e-4), geom="freqpoly")
summary(all)
qplot(meme, data=all, binwidth=1e-6, xlim=c(0,1e-4), geom="freqpoly", color=in.gs)
qplot(meme, ..density.., data=all, binwidth=1e-6, xlim=c(0,1e-4), geom="freqpoly", color=in.gs)
qplot(meme, ..density.., data=all, binwidth=1e-6, xlim=c(0,1e-4), geom="freqpoly", color=in.gs, log="y")
qplot(meme, ..density.., data=all, binwidth=1e-6, xlim=c(0,1e-4), geom="freqpoly", color=in.gs)
qplot(meme, ..density.., data=all, binwidth=1e-6, xlim=c(0,1e-4), color=in.gs)
qplot(meme, data=all, binwidth=1e-6, xlim=c(0,1e-4), color=in.gs)
qplot(meme, data=all, binwidth=1e-6, xlim=c(0,1e-4), fill=in.gs)
qplot(meme, data=all, binwidth=1e-6, xlim=c(0,1e-4), fill=ifelse(in.gs,"Validated","Predicted")
)
qplot(meme, data=all, binwidth=1e-6, xlim=c(0,1e-4), fill=ifelse(!in.gs,"Validated","Predicted")
)
qplot(meme, data=all, binwidth=1e-6, xlim=c(0,1e-4), fill=ifelse(!in.gs,"Validated","Predicted"), color="black")
qplot(meme, data=all, binwidth=1e-6, xlim=c(0,1e-4), fill=ifelse(!in.gs,"Validated","Predicted"), color=I("black"))
qplot(meme, data=all, binwidth=1e-6, xlim=c(0,1e-4), fill=ifelse(!in.gs,"Predicted","Validated"), color=I("black"))
qplot(meme, data=all, binwidth=1e-6, xlim=c(0,1e-4), fill=ifelse(!in.gs,"Validated","Predicted"), color=I("black"))
qplot(meme, data=all, binwidth=2e-6, xlim=c(0,1e-4), fill=ifelse(!in.gs,"Validated","Predicted"), color=I("black"))
qplot(meme, data=all, binwidth=2e-6, xlim=c(0,1e-4), fill=ifelse(!in.gs,"Validated","Predicted"), color=I("black")) + labs("labs")
qplot(meme, data=all, binwidth=2e-6, xlim=c(0,1e-4), fill=ifelse(!in.gs,"Validated","Predicted"), color=I("black")) + theme(legend.title="")
qplot(meme, data=all, binwidth=2e-6, xlim=c(0,1e-4), fill=ifelse(!in.gs,"Validated","Predicted"), color=I("black")) + theme(legend.title=NULL)
qplot(meme, data=all, binwidth=2e-6, xlim=c(0,1e-4), fill=ifelse(!in.gs,"Validated","Predicted"), color=I("black")) + theme(legend.title=element_text(""))
qplot(meme, data=all, binwidth=2e-6, xlim=c(0,1e-4), fill=ifelse(!in.gs,"Validated","Predicted"), color=I("black")) + theme(legend.title=element_text("sss"))
qplot(meme, data=all, binwidth=2e-6, xlim=c(0,1e-4), fill=ifelse(!in.gs,"Validated","Predicted"), color=I("black")) + theme_bw()
all$t=factor(ifelse(!in.gs,"Validated","Predicted"), levels=c("Validated","Predicted"))
all$t=factor(ifelse(!all$in.gs,"Validated","Predicted"), levels=c("Validated","Predicted"))
qplot(meme, data=all, binwidth=2e-6, xlim=c(0,1e-4), fill=t, color=I("black")) + theme_bw()
all$t=factor(ifelse(!all$in.gs,"Validated","Predicted"), levels=c("Predicted","Validated"))
qplot(meme, data=all, binwidth=2e-6, xlim=c(0,1e-4), fill=t, color=I("black")) + theme_bw()
all$t=factor(ifelse(all$in.gs,"Validated","Predicted"), levels=c("Predicted","Validated"))
qplot(meme, data=all, binwidth=2e-6, xlim=c(0,1e-4), fill=t, color=I("black")) + theme_bw()

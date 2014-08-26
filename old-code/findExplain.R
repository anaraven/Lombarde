PDF()
library(igraph)
library(parallel)
library(minet)
source("funcExplain.R")

w <- 1500 - oo2$weight
L <- list(3,5) # list(2,3,5)
thr <- lapply(L,  function(n) quantile(w, (0:n)/n))
wgt  <- lapply(thr, function(t) 10^(cut(w, breaks=t, include.lowest=T, labels=F)-1))
wgt[[length(L)+1]] <- w
oo3  <- lapply(wgt, function(ww) {o <-oo2; o$weight <- ww; o})
# oo3[[length(L)+1]] <- oo
nets <- lapply(oo3, graph.data.frame, directed=TRUE)

N <- length(nets)
W <- lapply(nets, shortest.paths, mode="out")
# filename <- paste("MutualInfo", "mi-pearson-mrnet-133", "coexp-operons.txt", sep="/")
# v.obs <- read.obs(filename, rownames(W[[1]]))
filename <- "Expression/mi-pearson-aracne-133.txt"
v.obs <- read.obs(filename, names(ops))
colnames(v.obs) <- c("from", "to", "weight")
v.obs <- v.obs[v.obs$from < v.obs$to, ]
v.obs <- contract.data.frame(v.obs, ops, max)
is.valid.obs <- v.obs$from %in% rownames(W[[1]]) & v.obs$to %in% rownames(W[[1]])
v.obs <- v.obs[is.valid.obs,]
v.obs <- v.obs[order(v.obs$weight, decreasing=T),]

all.paths <- mclapply(1:N, function(n) find.explanations(v.obs, nets[[n]], W[[n]]))
found <- mclapply(all.paths, count.found, rownames(ref.oo))
perf <- lapply(found, table.found)
cat(sapply(perf, auc.pr)*100,"\n")

plot(0,0, xlab = "recall", ylab = "precision", main = "PR-Curve", xlim = c(0,0.4), ylim=c(0,0.4),
     type="n")
# show.pr(perf[[1]], device=2, col=1, cex=0.4, type="o")
sapply(1:length(perf), function(n) show.pr(perf[[n]], device=2, col=n, type="o", cex=0.4))
legend("topright",leg=c("3 levels","continuous"),fill=c("black","red"))

#### fixes W, iterates over several v.obs
W <- W[[1]]
net <- nets[[1]]
method <- list("corr-pearson",   "mi-pearson-aracne","mi-pearson-c3",
               "mi-pearson-clr", "mi-pearson-mrnet")

filename <- paste("Expression/", method, "-133.txt", sep="")
v.obs <- mclapply(filename, read.obs, names(ops), colnames=c("from", "to", "weight"))
v.obs <- mclapply(v.obs, contract.data.frame, ops, max)
v.obs <- mclapply(v.obs, function(v) {
  is.valid.obs <- v$from %in% rownames(W) & v$to %in% rownames(W)
  v <- v[is.valid.obs,]
  v[order(v$weight, decreasing=T),]
  })

all.paths <- mclapply(v.obs, find.explanations, net, W)
found <- mclapply(all.paths, count.found, rownames(ref.oo))
perf <- lapply(found, table.found)
cat(sapply(perf, auc.pr)*100,"\n")

plot(0,0, xlab = "recall", ylab = "precision", main = "PR-Curve", xlim = c(0,0.4), ylim=c(0,0.4),
     type="n")
show.pr(perf[[1]], device=2, col=1, cex=0.4, type="o")
sapply(2:length(perf), function(n) show.pr(perf[[n]], device=2, col=n, type="o", cex=0.4))
dev.off()

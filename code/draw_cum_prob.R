library(igraph)
basedir <- "~/RLombarde/Out1/EvalBase"

NETWORKS = c('RegulonDB','Prodoric')
LEVELS = c(3, 9)
BASES = c(1.2, 2, 10)
MEMES = c('MEME', 'MEME-GS-N', 'MEME-GS-0')
COEXPS = c('M3D', 'RND1', 'RND2', 'RND3')


analize.graph <- function(g, gs.edges, induce=FALSE) {
  g.v.names <- V(g)$name
  if(induce){
    g <- delete.vertices(g, g.v.names[!g.v.names %in% gs.v.names])
    g.v.names <- V(g)$name  
  }
  
  g.edges <- get.edges(g, E(g))
  g.edges <- split(g.v.names[g.edges[,2]], g.v.names[g.edges[,1]])
  
  in.both <- sum(sapply(names(g.edges),
                        function(x) length(intersect(gs.edges[[x]], g.edges[[x]]))) )
  list(vertices=vcount(g), edges=ecount(g), in.both=in.both)
}

intersect.graphs <- function(gs.file, net.file) {
  gs <- read.graph(gs.file, format="ncol", names=T, direct=T, weights="no")
  # we store the gold-standard vertices names in a vector for easy access 
  gs.v.names <- V(gs)$name
  # and the edges as a two column matrix...
  gs.edges <- get.edges(gs, E(gs))
  # ...trasnformed in a two-level tree of node names
  gs.edges <- split(gs.v.names[gs.edges[,2]], gs.v.names[gs.edges[,1]])
  net <- read.graph(net.file, format="ncol", names=TRUE, direct=TRUE, weights="yes")
  analize.graph(net, gs.edges)
}

png(filename="Prob%03d.png", width=8.5, height=11, units='in', res=200)
par(mfrow=c(3,2))
g=list()
pr=list()
for(network in NETWORKS) {
  for(meme in MEMES) {
    l <- intersect.graphs("~/RLombarde/Input/RDB8.1/gold-std.ncol",
              paste("~/RLombarde/Input/RDB8.1", network, meme, "g0-3.ncol", sep="/"))
    VV <- l$in.both
    NN <- l$edges
    for(base in BASES) {
      for(level in LEVELS) {
        for(coexp in COEXPS) {
          g[[coexp]] <- read.table(
            paste(basedir,coexp,network,meme,level,base,"gl.txt",sep="/"),
            header=TRUE)
        }
        for(step in c(20, 50, 100, 200, 500, 1000)) {
          for(coexp in COEXPS) {
            pr[[coexp]] <- NULL
            TT <- floor(min(nrow(g[[coexp]]),30000)/step)
            cat(network, level, base, coexp, meme, step, VV, NN, TT, "\n")
            t <- 1:TT*step
            V <- c(0, g[[coexp]][t, "n_valid"])[1:TT]
            N <- c(0, g[[coexp]][t, "n_arcs" ])[1:TT]
            v <- diff(V)
            n <- diff(N)
            pr[[coexp]] <- phyper(q=v, m=VV-V, n=NN-VV-N+V, k=n, lower.tail=FALSE)
          }
          plot(pr[["M3D"]]~t, log="y", ylab="prob", pch=16, cex=1)
          abline(h=0.01)
          title(paste(network, meme, level, base, "step", step))
          lines(pr[["RND1"]]~t, col="red")
          
        }
      }
    }    
  }
}


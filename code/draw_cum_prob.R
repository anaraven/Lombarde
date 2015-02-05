library(igraph)
basedir <- "~/Documents/RLombarde/Out1/EvalBase"

METHODS = c('aracne','c3','clr','mrnetb')
NETWORKS = c('RegulonDB','Prodoric')
LEVELS = c(3, 5, 7, 9)
BASES = c(1.2, 2, 5, 10)
MEMES = c('MEME', 'MEME-GS-N', 'MEME-GS-0')
COEXPS = c('M3D', 'RND1', 'RND2', 'RND3')

intersect.graphs <- function(gs.file, net.file) {
  gs <- read.graph(gs.file, format="ncol", names=T, direct=T, weights="no")
  # we store the gold-standard vertices names in a vector for easy access 
  gs.v.names <- V(gs)$name
  # and the edges as a two column matrix...
  gs.edges <- get.edges(gs, E(gs))
  # ...trasnformed in a two-level tree of node names
  gs.edges <- split(gs.v.names[gs.edges[,2]], gs.v.names[gs.edges[,1]])
  
  g <- read.graph(net.file, format="ncol", names=TRUE, direct=TRUE, weights="yes")
  # analize graph g versus  gs.edges
  g.v.names <- V(g)$name
  # if(induce){
  #  g <- delete.vertices(g, g.v.names[!g.v.names %in% gs.v.names])
  #  g.v.names <- V(g)$name  
  #}
  
  g.edges <- get.edges(g, E(g))
  g.edges <- split(g.v.names[g.edges[,2]], g.v.names[g.edges[,1]])
  
  in.both <- sum(sapply(names(g.edges),
                        function(x) length(intersect(gs.edges[[x]], g.edges[[x]]))) )
  list(vertices=vcount(g), edges=ecount(g), in.both=in.both)
}

calc.prob <- function(g, lag, N.tot, V.tot) {
  V <- c(0, g[, "n_valid"])
  N <- c(0, g[, "n_arcs" ])
  v <- diff(V, lag=lag)
  n <- diff(N, lag=lag)
  rng <- 1:length(n)
  phyper(q=v[rng], m=V.tot-V[rng], n=N.tot-V.tot-N[rng]+V[rng], k=n[rng], lower.tail=FALSE)
}

png(filename="Prob%03d.png", width=8.5, height=11, units='in', res=200)
par(mfrow=c(3,2))
g=list()
pr=list()
for(method in METHODS) {
  for(network in NETWORKS) {
    for(meme in MEMES) {
      l <- intersect.graphs("~/RLombarde/Input/RDB8.1/gold-std.ncol",
                            paste("~/RLombarde/Input/RDB8.1", network, meme, "g0-3.ncol", sep="/"))
      V.tot <- l$in.both
      N.tot <- l$edges
      for(base in BASES) {
        for(level in LEVELS) {
          for(coexp in COEXPS) {
            g[[coexp]] <- read.table(
              paste(basedir, coexp, network, meme, level, base, "gl.txt", sep="/"),
              header=TRUE)
          }
          for(lag in c(20, 50, 100, 200, 500, 1000)) {
            for(coexp in COEXPS) {
              pr[[coexp]] <- calc.prob(g[[coexp]], lag, N.tot, V.tot)
              pr[[coexp]][pr[[coexp]]==0] <- NA
            }
            mi <- min(sapply(pr, min, na.rm=T))
            ma <- max(sapply(pr, max, na.rm=T))
            plot( pr[["M3D" ]], log="y", ylab="prob", ylim=c(mi,ma), type="n")
            lines(pr[["RND1"]], col="red") 
            lines(pr[["RND2"]], col="red") 
            lines(pr[["RND3"]], col="red") 
            points(pr[["M3D"]], cex=0.6)
            abline(h=0.01)
            title(paste(network, meme, level, base, "lag", lag))
          }
        }
      }    
    }
  }  
}


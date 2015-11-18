intersection.by.name <- function (g1, g2) 
{
  dv1 = igraph::get.vertices.as.data.frame(g1)
  dv2 = igraph::get.vertices.as.data.frame(g2)
  de1 = igraph::get.edges.as.data.frame(g1)
  de2 = igraph::get.edges.as.data.frame(g2)
  dv = igraph::safer.merge(dv1, dv2)
  de = igraph::safer.merge(de1, de2)
  g = igraph::graph.data.frame(de, directed = FALSE, vertices = dv)
  return(g)
}

valenteForeman <- function(g, what="both") {
  sp <- shortest.paths(g, mode="out", weight=NA)
  rd <- max(sp[is.finite(sp)])+1-sp
  rd[is.infinite(rd)] <- 0
  diag(rd) <- 0
  what <- match.arg(what, c("integration","radiality","both"))
  if(what=="both"){
    return( list(integration=colSums(rd)/(nrow(rd)-1), radiality=rowSums(rd)/(ncol(rd)-1)) )
  } else if(what=="radiality") {
    return( rowSums(rd)/(ncol(rd)-1) )
  } else {
    return( colSums(rd)/(ncol(rd)-1) )
  }
}

topNodes <- function(x,n=length(x)) names(x)[ order(x, decreasing=T)[1:n] ]
infoList <- function(x,field="note",genes=ops) lapply(genes[x], function(y) nom[y,field])

showTop <- function(g, n=5, e=1) {
  gg  <- delete.vertices(g, degree(g, mode="out")==0)
  tg  <- t(g[,])
  tg  <- graph.adjacency(tg)
  ec  <- evcent(g,  directed=T,scale=T)
  ecc <- evcent(gg, directed=T,scale=T)
  tec <- evcent(tg, directed=T,scale=T)  
  ir  <- valenteForeman(g )
  irr <- valenteForeman(gg)
  b   <- bonpow(g, exponent=e/ec$value )
  bb  <- bonpow(gg,exponent=e/ecc$value)
  
  list(ec=topNodes(ec$vector, n), ecc=topNodes(ecc$vector,n), tec=topNodes(tec$vector,n),
      integ=topNodes(ir$integration,n), radial=topNodes(ir$radiality,n),
      integ2=topNodes(irr$integration,n), radial2=topNodes(irr$radiality,n),
       b=topNodes(b,n), bb=topNodes(bb,n))
}


tabTF <- function(x) {
  data.frame(name=unlist(infoList(x,genes=TF,field="gene")), desc=unlist(infoList(x,genes=TF)))
}
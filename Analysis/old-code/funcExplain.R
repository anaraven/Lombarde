read.obs <- function(filename, nom, n.max=Inf, ...) {
  # `obs` elements are associated operons
  obs <- read.table(filename, as.is=T, ...)
  obs[,1] <- tolower(obs[,1])
  obs[,2] <- tolower(obs[,2])
  is.valid.obs <- obs[,1] %in% nom & obs[,2] %in% nom
  v.obs <- obs[is.valid.obs,]
  colnames(v.obs) <- c("from", "to", "weight")
  v.obs <- v.obs[v.obs$from < v.obs$to, ]
  if(nrow(v.obs)>n.max) {
    position <- order(abs(v.obs$weight), decreasing=T)
    v.obs <- v.obs[position[1:n.max], ]
  }
  return(v.obs)  
}

contract.data.frame <- function(m, ops, f) {
  m$from <- as.character(ops[m$from])
  m$to   <- as.character(ops[m$to])
  m <- m[!is.na(m$from) & !is.na(m$to),]
  m <- m[m$from!=m$to,]
  link <- paste(m$from, m$to)
  out <- data.frame(from  =sapply(split(m$from,   link), f),
                    to    =sapply(split(m$to,     link), f),
                    weight=sapply(split(m$weight, link), f))
  out$from <- as.character(out$from)
  out$to   <- as.character(out$to)
  out
}

discretize.weight <- function(w, n) {
  thr <- quantile(w, (0:n)/n)
  wgt <- 10^(cut(w, breaks=thr, include.lowest=T, labels=F)-1)
  wgt
}

library(igraph)
find.all.paths <- function(oo, v.obs) {
  net <- graph.data.frame(oo, directed=TRUE)
  W <- shortest.paths(net, mode="out")
  is.valid.obs <- v.obs$from %in% rownames(W) & v.obs$to %in% rownames(W)
  v.obs <- v.obs[is.valid.obs,]
  v.obs <- v.obs[order(v.obs$weight, decreasing=T),]
  find.explanations(v.obs, net, W)
}


list.paths <- function(i, exps, v.obs, nom, net) {
  sp <- lapply(names(exps[[i]]), 
               function(r) get.all.shortest.paths(net, r, to=unlist(v.obs[i,1:2]), mode="out")$res)
  lapply(sp, lapply, function(x) nom[x])
}

find.explanations <- function(v.obs, net, W) {
  N <- nrow(v.obs)
  exps <- lapply(1:N, function(i) {v <- rowSums(W[,unlist(v.obs[i,1:2])]); which(v==min(v))})
  lapply(1:N, list.paths, exps, v.obs, rownames(W), net)  
}

count.found <- function(all.paths, ref.names, oo, w.col="w1") {
  N <- length(all.paths)
  n.edges  <- vector("integer", N)
  n.in.ref <- vector("integer", N)
  out.list <- list()
  for(i in 1:N) {
    one.obs <- unlist(all.paths[[i]], recursive=FALSE)
    for(pth in one.obs) {
      n <- length(pth)
      if(n>1) {
        for(j in 2:n) {
          id <- paste(pth[j-1], pth[j])
          out.list[[id]] <- oo[id, w.col]  
        }  
      }
    }
    n.edges [i] <- length(out.list)
    n.in.ref[i] <- length(intersect(names(out.list), ref.names))
  }
  return(list(n.edges=n.edges, n.in.ref=n.in.ref, out.list=out.list, n.ref=length(ref.names)))
}

count.found.restricted <- function(all.paths, ref.names) {
  N <- length(all.paths)
  n.edges  <- vector("integer", N)
  n.in.ref <- vector("integer", N)
  src <- unique(sapply(strsplit(ref.names," "),`[[`,1))
  dst <- unique(sapply(strsplit(ref.names," "),`[[`,2))
  out.list <- list()
  for(i in 1:N) {
    one.obs <- unlist(all.paths[[i]], recursive=FALSE)
    for(pth in one.obs) {
      n <- length(pth)
      if(n>1) {
        inside <- pth[1:(n-1)] %in% src & pth[2:n] %in% dst
        for(j in 2:n) {
          if(inside[j-1]) {
            id <- paste(pth[j-1], pth[j])
#            out.list[[id]] <- oo[id,"weight"]            
            out.list[[id]] <- 1    
          }
        }
      }
    }
    n.edges [i] <- length(out.list)
    n.in.ref[i] <- length(intersect(names(out.list), ref.names))
  }
  return(list(n.edges=n.edges, n.in.ref=n.in.ref, out.list=out.list, n.ref=length(ref.names)))
}

table.found <- function(found) {
  with(found, data.frame(thrsh=1:length(n.edges), tp=n.in.ref, fp=n.edges-n.in.ref, fn=n.ref-n.in.ref)) 
}


plot.found <- function(found, n, ...) {
  plot(1-found$n.in.ref/n, found$n.in.ref/found$n.edges,
       xlab="Recall", ylab="1-Precision", ...)
}

lines.found <- function(found, n, ...) {
  lines(1-found$n.in.ref/n, found$n.in.ref/found$n.edges, ...)
}

edge.i.each.path <- function(all.paths, oo) {
  each.path <- unlist(unlist(all.paths, recursive=FALSE), recursive=FALSE)
  ll <- sapply(each.path, length)
  each.path <- each.path[ll>1]
  
  out.list0 <- list()
  for(pth in each.path) {
    n <- length(pth)
    for(j in 2:n) {
      id <- paste(pth[j-1], pth[j])
      out.list0[[id]] <- oo[id,"weight"]  
    }
  }
}

load.operons <- function() {
  oper <- read.delim("Operons/inOperon3.txt", header=F, as.is=T)
  ops <- oper$V1
  names(ops) <- oper$V2
  rm(oper)
  ops
}

# given gene-to-protein and protein-to-gene tables, build a list of gene-to-gene links
gene.to.gene.list <- function(gp, pg) {
  links <- list()
  k <- 1
  for(i in 1:nrow(gp)) {
    select <- pg$from == gp$to[i]
    for(j in (1:nrow(pg))[select]) {
      links[[k]] <- c(from=gp$from[i], to=pg$to[j], w1=gp$e[i]+pg$p[j], w2=gp$e[i]+pg$q[j])
      k <- k + 1
    }
  }
  links
}


# contract gene-to-gene links into operon-to-operon table
getLinkValues <- function(i, x) unlist(lapply(x, `[[`, i), use.names=FALSE)

o2o <- function(links, ops) {
  m <- data.frame(from=ops[getLinkValues("from", links)], to=ops[getLinkValues("to", links)],
                  w1=as.numeric(getLinkValues("w1", links)), w2=as.numeric(getLinkValues("w2", links)))
  m <- m[!is.na(m$to),]
  m$from <- as.character(m$from)
  m$to   <- as.character(m$to)
  m$link <- paste(m$from, m$to)
  m <- m[m$from!=m$to,]
  
  data.frame(from =sapply(split(m$from, m$link), min),
             to   =sapply(split(m$to,   m$link), min),
             w1   =sapply(split(m$w1,   m$link), min),
             w2   =sapply(split(m$w2,   m$link), min))
}

gold.o2o <- function(gold, ops) {
  m.ref <- data.frame(from=ops[getLinkValues("from", gold)], to=ops[getLinkValues("to", gold)])
  m.ref <- m.ref[!is.na(m.ref$to),]
  m.ref$from <- as.character(m.ref$from)
  m.ref$to   <- as.character(m.ref$to)
  m.ref$link <- paste(m.ref$from, m.ref$to)
  m.ref <- m.ref[m.ref$from!=m.ref$to,]
  
  data.frame(from  =sapply(split(m.ref$from, m.ref$link), `[[`, 1),
             to    =sapply(split(m.ref$to,   m.ref$link), `[[`, 1))
}


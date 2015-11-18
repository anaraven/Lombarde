jcv    <- read.delim("RegulonDB/Data/equiv-jcv.txt",   header=F, as.is=T)
geneTF <- read.delim("RegulonDB/Data/geneTF.txt",      header=F, as.is=T)
tfReg  <- read.delim("RegulonDB/Data/tfRegulGene.txt",      header=F, as.is=T)
equiv <- jcv$V1
names(equiv) <- jcv$V2

gp <- data.frame(from=equiv[geneTF$V1], to=geneTF$V2, e=0)
gp <- gp[!is.na(gp$from),]
gp$from <- as.character(gp$from)
gp$to   <- as.character(gp$to)

pg.ref <- data.frame(from=tfReg$V1, to=equiv[tfReg$V2], sign=tfReg$V3, p=0)
pg.ref <- pg.ref[!is.na(pg.ref$from),]
pg.ref$from <- as.character(pg.ref$from)
pg.ref$to   <- as.character(pg.ref$to)

gold <- list()
k <- 1
for(i in 1:nrow(gp)) {
  select <- pg.ref$from == gp$to[i]
  for(j in (1:nrow(pg.ref))[select]) {
    gold[[k]] <- c(from=gp$from[i], to=pg.ref$to[j])
    k <- k + 1
  }
}

rm(select, i, j, k)
rm(equiv, jcv, geneTF, gp, pg.ref)


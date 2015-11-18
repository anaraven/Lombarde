score <- read.table("rdb-vf.txt",row.names=1, col.names=c("id","integ","radial"))
tmp <- read.table("rdb-tfnames.txt",as.is=T)
n <- tmp$V2
names(n) <- tmp$V1
rm(tmp)

good <- score$integ>2 & score$radial>0.1
nom <- n[rownames(score)]

mc <- read.delim("~/Desktop/IRISA/branches/PWM_MI/Ecoli/martinez-collado.txt")
in.mc <- nom %in% mc$gene

plot(score[good,],type="n",xlab="Intergration",ylab="Radiality")
text(score[good,], nom[good], col=nom[good] %in% mc$gene+1)

rr <- length(score$radiality)-rank(score$radiality, ties.method="max")[in.mc]
names(rr) <- nom[in.mc]
mc$our <- rr[as.character(mc$gene)]+1


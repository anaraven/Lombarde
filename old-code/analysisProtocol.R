nom <- read.table("annotation.txt", sep="\t", as.is=T,quote="")
rownames(nom) <- tolower(rownames(nom))

chrom <- read.table("Operons/operons-raw.txt", skip=3, header=T,
            sep="\t", as.is=T, row.names=3)
ops <- split(rownames(chrom), sprintf("op%04d",chrom$Operon))

basedir <- "Out/RegulonDB/Meme"
# basedir <- "Out/Prodoric/Meme"
tf <- read.table(paste(basedir,"bipartite-genes/tfInOperon.txt",sep="/"), as.is=T)
TF <- split(tf$V1,tf$V2)

library(igraph)
source("analysisFunctions.R")

#g <- read.graph(paste(basedir,"graph-operons/full-p.txt",sep="/"),
#             format="ncol", names=T, direct=T, weights="no")
#
#rdb <- read.graph("RegulonDB/Data/graph-operons.csv",
#            format="ncol", names=T, direct=T, weights="no")

out <- read.graph(paste(basedir,"mi-pearson-mrnet-133/out_graph/full-p.txt",sep="/"),
                format="ncol", names=T, direct=T, weights="no")

#rdb10 <- showTop(rdb,20)
#g10   <- showTop(g,  20)
out10 <- showTop(out,20)

#sapply(names(rdb10), function(x) length(intersect(rdb10[[x]],g10[[x]])))
#sapply(names(rdb10), function(x) length(intersect(rdb10[[x]],out10[[x]])))
#sapply(names(out10), function(x) length(intersect(out10[[x]],g10[[x]])))

#unlist(infoList(topNodes(ecc$vector,5)), use.names=F)

vf <- valenteForeman(out)
write.table(as.data.frame(vf),"rdb-vf.txt",sep="\t",quote=FALSE, col.names=FALSE)
tfnames <- sapply(TF,function(y) paste(nom[y,"gene"],collapse="-"))
write.table(as.data.frame(tfnames),"rdb-tfnames.txt",sep="\t",quote=FALSE, col.names=FALSE)
gg <- delete.vertices(out, degree(out, mode="out")==0)
write.graph(gg,"rdb-core.ncol","ncol")

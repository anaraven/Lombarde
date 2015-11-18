read.node.score <- function(fname) {
  I <- read.table(fname,col.names=c("N","FN","TP","FP"),row.names=1)
  I$PR <- I$TP/(I$TP+I$FP)
  I$RC <- I$TP/(I$TP+I$FN)
  return(I)
}

# apply to
# RegulonDB/Meme/node-score.txt
# RegulonDB/Meme/cardinality-node-score.txt
# RegulonDB/Meme/wgt/node-score.txt

conserved.names <- function(I,C,W) {
  intersect(rownames(I),intersect(rownames(C),rownames(W)))
}

eval.matrix <- function(I,W,kept) {
  data.frame(TP=W[kept,"TP"]-I[kept,"TP"], FP=I[kept,"FP"]-W[kept,"FP"], FN=I[kept,"FN"]-W[kept,"FN"])
}
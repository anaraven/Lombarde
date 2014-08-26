a <- read.table("noinduced-summary.txt", header=T, sep="\t")
meme.pearson.initial <- a[(a$estimator=="pearson" | a$estimator=="Initial") & a$BS=="MEME", ]
meme.pearson <- a[a$estimator=="pearson" & a$BS=="MEME", ]
gs0.pearson  <- a[a$estimator=="pearson" & a$BS=="MEME-GS-0", ]


attach(a)
plot(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"), xaxp=c(3,9,3), ylim=c(0,6500), main="Ab initio")
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="aracne"),lty=2,pch=2)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="c3"),lty=3,pch=3)
points(Edges~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="clr"),lty=4,pch=4)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="mrnet"), col="blue")
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="aracne"), col="blue", lty=2, pch=2)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="c3"), col="blue", lty=3, pch=3)
points(In.gold~Levels, type="o", subset=(BS=="MEME" & estimator=="pearson" & MI=="clr"), col="blue", lty=4, pch=4)


plot(Edges~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="mrnet"), xaxp=c(3,9,3), ylim=c(0,6500), main="Extended")
points(Edges~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="aracne"),lty=2,pch=2)
points(Edges~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="c3"),lty=3,pch=3)
points(Edges~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="clr"),lty=4,pch=4)
points(In.gold~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="mrnet"), col="blue")
points(In.gold~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="aracne"), col="blue", lty=2, pch=2)
points(In.gold~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="c3"), col="blue", lty=3, pch=3)
points(In.gold~Levels, type="o", subset=(BS=="MEME-GS-0" & estimator=="pearson" & MI=="clr"), col="blue", lty=4, pch=4)

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


plot(Edges~factor(MI, levels=c("aracne","c3","clr","mrnet")), type="o", subset=(BS=="MEME" & estimator=="pearson"), ylim=c(2000,6500), main="Ab initio", 
     xlab="Co-expr. detect. meth.", ylab="Num. Arcs")

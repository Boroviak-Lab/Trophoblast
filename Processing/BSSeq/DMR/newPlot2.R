library(Seurat)
library(ggplot2)
library(caret)

#D1 <- read.table('AllDEInViv.mat',sep="\t",header = F)
#D2 <- read.table('AllDEInVit.mat',sep="\t",header = F)
#D3 <- read.table('AllInViv.mat',sep="\t",header = F)
#D4 <- read.table('AllInVit.mat',sep="\t",header = F)

D1 <- read.table('AllDEInVit.mat',sep="\t",header = F)
#D2 <- read.table('AllInAno.mat',sep="\t",header = F)


Mapping1 <- read.table('mapping3.txt',sep="\t",header = F)
#computeMatrix scale-regions -S Primed1.bw Primed2.bw Primed3.bw Niave1.bw Niave2.bw Niave3.bw TSC1.bw TSC2.bw TSC3.bw EmDisc.bw Tb.bw cMor.bw -R ../TSS_5kb.sh --binSize 1 --beforeRegionStartLength 0 --regionBodyLength 1 --afterRegionStartLength 0 --skipZeros -o AllIn.mat.gz


#X <- data.frame(x=(D2$V9+D2$V8+D2$V7),  y=(D2$V10+D2$V11+D2$V12), z=(D2$V13+D2$V14+D2$V15), u=D2$V16, v=D2$V17, w=D2$V18, genenames=D2$V23 )
#Do pairwise to keep code simpler i.e. alwaays x vs y #primed vs niave
#X <- data.frame(x=(D2$V9+D2$V8+D2$V7),  y=(D2$V10+D2$V11+D2$V12), genenames=D2$V23 )
X <- data.frame(x=log(D1$V16+1),  y=log(D1$V18+1), genenames=Mapping1$V2 )
X <- na.omit(X)
genenames <- X$genenames
lrfit <- train(y~poly(x,degree=3), data=X, method = "lm")
predictedValues<-predict(lrfit,newdata = data.frame(x=X$x) )
Xhat <- data.frame(x=X$x,y=X$y-predictedValues)
xnew <- seq(min(X$x), max(X$x), length.out= 500)
predictedValues2<-predict(lrfit,newdata = data.frame(x=xnew) )
vals <- quantile(Xhat$y,c(0.01,0.99))
ggplot(X, aes(x = x, y = y)) + geom_point(size=2.5)  + geom_point(color='blue') + theme_bw() + xlab("EmDisc") + ylab("cMor") +
geom_line(color='red',data = data.frame(x=xnew,y=predictedValues2), aes(x=x, y=y))
ggsave(filename=paste("Test9.pdf",sep=""),width = 20, height = 20) 
Xhat$Ext <- 0
Xhat$Ext[which(Xhat$y<vals[1])] <- 1
Xhat$Ext[which(Xhat$y>vals[2])] <- 2
Xhat$Ext <- as.factor(Xhat$Ext)
ggplot(Xhat, aes(x = x, y = y, color=Ext)) + geom_point(size=2.5)  + geom_point() + theme_bw() + xlab("EmDisc") + ylab("cMor") +
scale_color_manual(values=c("lightgrey", "#cc3761", "#457cb0"))+
geom_hline(yintercept = 0)+
geom_hline(yintercept = vals[1])
geom_hline(yintercept = vals[2])
ggsave(filename=paste("Test9_v2.pdf",sep=""),width = 20, height = 20)
X$Ext <- Xhat$Ext
ggplot(X, aes(x = x, y = y, color=Ext)) + geom_point(size=2.5)  + geom_point() + theme_bw() + xlab("EmDisc") + ylab("cMor") +
scale_color_manual(values=c("lightgrey", "#cc3761", "#457cb0")) +
geom_line(color='red',data = data.frame(x=xnew,y=predictedValues2), aes(x=x, y=y))
ggsave(filename=paste("Test9_v3.pdf",sep=""),width = 20, height = 20)
p1 <- ggplot(X, aes(x = x, y = y, color=Ext)) + geom_point(size=2.5)  + geom_point() + theme_bw() + xlab("EmDisc") + ylab("cMor") +
scale_color_manual(values=c("lightgrey", "#cc3761", "#457cb0")) +
geom_line(color='red',data = data.frame(x=xnew,y=predictedValues2), aes(x=x, y=y)) 
p1 <- p1 + geom_text(data=X[c(which(Xhat$y<vals[1]),which(Xhat$y>vals[2])),],aes(x,y,label=genenames))
ggsave(filename=paste("Test9_v4.pdf",sep=""),width = 40, height = 40, p1)
saveext <- "./"
write.csv(as.data.frame(genenames[which(Xhat$y<vals[1])]), file=paste(saveext,"/Run9_lower.csv",sep=""))
write.csv(as.data.frame(genenames[which(Xhat$y>vals[2])]), file=paste(saveext,"/Run9_upper.csv",sep=""))







X <- data.frame(x=log(D1$V16+1),  y=log(D1$V17+1), genenames=Mapping1$V2 )
X <- na.omit(X)
genenames <- X$genenames
lrfit <- train(y~poly(x,degree=3), data=X, method = "lm")
predictedValues<-predict(lrfit,newdata = data.frame(x=X$x) )
Xhat <- data.frame(x=X$x,y=X$y-predictedValues)
xnew <- seq(min(X$x), max(X$x), length.out= 500)
predictedValues2<-predict(lrfit,newdata = data.frame(x=xnew) )
vals <- quantile(Xhat$y,c(0.01,0.99))
ggplot(X, aes(x = x, y = y)) + geom_point(size=2.5)  + geom_point(color='blue') + theme_bw() + xlab("EmDisc") + ylab("Tb") +
geom_line(color='red',data = data.frame(x=xnew,y=predictedValues2), aes(x=x, y=y))
ggsave(filename=paste("Test10.pdf",sep=""),width = 20, height = 20)
Xhat$Ext <- 0
Xhat$Ext[which(Xhat$y<vals[1])] <- 1
Xhat$Ext[which(Xhat$y>vals[2])] <- 2
Xhat$Ext <- as.factor(Xhat$Ext)
ggplot(Xhat, aes(x = x, y = y, color=Ext)) + geom_point(size=2.5)  + geom_point() + theme_bw() + xlab("EmDisc") + ylab("Tb") +
scale_color_manual(values=c("lightgrey", "#cc3761", "#457cb0"))+
geom_hline(yintercept = 0)+
geom_hline(yintercept = vals[1])
geom_hline(yintercept = vals[2])
ggsave(filename=paste("Test10_v2.pdf",sep=""),width = 20, height = 20)
X$Ext <- Xhat$Ext
ggplot(X, aes(x = x, y = y, color=Ext)) + geom_point(size=2.5)  + geom_point() + theme_bw() + xlab("EmDisc") + ylab("Tb") +
scale_color_manual(values=c("lightgrey", "#cc3761", "#457cb0")) +
geom_line(color='red',data = data.frame(x=xnew,y=predictedValues2), aes(x=x, y=y))
ggsave(filename=paste("Test10_v3.pdf",sep=""),width = 20, height = 20)
p1 <- ggplot(X, aes(x = x, y = y, color=Ext)) + geom_point(size=2.5)  + geom_point() + theme_bw() + xlab("EmDisc") + ylab("Tb") +
scale_color_manual(values=c("lightgrey", "#cc3761", "#457cb0")) +
geom_line(color='red',data = data.frame(x=xnew,y=predictedValues2), aes(x=x, y=y))
p1 <- p1 + geom_text(data=X[c(which(Xhat$y<vals[1]),which(Xhat$y>vals[2])),],aes(x,y,label=genenames))
#p1 <- p1 + LabelPoints(plot = p1, points=c(which(Xhat$y<vals[1]),which(Xhat$y>vals[2])), labels = c(genenames[which(Xhat$y<vals[1])],genenames[which(Xhat$y>vals[2])]), repel = TRUE, color = 'red')
ggsave(filename=paste("Test10_v4.pdf",sep=""),width = 40, height = 40)
write.csv(as.data.frame(genenames[which(Xhat$y<vals[1])]), file=paste(saveext,"/Run10_lower.csv",sep=""))
write.csv(as.data.frame(genenames[which(Xhat$y>vals[2])]), file=paste(saveext,"/Run10_upper.csv",sep=""))


sadsadsa

X <- data.frame(x=(D1$V10+D1$V11+D1$V12), y=(D1$V13+D1$V14+D1$V15), genenames=Mapping1$V2 )
X <- na.omit(X)
genenames <- X$genenames
lrfit <- train(y~poly(x,degree=3), data=X, method = "lm")
predictedValues<-predict(lrfit,newdata = data.frame(x=X$x) )
Xhat <- data.frame(x=X$x,y=X$y-predictedValues)
xnew <- seq(min(X$x), max(X$x), length.out= 500)
predictedValues2<-predict(lrfit,newdata = data.frame(x=xnew) )
vals <- quantile(Xhat$y,c(0.01,0.99))
ggplot(X, aes(x = x, y = y)) + geom_point(size=2.5)  + geom_point(color='blue') + theme_bw() + xlab("Primed") + ylab("Naive") +
geom_line(color='red',data = data.frame(x=xnew,y=predictedValues2), aes(x=x, y=y))
ggsave(filename=paste("Test6.pdf",sep=""),width = 20, height = 20)
Xhat$Ext <- 0
Xhat$Ext[which(Xhat$y<vals[1])] <- 1
Xhat$Ext[which(Xhat$y>vals[2])] <- 2
Xhat$Ext <- as.factor(Xhat$Ext)
ggplot(Xhat, aes(x = x, y = y, color=Ext)) + geom_point(size=2.5)  + geom_point() + theme_bw() + xlab("Primed") + ylab("Naive") +
scale_color_manual(values=c("lightgrey", "#cc3761", "#457cb0"))+
geom_hline(yintercept = 0)+
geom_hline(yintercept = vals[1])+
geom_hline(yintercept = vals[2])
ggsave(filename=paste("Test6_v2.pdf",sep=""),width = 20, height = 20)
X$Ext <- Xhat$Ext
ggplot(X, aes(x = x, y = y, color=Ext)) + geom_point(size=2.5)  + geom_point() + theme_bw() + xlab("Primed") + ylab("Naive") +
scale_color_manual(values=c("lightgrey", "#cc3761", "#457cb0")) +
geom_line(color='red',data = data.frame(x=xnew,y=predictedValues2), aes(x=x, y=y))
ggsave(filename=paste("Test6_v3.pdf",sep=""),width = 20, height = 20)

p1 <- ggplot(X, aes(x = x, y = y, color=Ext)) + geom_point(size=2.5)  + geom_point() + theme_bw() + xlab("Primed") + ylab("Naive") +
scale_color_manual(values=c("lightgrey", "#cc3761", "#457cb0"))
p1 <- p1 + geom_text(data=X[c(which(Xhat$y<vals[1]),which(Xhat$y>vals[2])),],aes(x,y,label=genenames))
#geom_line(color='red',data = data.frame(x=xnew,y=predictedValues2), aes(x=x, y=y))
#p1 <- p1 + LabelPoints(plot = p1, points=c(which(Xhat$y<vals[1]),which(Xhat$y>vals[2])), labels = c(genenames[which(Xhat$y<vals[1])],genenames[which(Xhat$y>vals[2])]), repel = TRUE, color = 'red')
ggsave(filename=paste("Test6_v4.pdf",sep=""),width = 40, height = 40,p1)
write.csv(as.data.frame(genenames[which(Xhat$y<vals[1])]), file=paste(saveext,"/Run6_lower.csv",sep=""))
write.csv(as.data.frame(genenames[which(Xhat$y>vals[2])]), file=paste(saveext,"/Run6_upper.csv",sep=""))

asdsadsadda

ggplot(X, aes(x = x, y = y-predictedValues)) + geom_point(size=2.5)  + geom_point(color='blue') + theme_bw() + xlab("Primed") + ylab("Naive") +
geom_abline(intercept = 0, slope = 1, size = 0.5) +
geom_abline(intercept = 50, slope = 1, size = 0.5) +
geom_abline(intercept = -50, slope = 1, size = 0.5)
ggsave(filename=paste("Test1.pdf",sep=""),width = 20, height = 20)ggplot(X, aes(x = x, y = y)) + geom_point(size=2.5)  + geom_point(color='blue') + theme_bw() + xlab("Primed") + ylab("Naive") +

ggplot(X, aes(x = x, y = z)) + geom_point(size=2.5)  + geom_point(color='blue') + theme_bw() + xlab("Primed") + ylab("TSC") +
geom_abline(intercept = 0, slope = 1, size = 0.5) +
geom_abline(intercept = 50, slope = 1, size = 0.5) +
geom_abline(intercept = -50, slope = 1, size = 0.5)
ggsave(filename=paste("Test2.pdf",sep=""),width = 20, height = 20)

ggplot(X, aes(x = y, y = z)) + geom_point(size=2.5)  + geom_point(color='blue') + theme_bw() + xlab("Naive") + ylab("TSC") +
geom_abline(intercept = 0, slope = 1, size = 0.5) +
geom_abline(intercept = 50, slope = 1, size = 0.5) +
geom_abline(intercept = -50, slope = 1, size = 0.5)
ggsave(filename=paste("Test3.pdf",sep=""),width = 20, height = 20)

ggplot(X, aes(x = y, y = z)) + geom_point(size=2.5)  + geom_point(color='blue') + theme_bw() + xlab("Naive") + ylab("TSC") +
geom_abline(intercept = 0, slope = 1, size = 0.5)+
geom_abline(intercept = 50, slope = 1, size = 0.5) +
geom_abline(intercept = -50, slope = 1, size = 0.5)
ggsave(filename=paste("Test3.pdf",sep=""),width = 20, height = 20)

#X <- data.frame(x=D1$V7,y=D1$V11)

ggplot(X, aes(x = u, y = w)) + geom_point(size=2.5)  + geom_point(color='blue') + theme_bw() + xlab("EmDisc") + ylab("Tb") +
geom_abline(intercept = 0, slope = 1, size = 0.5) +
geom_abline(intercept = 50, slope = 1, size = 0.5) +
geom_abline(intercept = -50, slope = 1, size = 0.5)
ggsave(filename=paste("Test4.pdf",sep=""),width = 20, height = 20) 


#X <- data.frame(x=(D4$V9+D4$V8+D4$V7),  y=(D4$V10+D4$V11+D4$V12), z=(D4$V13+D4$V14+D4$V15) )

ggplot(X, aes(x = u, y = v)) + geom_point(size=2.5)  + geom_point(color='blue') + theme_bw() + xlab("Primed") + ylab("Naive") +
geom_abline(intercept = 0, slope = 1, size = 0.5) +
geom_abline(intercept = 50, slope = 1, size = 0.5) +
geom_abline(intercept = -50, slope = 1, size = 0.5)
ggsave(filename=paste("Test5.pdf",sep=""),width = 20, height = 20)

ggplot(X, aes(x = v, y = w)) + geom_point(size=2.5)  + geom_point(color='blue') + theme_bw() + xlab("Primed") + ylab("TSC") +
geom_abline(intercept = 0, slope = 1, size = 0.5) +
geom_abline(intercept = 50, slope = 1, size = 0.5) +
geom_abline(intercept = -50, slope = 1, size = 0.5)
ggsave(filename=paste("Test6.pdf",sep=""),width = 20, height = 20)

#ggplot(X, aes(x = y, y = z)) + geom_point(size=2.5)  + geom_point(color='blue') + theme_bw() + xlab("Naive") + ylab("TSC")
#ggsave(filename=paste("Test7.pdf",sep=""),width = 20, height = 20)



#computeMatrix scale-regions -S Primed1.bw Primed2.bw Primed3.bw Niave1.bw Niave2.bw Niave3.bw TSC1.bw TSC2.bw TSC3.bw -R ../TSS_5kb.sh --binSize 1 --beforeRegionStartLength 0 --regionBodyLength 1 --afterRegionStartLength 0 --skipZeros -o AllInVit.mat.gz
#computeMatrix scale-regions -S EmDisc.bw Am.bw VE.bw ExMes.bw Tb.bw Gland.bw -R ../TSS_5kb.sh --binSize 1 --beforeRegionStartLength 0 --regionBodyLength 1 --afterRegionStartLength 0 --skipZeros -o AllInViv.mat.gz

#bedtools merge -d 1000 -i allPeaks.bed > allMerge.bed


#computeMatrix scale-regions -S Primed1.bw Primed2.bw Primed3.bw Niave1.bw Niave2.bw Niave3.bw TSC1.bw TSC2.bw TSC3.bw -R allMerge.bed --binSize 1 --beforeRegionStartLength 0 --regionBodyLength 1 --afterRegionStartLength 0 --skipZeros -o AllDEInVit.mat.gz
#computeMatrix scale-regions -S EmDisc.bw Am.bw VE.bw ExMes.bw Tb.bw Gland.bw -R allMerge.bed --binSize 1 --beforeRegionStartLength 0 --regionBodyLength 1 --afterRegionStartLength 0 --skipZeros -o AllDEInViv.mat.gz



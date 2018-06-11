data<-read.xlsx("Monocytes in EAE_fem_apoefloxed.xlsx", rowNames=TRUE, sheet=3)
rownames(data)<-data[,1]
Match<-data[,1:2]
Match[,1]<-NULL
data[,1:2]<-NULL

source("http://www.bioconductor.org/biocLite.R")
biocLite("edgeR")

library(edgeR)
data[,8:15]<-NULL
group<-c(rep("AWT",4),rep("BFlox",3))
group<-factor(group)
cds<-DGEList(data,group=group)
cds$genes<-data.frame(Symbol=Match[,1])
keep<-filterByExpr(cds)
summary(keep)
cds<-cds[keep, , keep.lib.sizes=FALSE]
cds<-calcNormFactors(cds)
cds$samples
plotMD(cpm(cds,log=TRUE),column=1)
abline(h=0, col="red", lty=2, lwd=2)

points <- c(rep(0,4),rep(1,3))
colors <- c(rep("black",4),rep("red",3))
plotMDS(cds, col=colors, pch=points)
legend("topleft", legend=((group)), pch=points, col=colors)

design <- model.matrix(~0+group, data=cds$samples)
colnames(design)<-levels(cds$samples$group)
cds<-estimateDisp(cds,design,robust=TRUE)
fit<-glmQLFit(cds,design,robust=TRUE)
con<-makeContrasts(BFlox-AWT, levels=design)
qlf<-glmQLFTest(fit,contrast=con)

#####################trying to figure out why I'm getting so few sig genes
ex<-exactTest(cds,pair=c("AWT","BFlox"))
topTags(ex)
summary(decideTests(ex))




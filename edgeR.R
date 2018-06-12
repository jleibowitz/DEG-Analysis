setwd("/Users/jleibowitz/Desktop/Github/DEG-Analysis/Raw counts")
data<- read.csv(file = "mergeddata.csv", header = TRUE, sep = ",", quote= "'", dec = ".")
rownames(data)<-data[,2]
my_genes_ann3 <- read.delim("my_genes_ann2.csv", header = TRUE, sep = ",", quote= "'", dec = ".")
Match<-data.frame(my_genes_ann3[,3],my_genes_ann3[,4])
data[,1:2]<-NULL
sampleID <- read.delim("Metadata.csv", header = TRUE, sep = ",", quote= "'", dec = ".")

source("http://www.bioconductor.org/biocLite.R")
biocLite("edgeR")

library(edgeR)
group<-sampleID$cell_type##group<-c(rep("AWT",4),rep("BFlox",3))
group<-factor(group)
cds<-DGEList(data,group=group)
cds$genes<-data.frame(Symbol=Match[,2])
cpm(10, mean(cds$samples$lib.size))
keep <- rowSums(cpm(cds) > 0.81) >= 4 ##keep<-filterByExpr(cds)
summary(keep)
cds<-cds[keep, ]
cds<-calcNormFactors(cds)
cds$samples
plotMD(cpm(cds,log=TRUE),column=1)
abline(h=0, col="red", lty=2, lwd=2)

points <- c(0,0,0,0,1,1,1,1)
colors <- c("black","black","black","black","red","red","red","red")
plotMDS(cds, col=colors, pch=points)
legend("topleft", legend=((group)), pch=points, col=colors)

design <- model.matrix(~0+group, data=cds$samples)
colnames(design)<-levels(cds$samples$group)
colnames(design)<-c("Mo_R","Mo_sp")
cds<-estimateDisp(cds,design,robust=TRUE)
fit<-glmQLFit(cds,design,robust=TRUE)
con<-makeContrasts(Mo_Sp-Mo_R, levels=design)
qlf<-glmQLFTest(fit,contrast=con)

#####################trying to figure out why I'm getting so few sig genes
ex<-exactTest(cds,pair=c("Mo_R","Mo_Sp"))
t<-topTags(ex,n=1482)
summary(decideTests(ex))





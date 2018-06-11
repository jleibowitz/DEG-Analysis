setwd("/Users/jleibowitz/Desktop/Github/DEG-Analysis/Raw counts")
data1<-read.table(file="000008545469_A07.counts", header=TRUE)
filenames<-list.files(pattern="*.counts")
for(i in 1:length(list.files(pattern="*.counts"))){
  assign(paste0("data",i),read.table(file=filenames[i],header=TRUE))
}


mergeddata<-data1
for(i in 2:length(filenames)){
  mergeddata<-merge(x=mergeddata, y=get(paste0("data",i)),by.x="gene_id", by.y="gene_id",all.x=T,all.y=F)
}


mergeddata <- mergeddata[, !duplicated(colnames(mergeddata))]

mergeddata[,2]<-NULL
mergeddata[,3]<-NULL

#######ANNOTATING
source("https://bioconductor.org/biocLite.R")
biocLite("annotate")
biocLite("AnnotationDbi")
biocLite("org.Mm.eg.db")

geneLists<-mergeddata

biocLite("biomaRt")
library(biomaRt)
listEnsembl()
ensembl = useEnsembl(biomart="ensembl")
head(listDatasets(ensembl))
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")

listDatasets() 

library(biomaRt)

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

listFilters(mouse)

getBM( attributes=c("ensembl_gene_id", "mgi_symbol"), filters= "mgi_symbol", mart=mouse)
res <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol","chromosome_name",'strand','transcript_start','transcript_end'), mart = mouse)

annot.table <- data.frame("geneLists")
gene_ids <- character()
genes.table = NULL 
mart <- useMart("ensembl") 
mart <- useDataset("mmusculus_gene_ensembl", mart = mart) 
genes.table <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "description"), values= mergeddata[,1], mart= mart) 

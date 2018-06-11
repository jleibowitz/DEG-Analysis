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
source("http://www.bioconductor.org/biocLite.R")
biocLite("biomaRt")
require(biomaRt)
ensMart<-useMart("ensembl")
ensembl_ms_mart<-useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensembl_df<-getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version", "mgi_symbol","chromosome_name",'strand','transcript_start','transcript_end'), mart = ensembl_ms_mart)
my_genes=mergeddata[,1]
my_genes_ann=ensembl_df[match(my_genes, ensembl_df$ensembl_gene_id_version),]



my_genes_ann2=ensembl_df[match(test,ensembl_df$ensembl_gene_id),]

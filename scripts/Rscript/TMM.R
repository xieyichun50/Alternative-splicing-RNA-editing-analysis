####Load package
setwd("expression/")
library(edgeR)

####Read in TMM expression matrix
rawcount<-read.delim("gene_count_matrix.csv", header = TRUE)
row.names(rawcount)<-rawcount$gene_id
TMM<-rawcount[,2:ncol(rawcount)]
sampleorder<-c("BS1","BS2","BS3",
               "BS12h1","BS12h2","BS12h3",
               "BS24h1","BS24h2","BS24h3",
               "M1","M2","M3",
               "O1","O2","O3",
               "S1","S2","S3",
               "K1","K2","K3",
               "P1","P2","P3",
               "YFB1","YFB2","YFB3")
TMM<-TMM[, sampleorder]
TMM = TMM[rowSums(cpm(TMM) > 1) >= 2,]
####Read in sample group
trait<-read.delim("samples_n_reads_decribed.txt", header = FALSE)
trait<-trait[trait$V2 %in% sampleorder,]
#TMM<-TMM[,names(TMM) %in% trait$V2]
group <- factor(trait[,2])

data=as.matrix(TMM) 
y<-DGEList(counts=data,group=group)
y <- calcNormFactors(y, method = "TMM")
TMM.norm<-cpm(y)

write.table(TMM.norm, file = "gene_count_TMMmatrix.txt", sep = "\t", quote = F)

TMM.norm.log2<-log2(TMM.norm+1)
write.table(TMM.norm.log2, file = "gene_count_log2TMMmatrix.txt", sep = "\t", quote = F)

##mean
treatment<-as.data.frame(unique(trait[,1]))
names(treatment)[1]<-"treatment"
TMM.norm<-as.data.frame(TMM.norm)
TMM.norm$X<-row.names(TMM.norm)
TMM.mean<-matrix(NA, nrow = nrow(TMM.norm), ncol = nrow(treatment),
                 dimnames = list(TMM.norm$X,treatment$treatment))
TMM.mean<-as.data.frame(TMM.mean)

for (i in 1:nrow(treatment)) {
  TMM.mean[i] <-(TMM.norm[(i*3)-1]+TMM.norm[i*3]+TMM.norm[(i*3)-2])/3
  
}

write.csv(TMM.mean, file = "gene_count_matrix_TMMmean.txt", quote = FALSE)

TMM.mean.log2<-log2(TMM.mean+1)
write.table(TMM.mean.log2, file = "gene_count_log2TMMmean.txt", sep = "\t", quote = F)

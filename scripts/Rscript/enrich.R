setwd("expression/")
save.image("DEG.RData")
load("DEG.RData")
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(eoffice)
library(pheatmap)
library(pathview)
library(clusterProfiler)
filelist<-read.delim("edgeR_gene.min_reps2.min_cpm1/filelist.txt", header = TRUE)

#DEG.bk<-read.delim("DEG.allpairs.txt", header = T)
#DEG.bk<-DEG.bk[DEG.bk$change != "ns",]

DEG<-as.data.frame(matrix(NA, ncol = 7, nrow = 0))
names(DEG)<-c("Genes","sampleA","sampleB","logFC","logCPM","PValue","FDR")

for (i in 1:nrow(filelist)) {
  DEG.sub<-read.delim(file = paste0("edgeR_gene.min_reps2.min_cpm1/",filelist$file[i]), 
                      header = TRUE)
  DEG.sub$Genes<-row.names(DEG.sub)
  row.names(DEG.sub)<-1:nrow(DEG.sub)
  DEG<-rbind(DEG, DEG.sub)
}

rm(DEG.sub)
rm(i, filelist)

DEG.BK<-DEG

DEG$change<-NA
DEG$change="ns"
DEG$change[DEG$logFC< (-2) & DEG$FDR<0.05]<-"Down"
DEG$change[DEG$logFC> (2) & DEG$FDR<0.05]<-"Up"
DEG$Groups<-paste0(DEG$sampleB," > ", DEG$sampleA)
DEG.BK<-DEG

write.table(DEG, file = "DEG.allpairs.txt",
            sep = "\t", quote = F, row.names = F)

DEG<-read.delim("DEG.allpairs.txt", header = T)
##Expression table merged

##Matrix
DEG.matrix<-as.data.frame(unique(DEG$Genes))

group.order<-c("BS > BS12h","BS > BS24h","BS12h > BS24h",
               "Knot > Pri","Knot > YFB","Knot > Scl",
               "Myc > BS","Myc > BS12h","Myc > BS24h",  
               "Myc > Oidia","Myc > Scl", "Myc > Knot","Myc > Pri","Myc > YFB",
               "Oidia > BS","Oidia > Scl","Oidia > Knot",
               "Pri > YFB","YFB > BS",
               "Scl > BS","Scl > Pri","Scl > YFB")
group.order<-unique(DEG$Groups)
group.order<-group.order[order(group.order)]

treat.order<-c("BS",
               "BS12h",
               "BS24h",
               "Myc",
               "Oidia",
               "Scl",
               "Knot",
               "Pri",
               "YFB")

names(DEG.matrix)<-"Genes"

for (i in 1:length(group.order)) {
  DEG.matrix.sub<-DEG[DEG$Groups == group.order[i],c("Genes", "logFC", "change")]
  names(DEG.matrix.sub)=c("Genes", paste0("logFC.",group.order[i]), paste0("change.",group.order[i]))
  DEG.matrix<-merge(DEG.matrix, DEG.matrix.sub, by = "Genes", all.x = TRUE)
}

write.table(DEG.matrix, file = "DEG.allpairs.matrix.txt",
            sep = "\t", quote = F, row.names = F)
DEG.matrix.bk<-DEG.matrix
rm(DEG.matrix.sub)

##Filter DEGs →
DEG<-DEG[DEG$change != "ns",]

DEG.summary<-as.data.frame(xtabs(~change+sampleA+sampleB,DEG.BK))

piep<-DEG.summary %>%
  mutate(sampleA = factor(sampleA, levels = treat.order),
         sampleB = factor(sampleB, levels = treat.order)) %>%
  ggplot(aes(x="", y=Freq, fill=change))+
  geom_bar(stat = "identity", position = position_fill())+
  labs(fill = "Change")+
  scale_fill_manual(values=c("skyblue", "gray", "plum2"))+
  geom_text(aes(label = Freq), position = position_fill(vjust = 0.5), color = "white")+
  coord_polar(theta = "y", start=0)+
  facet_wrap(sampleB~sampleA, ncol = 9)+
  theme_void()

piep

f = "DEG.pie.pptx"
topptx(piep,f, width = 6.5, height = 7.8, units = "in")

write.table(DEG.summary, "DEG.summary.txt", row.names = F, quote = F, sep = "\t")

##Volcano plot
DEG.bk<-DEG.BK
DEG.BK$FDR[DEG.BK$FDR<10^(-10)]<-10^(-10)
a<-DEG.BK %>%
  mutate(sampleA = factor(sampleA, levels = treat.order),
         sampleB = factor(sampleB, levels = treat.order)) %>%
  ggplot(aes(x = logFC, y = -log10(FDR),
             colour = change))+
  geom_point(shape = 20, size = 1)+
  scale_colour_manual(values = c("skyblue", "gray", "plum2"))+
  scale_x_continuous(limits = c(-15,15), breaks = seq(-15,15,10))+
  scale_y_continuous(limits = c(0,10), breaks = seq(0,10,5))+
  labs(x = "log2(Fold change)", y = "-log10(FDR)", legend = "")+
  geom_hline(yintercept = -log10(0.05), color = "grey30", linetype = "dashed")+
  geom_vline(xintercept = 2, color = "grey30", linetype = "dashed")+
  geom_vline(xintercept = -2, color = "grey30", linetype = "dashed")+
  facet_grid(sampleB~sampleA, scales = "free", space = "free")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 0),
        legend.position = "none",
        panel.grid.minor.x = element_line(linetype = "blank"),
        strip.text = element_text(size = 12))
a

ggsave("DEG.Volcano.png", width = 6.5, height = 6.5, units = "in", dpi = 300)
ggsave("DEG.Volcano.tiff", width = 6.5, height = 6.5, units = "in", dpi = 300)
ggsave("DEG.VolcanoL.png", width = 20, height = 20, units = "in", dpi = 300)
ggsave("DEG.VolcanoL.tiff", width = 20, height = 20, units = "in", dpi = 300)

rm(piep, a, f, i)
##Gene name
gene2name<-read.delim("D:/scRNA/eggnog/gene2name.txt", header = TRUE)
##DEG freq
DEG.freq<-as.data.frame(xtabs(~Genes, DEG))
DEG.freq<-merge(DEG.freq, gene2name, by = "Genes", all.x = T)

DEG.freq.ud<-as.data.frame(xtabs(~Genes+change, DEG))
DEG.freq.ud<-merge(DEG.freq.ud, gene2name, by = "Genes", all.x = T)

#continuous changes
TMM.norm.log2<-read.delim("gene_count_log2TMMmatrix.txt",header = T)
heatmap.input<-TMM.norm.log2[row.names(TMM.norm.log2) %in% unique(DEG$Genes),]

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
heatmap.input<-as.data.frame(t(apply(heatmap.input, 1, normalize)), row.names = row.names(heatmap.input))
names(heatmap.input)<-c("BS1","BS2","BS3",
                        "BS12h1","BS12h2","BS12h3",
                        "BS24h1","BS24h2","BS24h3",
                        "Myc1","Myc2","Myc3",
                        "Oidia1","Oidia2","Oidia3",
                        "Scl1","Scl2","Scl3",
                        "Knot1","Knot2","Knot3",
                        "Pri1","Pri2","Pri3",
                        "YFB1","YFB2","YFB3")

heatp<-pheatmap(heatmap.input,
                legend = TRUE,
                scale = "none",
                #scale = "row",
                cluster_rows = TRUE, cluster_cols = FALSE,
                cutree_rows = 6, treeheight_row = 100, 
                border_color = NA,
                cellwidth = 30, cellheight = 0.1, fontsize = 24,
                angle_col = c("90"),
                width = 16, height = 16,
                show_colnames = TRUE,
                show_rownames = T,
                filename = "heatmap6.png")
#filename = "heatmap.DEG.all.png")
#filename = "heatmap.DEG.mean.png")

heatp.clust<-cbind(heatmap.input, 
                   order = as.list(heatp$tree_row["order"]),
                   cluster = cutree(heatp$tree_row, k=6))
write.table(heatp.clust, "heatmap6.clust.txt", sep = "\t", row.names = TRUE, quote = FALSE)

heatp.clust<-read.delim("heatmap6.clust.txt", header = TRUE)
heatp.clust$Genes<-row.names(heatp.clust)
DEG<-heatp.clust[,c("Genes", "cluster")]
names(DEG)[2]="Groups"
DEG$Groups<-paste0("cluster", DEG$Groups)

rm(heatp.clust)

heatp<-pheatmap(heatmap.input[row.names(heatmap.input) %in% DEG$Genes[DEG$Groups =="cluster6"],],
                legend = TRUE,
                scale = "none",
                #scale = "row",
                cluster_rows = TRUE, cluster_cols = FALSE,
                cutree_rows = 1, treeheight_row = 100, 
                border_color = NA,
                cellwidth = 30, cellheight = 0.1, fontsize = 24,
                angle_col = c("90"),
                width = 16, height = 16,
                show_colnames = TRUE,
                show_rownames = T,
                filename = "heatmap.cluster6.png")

ref= "heatmap6."

##Functional annotation
addline_format <- function(x,...){ 
  gsub('_','\n',x) 
} 

GenesKOGpair.1v1<-read.delim("C:/coprinopsis/eggnog/KOG.1v1.txt", header = TRUE)
GenesKEGGpair.1v1<-read.delim("C:/coprinopsis/eggnog/KEGG.1v1.txt", header = TRUE)
GenesGOpair.1v1<-read.delim("C:/coprinopsis/eggnog/GO.1v1.txt", header = TRUE)
Geneskopair.1v1<-read.delim("C:/coprinopsis/eggnog/KO.1v1.txt", header = TRUE)

kog2name<-read.delim("D:/3enrichment/kog2name.txt", header = TRUE)
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage\nand processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes\nand signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

ko2name<-read.delim("D:/3enrichment/ko2name.txt")
kegg2name<-read.delim("D:/3enrichment/kegg2name.txt")

kegg2ont<- read.delim("D:/3enrichment/kegglevel.AC.txt", 
                      sep = "\t", colClasses = "character")
names(kegg2ont)[2]="ID"

go2name<-read.delim("D:/3enrichment/go2name.txt", 
                    sep = "\t", colClasses = "character",
                    header = FALSE)
names(go2name)[1]="goClass"
names(go2name)[2]="goName"
names(go2name)[3]="ONTOLOGY"
go2name$ONTOLOGY<-gsub("biological_process", "Biological Process", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("molecular_function", "Molecular Function", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("cellular_component", "Cellular Component", go2name$ONTOLOGY)

DEG.bk<-DEG

#compare with start

####KOG enrich
{
  KOG.all.1<-compareCluster(Genes ~ Groups, 
                            data = DEG, 
                            fun = 'enricher',
                            TERM2GENE = GenesKOGpair.1v1,
                            TERM2NAME = kog2name,
                            pvalueCutoff = 1,
                            pAdjustMethod = "BH",
                            qvalueCutoff = 1,
                            minGSSize = 1,
                            maxGSSize = 200000)
  
  ##Plot
  plotin<-as.data.frame(KOG.all.1)
  View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  colnames(plotinsep)[names(plotinsep)=="ID"]<-"kogClass"
  
  plotdata<-merge(plotinsep, kog2name, by = c("kogClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  
  write.table(plotdata, file = paste0(ref,"KOGenrich.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  list(plotdata$Group)
  a<-ggplot(plotdata, aes(x = Groups, y = paste0(Description," (",BGnumerator,")"), size = ratio1, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave(paste0(ref,"KOG.tiff"), width = 12, height = 8, units = "in", dpi = 300)
  ggsave(paste0(ref,"KOG.png"), width = 12, height = 8, units = "in", dpi = 300)
}

library("AnnotationHub")
annotationhub<-AnnotationHub()
query(annotationhub, "Coprinopsis cinerea")
ah.ccin<-AnnotationHub()[["AH95203"]]
IDmatch<-read.delim("IDmatch.txt", header = TRUE)
DEG<-merge(DEG, IDmatch, by = "Genes", all.x = TRUE)

####KEGG enrich
{
  KEGG.all.1<-compareCluster(Genes ~ Groups, 
                             data = DEG, 
                             fun = 'enricher',
                             TERM2GENE = GenesKEGGpair.1v1,
                             TERM2NAME = kegg2name,
                             pvalueCutoff = 1,
                             pAdjustMethod = "BH",
                             qvalueCutoff = 1,
                             minGSSize = 1,
                             maxGSSize = 200000)
  
  KEGG.all.1<-compareCluster(CC130~Groups, fun = "enrichKEGG", data = DEG, 
                             organism = "cci", keyType = "kegg", pvalueCutoff = 1, 
                             pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,maxGSSize = 20000)
  ##Plot
  plotin<-as.data.frame(KEGG.all.1)
  View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotdata<-subset(plotinsep, is.na(Description) == FALSE & p.adjust < 1)
  
  kegg2ont$ID<-gsub("ko", "cci", kegg2ont$ID)
  plotdata<-merge(plotdata, kegg2ont, by = "ID", all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  
  write.table(plotdata, file = paste0(ref,"KEGG.enrich.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  plotdata$ONTOLOGY<-gsub(" ", "\n", plotdata$ONTOLOGY)
  list(plotdata$Group)
  
  plotdata1<-subset(plotdata, ratio1>0 & p.adjust< 0.2 & ONTOLOGY != "Human\nDiseases")
  plotdata1$ONTOLOGY<-gsub("Environmental\nInformation\nProcessing", 
                           "Environmental\nInformation Processing",
                           plotdata1$ONTOLOGY)
  a<-ggplot(plotdata1, aes(x = Groups, y = paste0(Description," (",BGnumerator,")"), size = ratio1, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12),
          strip.text.x = element_text(size = 11, angle = 0),
          strip.text.y = element_text(size = 9, angle = 0))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")
  
  a
  ggsave(paste0(ref,"KEGGr1.tiff"), width = 10, height = 10, units = "in", dpi = 300)
  ggsave(paste0(ref,"KEGGr1.png"), width = 10, height = 10, units = "in", dpi = 300)
  
}

####ko enrich
{
  ko.all.1<-compareCluster(Genes ~ Groups, 
                           data = DEG, 
                           fun = 'enricher',
                           TERM2GENE = Geneskopair.1v1,
                           TERM2NAME = ko2name,
                           pvalueCutoff = 1,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 1,
                           minGSSize = 1,
                           maxGSSize = 200000)
  
  ##Plot
  plotin<-as.data.frame(ko.all.1)
  View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotdata<-subset(plotinsep, is.na(Description) == FALSE & p.adjust < 1)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  
  write.table(plotdata, file = paste0(ref,"ko.enrich.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  list(plotdata$Group)
  
  plotdata1<-subset(plotdata, ratio2>0.01 & p.adjust< 0.2 )
  a<-ggplot(plotdata1, aes(x = Groups, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12),
          strip.text.x = element_text(size = 11, angle = 0),
          strip.text.y = element_text(size = 9, angle = 0))
  
  a
  ggsave(paste0(ref,"ko.tiff"), width = 18, height = 24, units = "in", dpi = 300)
  ggsave(paste0(ref,"ko.png"), width = 18, height = 24, units = "in", dpi = 300)
  
}

####GO enrich
{
  GO.all.1<-compareCluster(Genes ~ Groups, 
                           data = DEG, 
                           fun = 'enricher',
                           TERM2GENE = GenesGOpair.1v1,
                           TERM2NAME = go2name,
                           pvalueCutoff = 1,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 1,
                           minGSSize = 1,
                           maxGSSize = 2000000)
  
  GO.all.1<-compareCluster(CC130 ~ Groups, fun = "enrichGO", data = DEG, 
                           OrgDb = ah.ccin, keyType = "SYMBOL", ont = "ALL", pvalueCutoff = 1, 
                           pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,maxGSSize = 200000, 
                           readable = FALSE, pool = FALSE)
  GO.all.1.MF<-compareCluster(CC130 ~ Groups, fun = "enrichGO", data = DEG, 
                           OrgDb = ah.ccin, keyType = "SYMBOL", ont = "MF", pvalueCutoff = 1, 
                           pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,maxGSSize = 200000, 
                           readable = FALSE, pool = FALSE)
  GO.all.1.MF<-gofilter(GO.all.1.MF, level = 5)
  GO.all.1.MF<-as.data.frame(GO.all.1.MF)
  GO.all.1.MF$ONTOLOGY<-"MF"
  GO.all.1.CC<-compareCluster(CC130 ~ Groups, fun = "enrichGO", data = DEG, 
                           OrgDb = ah.ccin, keyType = "SYMBOL", ont = "CC", pvalueCutoff = 1, 
                           pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,maxGSSize = 200000, 
                           readable = FALSE, pool = FALSE)
  GO.all.1.CC<-gofilter(GO.all.1.CC, level = 5)
  GO.all.1.CC<-as.data.frame(GO.all.1.CC)
  GO.all.1.CC$ONTOLOGY<-"CC"
  GO.all.1.BP<-compareCluster(CC130 ~ Groups, fun = "enrichGO", data = DEG, 
                           OrgDb = ah.ccin, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, 
                           pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,maxGSSize = 200000, 
                           readable = FALSE, pool = FALSE)
  GO.all.1.BP<-gofilter(GO.all.1.BP, level = 5)
  GO.all.1.BP<-as.data.frame(GO.all.1.BP)
  GO.all.1.BP$ONTOLOGY<-"BP"
  
  #Plot
  plotin<-as.data.frame(GO.all.1)
  plotin<-rbind(GO.all.1.MF, GO.all.1.CC, GO.all.1.BP)
  View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotdata<-plotinsep
  plotdata$ONTOLOGY<-gsub("BP", "Biological\nProcess", plotdata$ONTOLOGY)
  plotdata$ONTOLOGY<-gsub("CC", "Cellular\nComponent", plotdata$ONTOLOGY)
  plotdata$ONTOLOGY<-gsub("MF", "Molecular\nFunction", plotdata$ONTOLOGY)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  
  write.table(plotdata, file = paste0(ref,"GOenrich.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  ##all GO
  plotdata1<-subset(plotdata, ratio2>0 & p.adjust< 0.2 )
  a<-ggplot(plotdata1,aes(x = Groups, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  
  ggsave(paste0(ref,"GOr2.tiff"), width = 15, height = 20, units = "in", dpi = 300)
  ggsave(paste0(ref,"GOr2.png"), width = 15, height = 20, units = "in", dpi = 300)
  
  plotdata1<-subset(plotdata, ratio1>0 & p.adjust< 0.2 )
  a<-ggplot(plotdata1,aes(x = Groups, y = paste0(Description," (",BGnumerator,")"), size = ratio1, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  
  ggsave(paste0(ref,"GOr1.tiff"), width = 15, height = 20, units = "in", dpi = 300)
  ggsave(paste0(ref,"GOr1.png"), width = 15, height = 20, units = "in", dpi = 300)
  
  #top 5 function
  plotdata.top5<-plotdata1 %>% group_by(Groups) %>% top_n(-10, p.adjust)
  #plotdata.top5<-plotdata1 %>% group_by(Cluster) %>% arrange(-ratio2/p.adjust) %>% slice_head(n=5)
  write.table(plotdata.top5, 
              paste0(ref,"GO.topbypadj.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
}


##pathview mapping
GenesKOpair.1v1<-read.delim("D:/scRNA/eggnog/KO.1v1.txt", header = TRUE)
DEG.sex<-merge(DEG.sex, GenesKOpair.1v1, all.x = TRUE, by = "Genes")
DEG.sex<-DEG.sex[is.na(DEG.sex$ko)==F,]
KEGG2koID<-read.delim("D:/3enrichment/kegg2koID.txt",header = TRUE)
DEG.sex<-merge(DEG.sex, KEGG2koID, by = "ko", all.x = TRUE)
pathways<-as.data.frame(unique(DEG.sex$KEGG))
names(pathways)[1]="KEGG"
pathways<-merge(pathways, kegg2name, by = "KEGG", all.x = TRUE)
names(pathways)[1]="pathways"
write.table(pathways, "pathways.txt", row.names = F, quote = F, sep = "\t")

DEG.matrix<-as.data.frame(unique(DEG$Genes))
names(DEG.matrix)[1]="Genes"

for (i in 1:nrow(DEG.sum.sub)) { 
  DEG.sub<-subset(DEG, sampleA == DEG.sum.sub$sampleA[i] & sampleB == DEG.sum.sub$sampleB[i],
                  select = c("Genes", "logFC"))
  names(DEG.sub)[2]=paste0(DEG.sum.sub$sampleA[i],".", DEG.sum.sub$sampleB[i])
  DEG.sub<-unique(DEG.sub)
  DEG.matrix<-merge(DEG.matrix, DEG.sub, all.x = TRUE, by = "Genes")
}

DEG.matrix<-merge(DEG.matrix, GenesKOpair.1v1, all.x = TRUE, by = "Genes")
DEG.matrix<-DEG.matrix[,-1]
DEG.ko.matrix<-as.data.frame(unique(DEG.matrix$ko))
names(DEG.ko.matrix)[1]="ko"

for (i in 1:(ncol(DEG.matrix)-2)){
  DEG.matrix.sub<-DEG.matrix[,c(i,9)]
  names(DEG.matrix.sub)[1]="FC"
  DEG.ko.matrix.sub<-DEG.matrix.sub %>% group_by(ko) %>% summarise(sum = sum(FC))
  names(DEG.ko.matrix.sub)[2]=names(DEG.matrix)[i]
  DEG.ko.matrix<-merge(DEG.ko.matrix, DEG.ko.matrix.sub, all.x = TRUE, by = "ko")
}
rm(DEG.matrix.sub, DEG.ko.matrix.sub)

DEG.ko.matrix<-DEG.ko.matrix[is.na(DEG.ko.matrix$ko)==F,]
row.names(DEG.ko.matrix)<-DEG.ko.matrix$ko
DEG.ko.matrix<-DEG.ko.matrix[,-1]
DEG.ko.matrix<-DEG.ko.matrix[,c(1,2,3,6,4,7,5,8)]
write.table(DEG.ko.matrix, "DEG.ko.matrix.sex.txt",
            sep = "\t", quote = F, row.names = T)
setwd("pv/")
for (i in 1:nrow(pathways)) {
  pathwayid<-pathways$pathways[i]
  pathwayid<-"ko04150"
  pathview(gene.data = DEG.ko.matrix,
           pathway.id = pathwayid, 
           species = "ko", 
           gene.idtype = "KEGG", 
           limit = list(gene = 1), 
           bins = list(gene=10), 
           multi.state = TRUE, 
           na.col="transparent", 
           out.suffix = ref)
}

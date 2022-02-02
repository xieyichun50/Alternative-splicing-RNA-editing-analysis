library(agricolae)
library(ggplot2)
library(scales)
library(dplyr)
library(eoffice)
setwd("qPCR/")
genelist<-c("CC2G_010905",
            "CC2G_005289",
            "CC2G_012760",
            "CC2G_012607",
            "CC2G_001628",
            "CC2G_012163",
            "CC2G_005434",
            "CC2G_011103"
)

treat.order<-c("Myc",
               "Oidia",
               "Scl",
               "Knot",
               "Pri",
               "YFB")

taxoncolor<-c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
              "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
              "#BC80BD", "#CCEBC5", "#FFED6F", "#D9D9D9",  "white")

##Read in RNA-seq
TMM<-read.delim("/home/yichun/RNAmodification/expression/gene_count_TMMmatrix.txt", header = T)

TMM<-TMM[row.names(TMM) %in% genelist, 10:27]
TMMfull<-as.data.frame(matrix(NA, nrow = 0, ncol = 3))
names(TMMfull)<-c("Genes","Sample","exp")

for (i in 1:nrow(TMM)) {
  TMM.sub<-as.data.frame(t(as.data.frame(as.list(TMM[i,1:ncol(TMM)]))))
  names(TMM.sub)<-"exp"
  TMM.sub$Sample<-row.names(TMM.sub)
  TMM.sub$Genes<-row.names(TMM)[i]
  TMMfull<-rbind(TMMfull,TMM.sub)
}

row.names(TMMfull)<-1:nrow(TMMfull)
rm(TMM.sub)

TMMfull$experiment<-"RNA-seq"
rm(TMM,i)
TMMfull$Sample<-gsub("YFB","Y", TMMfull$Sample)

##Read in qPCR
qpcr<-read.delim("qPCR_exp.txt", header = T)
qpcr<-unique(qpcr)
qpcr<-qpcr[is.na(qpcr$exp) == F, c("Genes","Sample","exp")]

qpcr$experiment<-"qPCR"

##merge RNA-seq and qPCR
expall<-rbind(TMMfull, qpcr)

expall$Stage<-expall$Sample
expall$Stage<-gsub("1","", expall$Stage)
expall$Stage<-gsub("2","", expall$Stage)
expall$Stage<-gsub("3","", expall$Stage)
expall$Stage<-gsub("4","", expall$Stage)
expall$Stage<-gsub("5","", expall$Stage)

expall$Stage<-gsub("M","Myc", expall$Stage)
expall$Stage<-gsub("O","Oidia", expall$Stage)
expall$Stage<-gsub("S","Scl", expall$Stage)
expall$Stage<-gsub("K","Knot", expall$Stage)
expall$Stage<-gsub("P","Pri", expall$Stage)
expall$Stage<-gsub("Y","YFB", expall$Stage)

exp.myc<-expall[expall$Stage == "Myc",]

exp.mycmean<-exp.myc %>% 
  group_by(experiment,Genes) %>% 
  summarise_at(vars(exp), list(mycmean = mean))

expall<-merge(expall, exp.mycmean, by = c("experiment", "Genes"), all.x = TRUE)
expall$exp.norm<-expall$exp/expall$mycmean

exp.summary<-expall %>% 
  group_by(experiment,Genes,Stage) %>% 
  summarise(mean = mean(exp.norm), std = sd(exp.norm), n = n())
exp.summary$n<-as.numeric(exp.summary$n)
exp.summary$se<-exp.summary$std/sqrt(exp.summary$n)
exp.summary$cv<-exp.summary$std/exp.summary$mean

##Plot
gene.name<-read.delim("genelist.txt", header = T)
gene.name$name<-gsub("cytosine deaminase-uracil phosphoribosyltransferase fusion protein",
                       paste0("cytosine deaminase-uracil","\n","phosphoribosyltransferase fusion protein"),
                     gene.name$name)
gene.name$head<-paste0(gene.name$name, "\n(", gene.name$Genes, ")")
exp.summary<-merge(exp.summary, gene.name, by = "Genes", all.x = TRUE)

a<-exp.summary %>%
  mutate(head = factor(head, levels = gene.name$head),
         Genes = factor(Genes, levels = genelist),
         Stage = factor(Stage, levels = treat.order),
         experiment = factor(experiment, levels = c("RNA-seq", "qPCR"))) %>%
  ggplot(aes(x = Stage, y = log2(mean), fill = experiment))+
  geom_bar(stat = "identity",  position = "dodge", colour = "black", width = 0.75)+
  geom_errorbar(aes(ymin=log2(mean-se), ymax=log2(mean+se), width = 0.15), position = position_dodge(0.75))+
  scale_fill_manual(values = c("grey50", "white"))+
  labs(title = "", y = "log2(Fold change)", x = NULL, colour = NULL)+
  #scale_y_continuous(breaks = seq(0, 6, 1.5))+
  geom_hline(yintercept = -1, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 0, color = "black", linetype = "solid")+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed")+
  #facet_wrap(~Genes, ncol = 3)+
  facet_wrap(~head, ncol = 2, scales = "free_y")+
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks.y = element_line(colour = "black", size = 0.5),
        axis.ticks.x = element_line(colour = "black", size = 0.5),
        axis.text.x = element_text(size = 10, colour = "black", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 10, hjust = 0.5, face = "plain"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 0),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "bottom",
        strip.text = element_text(size = 10))
a

ggsave("deaminase_exp.png", width = 6, height = 8.5, units = "in", dpi = 300)
ggsave("deaminase_exp.tiff", width = 6, height = 8.5, units = "in", dpi = 300)

f = "deaminase_exp.pptx"
topptx(a,f, width = 6, height = 9, units = "in")

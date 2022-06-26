setwd("RE_statistics")
library(dplyr)
library(tidyr)
library(ggplot2)
library(eoffice)

##general arguements
treat.order<-c("BS",
               "BS12h",
               "BS24h",
               "Myc",
               "Oidia",
               "Scl",
               "Knot",
               "Pri",
               "YFB")
taxoncolor<-c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
              "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
              "#BC80BD", "#CCEBC5", "#FFED6F", "#D9D9D9",  "white")

REcolor<-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
           "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
           "#cab2d6", "#6a3d9a", "gray75", "#b15928")

cutoff.RE=0.03 ##At least 3 % reads edited but no more than 97 % edited
cutoff.rep=2 ##At least two replicates in one stage edited
cutoff.cov=10 ##At least ten reads of the position
filelist<-read.delim("filelist.txt", header = T)

###########################################################
##Gather all editing levels and editing types
###########################################################

##Create matrix of stage.list.bk
stage.list.bk<-as.data.frame(matrix(NA, nrow = 0, ncol = 10))
names(stage.list.bk)<-c("Region","Position","Reference","Strand",
                        "AllSubs","gCoverage.q25","gMeanQ",
                        "gBaseCount.A.C.G.T.",
                        "gAllSubs","gFrequency")

RE.list<-as.data.frame(matrix(NA, nrow = 0, ncol = 11))
names(RE.list)<-c("Region","Position","Reference","Strand",
                  "AllSubs","gCoverage.q25","gMeanQ",
                  "gBaseCount.A.C.G.T.",
                  "gAllSubs","gFrequency", "Stage")

for (i in 1:length(treat.order)) {
  stage.list<-stage.list.bk
  for (j in 1:3) {
    raw.table<-read.delim(paste0("/home/yichun/RNAmodification/REDIresult", 
                               filelist[filelist$Sample==paste0(treat.order[i],j),"filelist"]), 
                          header = T)
    raw.table.filter<-raw.table[raw.table$Coverage.q25 >= cutoff.cov &
                                  raw.table$Frequency >= cutoff.RE &
                                  raw.table$Frequency <= (1-cutoff.RE) & 
                                  raw.table$gAllSubs == "-" &
                                  raw.table$gMeanQ != "-",
                                c("Region","Position","Reference","Strand",
                                  "AllSubs","gCoverage.q25","gMeanQ",
                                  "gBaseCount.A.C.G.T.",
                                  "gAllSubs","gFrequency")]
    stage.list<-rbind(stage.list, raw.table.filter)
    rm(raw.table, raw.table.filter)
  }
  stage.list.freq<-as.data.frame(xtabs(~Region+Position+AllSubs, stage.list))
  stage.list.freq<-stage.list.freq[stage.list.freq$Freq >= cutoff.rep,]
  stage.list<-unique(stage.list)
  stage.list.filter<-merge(stage.list, stage.list.freq[,c("Region","Position","AllSubs")],
                           by = c("Region","Position","AllSubs"), all.y = T)
  
  stage.list.filter$Stage<-treat.order[i]
  RE.list<-rbind(RE.list,stage.list.filter)
  rm(stage.list,stage.list.freq,stage.list.filter)
}

RE.list<-RE.list[RE.list$gMeanQ != "-" & is.na(RE.list$gMeanQ)==F,]
write.table(RE.list, "RE.feature.stages.txt",
            row.names = F, quote = F, sep = "\t")
rm(stage.list.bk)
##Obtained filtered RNA editing position list

##With the list, obtain all RE events
RE.list<-read.delim("RE.feature.stages.txt", header = T)

RE.matrix<-RE.list[,c("Region","Position","Reference","AllSubs","Strand",
                      "gCoverage.q25","gMeanQ",
                      "gBaseCount.A.C.G.T.",
                      "gAllSubs","gFrequency")]
RE.matrix<-unique(RE.matrix)
RE.matrix.bk<-RE.matrix

RE.event<-as.data.frame(matrix(NA, nrow =0, ncol = 8))
names(RE.event)<-c("Region","Position","Reference","AllSubs",
                   "Strand","Frequency","Stage","Sample")

for (i in 1:length(treat.order)) {
  for (j in 1:3) {
    raw.table<-read.delim(paste0("/home/yichun/RNAmodification/REDIresult", 
                                 filelist[filelist$Sample==paste0(treat.order[i],j),"filelist"]), 
                          header = T)
    raw.table.filter<-merge(RE.matrix.bk, 
                            raw.table[,c("Region","Position","Reference","Strand","AllSubs","Frequency")],
                            by = c("Region","Position","Reference","Strand","AllSubs"),
                            all.x = TRUE)
    ##Event long list
    raw.table.filter$Sample<-paste0(treat.order[i],j)
    raw.table.filter$Stage<-treat.order[i]
    RE.event<-rbind(RE.event, raw.table.filter[raw.table.filter$Frequency >= cutoff.RE &
                                                 raw.table.filter$Frequency <= (1-cutoff.RE) &
                                                 is.na(raw.table.filter$Frequency)==F,])
    
    ##Editing level matrix
    raw.table.filter<-raw.table.filter[,c("Region","Position","Frequency")]
    names(raw.table.filter)[names(raw.table.filter)=="Frequency"]<-paste0(treat.order[i],j)
    RE.matrix<-merge(RE.matrix, raw.table.filter,
                     by = c("Region","Position"), all.x = T)
    
    rm(raw.table, raw.table.filter)
  }

}

rm(RE.matrix.bk, filelist)
write.table(RE.event, "RE.event.txt",
            row.names = F, quote = F, sep = "\t")

write.table(RE.matrix, "RE.matrix.txt",
            row.names = F, quote = F, sep = "\t")

###########################################################
##Summarise the RNA editing
###########################################################
RE.list<-read.delim("RE.feature.stages.txt", header = T)
RE.event<-read.delim("RE.event.txt", header = T)
RE.matrix<-read.delim("RE.matrix.txt", header = T)

RE.event<-merge(RE.list, RE.event[,c("Region","Position","Stage","Sample","Frequency")], 
                by = c("Region","Position","Stage"),
                all = FALSE)
write.table(RE.event, "RE.event.T.txt",
            row.names = F, quote = F, sep = "\t")
RE.event<-read.delim("RE.event.T.txt", header = T)

##Number of RNA editing events by samples and stages
RE.event.freq.sample<-as.data.frame(xtabs(~Sample, RE.event))
names(RE.event.freq.sample)[2]="Freq.Sample"
RE.event.freq.stage<-as.data.frame(xtabs(~Stage, unique(RE.event[,c("Region","Position","Stage")])))
names(RE.event.freq.stage)[2]="Freq.Stage"
RE.event.freq<-merge(RE.event.freq.sample,
                     filelist[,c("Stage","Sample")],
                     by = c("Sample"), all.x = T)
RE.event.freq<-merge(RE.event.freq,
                     RE.event.freq.stage,
                     by = c("Stage"), all.x = T)
write.table(RE.event.freq, "RE.event.freq.txt",
            row.names = F, quote = F, sep = "\t")
rm(RE.event.freq.sample,RE.event.freq.stage)

##Barchart, Types to per million mapped reads
mappedread<-read.delim("mappedreads.txt", header = T)

RE.event2<-RE.event
RE.event2$Stage<-"All"
RE.event2<-rbind(RE.event, RE.event2)
RE.event.type.freq<-as.data.frame(xtabs(~AllSubs+Stage,RE.event2))
RE.event.type.freq<-merge(RE.event.type.freq, mappedread,
                          by = "Stage", all.x = T)
RE.event.type.freq$Freq.MR<-RE.event.type.freq$Freq/RE.event.type.freq$MR*1000000

subtype.stage.summary.plot<-RE.event.type.freq %>%
  mutate(Stage = factor(Stage, levels = c(treat.order,"All"))) %>%
  ggplot(aes(x = Stage, y = Freq.MR, 
             fill = factor(AllSubs, 
                           levels = c("AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"))))+
  geom_bar(stat = "identity", position = "stack", width = 0.5)+
  scale_y_continuous(breaks = seq(0,6,1))+
  scale_fill_manual(values = REcolor,
                    labels = c("A?C", "A?G", "A?T", "C?A", "C?G", "C?T", "G?A", "G?C", "G?T", "T?A", "T?C", "T?G"))+
  labs(y = "Number of RNA editing events \n per million mapped reads", 
       x = "",
       title = "",
       fill = "")+
  theme(axis.title = element_text(size = 9), 
        axis.text = element_text(size = 9, colour = "black"), 
        axis.text.x = element_text(size = 9, angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(size = 9), 
        plot.title = element_text(size = 0, hjust = 0.5), 
        legend.text = element_text(size = 9), 
        legend.title = element_text(size = 0),
        legend.key.size = unit(0.15, "in"),
        axis.line = element_line(size = 0.5),
        #panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        panel.background = element_rect(fill = NA))
subtype.stage.summary.plot
ggsave("RE.subtype.stage.mmapped.png", width = 3, height = 3, units = "in", dpi = 300)
ggsave("RE.subtype.stage.mmapped.tiff", width = 3, height = 3, units = "in", dpi = 300)
f = "RE.subtype.stage.mmapped.pptx"
topptx(subtype.stage.summary.plot, f, width = 3.2, height = 3, units = "in")
rm(subtype.stage.summary.plot, RE.event.type.freq, mappedread, RE.event.freq )

##Violin plot, stage 
library(ggpubr)
#violin<-RE.event2[RE.event2$Stage !="All",c("Stage","Frequency")] %>% 
violin<-RE.event2[,c("Stage","Frequency")] %>% 
  mutate(Stage = factor(Stage, levels = c(treat.order,"All"))) %>%
  ggplot(aes(x=Stage, y = Frequency, fill = Stage))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(position = position_dodge(width = 0.9), alpha = 0.5, width = 0.5, fill = "white")+
  scale_fill_manual(values=taxoncolor)+
  labs(x = "", y = "Editing level")+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  stat_compare_means(label.y = 1.05)+
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 9, angle = 45,  hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 0, hjust = 0.5),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        axis.line = element_line(size = 0.5),
        panel.background = element_rect(fill = NA),
        legend.position = "none")
violin
my_comparisons<-compare_means(Frequency ~ Stage, data = RE.event2[RE.event2$Stage != "All",])
write.table(my_comparisons, "RE.stage.violin.KWcomparison.txt",
            row.names = F, quote = F, sep = '\t')
f = "RE.stage.violin.pptx"
topptx(violin, f, width = 3, height = 3, units = "in")

RE.event2 %>% group_by(Stage) %>% summarise(mean = mean(Frequency))
rm(RE.event2,f,violin,subtype.stage.summary.plot,my_comparisons)

##Violin plot, type 
violin<-RE.event[,c("AllSubs","Frequency")] %>% 
  ggplot(aes(x=AllSubs, y = Frequency,
             fill = factor(AllSubs, 
                           levels = c("AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"))))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(position = position_dodge(width = 0.9), alpha = 0.5, width = 0.5, fill = "white")+
  labs(x = "", y = "Editing level")+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  scale_x_discrete(labels = c("A?C", "A?G", "A?T", "C?A", "C?G", "C?T", "G?A", "G?C", "G?T", "T?A", "T?C", "T?G"))+
  scale_fill_manual(values = REcolor,
                    labels = c("A?C", "A?G", "A?T", "C?A", "C?G", "C?T", "G?A", "G?C", "G?T", "T?A", "T?C", "T?G"))+
  stat_compare_means(label.y = 1.05)+
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 9, angle = 45,  hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 0, hjust = 0.5),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        axis.line = element_line(size = 0.5),
        panel.background = element_rect(fill = NA),
        legend.position = "none")
violin
my_comparisons<-compare_means(Frequency ~ AllSubs, data = RE.event)
write.table(my_comparisons, "RE.type.violin.KWcomparison.txt",
            row.names = F, quote = F, sep = '\t')
f = "RE.type.violin.pptx"
topptx(violin, f, width = 3, height = 3, units = "in")
rm(RE.event2,f,violin,subtype.stage.summary.plot,my_comparisons)

##histogram
freq.hist<-as.data.frame(xtabs(~Frequency,RE.event))
freq.hist$Frequency<-as.numeric(as.character(freq.hist$Frequency))
freq.hist$percent<-freq.hist$Freq/sum(freq.hist$Freq)

freq.hist.plot<-ggplot(freq.hist, 
                       aes(x = Frequency, y = percent))+
  geom_bar(stat = "identity", position = "dodge", width = 0.01)+
  scale_fill_manual(values = c("black"))+
  #scale_x_discrete(limits = c("3-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"))+
  scale_y_continuous(labels = scales::percent, breaks = seq(0,0.25,0.05))+
  labs(y = "Percentage", x = "Editing levels")+
  theme(axis.title = element_text(size = 9), 
        axis.text = element_text(size = 9, colour = "black"), 
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size = 9), 
        plot.title = element_text(size = 0, hjust = 0.5), 
        legend.text = element_text(size = 9), 
        legend.title = element_text(size = 9),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5),
        #panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        panel.background = element_rect(fill = NA))
freq.hist.plot
ggsave("RE.hist.tiff", width = 3, height = 3, units = "in", dpi = 300)
ggsave("RE.hist.png", width = 3, height = 3, units = "in", dpi = 300)

f = "RE.hist.pptx"
topptx(freq.hist.plot, f, width = 3, height = 3, units = "in")
rm(freq.hist.plot, f, freq.hist)

##Feed to snpEff
library(stringr)
snpEff.vcf<-RE.matrix[,c("Region","Position","AllSubs")]
names(snpEff.vcf)<-c("CHROM","POS","AllSubs")
snpEff.vcf$REF<-str_sub(snpEff.vcf$AllSubs,1,1)
snpEff.vcf$ALT<-str_sub(snpEff.vcf$AllSubs,2,2)
snpEff.vcf$AllSubs<-"."
names(snpEff.vcf)[names(snpEff.vcf)=="AllSubs"]<-"ID"
snpEff.vcf$QUAL<-""
snpEff.vcf$FILTER<-""
snpEff.vcf$INFO<-""
write.table(snpEff.vcf, "RE.snpEff.vcf",
            row.names = F, quote = F, sep = '\t')

##Get snpEff results
snpEff.vcf<-read.delim("RE.snpEff.anno.filter.vcf", header = T)
snpEff.vcf<-snpEff.vcf[,c(1,2,8)]
names(snpEff.vcf)<-c("Region","Position","INFO")
snpEff.vcf<-separate(snpEff.vcf, INFO, 
                     c("ALT","Impact","Effect",
                       "Genes","Genes1","stype",
                       "Transcript","ptype",
                       "Exon","cds.change",
                       "aa.change","cDNA.position","cds.position","aa.position",
                       "Distance"), sep = "\\|")
snpEff.vcf<-snpEff.vcf[,c("Region","Position",
                          "Impact","Effect",
                          "Genes",
                          "Transcript",
                          "Exon",
                          "cds.change","cDNA.position","cds.position",
                          "aa.change","aa.position",
                          "Distance")]
snpEff.vcf$Gene.region<-NA
snpEff.vcf$Gene.region[snpEff.vcf$aa.position != ""]<-"CDS"
snpEff.vcf$Gene.region[grep("3_prime_UTR",snpEff.vcf$Impact)]<-"3'-UTR"
snpEff.vcf$Gene.region[grep("5_prime_UTR",snpEff.vcf$Impact)]<-"5'-UTR"
#snpEff.vcf$Gene.region[grep("upstream",snpEff.vcf$Impact)]<-"Upstream"
#snpEff.vcf$Gene.region[grep("downstream",snpEff.vcf$Impact)]<-"Downstream"
snpEff.vcf$Gene.region[grep("stream",snpEff.vcf$Impact)]<-"Intergenic"

Impact.freq<-as.data.frame(xtabs(~Impact, snpEff.vcf))
write.table(Impact.freq, "Impact.freq.txt",
            row.names = F, quote = F, sep = '\t')

genelist<-read.delim("/home/yichun/RNAmodification/genome/genelist", header = F)
genelist<-unique(genelist[,c(1,3)])
names(genelist)<-c("Transcript", "Genename")
genelist$Transcript<-gsub("ID=", "", genelist$Transcript)
genelist$Genename<-gsub("product=", "", genelist$Genename)
snpEff.vcf<-merge(snpEff.vcf, genelist,
                  by ="Transcript", all.x = T)
write.table(snpEff.vcf, "RE.feature.txt", 
            quote = F, row.names = F, sep = "\t")
snpEff.vcf<-read.delim("RE.feature.txt", header = T)

##Gene region and RE type/RE level
Generegion.order<-c(#"Upstream",
                    "5'-UTR",
                    "CDS",
                    "3'-UTR",
                    #"Downstream",
                    "Intergenic")
##Bar chart
RE.generegion.type.freq<-merge(RE.matrix[,c("Region","Position","AllSubs")],
                               snpEff.vcf[,c("Region","Position","Gene.region")],
                               by = c("Region","Position"),
                               all = T)
RE.generegion.type.freq<-as.data.frame(xtabs(~AllSubs+Gene.region, RE.generegion.type.freq))
RE.generegion.type.freq %>% group_by(Gene.region) %>% summarise(sum = sum(Freq))
barchart<-RE.generegion.type.freq %>%
  mutate(Gene.region = factor(Gene.region, levels = Generegion.order)) %>%
  ggplot(aes(x = Gene.region, y = Freq, 
             fill = factor(AllSubs, 
                           levels = c("AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"))))+
  geom_bar(stat = "identity", position = "stack", width = 0.5)+
  scale_y_continuous(breaks = seq(0,120,20))+
  scale_fill_manual(values = REcolor,
                    labels = c("A?C", "A?G", "A?T", "C?A", "C?G", "C?T", "G?A", "G?C", "G?T", "T?A", "T?C", "T?G"))+
  labs(y = "Number of RNA editing site", 
       x = "",
       title = "",
       fill = "")+
  theme(axis.title = element_text(size = 9), 
        axis.text = element_text(size = 9, colour = "black"), 
        axis.text.x = element_text(size = 9, angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(size = 9), 
        plot.title = element_text(size = 0, hjust = 0.5), 
        legend.text = element_text(size = 9), 
        legend.title = element_text(size = 0),
        legend.key.size = unit(0.15, "in"),
        axis.line = element_line(size = 0.5),
        #panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        panel.background = element_rect(fill = NA))
barchart
f = "RE.generegion.subtype.pptx"
topptx(barchart, f, width = 2.5, height = 3, units = "in")
rm(barchart, RE.generegion.type.freq)

##stage and gene region
RE.generegion.stage.freq<-merge(RE.list[,c("Region","Position","Stage")],
                               snpEff.vcf[,c("Region","Position","Gene.region")],
                               by = c("Region","Position"),
                               all = T)
RE.generegion.stage.freq<-as.data.frame(xtabs(~Stage+Gene.region, RE.generegion.stage.freq))
RE.generegion.stage.freq %>% group_by(Stage) %>% summarise(sum = sum(Freq))
barchart<-RE.generegion.stage.freq %>%
  mutate(Gene.region = factor(Gene.region, levels = Generegion.order),
         Stage = factor(Stage, levels = treat.order)) %>%
  ggplot(aes(x = Stage, y = Freq, 
             fill = Gene.region))+
  geom_bar(stat = "identity", position = "stack", width = 0.5)+
  scale_y_continuous(breaks = seq(0,180,30))+
  scale_fill_manual(values = taxoncolor)+
  labs(y = "Number of RNA editing site", 
       x = "",
       title = "",
       fill = "")+
  theme(axis.title = element_text(size = 9), 
        axis.text = element_text(size = 9, colour = "black"), 
        axis.text.x = element_text(size = 9, angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(size = 9), 
        plot.title = element_text(size = 0, hjust = 0.5), 
        legend.text = element_text(size = 9), 
        legend.title = element_text(size = 0),
        legend.key.size = unit(0.15, "in"),
        axis.line = element_line(size = 0.5),
        #panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        panel.background = element_rect(fill = NA))
barchart
f = "RE.generegion.stage.pptx"
topptx(barchart, f, width = 3, height = 2.8, units = "in")
rm(barchart,RE.generegion.stage.freq)

##Violin plot 
RE.event2<-merge(RE.event[,c("Region","Position","Frequency")], 
                 snpEff.vcf[,c("Region","Position","Gene.region")],
                 by = c("Region","Position"), all = T)
violin<-RE.event2[,c("Gene.region","Frequency")] %>% 
  mutate(Gene.region = factor(Gene.region, levels = Generegion.order)) %>%
  ggplot(aes(x=Gene.region, y = Frequency, fill = Gene.region))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(position = position_dodge(width = 0.9), alpha = 0.5, width = 0.5, fill = "white")+
  scale_fill_manual(values=taxoncolor)+
  labs(x = "", y = "Editing level")+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  stat_compare_means(label.y = 1.05)+
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 9, angle = 45,  hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 0, hjust = 0.5),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        axis.line = element_line(size = 0.5),
        panel.background = element_rect(fill = NA),
        legend.position = "none")
violin
my_comparisons<-compare_means(Frequency ~ Gene.region, data = RE.event2)
write.table(my_comparisons, "RE.Generegion.violin.KWcomparison.txt",
            row.names = F, quote = F, sep = '\t')
f = "RE.generegion.violin.pptx"
topptx(violin, f, width = 2, height = 3, units = "in")

RE.event2 %>% group_by(Gene.region) %>% summarise(mean = mean(Frequency))
rm(RE.event2,f,violin)

###########################################################
##Functional annotation
###########################################################
addline_format <- function(x,...){ 
  gsub('_','\n',x) 
} 

GenesKOGpair.1v1<-read.delim("/home/yichun/RNAmodification/genome/eggnog/KOG.1v1.txt", header = TRUE)
GenesKEGGpair.1v1<-read.delim("/home/yichun/RNAmodification/genome/eggnog/KEGG.1v1.txt", header = TRUE)
GenesGOpair.1v1<-read.delim("/home/yichun/RNAmodification/genome/eggnog/GO.1v1.txt", header = TRUE)
Geneskopair.1v1<-read.delim("/home/yichun/RNAmodification/genome/eggnog/ko.1v1.txt", header = TRUE)

kog2name<-read.delim("/home/yichun/3enrichment/kog2name.txt", header = TRUE)
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage\nand processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes\nand signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

ko2name<-read.delim("/home/yichun/3enrichment/ko2name.txt")
kegg2name<-read.delim("/home/yichun/3enrichment/kegg2name.txt")

kegg2ont<- read.delim("/home/yichun/3enrichment/kegglevel.AC.txt", 
                      sep = "\t", colClasses = "character")
names(kegg2ont)[2]="ID"

go2name<-read.delim("/home/yichun/3enrichment/go2name.txt", 
                    sep = "\t", colClasses = "character",
                    header = FALSE)
names(go2name)[1]="goClass"
names(go2name)[2]="goName"
names(go2name)[3]="ONTOLOGY"
go2name$ONTOLOGY<-gsub("biological_process", "Biological Process", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("molecular_function", "Molecular Function", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("cellular_component", "Cellular Component", go2name$ONTOLOGY)

##clusterprofiler
library(clusterProfiler)
library(ggplot2)
library(tidyr)
library(dplyr)

ref="RE.stage."

RE.event2<-merge(RE.list[,c("Region","Position","Stage")], 
                 snpEff.vcf[snpEff.vcf$Gene.region != "Intergenic",
                            c("Region","Position","Genes")],
                 by = c("Region","Position"), all.y = T)
names(RE.event2)[names(RE.event2)=="Stage"]="Groups"
####KOG enrich
{
  KOG.all.1<-compareCluster(Genes ~ Groups, 
                            data = RE.event2, 
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
  a<-plotdata %>%
    mutate(Groups = factor(Groups, levels = treat.order)) %>%
    ggplot(aes(x = Groups,
               y = paste0(Description," (",BGnumerator,")"),
               size = ratio2, colour = Count))+
    labs(title = "", size = "Gene ratio", colour = "Count", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,max(plotdata$Count)))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 9), 
          axis.text.x = element_text(size = 9), 
          axis.text.y = element_text(size = 9),
          legend.title = element_text(size = 9), 
          plot.title = element_text(size = 9))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 9))
  a
  ggsave(paste0(ref,"KOG.tiff"), width = 10, height = 6, units = "in", dpi = 300)
  ggsave(paste0(ref,"KOG.png"), width = 10, height = 6, units = "in", dpi = 300)
}

####KEGG enrich
IDmatch<-read.delim("IDmatch.txt", header = TRUE)
RE.event2<-merge(RE.event2, IDmatch[,c("Genes","CC130")], by = "Genes", all.x = TRUE)

{
  KEGG.all.1<-compareCluster(Genes ~ Groups, 
                             data = RE.event2, 
                             fun = 'enricher',
                             TERM2GENE = GenesKEGGpair.1v1,
                             TERM2NAME = kegg2name,
                             pvalueCutoff = 1,
                             pAdjustMethod = "BH",
                             qvalueCutoff = 1,
                             minGSSize = 1,
                             maxGSSize = 200000)
  
  KEGG.all.1<-compareCluster(CC130~Groups, fun = "enrichKEGG", data = RE.event2, 
                             organism = "cci", keyType = "kegg", pvalueCutoff = 1, 
                             pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 10,maxGSSize = 20000)
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
  
  plotdata1<-subset(plotdata, ratio1>0 & p.adjust< 1 & ONTOLOGY != "Human\nDiseases" & is.na(ONTOLOGY)==F)
  plotdata1$ONTOLOGY<-gsub("Environmental\nInformation\nProcessing", 
                           "Environmental\nInformation Processing",
                           plotdata1$ONTOLOGY)
  a<-plotdata1 %>%
    mutate(Groups = factor(Groups, levels = treat.order)) %>%
    ggplot(aes(x = Groups, 
               y = paste0(Description," (",BGnumerator,")"), 
               size = ratio2, colour = Count))+
    labs(title = "", size = "Gene ratio", colour = "Count", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,max(plotdata1$Count)))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 9), 
          axis.text.y = element_text(size = 9),
          legend.title = element_text(size = 9), 
          plot.title = element_text(size = 9),
          strip.text.x = element_text(size = 9, angle = 0),
          strip.text.y = element_text(size = 9, angle = 0))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")
  
  a
  ggsave(paste0(ref,"KEGGr1.tiff"), width = 10, height = 8, units = "in", dpi = 300)
  ggsave(paste0(ref,"KEGGr1.png"), width = 10, height = 8, units = "in", dpi = 300)
  
}

####ko enrich
{
  ko.all.1<-compareCluster(Genes ~ Groups, 
                           data = RE.event2, 
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
  
  plotdata1<-subset(plotdata, ratio1>0 & p.adjust< 1)
  a<-plotdata1 %>%
    mutate(Groups = factor(Groups, levels = treat.order)) %>%
    ggplot(aes(x = Groups, y = Description, size = ratio1, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 9), 
          axis.text.y = element_text(size = 9),
          legend.title = element_text(size = 9), 
          plot.title = element_text(size = 9),
          strip.text.x = element_text(size = 9, angle = 0),
          strip.text.y = element_text(size = 9, angle = 0))
  
  a
  ggsave(paste0(ref,"ko.tiff"), width = 9, height = 12, units = "in", dpi = 300)
  ggsave(paste0(ref,"ko.png"), width = 9, height = 12, units = "in", dpi = 300)
  
}

###########################################################
##Other calculation
###########################################################
path.order<-c("Spore\ngermination",
              "Oidiation",
              "Sclerotia\nformation",
              "Fruiting",
              "Sporulation")
RE.list<-read.delim("RE.feature.stages.txt", header = T)
RE.matrix<-read.delim("RE.matrix.txt", header = T)
RE.event<-read.delim("RE.event.T.txt", header = T)
snpEff.vcf<-read.delim("RE.feature.txt", header = T)

RE.list.freq<-as.data.frame(xtabs(~Region+Position, RE.list))
RE.list.freq<-RE.list.freq[RE.list.freq$Freq > 0,]

##Stage specific
RE.stage.spc<-merge(snpEff.vcf, RE.list.freq[RE.list.freq$Freq ==1,],
                    by = c("Region","Position"), all = F)
nrow(RE.stage.spc) #207/286
RE.stage.spc<-merge(RE.list, RE.list.freq[RE.list.freq$Freq ==1,],
                    by = c("Region","Position"), all = F)
nrow(RE.stage.spc) #207/286
xtabs(~Stage, RE.stage.spc)

##DEG
DEG<-read.delim("/home/yichun/RNAmodification/expression/DEG.allpairs.txt", header = T)
DEG<-DEG[DEG$change != "ns",]
RE.DEG<-snpEff.vcf[snpEff.vcf$Genes %in% DEG$Genes & snpEff.vcf$Gene.region != "Intergenic",]
nrow(RE.DEG) #125 DEG gene
nrow(snpEff.vcf[snpEff.vcf$Gene.region != "Intergenic",]) #169 Gene region

##codon position  nrow = 82
RE.coding<-snpEff.vcf[snpEff.vcf$Gene.region == "CDS",]
RE.coding<-separate(RE.coding, cds.position, 
                    c("cds.position", "transcript.length"),
                    sep = "/", convert = T)
RE.coding$codon.position<-RE.coding$cds.position%%3
RE.coding$codon.position[RE.coding$codon.position == 0] <- 3
RE.codon.freq<-as.data.frame(xtabs(~codon.position, RE.coding))
RE.codon.freq$percent<-RE.codon.freq$Freq/sum(RE.codon.freq$Freq)
freq.hist.plot<-ggplot(RE.codon.freq, 
                       aes(x = codon.position, y = percent))+
  geom_bar(stat = "identity", position = "dodge", width = 0.75, colour = "black")+
  scale_fill_manual(values = c("gray50"))+
  scale_y_continuous(labels = scales::percent, breaks = seq(0,0.6,0.1))+
  labs(y = "Percentage", x = "Codon position")+
  theme(axis.title = element_text(size = 9), 
        axis.text = element_text(size = 9, colour = "black"), 
        axis.text.x = element_text(size = 9), 
        axis.text.y = element_text(size = 9), 
        plot.title = element_text(size = 0, hjust = 0.5), 
        legend.text = element_text(size = 9), 
        legend.title = element_text(size = 9),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5),
        #panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        panel.background = element_rect(fill = NA))
freq.hist.plot

f = "RE.codon.hist.pptx"
topptx(freq.hist.plot, f, width = 2, height = 2, units = "in")
rm(freq.hist.plot, f, RE.codon.freq)

#amino acid change
RE.coding$aa.change1<-gsub("\\d","", RE.coding$aa.change)
RE.coding$aa1<-str_sub(RE.coding$aa.change1,3,5)
RE.coding$aa2<-str_sub(RE.coding$aa.change1,6,8)
RE.coding$aa1<-gsub("Ter","*", RE.coding$aa1)
nrow(RE.coding[RE.coding$aa1 == RE.coding$aa2,]) #35 synonymous_variant
RE.aa.freq<-as.data.frame(xtabs(~aa1+aa2, RE.coding[RE.coding$aa1 != RE.coding$aa2,]))
RE.aa.freq<-RE.aa.freq[RE.aa.freq$Freq>0,]
RE.aa.freq$change<-paste0(RE.aa.freq$aa1, "?",RE.aa.freq$aa2)
RE.aa.freq<-RE.aa.freq[order(RE.aa.freq$Freq, decreasing = T),]
aa.order<-RE.aa.freq$change

freq.hist.plot<-RE.aa.freq %>%
  mutate(change = factor(change, levels = aa.order)) %>% 
  ggplot(aes(x = change, y = Freq))+
  geom_bar(stat = "identity", position = "dodge", width = 0.75, colour = "black")+
  scale_fill_manual(values = c("gray50"))+
  scale_y_continuous(breaks = seq(0,5,1))+
  labs(y = "Number of editing sites", x = "")+
  theme(axis.title = element_text(size = 9), 
        axis.text = element_text(size = 9, colour = "black"), 
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size = 9), 
        plot.title = element_text(size = 0, hjust = 0.5), 
        legend.text = element_text(size = 9), 
        legend.title = element_text(size = 9),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5),
        #panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        panel.background = element_rect(fill = NA))

freq.hist.plot

f = "RE.aa.hist.pptx"
topptx(freq.hist.plot, f, width = 5.6, height = 2.3, units = "in")
rm(freq.hist.plot, f, RE.aa.freq)

##RE similarity
library(pheatmap)
heatmap.input<-RE.matrix
row.names(heatmap.input)<-paste0(heatmap.input$Region, "_", heatmap.input$Position)
heatmap.input<-heatmap.input[,c(11:37)]
heatmap.input[is.na(heatmap.input)==F]<-1
heatmap.input[is.na(heatmap.input)]<-0

heatp<-pheatmap(heatmap.input,
                legend = TRUE,
                scale = "none",
                #scale = "row",
                cluster_rows = TRUE, cluster_cols = TRUE,
                treeheight_row = 0, treeheight_col = 50,
                border_color = NA,
                cellwidth = 12, cellheight = 0.8, fontsize = 9,
                angle_col = c("45"),
                width = 6, height = 6,
                show_colnames = TRUE,
                show_rownames = F,
                na_col = "gray",
                filename = "RE.level.heatmap.png")

##venn diagram
library("VennDiagram")
RE.venn<-RE.list
RE.venn$location<-paste0(RE.venn$Region,"_",RE.venn$Position)

venn.diagram(x = list(`Sclerotia\nformation` = RE.venn[RE.venn$Stage %in% c("Scl"),"location"],
                      `Spore\ngermination` = RE.venn[RE.venn$Stage %in% c("BS", "BS12h", "BS24h"),"location"], 
                      `Mycelia-type` = RE.venn[RE.venn$Stage %in% c("Myc","Oidia","Knot"),"location"],
                      `FB-type` = RE.venn[RE.venn$Stage %in% c("Pri", "YFB"),"location"]),
             margin = 0.15, cat.dist = 0.25,
             fontfamily = "arial", sub.fontfamily = "arial", main.fontfamily = "arial", cat.fonfamily = "arial",
             filename = "RE.share.venn1.png", imagetype = "png" , 
             units = "in", height = 4, width = 4, resolution = 300,
             fill = c("skyblue", "plum2", "burlywood2", "darkgray"), na = "remove")

venn.diagram(x = list(`BS` = RE.venn[RE.venn$Stage %in% c("BS"),"location"], 
                      `BS12h` = RE.venn[RE.venn$Stage %in% c("BS12h"),"location"],
                      `BS24h` = RE.venn[RE.venn$Stage %in% c("BS24h"),"location"],
                      `Pri` = RE.venn[RE.venn$Stage %in% c("Pri"),"location"],
                      `YFB` = RE.venn[RE.venn$Stage %in% c("YFB"),"location"]),
             margin = 0.15, cat.dist = 0.25,
             fontfamily = "arial", sub.fontfamily = "arial", main.fontfamily = "arial", cat.fonfamily = "arial",
             filename = "RE.share.venn2.png", imagetype = "png" , 
             units = "in", height = 4, width = 4, resolution = 300,
             fill = c("skyblue", "plum2", "burlywood2", "darkgray", "tan1"), na = "remove")

venn.diagram(x = list(`BS` = RE.venn[RE.venn$Stage %in% c("BS"),"location"], 
                      `BS12h` = RE.venn[RE.venn$Stage %in% c("BS12h"),"location"],
                      `BS24h` = RE.venn[RE.venn$Stage %in% c("BS24h"),"location"]),
             margin = 0.15, cat.dist = 0.25,
             fontfamily = "arial", sub.fontfamily = "arial", main.fontfamily = "arial", cat.fonfamily = "arial",
             filename = "RE.share.venn5.png", imagetype = "png" , 
             units = "in", height = 4, width = 4, resolution = 300,
             fill = c("skyblue", "plum2", "burlywood2"), na = "remove")

venn.diagram(x = list(`Pri` = RE.venn[RE.venn$Stage %in% c("Pri"),"location"],
                      `YFB` = RE.venn[RE.venn$Stage %in% c("YFB"),"location"]),
             margin = 0.15, cat.dist = 0.25,
             fontfamily = "arial", sub.fontfamily = "arial", main.fontfamily = "arial", cat.fonfamily = "arial",
             filename = "RE.share.venn4.png", imagetype = "png" , 
             units = "in", height = 4, width = 4, resolution = 300,
             fill = c("skyblue", "plum2"), na = "remove")

venn.diagram(x = list(`Myc` = RE.venn[RE.venn$Stage %in% c("Myc"),"location"],
                      `Oidia` = RE.venn[RE.venn$Stage %in% c("Oidia"),"location"],
                      `Knot` = RE.venn[RE.venn$Stage %in% c("Knot"),"location"]),
             margin = 0.15, cat.dist = 0.1,
             fontfamily = "arial", sub.fontfamily = "arial", main.fontfamily = "arial", cat.fonfamily = "arial",
             filename = "RE.share.venn3.png", imagetype = "png" , 
             units = "in", height = 4, width = 4, resolution = 300,
             fill = c("skyblue", "plum2", "burlywood2"), na = "remove")

##RE level + gene age + dNdS + stage
setwd("/home/yichun/RNAmodification/RE_statistics")
treat.order<-c("Scl",
               "Oidia",
               "BS",
               "BS12h",
               "BS24h",
               "Myc",
               "Knot",
               "Pri",
               "YFB")

RE.event<-read.delim("RE.event.T.txt", header = T)
snpEff.vcf<-read.delim("RE.feature.txt", header = T)
species="Coprinopsis_cinerea_A43mutB43mut_pab1-1_326"
refspecies="Coprinopsis_sclerotiger"
genes.PS<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/blast/",
                            species,"/e10-5aa30/", 
                            species,".proteins.fa.PSfinal.txt"), header = T)
names(genes.PS)[names(genes.PS)=="Gene"]="Genes"
genes.NS<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/blast/",
                            species,
                            "/dnds/Coprinopsis_cinerea_VS_",
                            refspecies,".NSfinal.txt"))
names(genes.NS)[names(genes.NS)=="dNdS"]<-"NS"
genes.NS<-separate(genes.NS, "Genes", c("Genes"), sep = "T")
genes.hour<-genes.NS[is.na(genes.NS$NS)==F & genes.NS$NS<=2,c("Genes","NS")]

RE.phylo<-merge(RE.event[,c("Region","Position","Stage","Frequency")],
                snpEff.vcf[,c("Region","Position","Impact","Genes", "Gene.region")],
                by = c("Region","Position"),
                all.x = T)
RE.phylo<-merge(RE.phylo, genes.PS, by="Genes", all.x = T)
RE.phylo<-merge(RE.phylo, genes.NS[,c("Genes","NS")], by="Genes", all.x = T)
RE.phylo<-RE.phylo[is.na(RE.phylo$PS)==F,]
RE.phylo$PS<-paste0("PS",RE.phylo$PS)
RE.phylo$Mtype="RNA editing"

tab.size<-as.data.frame(xtabs(~PS, unique(RE.phylo[,c("Genes","PS")])))
genes.PS.stat<-as.data.frame(xtabs(~PS, genes.PS))
genes.PS.stat$PS<-paste0("PS",genes.PS.stat$PS)
tab.size<-merge(tab.size,genes.PS.stat, by = "PS", all.x = T)
tab.size$dotsize<-tab.size$Freq.x/tab.size$Freq.y
RE.phylo<-merge(RE.phylo, tab.size[,c("PS","dotsize")], by = "PS", all.x = T)

RE.phylo$Gene.region[RE.phylo$Impact == "synonymous_variant"]<-"Synonymous"
RE.phylo$Gene.region[RE.phylo$Impact != "synonymous_variant" & RE.phylo$Gene.region == "CDS"]<-"Non-synonymous"

write.table(RE.phylo, "RE.phylo.txt", sep = "\t", quote = F, row.names = F)

a<-RE.phylo %>%
  mutate(Stage = factor(Stage, levels = treat.order),
         PS = factor(PS, levels = paste0("PS",1:12))) %>%
  ggplot(aes(x = Frequency, y = NS, size = dotsize, alpha = dotsize))+
  geom_point(shape = 1)+
  scale_colour_manual(values = "black")+
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,1))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,1))+
  labs(x = "Editing levels", y = "dN/dS", legend = "")+
  facet_grid(PS~., scales = "free", space = "free")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 0),
        legend.position = "none",
        panel.grid.minor.x = element_line(linetype = "blank"),
        strip.text = element_text(size = 12))
a

ggsave("RE.geneage.png", width = 3, height = 6.5, units = "in", dpi = 300)
ggsave("RE.geneage.tiff", width = 3, height = 6.5, units = "in", dpi = 300)

a<-RE.phylo %>%
  mutate(Stage = factor(Stage, levels = treat.order),
         PS = factor(PS, levels = paste0("PS",1:12)),
         Gene.region = factor(Gene.region, 
                              levels = c("5'-UTR","Synonymous","Non-synonymous","3'-UTR","Intergenic"))) %>%
  ggplot(aes(x = Frequency, y = NS, size = dotsize, alpha = dotsize))+
  geom_point(shape = 16)+
  scale_colour_manual(values = "black")+
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,1))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,1))+
  labs(x = "Editing levels", y = "dN/dS", legend = "")+
  facet_grid(PS~Gene.region, scales = "free", space = "free")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 0),
        legend.position = "none",
        panel.grid.minor.x = element_line(linetype = "blank"),
        strip.text = element_text(size = 12))
a

ggsave("RE.geneage.type.png", width = 6.5, height = 6.5, units = "in", dpi = 300)
ggsave("RE.geneage.type.tiff", width = 6.5, height = 6.5, units = "in", dpi = 300)


##RE+AS, hourglass
all.phylo<-rbind(AS.phylo[,c("PS","NS","Genes","Frequency","Mtype","dotsize")],
                 RE.phylo[,c("PS","NS","Genes","Frequency","Mtype","dotsize")])
a<-all.phylo %>%
  mutate(PS = factor(PS, levels = paste0("PS",1:12))) %>%
  ggplot(aes(x = Frequency, y = NS, size = dotsize, alpha = dotsize))+
  geom_point(shape = 1)+
  scale_colour_manual(values = "black")+
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,1))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,1))+
  labs(x = "Relative modification intensity", y = "dN/dS", legend = "")+
  facet_grid(PS~Mtype, scales = "free", space = "free")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 0),
        legend.position = "none",
        panel.grid.minor.x = element_line(linetype = "blank"),
        strip.text = element_text(size = 12))
a

ggsave("all.geneage.png", width = 4, height = 8, units = "in", dpi = 300)
ggsave("all.geneage.tiff", width = 4, height = 8, units = "in", dpi = 300)

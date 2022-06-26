setwd("AS_cash/")
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
jad=6 ##minimum distance to mismatch
cutoff.PSI=0.05 ##At least one stage showed minor isoform > 0.05 
cutoff.expression=5 ##At least one stage expression > 5

###########################################################
##Gather the PSI and expression of all AS
###########################################################
setwd("alldiff/")

##Create blank matrix
PSI.matrix<-read.delim("sitelist.txt", header = TRUE)
exp.matrix<-read.delim("sitelist.txt", header = TRUE)
event.alllist<-as.data.frame(matrix(NA, ncol = 3, nrow = 0))
names(event.alllist)<-c("Location", "SplicingType","Stage")
##Create stage list, PSI
PSI.stage.bk<-as.data.frame(matrix(NA, ncol = 3, nrow = 0))
names(PSI.stage.bk)<-c("Location", "SplicingType","PSI")
##Create stage list, expression
exp.stage.bk<-as.data.frame(matrix(NA, ncol = 3, nrow = 0))
names(exp.stage.bk)<-c("Location", "SplicingType","expression")
##Create stage list, stage positive list
event.stage.bk<-as.data.frame(matrix(NA, ncol = 2, nrow = 0))
names(event.stage.bk)<-c("Location", "SplicingType")

##Merge PSI matrix and expression matrix of stages
for (i in 1:length(treat.order)) {
  PSI.stage<-PSI.stage.bk
  exp.stage<-exp.stage.bk
  event.stage<-event.stage.bk
  cont<-treat.order[i]
  case<-treat.order[treat.order != treat.order[i]]
  for (j in 1:length(case)) {
    raw.sub<-read.delim(paste0(cont,"_",case[j],".",case[j],"vs",cont,".alldiff.txt"), header = TRUE)
    raw.sub<-raw.sub[, c("Location", "SplicingType",
                         paste0(cont,"_Junc_Inclusive..Exclusive"))]
    raw.sub<-separate(raw.sub, 
                      paste0(cont,"_Junc_Inclusive..Exclusive"), 
                      c("Inclusive","Exclusive"), 
                      sep = "::", convert = T)
    raw.sub$expression=raw.sub$Inclusive+raw.sub$Exclusive
    raw.sub$PSI<-raw.sub$Inclusive/raw.sub$expression
    
    PSI.sub<-raw.sub[,c("Location", "SplicingType","PSI")]
    exp.sub<-raw.sub[,c("Location", "SplicingType","expression")]
    event.sub<-raw.sub[raw.sub$PSI>0.05 & raw.sub$PSI<0.95 & raw.sub$expression>5,c("Location", "SplicingType")]
    PSI.stage<-rbind(PSI.stage, PSI.sub)
    exp.stage<-rbind(exp.stage, exp.sub)
    event.stage<-rbind(event.stage, event.sub)
  }
  ##PSI matrix 
  PSI.stage<-unique(PSI.stage)
  PSI.stage<-PSI.stage[PSI.stage$PSI != "NaN",]
  PSI.stage<-PSI.stage[!duplicated(PSI.stage[,c("Location", "SplicingType")]),]
  names(PSI.stage)[names(PSI.stage)=="PSI"]<-cont
  PSI.matrix<-merge(PSI.matrix, PSI.stage, by = c("Location", "SplicingType"), all = TRUE)
  
  ##expression matrix 
  exp.stage<-unique(exp.stage)
  exp.stage<-exp.stage[exp.stage$expression != "NaN",]
  exp.stage<-exp.stage[is.na(exp.stage$expression) == F,]
  exp.stage<-exp.stage[!duplicated(exp.stage[,c("Location", "SplicingType")]),]
  names(exp.stage)[names(exp.stage)=="expression"]<-cont
  exp.matrix<-merge(exp.matrix, exp.stage, by = c("Location", "SplicingType"), all = TRUE)
  
  ##event list
  event.stage<-unique(event.stage)
  event.stage$Stage<-cont
  event.alllist<-rbind(event.alllist, event.stage)
}
rm(exp.stage, exp.stage.bk, exp.sub, 
   PSI.stage, PSI.stage.bk, PSI.sub,
   event.stage, event.stage.bk, event.sub,
   raw.sub, cont, i, j)

setwd("/home/yichun/RNAmodification/AS_cash/")

##At least one stage showed minor isoform > 0.05 & at least one stage expression > 5
PSI.matrix$max<-apply(PSI.matrix[,3:11], 1, max, na.rm = T)
PSI.matrix$min<-apply(PSI.matrix[,3:11], 1, min, na.rm = T)
exp.matrix$min<-apply(exp.matrix[,3:11], 1, min, na.rm = T)

PSI.matrix<-PSI.matrix[PSI.matrix$max > cutoff.PSI & PSI.matrix$min < (1-cutoff.PSI), ]
exp.matrix<-exp.matrix[exp.matrix$min > cutoff.expression, ]
PSI.matrix<-merge(exp.matrix[,c("Location", "SplicingType")], PSI.matrix,
                  by = c("Location", "SplicingType"), all = FALSE)
rm(exp.matrix)
PSI.matrix<-PSI.matrix[,c("Location", "SplicingType", treat.order)]
row.names(PSI.matrix)<-1:nrow(PSI.matrix)

##Remove AS with mismatches around, jad = 6

##Open a new data.frame
AS<-PSI.matrix[,c("Location", "SplicingType")]
AS<-separate(AS, Location, c("CHROM","POS"), sep = ":", remove = F)
AS<-separate(AS, POS, c("Start","End"), sep = "-", remove = T, convert = TRUE)
AS$Start<-as.numeric(AS$Start)
AS$End<-as.numeric(AS$End)

##Read in SNP and INDEL
vcf<-as.data.frame(matrix(NA, ncol = 2, nrow = 0))
names(vcf)<-c("CHROM","POS")
for (i in 1:length(treat.order)) {
  vcf.sub<-read.delim(paste0("/home/yichun/RNAmodification/mergebam/",treat.order[i],".snpindel.filter.strict.txt"), header = TRUE)
  vcf<-rbind(vcf, vcf.sub[,c("CHROM","POS")])
}
rm(vcf.sub, i)
vcf<-vcf[order(vcf$CHROM, vcf$POS),c("CHROM","POS")]
vcf<-unique(vcf)

##Read in gff3 annotation
gff<-read.delim("/home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gff3", 
                skip = 1, header = FALSE)
names(gff)[1]="CHROM"
names(gff)[2]="POS"
names(gff)[7]="Strand"
gff$Start<-apply(gff[,c("V4","V5")],1, min)
gff$End<-apply(gff[,c("V4","V5")],1, max)
gff$Start<-as.numeric(gff$Start)
gff$End<-as.numeric(gff$End)
gff<-gff[,c("CHROM","POS","V3",
            "Strand","Start","End",
            "V9")]

##Get Gene ID/Gene Names for later matching with AS$checkgeneID
gff.genelist<-gff[gff$V3=="gene",]
gff.genelist<-separate(gff.genelist, V9, c("ID","Name"), sep = ";Name=")
gff.genelist$ID<-gsub("ID=","", gff.genelist$ID)
gff.genelist$ID<-gsub(";","", gff.genelist$ID)
gff.genelist$Name<-gsub(";","", gff.genelist$Name)

##Get exon region of genes
gff.exon<-gff[gff$V3=="exon",]
gff.exon<-separate(gff.exon, V9, c("V9"), sep = ";")
gff.exon<-separate(gff.exon, V9, c("ID","Transcript"), sep = "-")
gff.exon<-separate(gff.exon, Transcript, c("Transcript", "exon"), sep = "\\.")
gff.exon$ID<-gsub("ID=","",gff.exon$ID)

##Correct gene location
AS$checkexon<-NA
AS$checkgeneID<-NA
AS$checkstrand<-NA

for (i in 1:nrow(AS)) {
  exon.tmp.chr<-gff.exon[gff.exon$CHROM==AS$CHROM[i],]
  exon.tmp<-exon.tmp.chr[exon.tmp.chr$Start<=AS$Start[i] & exon.tmp.chr$End>=AS$End[i],]
  exon.tmp<-unique(exon.tmp[,-8])
  if (nrow(exon.tmp)==1) {
    ##Known exon
    AS$checkgeneID[i]<-exon.tmp$ID
    AS$checkexon[i]<-gsub("exon","",exon.tmp$exon)
    AS$checkstrand[i]<-exon.tmp$Strand
  } else {
    exon.tmp<-exon.tmp.chr[exon.tmp.chr$Start == min(exon.tmp.chr$Start[exon.tmp.chr$Start >= AS$End[i]]) |
                             exon.tmp.chr$End == max(exon.tmp.chr$End[exon.tmp.chr$End <= AS$Start[i]]),]
    exon.tmp$exon<-gsub("exon","",exon.tmp$exon)
    exon.tmp$exon<-as.numeric(exon.tmp$exon)
    ##test if within one gene, novel exon
    if (length(unique(exon.tmp$ID))==1) {
      AS$checkgeneID[i]<-unique(exon.tmp$ID)
      AS$checkstrand[i]<-unique(exon.tmp$Strand)
      ##Whether more than one transcript
      if (length(unique(exon.tmp$Transcript))>1) {
        transcript<-as.data.frame(xtabs(~Transcript, exon.tmp))
        transcript<-transcript[max(transcript$Freq),1]
        if (min(exon.tmp$exon[exon.tmp$Transcript==transcript]) != max(exon.tmp$exon[exon.tmp$Transcript==transcript])) {
          AS$checkexon[i]<-paste0(min(exon.tmp$exon[exon.tmp$Transcript==transcript]),
                                  "-",
                                  max(exon.tmp$exon[exon.tmp$Transcript==transcript]))
        } else {
          ##Avoid reprinting the same exon number
          AS$checkexon[i]<-min(exon.tmp$exon[exon.tmp$Transcript==transcript])
        }
        
      } else {
        if (min(exon.tmp$exon) != min(exon.tmp$exon)) {
          AS$checkexon[i]<-paste0(min(exon.tmp$exon),
                                  "-",
                                  max(exon.tmp$exon))
        } else {
          ##Avoid reprinting the same exon number
          AS$checkexon[i]<-min(exon.tmp$exon)
        }
      }
      
    } else {
      ##Intergenic or overlap
      AS$checkgeneID[i]<-"Intergenic"
      AS$checkstrand[i]<-"unknown"
      AS$checkexon[i]<-"Intergenic"
    }
  }
}

rm(exon.tmp,exon.tmp.chr,transcript, i)

##Remove intergenic and overlapping
AS<-AS[AS$checkgeneID != "Intergenic",]

##Start detecting mismatch
row.names(AS)<-1:nrow(AS)
AS$SNPINDEL<-NA
for (i in 1:nrow(AS)) {  
  gff.geneID<-gff.genelist[which(gff.genelist$ID == AS$checkgeneID[i]),]
  gff.genesub<-gff.exon[gff.exon$ID==gff.geneID$ID[1],]
  
  testpos<-NA
  posA<-NA
  posB<-NA
  posx<-NA
  posy<-NA
  
  if (gff.geneID$Strand == "+") {
    if (AS$checkexon[i]=="before first") {
      ##Before exon1
      posA=gff.genesub[gff.genesub$exon=="exon1","Start"]
      posA=c((posA-jad):(posA+jad))
      ##AS upstream
      posx=AS$Start[i]
      posx=c((posx-jad):(posx+jad))
      ##AS downstream
      posy=AS$End[i]
      posy=c((posy-jad):(posy+jad))
    } else if(grepl("-", AS$checkexon[i])==TRUE) {
      exonsub<-as.numeric(unlist(strsplit(AS$checkexon[i], "-")))
      ##In case of cross multiple exon, generate an exon list
      exonsub<-exonsub[1]:exonsub[2]
      ##Upstream exon end
      posA=gff.genesub[gff.genesub$exon==paste0("exon",exonsub[1]),"End"]
      posA=c((min(posA)-jad):(max(posA)+jad))
      ##Downstream exon start
      posB=gff.genesub[gff.genesub$exon==paste0("exon",exonsub[2]),"Start"]
      posB=c((min(posB)-jad):(max(posB)+jad))
      ##AS upstream
      posx=AS$Start[i]
      posx=c((posx-jad):(posx+jad))
      ##AS downstream
      posy=AS$End[i]
      posy=c((posy-jad):(posy+jad))
    } else if (grepl("-", AS$checkexon[i])==FALSE) {
      ##Upstream closest end
      posA<-gff.genesub[which(gff.genesub$End <= AS$Start[i]),"End"]
      if (length(posA)>0) {
        posA<-AS$Start[i]-min(abs(posA-AS$Start[i]))
        posA=c((posA-jad):(posA+jad))
      }
      
      ##Downstream closest start
      posB<-gff.genesub[which(gff.genesub$Start >= AS$End[i]),"Start"]
      if (length(posB)>0) {
        posB<-min(abs(posB-AS$End[i]))+AS$End[i]
        posB=c((posB-jad):(posB+jad))
      }
      
      ##AS upstream
      posx=AS$Start[i]
      posx=c((posx-jad):(posx+jad))
      ##AS downstream
      posy=AS$End[i]
      posy=c((posy-jad):(posy+jad))
    }
  } else if (gff.geneID$Strand == "-") {
    if (AS$checkexon[i]=="before first") {
      ##Before exon1
      posA=gff.genesub[gff.genesub$exon=="exon1","End"]
      posA=c((posA-jad):(posA+jad))
      ##AS upstream
      posB=AS$End[i]
      posB=c((posB-jad):(posB+jad))
      ##AS downstream
      posC=AS$Start[i]
      posC=c((posC-jad):(posC+jad))
    } else if(grepl("-", AS$checkexon[i])==TRUE) {
      exonsub<-as.numeric(unlist(strsplit(AS$checkexon[i], "-")))
      ##In case of cross multiple exon, generate an exon list
      exonsub<-exonsub[1]:exonsub[2]
      ##Upstream exon end
      posA=gff.genesub[gff.genesub$exon==paste0("exon",exonsub[1]),"Start"]
      posA=c((min(posA)-jad):(max(posA)+jad))
      ##Downstream exon start
      posB=gff.genesub[gff.genesub$exon==paste0("exon",exonsub[2]),"End"]
      posB=c((min(posB)-jad):(max(posB)+jad))
      ##AS upstream
      posx=AS$End[i]
      posx=c((posx-jad):(posx+jad))
      ##AS downstream
      posy=AS$Start[i]
      posy=c((posy-jad):(posy+jad))
    } else if (grepl("-", AS$checkexon[i])==FALSE) {
      ##Upstream closest end
      
      posA<-gff.genesub[which(gff.genesub$Start >= AS$End[i]),"Start"]
      if (length(posA)>0) {
        posA<-min(abs(posA-AS$End[i]))+AS$End[i]
        posA=c((posA-jad):(posA+jad))
      }
      
      ##Downstream closest start
      posB<-gff.genesub[which(gff.genesub$End <= AS$Start[i]),"End"]
      if (length(posB)>0) {
        posB<-AS$Start[i]-min(abs(posB-AS$Start[i]))
        posB=c((posB-jad):(posB+jad))
      }
      
      ##AS upstream
      posx=AS$End[i]
      posx=c((posx-jad):(posx+jad))
      ##AS downstream
      posy=AS$Start[i]
      posy=c((posy-jad):(posy+jad))
    }
    
  }
  
  ##Gather all position pending for test
  POS<-c(posA, posB, posx, posy)
  POS<-unique(POS)
  POS<-POS[which(is.na(POS)==F)]
  testpos<-as.data.frame(POS)
  testpos$CHROM<-NA
  testpos$CHROM<-gff.geneID$CHROM[1]
  
  testpos<-merge(testpos, vcf, all = FALSE, by = c("CHROM","POS"))
  if (nrow(testpos)>0) {
    AS$SNPINDEL[i]<-"YES"
  } else {
    AS$SNPINDEL[i]<-"NO"
  }
}

rm(testpos, exonsub, i, POS,
   posA, posB,posx,posy,
   vcf, gff.geneID, gff.genesub,
   gff.genelist, gff.exon)
AS<-AS[AS$SNPINDEL=="NO",]
AS.new<-merge(AS, PSI.matrix, by = c("Location", "SplicingType"), all = FALSE)
PSI.matrix<-merge(AS[,c("Location", "SplicingType")], PSI.matrix, by = c("Location", "SplicingType"), all = FALSE)

write.table(AS.new, 
            "Alternative_splicing.exp5.minor05.jad6.txt",
            row.names = F, quote = F, sep = "\t")

write.table(PSI.matrix, 
            "Alternative_splicing.exp5.minor05.jad6.psi.txt",
            row.names = F, quote = F, sep = "\t")

write.table(AS, 
            "Alternative_splicing.exp5.minor05.jad6.feature.txt",
            row.names = F, quote = F, sep = "\t")

write.table(event.alllist,
            "Alternative_splicing.positive_by_stage.txt",
            row.names = F, quote = F, sep = "\t")

###########################################################
##Finish gathering and filtering
###########################################################

###########################################################
##Start summarising the results
###########################################################
AS.new<-read.delim("Alternative_splicing.exp5.minor05.jad6.txt", header = T)
PSI.matrix<-read.delim("Alternative_splicing.exp5.minor05.jad6.psi.txt", header = T)
AS.feature<-read.delim("Alternative_splicing.exp5.minor05.jad6.feature.txt", header = T)
event.alllist<-read.delim("Alternative_splicing.positive_by_stage.txt", header = T)

##format to long list, all expressed AS site, allow minor ratio < 0.05 | > 0.95, remove NaN
PSI.list<-as.data.frame(matrix(NA, ncol = 4, nrow = 0))
names(PSI.list)<-c("Location", "SplicingType", "Stage", "PSI")
for (i in 3:ncol(PSI.matrix)) {
  PSI.sub<-PSI.matrix[is.na(PSI.matrix[i])==F,c(1:2,i)]
  names(PSI.sub)[3]<-"PSI"
  PSI.sub$Stage<-names(PSI.matrix)[i]
  PSI.list<-rbind(PSI.list, PSI.sub)
}

write.table(PSI.list, "AS.PSI.list.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
##get each stage positive sites
PSI.list.stage<-merge(PSI.list, event.alllist, by = c("Location", "SplicingType", "Stage"), all = FALSE)
PSI.list.stage<-PSI.list.stage[PSI.list.stage$PSI > 0.05 & PSI.list.stage$PSI < 0.95,]

write.table(PSI.list.stage, "AS.PSI.list-stage.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

PSI.list<-PSI.list.stage
PSI.sub<-PSI.list
PSI.sub$Stage<-"All"
PSI.list<-rbind(PSI.list, PSI.sub)
rm(PSI.sub,i)

PSI.list.freq<-as.data.frame(xtabs(~SplicingType+Stage, unique(PSI.list[,c("Location", "SplicingType", "Stage")])))
PSI.list.freq.sum=as.data.frame(PSI.list.freq %>% group_by(Stage) %>% summarise(sum = sum(Freq)))
PSI.list.freq.sum$Stagesum<-paste0("(",PSI.list.freq.sum$sum, ") ", PSI.list.freq.sum$Stage)

PSI.list.freq<-merge(PSI.list.freq, PSI.list.freq.sum, by = c("Stage"))
PSI.list.freq$percentage<-PSI.list.freq$Freq/PSI.list.freq$sum
rm(PSI.list.freq.sum)

##Barchart
label.order<-c("(3239) BS",
               "(1972) BS12h",
               "(1892) BS24h",
               "(2575) Myc",
               "(2412) Oidia",
               "(2619) Scl",
               "(2829) Knot",
               "(2710) Pri",
               "(2865) YFB")

bar<-PSI.list.freq %>%
  mutate(Stagesum = factor(Stagesum, levels = c(label.order,"(6839) All"))) %>%
  ggplot(aes(x=Stagesum, y = percentage, fill=SplicingType))+
  geom_bar(stat = "identity", position = "fill", width = 0.6)+
  scale_fill_manual(values=taxoncolor)+
  labs(x = "", y = "Percentage", fill = "Splicing type")+
  scale_y_continuous(labels = scales::percent, breaks = seq(0,1,0.1))+
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 0, hjust = 0.5),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        axis.line = element_line(size = 0.5),
        panel.background = element_rect(fill = NA),
        legend.position = "none")
bar

f = "AS.type.stage.bar.pptx"
topptx(bar, f, width = 3.2, height = 3.2, units = "in")
write.table(PSI.list.freq, "AS.type.stage.txt", row.names = F, quote = F, sep = '\t')
rm(PSI.list.freq, bar, label.order)

##Violin plot
library(ggpubr)

##By AS type
PSI.list$SplicingType<-gsub("Cassette_multi","Cassette\nMulti",PSI.list$SplicingType)
violin<-PSI.list[PSI.list$Stage != "All",c("SplicingType","PSI")] %>% 
  ggplot(aes(x=SplicingType, y = PSI, fill=SplicingType))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(position = position_dodge(width = 0.9), alpha = 0.5, width = 0.5, fill = "white")+
  scale_fill_manual(values=taxoncolor)+
  labs(x = "", y = "PSI", fill = "Splicing type")+
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
PSI.list$SplicingType<-gsub("Cassette\nMulti","Cassette_multi",PSI.list$SplicingType)
my_comparisons<-compare_means(PSI ~ SplicingType, data = PSI.list[PSI.list$Stage != "All",])
write.table(my_comparisons, "AS.type.violin.KWcomparison.txt",
            row.names = F, quote = F, sep = '\t')
f = "AS.type.violin.pptx"
topptx(violin, f, width = 3, height = 3, units = "in")

##By AS stage
PSI.list$SplicingType<-gsub("Cassette_multi","Cassette\nMulti",PSI.list$SplicingType)
violin<-PSI.list[,c("Stage","PSI")] %>% 
  mutate(Stage = factor(Stage, levels = c(treat.order,"All"))) %>%
  ggplot(aes(x=Stage, y = PSI, fill = Stage))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(position = position_dodge(width = 0.9), alpha = 0.5, width = 0.5, fill = "white")+
  scale_fill_manual(values=taxoncolor)+
  labs(x = "", y = "PSI")+
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
PSI.list$SplicingType<-gsub("Cassette\nMulti","Cassette_multi",PSI.list$SplicingType)
my_comparisons<-compare_means(PSI ~ Stage, data = PSI.list[PSI.list$Stage != "All",])
write.table(my_comparisons, "AS.stage.violin.KWcomparison.txt",
            row.names = F, quote = F, sep = '\t')
f = "AS.stage.violin.pptx"
topptx(violin, f, width = 3, height = 3, units = "in")

##By AS type + stage
PSI.list$SplicingType<-gsub("Cassette_multi","Cassette\nMulti",PSI.list$SplicingType)
violin<-PSI.list[PSI.list$Stage != "All", c("SplicingType","Stage","PSI")] %>% 
  mutate(Stage = factor(Stage, levels = c(treat.order,"All"))) %>%
  ggplot(aes(x=SplicingType, y = PSI, fill = SplicingType))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(position = position_dodge(width = 0.9), alpha = 0.5, width = 0.5, fill = "white")+
  scale_fill_manual(values=taxoncolor)+
  labs(x = "", y = "PSI")+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  stat_compare_means(label.y = 1.05)+
  facet_wrap(~Stage, ncol = 3,)+
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
PSI.list$SplicingType<-gsub("Cassette\nMulti","Cassette_multi",PSI.list$SplicingType)
my_comparisons<-compare_means(PSI ~ SplicingType, data = PSI.list[PSI.list$Stage != "All",],
                              group.by = "Stage")
write.table(my_comparisons, "AS.type-stage.violin.KWcomparison.txt",
            row.names = F, quote = F, sep = '\t')
f = "AS.type-stage.violin.pptx"
topptx(violin, f, width = 6, height = 6, units = "in")

##By stage + AS type
PSI.list$SplicingType<-gsub("Cassette_multi","Cassette\nMulti",PSI.list$SplicingType)
violin<-PSI.list[PSI.list$Stage != "All", c("SplicingType","Stage","PSI")] %>% 
  mutate(Stage = factor(Stage, levels = c(treat.order))) %>%
  ggplot(aes(x=Stage, y = PSI, fill = Stage))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(position = position_dodge(width = 0.9), alpha = 0.5, width = 0.5, fill = "white")+
  scale_fill_manual(values=taxoncolor)+
  labs(x = "", y = "PSI")+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  stat_compare_means(label.y = 1.05)+
  facet_wrap(~SplicingType, ncol = 3,)+
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
PSI.list$SplicingType<-gsub("Cassette\nMulti","Cassette_multi",PSI.list$SplicingType)
my_comparisons<-compare_means(PSI ~ Stage, data = PSI.list[PSI.list$Stage != "All",],
                              group.by = "SplicingType")
write.table(my_comparisons, "AS.stage-type.violin.KWcomparison.txt",
            row.names = F, quote = F, sep = '\t')
f = "AS.stage-type.violin.pptx"
topptx(violin, f, width = 6, height = 6, units = "in")
rm(my_comparisons,f,violin)

###########################################################
##AS and expression
###########################################################

PSI.list<-read.delim("AS.PSI.list-stage.txt", header = TRUE)
TMM<-read.delim("/home/yichun/RNAmodification/expression/gene_count_log2TMMmatrix.txt", header = TRUE)
TMM.mean<-as.data.frame(matrix(NA, nrow = nrow(TMM), ncol = length(treat.order)))
row.names(TMM.mean)<-row.names(TMM)
names(TMM.mean)<-treat.order

for (i in 1:ncol(TMM.mean)) {
  TMM.mean[,i]<-(TMM[,names(TMM)==paste0(names(TMM.mean)[i],"1")]+TMM[,names(TMM)==paste0(names(TMM.mean)[i],"2")]+TMM[,names(TMM)==paste0(names(TMM.mean)[i],"3")])/3
}
TMM.mean$Gene<-row.names(TMM.mean)
row.names(TMM.mean)<-1:nrow(TMM.mean)

TMM.list<-as.data.frame(matrix(NA, ncol = 3, nrow = 0))
names(TMM.list)<-c("Gene", "Stage", "TMM")
for (i in 1:(ncol(TMM.mean)-1)) {
  TMM.list.sub<-TMM.mean[,c(10,i)]
  names(TMM.list.sub)[2]<-"TMM"
  TMM.list.sub$Stage<-names(TMM.mean)[i]
  TMM.list<-rbind(TMM.list, TMM.list.sub)
}
rm(i, TMM.list.sub)

AS.PSI.TMM<-merge(PSI.list, AS.feature[,c("Location","SplicingType","checkgeneID")],
                  by = c("Location","SplicingType"), all.x = TRUE)
names(AS.PSI.TMM)[names(AS.PSI.TMM)=="checkgeneID"]<-"Gene"
AS.PSI.TMM<-merge(AS.PSI.TMM, TMM.list, by = c("Gene","Stage"), all.x = TRUE)

write.table(AS.PSI.TMM, "AS.PSI.TMM.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

rm(PSI.sub, TMM.list, TMM.mean, TMM)

##scatter plot
AS.PSI.TMM<-read.delim("AS.PSI.TMM.txt", header = T)
#summary(lm(log2(PSI)~log2(TMM),AS.PSI.TMM))
a<-ggplot(data = AS.PSI.TMM, aes(x = TMM, y = PSI, 
                                 colour = SplicingType))+
  geom_point(shape = 16, size = 0.05, alpha = 0.8)+
  scale_colour_manual(values = taxoncolor)+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1))+
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1))+
  labs(y = "Alternative splicing (PSI)", x = "Gene expression (Log2TMM)", legend = "")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 9, colour = "black"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 0),
        legend.position = "none",
        panel.grid.minor.x = element_line(linetype = "blank"),
        panel.background = element_rect(fill = NA),
        strip.text = element_text(size = 9))
a
ggsave("AS.PSI.TMM.png", width = 4, height = 4, units = "in", dpi = 300)
ggsave("AS.PSI.TMM.tiff", width = 4, height = 4, units = "in", dpi = 300)


###########################################################
##Differential AS
###########################################################

##Gather all differential AS record

filelist<-read.delim("alldiff/filelist.txt", header = TRUE)
AS.diflist<-as.data.frame(matrix(NA, ncol = 10, nrow = 0))
names(AS.diflist)<-c("Location",
                     "Mut_Junc_Inclusive..Exclusive", 
                     "Ref_Junc_Inclusive..Exclusive",
                     "Mut_Exp_Inclusive..Exclusive", 
                     "Ref_Exp_Inclusive..Exclusive",
                     "delta_PSI","P.Value","FDR",
                     "SplicingType",
                     "samplepair")

for (i in 1:nrow(filelist)) {
  AS.list.sub<-NA
  AS.list.sub<-read.delim(paste0("alldiff/",filelist$file[i]), header = TRUE)
  names(AS.list.sub)[4:7]<-c("SampleA_Junc_Inclusive..Exclusive", 
                             "SampleB_Junc_Inclusive..Exclusive",
                             "SampleA_Exp_Inclusive..Exclusive", 
                             "SampleB_Exp_Inclusive..Exclusive")
  AS.list.sub<-AS.list.sub[,
                           c("Location","SplicingType",
                             "delta_PSI","P.Value","FDR",
                             "SampleA_Junc_Inclusive..Exclusive", 
                             "SampleB_Junc_Inclusive..Exclusive",
                             "SampleA_Exp_Inclusive..Exclusive", 
                             "SampleB_Exp_Inclusive..Exclusive")]
  if (nrow(AS.list.sub)>0) {
    AS.list.sub$samplepair<-filelist$file[i]
    AS.diflist<-rbind(AS.diflist, AS.list.sub)
  }
}
rm(AS.list.sub, i)

##Get only filtered pairs
AS.diflist<-merge(AS.feature[,c("Location","SplicingType")], AS.diflist, 
                  by = c("Location","SplicingType"),
                  all.x = TRUE)

AS.diflist<-separate(AS.diflist, samplepair, c("samplepair"), sep = "\\.")
AS.diflist<-separate(AS.diflist, samplepair, c("sampleB","sampleA"), sep = "_")
write.table(AS.diflist, "AS.cashtest.all.txt",
            row.names = F, quote = F, sep = "\t")

AS.diflist<-read.delim("AS.cashtest.all.txt", header = T)
AS.diflist<-AS.diflist[AS.diflist$FDR<0.05 & abs(AS.diflist$delta_PSI) > 0.1,]

##PSI heatmap
library(pheatmap)
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
##All
heatmap.input<-PSI.matrix
row.names(heatmap.input)<-paste0(heatmap.input$SplicingType, "_", heatmap.input$Location)
heatmap.input<-heatmap.input[,c(3:11)]
heatmap.input[is.na(heatmap.input)]<-0

heatp<-pheatmap(heatmap.input,
                legend = TRUE,
                scale = "none",
                #scale = "row",
                cluster_rows = TRUE, cluster_cols = FALSE,
                treeheight_row = 0, 
                border_color = NA,
                cellwidth = 30, cellheight = 0.1, fontsize =12,
                angle_col = c("45"),
                width = 6, height = 12,
                show_colnames = TRUE,
                show_rownames = F,
                na_col = "gray",
                filename = "AS.PSI.heatmap.png")

heatmap.input<-as.data.frame(t(apply(heatmap.input, 1, normalize)), row.names = row.names(heatmap.input))
heatp<-pheatmap(heatmap.input,
                legend = TRUE,
                scale = "none",
                #scale = "row",
                cluster_rows = TRUE, cluster_cols = FALSE,
                treeheight_row = 0, 
                border_color = NA,
                cellwidth = 30, cellheight = 0.1, fontsize =12,
                angle_col = c("45"),
                width = 6, height = 12,
                show_colnames = TRUE,
                show_rownames = F,
                na_col = "gray",
                filename = "AS.PSI.heatmap.norm.png")

##Differential PSI
heatmap.input<-merge(unique(AS.diflist[,c("Location","SplicingType")]), PSI.matrix,
                     by = c("Location","SplicingType"), all.x = TRUE)

row.names(heatmap.input)<-paste0(heatmap.input$SplicingType, "_", heatmap.input$Location)
heatmap.input<-heatmap.input[,c(3:11)]
heatmap.input[is.na(heatmap.input)]<-0

heatp<-pheatmap(heatmap.input,
                legend = TRUE,
                scale = "none",
                #scale = "row",
                cluster_rows = TRUE, cluster_cols = FALSE,
                treeheight_row = 0, 
                border_color = NA,
                cellwidth = 30, cellheight = 0.1, fontsize =12,
                angle_col = c("45"),
                width = 6, height = 6,
                show_colnames = TRUE,
                show_rownames = F,
                na_col = "gray",
                filename = "AS.PSIdif.heatmap.png")

heatmap.input<-as.data.frame(t(apply(heatmap.input, 1, normalize)), row.names = row.names(heatmap.input))
heatp<-pheatmap(heatmap.input,
                legend = TRUE,
                scale = "none",
                #scale = "row",
                cluster_rows = TRUE, cluster_cols = FALSE,
                treeheight_row = 0, 
                border_color = NA,
                cellwidth = 30, cellheight = 0.1, fontsize =12,
                angle_col = c("45"),
                width = 6, height = 6,
                show_colnames = TRUE,
                show_rownames = F,
                na_col = "gray",
                filename = "AS.PSIdif.heatmap.norm.png")

rm(heatmap.input, heatp, filelist)

##Combine delta-PSI and log2FC
FC<-read.delim("C:/coprinopsis/expression/DEG.allpairs.txt", header = TRUE)
FC<-FC[,c("Genes","sampleA","sampleB","logFC")]
AS.diflist<-read.delim("AS.cashtest.all.txt", header = TRUE)
AS.feature<-read.delim("Alternative_splicing.exp5.minor05.jad6.feature.txt", header = T)
AS.dPSI.FC<-merge(AS.diflist, AS.feature[,c("Location","SplicingType","checkgeneID")],
                  by = c("Location","SplicingType"), all.x = TRUE)
names(AS.dPSI.FC)[names(AS.dPSI.FC)=="checkgeneID"]<-"Genes"
AS.dPSI.FC<-merge(AS.dPSI.FC, FC, by = c("Genes","sampleA","sampleB"), all.x = TRUE)
AS.dPSI.FC<-AS.dPSI.FC[is.na(AS.dPSI.FC$logFC) ==F,]
write.table(AS.dPSI.FC, "AS.dPSI.FC.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

AS.dPSI.FC<-read.delim("AS.dPSI.FC.txt",header = TRUE)
#summary(lm(delta_PSI~logFC,AS.dPSI.FC))
a<-ggplot(data = AS.dPSI.FC, aes(y = delta_PSI, x = logFC, 
                              colour = SplicingType))+
  geom_point(shape = 16, size = 0.05, alpha = 0.8)+
  scale_colour_manual(values = taxoncolor)+
  scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,0.2))+
  scale_x_continuous(limits = c(-10,10), breaks = seq(-10,10,1))+
  labs(y = "delta_PSI", x = "-log2(Fold change)", legend = "")+
  geom_vline(xintercept = -1, color = "grey30", linetype = "dashed")+
  geom_vline(xintercept = 1, color = "grey30", linetype = "dashed")+
  geom_hline(yintercept = 0.1, color = "grey30", linetype = "dashed")+
  geom_hline(yintercept = -0.1, color = "grey30", linetype = "dashed")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 9, colour = "black"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 0),
        legend.position = "none",
        panel.grid.minor.x = element_line(linetype = "blank"),
        panel.background = element_rect(fill = NA),
        strip.text = element_text(size = 9))
a
ggsave("AS.dPSI.FC.png", width = 4, height = 4, units = "in", dpi = 300)
ggsave("AS.dPSI.FC.tiff", width = 4, height = 4, units = "in", dpi = 300)

##each stage
a<-AS.dPSI.FC %>% 
  mutate(sampleA = factor(sampleA, levels = treat.order),
         sampleB = factor(sampleB, levels = treat.order)) %>%
  ggplot(aes(y = delta_PSI, x = logFC, 
                                 colour = SplicingType))+
  geom_point(shape = 16, size = 0.05, alpha = 0.8)+
  scale_colour_manual(values = taxoncolor)+
  scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,0.5))+
  scale_x_continuous(limits = c(-10,10), breaks = seq(-10,10,5))+
  labs(y = "delta_PSI", x = "-log2(Fold change)", legend = "")+
  geom_vline(xintercept = -1, color = "grey30", linetype = "dashed")+
  geom_vline(xintercept = 1, color = "grey30", linetype = "dashed")+
  geom_hline(yintercept = 0.1, color = "grey30", linetype = "dashed")+
  geom_hline(yintercept = -0.1, color = "grey30", linetype = "dashed")+
  facet_grid(sampleB~sampleA, scales = "free", space = "free")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 9, colour = "black"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 0),
        legend.position = "none",
        panel.grid.minor.x = element_line(linetype = "blank"),
        panel.background = element_rect(fill = NA),
        strip.text = element_text(size = 9))
a
ggsave("AS.dPSI.FC.stage.png", width = 6.5, height = 6.5, units = "in", dpi = 300)
ggsave("AS.dPSI.FC.stage.tiff", width = 6.5, height = 6.5, units = "in", dpi = 300)

AS.dPSI.FC.summary<-as.data.frame(xtabs(~SplicingType+sampleA+sampleB, AS.dPSI.FC[abs(AS.dPSI.FC$delta_PSI)>0.1,]))
piep<-AS.dPSI.FC.summary %>%
  mutate(sampleA = factor(sampleA, levels = treat.order),
         sampleB = factor(sampleB, levels = treat.order)) %>%
  ggplot(aes(x="", y=Freq, fill=SplicingType))+
  geom_bar(stat = "identity", position = position_fill())+
  labs(fill = "SplicingType")+
  scale_fill_manual(values=taxoncolor)+
  coord_polar(theta = "y", start=0)+
  facet_wrap(sampleB~sampleA, ncol = 9)+
  theme_void()

piep

f = "AS.dPSI.FC.stage.pie.pptx"
topptx(piep,f, width = 6.5, height = 7.8, units = "in")

AS.dPSI.FC.summary1<-as.data.frame(xtabs(~sampleA+sampleB, AS.dPSI.FC[abs(AS.dPSI.FC$delta_PSI)>0.1,]))
piep<-AS.dPSI.FC.summary1 %>%
  mutate(sampleA = factor(sampleA, levels = treat.order),
         sampleB = factor(sampleB, levels = treat.order)) %>%
  ggplot(aes(x="", y=Freq))+
  geom_bar(stat = "identity", position = position_fill())+
  geom_text(aes(label = Freq), position = position_fill(vjust = 0.5), color = "black")+
  facet_wrap(sampleB~sampleA, ncol = 9)+
  theme_void()

piep

f = "AS.dPSI.FC.stage.count.pptx"
topptx(piep,f, width = 6.5, height = 7.8, units = "in")

##Functional enrichment
setwd("AS_cash/")
AS.new<-read.delim("Alternative_splicing.exp5.minor05.jad6.txt", header = T)
PSI.matrix<-read.delim("Alternative_splicing.exp5.minor05.jad6.psi.txt", header = T)
AS.feature<-read.delim("Alternative_splicing.exp5.minor05.jad6.feature.txt", header = T)
AS.diflist<-read.delim("AS.cashtest.all.txt", header = TRUE)
AS.diflist<-AS.diflist[AS.diflist$FDR<0.05 & abs(AS.diflist$delta_PSI) > 0.1,]
AS.diflist<-merge(AS.diflist, AS.feature[,c("Location","SplicingType","checkgeneID")],
                  by = c("Location","SplicingType"), all.x = TRUE)
names(AS.diflist)[names(AS.diflist)=="checkgeneID"]<-"Genes"

AS.diflist.en<-AS.diflist
AS.diflist.en<-AS.diflist.en[,c("sampleA","sampleB","Genes")]
AS.diflist.en$Groups<-NA
AS.diflist.en$Groups[AS.diflist.en$sampleB %in% c("BS", "BS12h") & 
                       AS.diflist.en$sampleA %in% c("BS12h","BS24h")]<-"Spore\ngermination"
AS.diflist.en$Groups[AS.diflist.en$sampleB %in% c("Myc") & 
                       AS.diflist.en$sampleA %in% c("Oidia")]<-"Oidiation"
AS.diflist.en$Groups[AS.diflist.en$sampleB %in% c("Myc") & 
                       AS.diflist.en$sampleA %in% c("Scl")]<-"Sclerotia\nformation"
AS.diflist.en$Groups[AS.diflist.en$sampleB %in% c("Myc", "Knot") & 
                       AS.diflist.en$sampleA %in% c("Knot","Pri")]<-"Fruiting"
AS.diflist.en$Groups[AS.diflist.en$sampleB %in% c("Pri", "YFB") & 
                       AS.diflist.en$sampleA %in% c("YFB","BS")]<-"Sporulation"

AS.diflist.en<-AS.diflist.en[is.na(AS.diflist.en$Groups)==F,]
unique(AS.diflist.en[,c("sampleA","sampleB","Groups")])

path.order<-c("Spore\ngermination",
              "Oidiation",
              "Sclerotia\nformation",
              "Fruiting",
              "Sporulation")

IDmatch<-read.delim("IDmatch.txt", header = TRUE)
AS.diflist.en<-merge(AS.diflist.en, IDmatch[,c("Genes","CC130")], by = "Genes", all.x = TRUE)

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
ref="AS.difPSI."
####KOG enrich
{
  KOG.all.1<-compareCluster(Genes ~ Groups, 
                            data = AS.diflist.en, 
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
    mutate(Groups = factor(Groups, levels = path.order)) %>%
    ggplot(aes(x = Groups,
               y = paste0(Description," (",BGnumerator,")"),
               size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
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

library("AnnotationHub")
annotationhub<-AnnotationHub()
query(annotationhub, "Coprinopsis cinerea")
ah.ccin<-AnnotationHub()[["AH95203"]]


####KEGG enrich
{
  KEGG.all.1<-compareCluster(Genes ~ Groups, 
                             data = AS.diflist.en, 
                             fun = 'enricher',
                             TERM2GENE = GenesKEGGpair.1v1,
                             TERM2NAME = kegg2name,
                             pvalueCutoff = 1,
                             pAdjustMethod = "BH",
                             qvalueCutoff = 1,
                             minGSSize = 1,
                             maxGSSize = 200000)
  
  KEGG.all.1<-compareCluster(CC130~Groups, fun = "enrichKEGG", data = AS.diflist.en, 
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
    mutate(Groups = factor(Groups, levels = path.order)) %>%
    ggplot(aes(x = Groups, 
               y = paste0(Description," (",BGnumerator,")"), 
               size = ratio1, colour = p.adjust))+
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
          strip.text.y = element_text(size = 9, angle = 0))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")
  
  a
  ggsave(paste0(ref,"KEGGr1.tiff"), width = 10, height = 12, units = "in", dpi = 300)
  ggsave(paste0(ref,"KEGGr1.png"), width = 10, height = 12, units = "in", dpi = 300)
  
}

####ko enrich
{
  ko.all.1<-compareCluster(Genes ~ Groups, 
                           data = AS.diflist.en, 
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
  
  plotdata1<-subset(plotdata, ratio1>0.8 & p.adjust< 0.2)
  a<-plotdata1 %>%
    mutate(Groups = factor(Groups, levels = path.order)) %>%
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
  ggsave(paste0(ref,"ko.tiff"), width = 18, height = 24, units = "in", dpi = 300)
  ggsave(paste0(ref,"ko.png"), width = 18, height = 24, units = "in", dpi = 300)
  
}

####GO enrich
{
  GO.all.1<-compareCluster(Genes ~ Groups, 
                           data = AS.diflist.en, 
                           fun = 'enricher',
                           TERM2GENE = GenesGOpair.1v1,
                           TERM2NAME = go2name,
                           pvalueCutoff = 1,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 1,
                           minGSSize = 1,
                           maxGSSize = 2000000)
  
  GO.all.1<-compareCluster(CC130 ~ Groups, fun = "enrichGO", data = AS.diflist.en, 
                           OrgDb = ah.ccin, keyType = "SYMBOL", ont = "ALL", pvalueCutoff = 1, 
                           pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 10,maxGSSize = 200000, 
                           readable = FALSE, pool = FALSE)
  
  #Plot
  plotin<-as.data.frame(GO.all.1)
  
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
  #top 5 function
  plotdata.top5<-plotdata %>% group_by(Groups) %>% top_n(-10, p.adjust)
  #plotdata.top5<-plotdata1 %>% group_by(Cluster) %>% arrange(-ratio2/p.adjust) %>% slice_head(n=5)
  write.table(plotdata.top5, 
              paste0(ref,"GO.topbypadj.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  plotdata1<-plotdata.top5[plotdata.top5$Count>2,]
  
  a<-plotdata1 %>%
    mutate(Groups = factor(Groups, levels = path.order)) %>%
    ggplot(aes(x = Groups, y = paste0(Description," (",BGnumerator,")"), size = ratio1, colour = p.adjust))+
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
  
  ggsave(paste0(ref,"GOr1.tiff"), width = 12, height = 10, units = "in", dpi = 300)
  ggsave(paste0(ref,"GOr1.png"), width = 12, height = 10, units = "in", dpi = 300)
  
  
  GO.all.1.MF<-compareCluster(CC130 ~ Groups, fun = "enrichGO", data = AS.diflist.en, 
                              OrgDb = ah.ccin, keyType = "SYMBOL", ont = "MF", pvalueCutoff = 1, 
                              pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 10,maxGSSize = 200000, 
                              readable = FALSE, pool = FALSE)
  GO.all.1.MF<-gofilter(GO.all.1.MF, level = 5)
  GO.all.1.MF<-as.data.frame(GO.all.1.MF)
  GO.all.1.MF$ONTOLOGY<-"MF"
  GO.all.1.CC<-compareCluster(CC130 ~ Groups, fun = "enrichGO", data = AS.diflist.en, 
                              OrgDb = ah.ccin, keyType = "SYMBOL", ont = "CC", pvalueCutoff = 1, 
                              pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 10,maxGSSize = 200000, 
                              readable = FALSE, pool = FALSE)
  GO.all.1.CC<-gofilter(GO.all.1.CC, level = 5)
  GO.all.1.CC<-as.data.frame(GO.all.1.CC)
  GO.all.1.CC$ONTOLOGY<-"CC"
  GO.all.1.BP<-compareCluster(CC130 ~ Groups, fun = "enrichGO", data = AS.diflist.en, 
                              OrgDb = ah.ccin, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, 
                              pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 10,maxGSSize = 200000, 
                              readable = FALSE, pool = FALSE)
  GO.all.1.BP<-gofilter(GO.all.1.BP, level = 5)
  GO.all.1.BP<-as.data.frame(GO.all.1.BP)
  GO.all.1.BP$ONTOLOGY<-"BP"
  plotin<-rbind(GO.all.1.MF, GO.all.1.CC, GO.all.1.BP)
  
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotdata<-plotinsep
  plotdata$ONTOLOGY<-gsub("BP", "Biological\nProcess", plotdata$ONTOLOGY)
  plotdata$ONTOLOGY<-gsub("CC", "Cellular\nComponent", plotdata$ONTOLOGY)
  plotdata$ONTOLOGY<-gsub("MF", "Molecular\nFunction", plotdata$ONTOLOGY)
  
  #colnames(plotinsep)[names(plotinsep)=="ID"]<-"goClass"
  #go2name$ONTOLOGY<-gsub("biological_process", "Biological Process", go2name$ONTOLOGY)
  #go2name$ONTOLOGY<-gsub("molecular_function", "Molecular Function", go2name$ONTOLOGY)
  #go2name$ONTOLOGY<-gsub("cellular_component", "Cellular Component", go2name$ONTOLOGY)
  
  #plotdata<-merge(plotinsep, go2name, by = c("goClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$enrichfold=plotdata$ratio2*plotdata$BGdenominator/plotdata$BGnumerator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  
  write.table(plotdata, file = paste0(ref,"GOenrichL5.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  plotdata1<-subset(plotdata, ratio2>0.1 & p.adjust< 1)
  plotdata1$Description<-gsub(", with", ",\nwith",plotdata1$Description)
  plotdata1$Description<-gsub(", in", ",\nin",plotdata1$Description)
  a<-plotdata1[plotdata1$Count>2,] %>%
    mutate(Groups = factor(Groups, levels = path.order)) %>%
    ggplot(aes(x = Groups, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,1))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 8), 
          axis.text.x = element_text(size = 8), 
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size = 8), 
          plot.title = element_text(size = 8))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  
  ggsave(paste0(ref,"GOr2L5f.tiff"), width = 8, height = 8, units = "in", dpi = 300)
  ggsave(paste0(ref,"GOr2L5f.png"), width = 8, height = 8, units = "in", dpi = 300)
  
  #Significant function
  plotdata.sig<-plotdata1 %>% group_by(Groups) %>% top_n(5, enrichfold)
  plotdata.sig<-plotdata.sig[,c("Groups","ID", "Description","ONTOLOGY")]
  plotdata.sig$ONTOLOGY<-gsub("\n"," ", plotdata.sig$ONTOLOGY)
  plotdata.sig$Groups<-gsub("\n"," ", plotdata.sig$Groups)
  write.table(plotdata.sig, 
              paste0(ref,"GO.topbyfold.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  plotdata.sig<-plotdata1 %>% group_by(Groups) %>% top_n(5, enrichfold)
  #plotdata.sig$Groups<-gsub("\n"," ", plotdata.sig$Groups)
  #plotdata.sig$ONTOLOGY<-gsub("\n"," ", plotdata.sig$ONTOLOGY)
  a<-plotdata.sig %>%
    mutate(Groups = factor(Groups, levels = path.order)) %>%
    ggplot(aes(x = Groups, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,1))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 18), 
          axis.text.x = element_text(size = 18), 
          axis.text.y = element_text(size = 18),
          legend.title = element_text(size = 18), 
          legend.text = element_text(size = 16), 
          #legend.position = "bottom",
          plot.title = element_text(size = 18))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 18))
  a
  
  f = paste0(ref,"GOr2L5sig.pptx")
  topptx(a, f, width = 16, height = 7, units = "in")
  
  ggsave(paste0(ref,"GOr2L5sig.tiff"), width = 16, height = 7, units = "in", dpi = 300)
  ggsave(paste0(ref,"GOr2L5sig.png"), width = 16, height = 7, units = "in", dpi = 300)
}

###########################################################
##Get path specific regulation
###########################################################
path.order<-c("Spore\ngermination",
              "Oidiation",
              "Sclerotia\nformation",
              "Fruiting",
              "Sporulation")
AS.feature<-read.delim("Alternative_splicing.exp5.minor05.jad6.feature.txt", header = T)
AS.diflist<-read.delim("AS.cashtest.all.txt", header = TRUE)
AS.diflist<-AS.diflist[AS.diflist$FDR<0.05 & abs(AS.diflist$delta_PSI) > 0.1,]
AS.diflist<-merge(AS.diflist, AS.feature[,c("Location","SplicingType","checkgeneID")],
                  by = c("Location","SplicingType"), all.x = TRUE)
names(AS.diflist)[names(AS.diflist)=="checkgeneID"]<-"Genes"

AS.diflist.en<-AS.diflist
AS.diflist.en<-AS.diflist.en[,c("sampleA","sampleB","Location","SplicingType")]
AS.diflist.en$Groups<-NA
AS.diflist.en$Groups[AS.diflist.en$sampleB %in% c("BS", "BS12h") & 
                       AS.diflist.en$sampleA %in% c("BS12h","BS24h")]<-"Spore\ngermination"
AS.diflist.en$Groups[AS.diflist.en$sampleB %in% c("Myc") & 
                       AS.diflist.en$sampleA %in% c("Oidia")]<-"Oidiation"
AS.diflist.en$Groups[AS.diflist.en$sampleB %in% c("Myc") & 
                       AS.diflist.en$sampleA %in% c("Scl")]<-"Sclerotia\nformation"
AS.diflist.en$Groups[AS.diflist.en$sampleB %in% c("Myc", "Knot") & 
                       AS.diflist.en$sampleA %in% c("Knot","Pri")]<-"Fruiting"
AS.diflist.en$Groups[AS.diflist.en$sampleB %in% c("Pri", "YFB") & 
                       AS.diflist.en$sampleA %in% c("YFB","BS")]<-"Sporulation"

AS.diflist.en<-AS.diflist.en[is.na(AS.diflist.en$Groups)==F, c("Location", "Groups","SplicingType")]
AS.diflist.en<-unique(AS.diflist.en)

##venn diagram
library("VennDiagram")
AS.diflist.en$region<-paste0(AS.diflist.en$Location,AS.diflist.en$SplicingType)

venn.diagram(x = list(`Spore\ngermination` = AS.diflist.en[AS.diflist.en$Groups=="Spore\ngermination","region"], 
                      `Oidiation` = AS.diflist.en[AS.diflist.en$Groups=="Oidiation","region"],
                      `Sclerotia\nformation` = AS.diflist.en[AS.diflist.en$Groups=="Sclerotia\nformation","region"],
                      `Fruiting` = AS.diflist.en[AS.diflist.en$Groups=="Fruiting","region"],
                      `Sporulation` = AS.diflist.en[AS.diflist.en$Groups=="Sporulation","region"]),
             margin = 0.15, cat.dist = 0.25,
             fontfamily = "arial", sub.fontfamily = "arial", main.fontfamily = "arial", cat.fonfamily = "arial",
             filename = "AS.dPSI.path.venn1.png", imagetype = "png" , 
             units = "in", height = 6, width = 6, resolution = 300,
             fill = c("skyblue", "plum2", "burlywood2", "darkgray", "tan1"), na = "remove")

##Get specific
AS.diflist.en.freq<-as.data.frame(xtabs(~Location+SplicingType, AS.diflist.en))
AS.diflist.en.freq<-AS.diflist.en.freq[AS.diflist.en.freq$Freq>0,]
AS.diflist.en.freq<-AS.diflist.en.freq[AS.diflist.en.freq$Freq==1,]

AS.diflist.en.spc<-merge(AS.diflist.en, AS.diflist.en.freq, by = c("Location","SplicingType"), all = FALSE)
AS.diflist.en.spc<-merge(AS.diflist.en.spc, 
                         AS.feature[,c("Location","SplicingType","checkgeneID")],
                         by = c("Location","SplicingType"), all.x = TRUE)
names(AS.diflist.en.spc)[names(AS.diflist.en.spc)=="checkgeneID"]<-"Genes"
gene2name<-read.delim("/home/yichun/RNAmodification/genome/eggnog/new/gene2name.txt", header = TRUE)
AS.diflist.en.spc<-merge(AS.diflist.en.spc, 
                         gene2name,
                         by = "Genes", all.x = TRUE)
write.table(AS.diflist.en.spc,"AS.dPSI.path.spclist.txt",
            row.names = F, quote = F, sep = '\t')

AS.diflist.en.spc.freq<-as.data.frame(xtabs(~SplicingType+Groups, AS.diflist.en.spc))
AS.diflist.en.spc.freq.sum=as.data.frame(AS.diflist.en.spc.freq %>% group_by(Groups) %>% summarise(sum = sum(Freq)))
AS.diflist.en.spc.freq.sum$Stagesum<-paste0(AS.diflist.en.spc.freq.sum$Stage, " (",AS.diflist.en.spc.freq.sum$sum, ")")

AS.diflist.en.spc.freq<-merge(AS.diflist.en.spc.freq, AS.diflist.en.spc.freq.sum, by = c("Groups"))
AS.diflist.en.spc.freq$percentage<-AS.diflist.en.spc.freq$Freq/AS.diflist.en.spc.freq$sum
rm(AS.diflist.en.spc.freq.sum)

##Barchart
bar<-AS.diflist.en.spc.freq %>%
  mutate(Groups = factor(Groups, levels = path.order)) %>%
  ggplot(aes(x=Groups, y = percentage, fill=SplicingType))+
  geom_bar(stat = "identity", position = "fill", width = 0.6)+
  scale_fill_manual(values=taxoncolor)+
  labs(x = "", y = "Percentage", fill = "Splicing type")+
  scale_y_continuous(labels = scales::percent, breaks = seq(0,1,0.1))+
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 9, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 0, hjust = 0.5),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        axis.line = element_line(size = 0.5),
        panel.background = element_rect(fill = NA),
        legend.position = "none")
bar

f = "AS.dPSI.path.spc.bar.pptx"
topptx(bar, f, width = 2.5, height = 3, units = "in")
write.table(AS.diflist.en.spc.freq, "AS.dPSI.path.spc.freq.txt", row.names = F, quote = F, sep = '\t')
rm(AS.diflist.en.spc.freq, bar)

##GO of path specific 
{
  AS.diflist<-merge(AS.diflist.en.spc, AS.feature[,c("Location","SplicingType","checkgeneID")],
                    by = c("Location","SplicingType"), all.x = TRUE)
  names(AS.diflist)[names(AS.diflist)=="checkgeneID"]<-"Genes"
  
  IDmatch<-read.delim("IDmatch.txt", header = TRUE)
  AS.diflist.en<-merge(AS.diflist.en.spc, IDmatch[,c("Genes","CC130")], by = "Genes", all.x = TRUE)
  ref="AS.dPSI.path.spc."
  
  library(clusterProfiler)
  library("AnnotationHub")
  annotationhub<-AnnotationHub()
  query(annotationhub, "Coprinopsis cinerea")
  #ah.ccin<-AnnotationHub()[["AH95203"]] #biocversion: 3.13
  ah.ccin<-AnnotationHub()[["AH97415"]] #biocversion: 3.14
  addline_format <- function(x,...){ 
    gsub('_','\n',x) 
  } 
  
  #GenesGOpair.1v1<-read.delim("/home/yichun/RNAmodification/genome/eggnog/GO.1v1.txt", header = TRUE)
  go2name<-read.delim("/home/yichun/3enrichment/go2name.txt", 
                      sep = "\t", colClasses = "character",
                      header = FALSE)
  names(go2name)[1]="goClass"
  names(go2name)[2]="goName"
  names(go2name)[3]="ONTOLOGY"
  go2name$ONTOLOGY<-gsub("biological_process", "Biological Process", go2name$ONTOLOGY)
  go2name$ONTOLOGY<-gsub("molecular_function", "Molecular Function", go2name$ONTOLOGY)
  go2name$ONTOLOGY<-gsub("cellular_component", "Cellular Component", go2name$ONTOLOGY)
  
  
  GO.all.1.MF<-compareCluster(CC130 ~ Groups, fun = "enrichGO", data = AS.diflist.en, 
                              OrgDb = ah.ccin, keyType = "SYMBOL", ont = "MF", pvalueCutoff = 1, 
                              pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 10,maxGSSize = 200000, 
                              readable = FALSE, pool = FALSE)
  GO.all.1.MF<-gofilter(GO.all.1.MF, level = 5)
  GO.all.1.MF<-as.data.frame(GO.all.1.MF)
  GO.all.1.MF$ONTOLOGY<-"MF"
  GO.all.1.CC<-compareCluster(CC130 ~ Groups, fun = "enrichGO", data = AS.diflist.en, 
                              OrgDb = ah.ccin, keyType = "SYMBOL", ont = "CC", pvalueCutoff = 1, 
                              pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 10,maxGSSize = 200000, 
                              readable = FALSE, pool = FALSE)
  GO.all.1.CC<-gofilter(GO.all.1.CC, level = 5)
  GO.all.1.CC<-as.data.frame(GO.all.1.CC)
  GO.all.1.CC$ONTOLOGY<-"CC"
  GO.all.1.BP<-compareCluster(CC130 ~ Groups, fun = "enrichGO", data = AS.diflist.en, 
                              OrgDb = ah.ccin, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, 
                              pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 10,maxGSSize = 200000, 
                              readable = FALSE, pool = FALSE)
  GO.all.1.BP<-gofilter(GO.all.1.BP, level = 5)
  GO.all.1.BP<-as.data.frame(GO.all.1.BP)
  GO.all.1.BP$ONTOLOGY<-"BP"
  plotin<-rbind(GO.all.1.MF, GO.all.1.CC, GO.all.1.BP)
  
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotdata<-plotinsep
  plotdata$ONTOLOGY<-gsub("BP", "Biological\nProcess", plotdata$ONTOLOGY)
  plotdata$ONTOLOGY<-gsub("CC", "Cellular\nComponent", plotdata$ONTOLOGY)
  plotdata$ONTOLOGY<-gsub("MF", "Molecular\nFunction", plotdata$ONTOLOGY)
  
  #colnames(plotinsep)[names(plotinsep)=="ID"]<-"goClass"
  #go2name$ONTOLOGY<-gsub("biological_process", "Biological Process", go2name$ONTOLOGY)
  #go2name$ONTOLOGY<-gsub("molecular_function", "Molecular Function", go2name$ONTOLOGY)
  #go2name$ONTOLOGY<-gsub("cellular_component", "Cellular Component", go2name$ONTOLOGY)
  
  #plotdata<-merge(plotinsep, go2name, by = c("goClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$enrichfold=plotdata$ratio2*plotdata$BGdenominator/plotdata$BGnumerator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  
  write.table(plotdata, file = paste0(ref,"GOenrichL5.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  plotdata1<-subset(plotdata, ratio2>0.1 & p.adjust< 1)
  plotdata1$Description<-gsub(", with", ",\nwith",plotdata1$Description)
  plotdata1$Description<-gsub(", in", ",\nin",plotdata1$Description)
  a<-plotdata1[plotdata1$Count>2,] %>%
    mutate(Groups = factor(Groups, levels = path.order)) %>%
    ggplot(aes(x = Groups, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,1))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 8), 
          axis.text.x = element_text(size = 8), 
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size = 8), 
          plot.title = element_text(size = 8))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  
  ggsave(paste0(ref,"GOr2L5f.tiff"), width = 8, height = 8, units = "in", dpi = 300)
  ggsave(paste0(ref,"GOr2L5f.png"), width = 8, height = 8, units = "in", dpi = 300)
  
  #Significant function
  plotdata1<-plotdata
  plotdata.sig<-plotdata1 %>% group_by(Groups) %>% top_n(5, enrichfold)
  plotdata.sig<-plotdata.sig[,c("Groups","ID", "Description","ONTOLOGY")]
  plotdata.sig$ONTOLOGY<-gsub("\n"," ", plotdata.sig$ONTOLOGY)
  plotdata.sig$Groups<-gsub("\n"," ", plotdata.sig$Groups)
  write.table(plotdata.sig, 
              paste0(ref,"GO.topbyfold.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  plotdata.sig<-plotdata1 %>% group_by(Groups) %>% top_n(5, enrichfold)
  plotdata.sig$Description<-gsub(", with incorporation or reduction of molecular oxygen", "", plotdata.sig$Description)
  plotdata.sig$Description<-gsub(", NAD\\(P\\)H as one donor, and incorporation of one atom of oxygen", "", plotdata.sig$Description)
  
  #plotdata.sig$Description<-gsub("phosphotransferase activity, ", "phosphotransferase activity,\n", plotdata.sig$Description)
  #plotdata.sig$Groups<-gsub("\n"," ", plotdata.sig$Groups)
  #plotdata.sig$ONTOLOGY<-gsub("\n"," ", plotdata.sig$ONTOLOGY)
  a<-plotdata.sig %>%
    mutate(Groups = factor(Groups, levels = path.order)) %>%
    ggplot(aes(x = Groups, y = Description, size = ratio2, colour = log10(enrichfold)))+
    labs(title = "", size = "Gene ratio", colour = "Enrich factor", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(high = "#FF0000", low = "#0000FF")+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 8), 
          axis.text.x = element_text(size = 8, angle = 45,hjust = 1,vjust = 1), 
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size = 8), 
          legend.text = element_text(size = 8), 
          legend.position = "none",
          plot.title = element_text(size = 8))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  
  f = paste0(ref,"GOr2L5sig.pptx")
  topptx(a, f, width = 6, height = 5, units = "in")
  ggsave(paste0(ref,"GOr2L5sig.tiff"), width = 6, height = 5, units = "in", dpi = 300)
  ggsave(paste0(ref,"GOr2L5sig.png"), width = 6, height = 5, units = "in", dpi = 300)
  
  KOG.all.1<-compareCluster(Genes ~ Groups, 
                            data = AS.diflist.en, 
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
    mutate(Groups = factor(Groups, levels = path.order)) %>%
    ggplot(aes(x = Groups,
               y = paste0(Description," (",BGnumerator,")"),
               size = ratio1, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,1))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 9), 
          axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1), 
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size = 9), 
          plot.title = element_text(size = 9))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave(paste0(ref,"KOG.tiff"), width = 8, height = 6, units = "in", dpi = 300)
  ggsave(paste0(ref,"KOG.png"), width = 8, height = 6, units = "in", dpi = 300)
  
}

##stage specific and differential expressed
FC<-read.delim("C:/coprinopsis/expression/DEG.allpairs.txt", header = TRUE)
FC<-FC[abs(FC$logFC)>1 & FC$FDR<0.05, ]
FC$Groups<-NA
FC$Groups[FC$sampleB %in% c("BS", "BS12h") & 
                       FC$sampleA %in% c("BS12h","BS24h")]<-"Spore germination"
FC$Groups[FC$sampleB %in% c("Myc") & 
                       FC$sampleA %in% c("Oidia")]<-"Oidiation"
FC$Groups[FC$sampleB %in% c("Myc") & 
                       FC$sampleA %in% c("Scl")]<-"Sclerotia formation"
FC$Groups[FC$sampleB %in% c("Myc", "Knot") & 
                       FC$sampleA %in% c("Knot","Pri")]<-"Fruiting"
FC$Groups[FC$sampleB %in% c("Pri", "YFB") & 
                       FC$sampleA %in% c("YFB","BS")]<-"Sporulation"
AS.diflist.en.spc<-merge(AS.diflist.en.spc, FC[is.na(FC$Groups)==F, c("Genes","Groups","logFC")],
                         by = c("Genes","Groups"), all = F)

##PSI level + gene age + dNdS + stage
setwd("/home/yichun/RNAmodification/AS_cash")
treat.order<-c("Scl",
               "Oidia",
               "BS",
               "BS12h",
               "BS24h",
               "Myc",
               "Knot",
               "Pri",
               "YFB")

AS.new<-read.delim("Alternative_splicing.exp5.minor05.jad6.txt", header = T)
PSI.matrix<-read.delim("Alternative_splicing.exp5.minor05.jad6.psi.txt", header = T)
AS.feature<-read.delim("Alternative_splicing.exp5.minor05.jad6.feature.txt", header = T)
AS.diflist<-read.delim("AS.cashtest.all.txt", header = TRUE)
AS.diflist<-AS.diflist[AS.diflist$FDR<0.05 & abs(AS.diflist$delta_PSI) > 0.1,]
AS.diflist<-merge(AS.diflist, AS.feature[,c("Location","SplicingType","checkgeneID")],
                  by = c("Location","SplicingType"), all.x = TRUE)
names(AS.diflist)[names(AS.diflist)=="checkgeneID"]<-"Genes"

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

AS.phylo<-AS.diflist[,c("Genes","delta_PSI","SplicingType")]
AS.phylo<-merge(AS.phylo, genes.PS, by="Genes", all.x = T)
AS.phylo<-merge(AS.phylo, genes.NS[,c("Genes","NS")], by="Genes", all.x = T)
AS.phylo<-AS.phylo[is.na(AS.phylo$PS)==F,]
AS.phylo$PS<-paste0("PS",AS.phylo$PS)
AS.phylo$Mtype="Alternative splicing"
names(AS.phylo)[names(AS.phylo)=="delta_PSI"]<-"Frequency"

tab.size<-as.data.frame(xtabs(~PS, unique(AS.phylo[,c("Genes","PS")])))
genes.PS.stat<-as.data.frame(xtabs(~PS, genes.PS))
genes.PS.stat$PS<-paste0("PS",genes.PS.stat$PS)
tab.size<-merge(tab.size,genes.PS.stat, by = "PS", all.x = T)
tab.size$dotsize<-tab.size$Freq.x/tab.size$Freq.y
AS.phylo<-merge(AS.phylo, tab.size[,c("PS","dotsize")], by = "PS", all.x = T)
write.table(AS.phylo, "AS.phylo.txt", sep = "\t", quote = F, row.names = F)

a<-AS.phylo %>%
  mutate(PS = factor(PS, levels = paste0("PS",1:12))) %>%
  ggplot(aes(x = Frequency, y = NS, size = dotsize, alpha = dotsize))+
  geom_point(shape = 1)+
  scale_colour_manual(values = "black")+
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,1))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,1))+
  labs(x = "Delta PSI", y = "dN/dS", legend = "")+
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

ggsave("AS.geneage.png", width = 3, height = 6.5, units = "in", dpi = 300)
ggsave("AS.geneage.tiff", width = 3, height = 6.5, units = "in", dpi = 300)

a<-AS.phylo %>%
  mutate(PS = factor(PS, levels = paste0("PS",1:12))) %>%
  ggplot(aes(x = Frequency, y = NS, size = dotsize, alpha = dotsize))+
  geom_point(shape = 16)+
  scale_colour_manual(values = "black")+
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,1))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,1))+
  labs(x = "Delta PSI", y = "dN/dS", legend = "")+
  facet_grid(PS~SplicingType, scales = "free", space = "free")+
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

ggsave("AS.geneage.type.png", width = 6.5, height = 6.5, units = "in", dpi = 300)
ggsave("AS.geneage.type.tiff", width = 6.5, height = 6.5, units = "in", dpi = 300)

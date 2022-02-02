#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plotly))

##Create parameters
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="cash.alldifff.txt [default %default]",
              dest="input_filename"),
  make_option(c("-r","--ref"), type="character", default=NULL,
              help="reference.gff3 [default %default]",
              dest="gff"),
  make_option(c("-o","--output"), type="character", default="filtered",
              help="output filename prefix [default %default]",
              dest="output_filename"),
  make_option(c("-a", "--xspecies"), type="character", default="x.species",
              help="xsnpindel.filter.strict.txt [default %default]",
              dest="vcf.x"),
  make_option(c("-b", "--yspecies"), type="character", default="y.species",
              help="ysnpindel.filter.strict.txt [default %default]",
              dest="vcf.y"),
  make_option(c("-m", "--min-jad"), type="numeric", default=6,
              help="junction alignment distance (JAD), minimum mismatch near [default %default]",
              dest="jad")
)
options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

jad=opt$jad

#opt$input_filename="AS_cash/Myc_Scl.SclvsMyc.alldiff.txt"
#opt$gff="genome/annotate_results/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gff3"
#opt$vcf.x="mergebam/Myc.snpindel.filter.strict.txt"
#opt$vcf.y="mergebam/Scl.snpindel.filter.strict.txt"

##Read in AS from CASH
AS<-read.delim(opt$input_filename, header = TRUE)
AS<-separate(AS, Location, c("CHROM","POS"), sep = ":")
AS<-separate(AS, POS, c("Start","End"), sep = "-", remove = FALSE, convert = TRUE)
AS$Start<-as.numeric(AS$Start)
AS$End<-as.numeric(AS$End)
AS<-AS[AS$P.Value< 0.05 & AS$FDR < 0.05,]
row.names(AS)<-1:nrow(AS)

##Read in SNP and INDEL
vcf.x<-read.delim(opt$vcf.x, header = TRUE)
vcf.y<-read.delim(opt$vcf.y, header = TRUE)
vcf<-rbind(vcf.x, vcf.y)
rm(vcf.x, vcf.y)
vcf<-vcf[order(vcf$CHROM, vcf$POS),c("CHROM","POS")]
vcf<-unique(vcf)

##Read in gff3 annotation
gff<-read.delim(opt$gff, skip = 1, header = FALSE)
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

rm(exon.tmp,exon.tmp.chr,transcript)

##Remove intergenic and overlapping
AS<-AS[AS$checkgeneID != "Intergenic",]

row.names(AS)<-1:nrow(AS)
##Filtering
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

al<-AS[AS$SNPINDEL=="NO",]
names(al)[7:10]<-c("Mut_Junc_Inclusive..Exclusive", "Ref_Junc_Inclusive..Exclusive",
                   "Mut_Exp_Inclusive..Exclusive", "Ref_Exp_Inclusive..Exclusive")
al<-separate(al, Mut_Junc_Inclusive..Exclusive, 
             c("Mut_Junc_Inclusive","Mut_Junc_Exclusive"),
             sep = "::", remove = T, convert = T)
al<-separate(al, Ref_Junc_Inclusive..Exclusive, 
             c("Ref_Junc_Inclusive","Ref_Junc_Exclusive"),
             sep = "::", remove = T, convert = T)
al<-separate(al, Mut_Exp_Inclusive..Exclusive, 
             c("Mut_Exp_Inclusive","Mut_Exp_Exclusive"),
             sep = "::", remove = T, convert = T)
al<-separate(al, Ref_Exp_Inclusive..Exclusive, 
             c("Ref_Exp_Inclusive","Ref_Exp_Exclusive"),
             sep = "::", remove = T, convert = T)
al$Mut_Junc_Inclusive_ratio = al$Mut_Junc_Exclusive/(al$Mut_Junc_Exclusive+al$Mut_Junc_Inclusive)
al$Ref_Junc_Inclusive_ratio = al$Ref_Junc_Exclusive/(al$Ref_Junc_Exclusive+al$Ref_Junc_Inclusive)


##At least 5 supporting the junction, either sample is fine
##Minor type > 5 %
al.filter<-al[(al$Mut_Junc_Inclusive > 5 | al$Ref_Junc_Inclusive > 5) & 
                (al$Mut_Junc_Exclusive >5 | al$Ref_Junc_Exclusive >5) &
                ((al$Mut_Junc_Inclusive_ratio > 0.05 & al$Mut_Junc_Inclusive_ratio < 0.95) |
                   (al$Ref_Junc_Inclusive_ratio > 0.05 & al$Ref_Junc_Inclusive_ratio < 0.95)),]

##Save results
write.table(al.filter, paste0(opt$output_filename,".", opt$input_filename), 
            row.names = FALSE, quote = FALSE, sep = "\t")
cat(paste0("After filtering, ", nrow(al.filter), " detected, results saved in ",
            paste0(opt$output_filename,".", opt$input_filename), ".\n"))

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plotly))

##Create parameters
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="SNPINDEL.filter.vcf [default %default]",
              dest="input_filename"),
  make_option(c("-o","--output"), type="character", default="strict",
              help="output filename prefix [default %default]",
              dest="output_filename")
)
options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

#opt$input_filename="mergebam/Oidia.snpindel.filter.vcf"

##read in SNP INDEL list and filter by 10 reads & > 10 %
vcf.x<-read.delim(opt$input_filename, stringsAsFactors = F, header = T)
vcf.x<-separate(vcf.x, INFO, c("Type", "I16"), sep = "I16=")
vcf.x<-separate(vcf.x, I16, c("ref.plus","ref.minor","sub.plus","sub.minor"), sep = ",", convert = TRUE)
vcf.x$type<-"SNP"
vcf.x$type[grep(pattern = "INDEL", vcf.x$Type)]<-"INDEL"
names(vcf.x)[1]="CHROM"
vcf.x$var.depth=(vcf.x$sub.plus+vcf.x$sub.minor)
vcf.x$var.ratio=(vcf.x$sub.plus+vcf.x$sub.minor)/(vcf.x$sub.plus+vcf.x$sub.minor+vcf.x$ref.plus+vcf.x$ref.minor)
vcf.filter.x<-vcf.x[vcf.x$var.depth>10 | (vcf.x$var.depth>5 & vcf.x$var.ratio > 0.1), 
                    c("CHROM","POS","REF","ALT",
                      "type","var.depth","var.ratio",
                      "ref.plus","ref.minor","sub.plus","sub.minor")]
row.names(vcf.filter.x)<-1:nrow(vcf.filter.x)
write.table(vcf.filter.x, 
            paste0(gsub(".vcf","", opt$input_filename),".",opt$output_filename,".txt"),
            row.names = F, sep = "\t", quote = F)
cat(paste0(opt$input_filename, " filtered, with ", nrow(vcf.filter.x), "records.\n"))
cat(paste0("result saved in ", gsub(".vcf","", opt$input_filename),".",opt$output_filename,".txt"))

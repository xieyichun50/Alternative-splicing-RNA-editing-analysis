#!/bin/bash
source /PARA/app/scripts/cn-module.sh
module load Bowtie2/2.2.6
module load hisat2/2.0.5
module load samtools/1.5
module load java/jdk1.8.0_11
module load picard/2.9.0
export LD_LIBRARY_PATH=/PARA/pp585/Yichun/lib:$LD_LIBRARY_PATH

#hisat2-build -f Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa -p 24

cat /PARA/pp585/Yichun/filelist_RNA | while read i;
do
echo "hisat2 -x /PARA/pp585/Yichun/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa -1 /PARA/pp585/Yichun/fastp/$i.r1.fastp.fq -2 /PARA/pp585/Yichun/fastp/$i.r2.fastp.fq -S /PARA/pp585/Yichun/bam/$i.sam --phred33 --dta-cufflinks --novel-splicesite-outfile /PARA/pp585/Yichun/bam/splice/$i.splicesite.txt  -p 24"
hisat2 -x /PARA/pp585/Yichun/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa -1 /PARA/pp585/Yichun/fastp/$i.r1.fastp.fq -2 /PARA/pp585/Yichun/fastp/$i.r2.fastp.fq -S /PARA/pp585/Yichun/bam/$i.sam --phred33 --dta-cufflinks --novel-splicesite-outfile /PARA/pp585/Yichun/bam/splice/$i.splicesite.txt  -p 24

echo "samtools view -b -S /PARA/pp585/Yichun/bam/$i.sam > /PARA/pp585/Yichun/bam/$i.bam"
samtools view --threads 24 -b -S /PARA/pp585/Yichun/bam/$i.sam > /PARA/pp585/Yichun/bam/$i.bam
echo "/PARA/pp585/Yichun/bam/$i.bam generated"

rm /PARA/pp585/Yichun/bam/$i.sam

echo "samtools sort /PARA/pp585/Yichun/bam/$i.bam -o /PARA/pp585/Yichun/bam/$i.sorted.bam"
samtools sort --threads 24 /PARA/pp585/Yichun/bam/$i.bam -o /PARA/pp585/Yichun/bam/$i.sorted.bam

####Remove duplicate reads (generate by PCR, etc.)

echo "java -Xmx40g -jar /WORK/app/picard/2.9.0/picard.jar MarkDuplicates INPUT=/PARA/pp585/Yichun/bam/$i.sorted.bam OUTPUT=/PARA/pp585/Yichun/redubam/$i.sorted.redu.bam METRICS_FILE=/PARA/pp585/Yichun/redubam/$i.sorted.redu.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true"

java -Xmx40g -jar /WORK/app/picard/2.9.0/picard.jar MarkDuplicates INPUT=/PARA/pp585/Yichun/bam/$i.sorted.bam OUTPUT=/PARA/pp585/Yichun/redubam/$i.sorted.redu.bam METRICS_FILE=/PARA/pp585/Yichun/redubam/$i.sorted.redu.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true

samtools index /PARA/pp585/Yichun/redubam/$i.sorted.redu.bam
echo "/PARA/pp585/Yichun/redubam/$i.sorted.redu.bam indexed"
done

wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff

vim snpEff.config
#Coprinopsis_cinerea_326
Cc326.genome:Coprinopsis cinerea 326

mkdir data
cd data/
mkdir Cc326
mkdir genomes
ln -s /home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.clean.gff3 Cc326/genes.gff
ln -s /home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.proteins.fa Cc326/protein.fa
ln -s /home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.cds-transcripts.fa Cc326/cds.fa
ln -s /home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa genomes/Cc326.fa
cd ..
java -jar snpEff.jar build -gff3 -v Cc326

java -jar /home/yichun/tools/snpEff/snpEff.jar Cc326 RE.snpEff.vcf > RE.snpEff.anno.vcf
gff=/home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gff3
dir=DEG

cat samples_n_reads_decribed.txt | while read id1 id2;
do  
	#hisat2 -p 8 --dta -x $ref -1 $fq1 -2 $fq2 --rf -S $id2.sam &>$id2.hisat2.log 
	#samtools sort -@ 8 -o $id2.bam $id2.sam 
	#stringtie -p 8 -G $ref.gtf -o $id2.stringtie.gtf -l $id2 $id2.bam
	stringtie -e -B -p 14 -G $gff -A ${id2}_gene_count.xls -o ballgown/${id2}_ballgown/$id2.gtf /home/yichun/RNAmodification/redubam/$id2.sorted.redu.bam
	#rm $id2.sam $id2.bam
done
perl /home/yichun/RNAmodification/stringtie_gene_count_matrix.pl *_gene_count.xls

prepDE.py -i ballgown/ -l 150
sed -i 's/,/\t/g;s/_ballgown//g' gene_count_matrix.csv transcript_count_matrix.csv

##file ready
##
Rscript /home/yichun/RNAmodification/edgeR_TMM_norm.R

cd ${dir};
/home/yichun/tools/trinityrnaseq/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix gene_count_matrix.csv --method edgeR --samples_file samples_n_reads_decribed.txt --output edgeR_gene.min_reps2.min_cpm1 --min_reps_min_cpm 2,1
cd ${dir}/edgeR_gene.min_reps2.min_cpm1;
/home/yichun/tools/trinityrnaseq/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../gene_count_matrix.csv -P 0.05 -C 1 --samples ../samples_n_reads_decribed.txt

find  ballgown/ -name '*.gtf' | while read i;
do 
	cp $i .;
done

perl /home/yichun/RNAmodification/stringtie_gtf_transcript_count_matrix.pl $( cut -f 2 samples_n_reads_decribed.txt | while read i; do echo -ne "$i.gtf ";done)
sed -i 's/.gtf//g' transcript_count_matrix.*


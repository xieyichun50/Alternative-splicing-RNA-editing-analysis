cat pairlist | while read x y;
do

Rscript /home/yichun/RNAmodification/AS_cash/AS_filtering.R --input ${x}_${y}.${y}vs${x}.alldiff.txt --ref /home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gff3 -a /home/yichun/RNAmodification/mergebam/${x}.snpindel.filter.strict.txt -b /home/yichun/RNAmodification/mergebam/${y}.snpindel.filter.strict.txt --min-jad 6

done

#!/bin/bash
source /PARA/app/scripts/cn-module.sh
module load python-tools/2.7.10
module load samtools/1.5
#module load bcftools/1.3.1
module load blat/36
export LD_LIBRARY_PATH=/PARA/pp585/Yichun/lib:$LD_LIBRARY_PATH

cat /PARA/pp585/Yichun/filelist_RNA | while read i ;
do
echo "python /PARA/pp585/Yichun/tools/REDItools-1.0.4/reditools/REDItoolBlatCorrection.py -i /PARA/pp585/Yichun/redubam/$i.sorted.redu.bam -f /PARA/pp585/Yichun/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa -F /PARA/pp585/Yichun/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa.2bit -o /PARA/pp585/Yichun/Blatcorrection/Blatcorrection_$i -V -T -t 24"

#python /PARA/pp585/Yichun/tools/REDItools-1.0.4/reditools/REDItoolBlatCorrection.py -i /PARA/pp585/Yichun/redubam/${i}1.sorted.redu.bam -f /PARA/pp585/Yichun/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa -F /PARA/pp585/Yichun/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa.2bit -o /PARA/pp585/Yichun/Blatcorrection/Blatcorrection_${i}1 -V -T -t 24 &
#python /PARA/pp585/Yichun/tools/REDItools-1.0.4/reditools/REDItoolBlatCorrection.py -i /PARA/pp585/Yichun/redubam/${i}2.sorted.redu.bam -f /PARA/pp585/Yichun/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa -F /PARA/pp585/Yichun/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa.2bit -o /PARA/pp585/Yichun/Blatcorrection/Blatcorrection_${i}2 -V -T -t 24 &
#python /PARA/pp585/Yichun/tools/REDItools-1.0.4/reditools/REDItoolBlatCorrection.py -i /PARA/pp585/Yichun/redubam/${i}3.sorted.redu.bam -f /PARA/pp585/Yichun/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa -F /PARA/pp585/Yichun/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa.2bit -o /PARA/pp585/Yichun/Blatcorrection/Blatcorrection_${i}3 -V -T -t 24 &


done
echo "Blatcorrection Finished"

#!/bin/bash
source /PARA/app/scripts/cn-module.sh

#module load python-tools/2.7.10
module load samtools/1.5
#module load bcftools/1.3.1
module load blat/36

export LD_LIBRARY_PATH=/WORK/paratera_share/python/2.7.10.bak/lib:$LD_LIBRARY_PATH
export PATH=/WORK/paratera_share/python/2.7.10.bak/bin:$PATH

export LD_LIBRARY_PATH=/PARA/pp585/Yichun/lib:$LD_LIBRARY_PATH

cat /PARA/pp585/Yichun/filelist_RNA | while read i ;
do
echo "python /PARA/pp585/Yichun/tools/REDItools-1.0.4/reditools/REDItoolDnaRna.py -i /PARA/pp585/Yichun/redubam/$i.sorted.redu.bam -b /PARA/pp585/Yichun/Blatcorrection/Blatcorrection_$i -j /PARA/pp585/Yichun/redubam/cc326DNA.sorted.redu.bam -B /PARA/pp585/Yichun/Blatcorrection/Blatcorrection_cc326DNA -f /PARA/pp585/Yichun/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa -t 24 -o /PARA/pp585/Yichun/REDIresult/REDI_$i -c 10,10 -Q 33,33 -q 25,25 -m 25,25 -G /PARA/pp585/Yichun/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gff3.gz -e -E -p -P -u -U -v 3 -n 0.01 -N 0 -z -Z "

python /PARA/pp585/Yichun/tools/REDItools-1.0.4/reditools/REDItoolDnaRna.py -i /PARA/pp585/Yichun/redubam/$i.sorted.redu.bam -b /PARA/pp585/Yichun/Blatcorrection/Blatcorrection_$i -j /PARA/pp585/Yichun/redubam/cc326DNA.sorted.redu.bam -B /PARA/pp585/Yichun/Blatcorrection/Blatcorrection_cc326DNA -f /PARA/pp585/Yichun/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa -t 24 -o /PARA/pp585/Yichun/REDIresult/REDI_$i -c 10,10 -Q 33,33 -q 25,25 -m 25,25 -G /PARA/pp585/Yichun/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gff3.gz -e -E -p -P -u -U -v 3 -n 0.01 -N 0 -z -Z 
done

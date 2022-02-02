dir=/home/yichun/RNAmodification/redubam
suffix=sorted.redu.bam
cat treatlist | while read treat sample;

do

echo "
java -jar -Xmx50g /home/yichun/tools/cash_v2.2.1/cash.jar --MergePval G --Combine True --DisplayAllEvent True --StrandSpecific NONE --SpliceCons True --JuncAllSample 25 --JuncOneGroup 10 --minAnchorLen 5 --minIntronLen 25 --minJuncReadsForNewIso 10 --runSepChr True --GTF /home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gff3 --Case:Myc ${dir}/M1.${suffix},${dir}/M2.${suffix},${dir}/M3.${suffix} --Control:${treat} ${dir}/${sample}1.${suffix},${dir}/${sample}2.${suffix},${dir}/${sample}3.${suffix} --Output ${treat}_Myc > ${treat}_Myc.log &

java -jar -Xmx50g /home/yichun/tools/cash_v2.2.1/cash.jar --MergePval G --Combine True --DisplayAllEvent True --StrandSpecific NONE --SpliceCons True --JuncAllSample 25 --JuncOneGroup 10 --minAnchorLen 5 --minIntronLen 25 --minJuncReadsForNewIso 10 --runSepChr True --GTF /home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gff3 --Case:Oidia ${dir}/O1.${suffix},${dir}/O2.${suffix},${dir}/O3.${suffix} --Control:${treat} ${dir}/${sample}1.${suffix},${dir}/${sample}2.${suffix},${dir}/${sample}3.${suffix} --Output ${treat}_Oidia > ${treat}_Oidia.log &

java -jar -Xmx50g /home/yichun/tools/cash_v2.2.1/cash.jar --MergePval G --Combine True --DisplayAllEvent True --StrandSpecific NONE --SpliceCons True --JuncAllSample 25 --JuncOneGroup 10 --minAnchorLen 5 --minIntronLen 25 --minJuncReadsForNewIso 10 --runSepChr True --GTF /home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gff3 --Case:Scl ${dir}/S1.${suffix},${dir}/S2.${suffix},${dir}/S3.${suffix} --Control:${treat} ${dir}/${sample}1.${suffix},${dir}/${sample}2.${suffix},${dir}/${sample}3.${suffix} --Output ${treat}_Scl > ${treat}_Scl.log &

java -jar -Xmx50g /home/yichun/tools/cash_v2.2.1/cash.jar --MergePval G --Combine True --DisplayAllEvent True --StrandSpecific NONE --SpliceCons True --JuncAllSample 25 --JuncOneGroup 10 --minAnchorLen 5 --minIntronLen 25 --minJuncReadsForNewIso 10 --runSepChr True --GTF /home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gff3 --Case:Knot ${dir}/K1.${suffix},${dir}/K2.${suffix},${dir}/K3.${suffix} --Control:${treat} ${dir}/${sample}1.${suffix},${dir}/${sample}2.${suffix},${dir}/${sample}3.${suffix} --Output ${treat}_Knot > ${treat}_Knot.log &

java -jar -Xmx50g /home/yichun/tools/cash_v2.2.1/cash.jar --MergePval G --Combine True --DisplayAllEvent True --StrandSpecific NONE --SpliceCons True --JuncAllSample 25 --JuncOneGroup 10 --minAnchorLen 5 --minIntronLen 25 --minJuncReadsForNewIso 10 --runSepChr True --GTF /home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gff3 --Case:Pri ${dir}/P1.${suffix},${dir}/P2.${suffix},${dir}/P3.${suffix} --Control:${treat} ${dir}/${sample}1.${suffix},${dir}/${sample}2.${suffix},${dir}/${sample}3.${suffix} --Output ${treat}_Pri > ${treat}_Pri.log &

java -jar -Xmx50g /home/yichun/tools/cash_v2.2.1/cash.jar --MergePval G --Combine True --DisplayAllEvent True --StrandSpecific NONE --SpliceCons True --JuncAllSample 25 --JuncOneGroup 10 --minAnchorLen 5 --minIntronLen 25 --minJuncReadsForNewIso 10 --runSepChr True --GTF /home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gff3 --Case:YFB ${dir}/YFB1.${suffix},${dir}/YFB2.${suffix},${dir}/YFB3.${suffix} --Control:${treat} ${dir}/${sample}1.${suffix},${dir}/${sample}2.${suffix},${dir}/${sample}3.${suffix} --Output ${treat}_YFB > ${treat}_YFB.log &

java -jar -Xmx50g /home/yichun/tools/cash_v2.2.1/cash.jar --MergePval G --Combine True --DisplayAllEvent True --StrandSpecific NONE --SpliceCons True --JuncAllSample 25 --JuncOneGroup 10 --minAnchorLen 5 --minIntronLen 25 --minJuncReadsForNewIso 10 --runSepChr True --GTF /home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gff3 --Case:BS ${dir}/BS1.${suffix},${dir}/BS2.${suffix},${dir}/BS3.${suffix} --Control:${treat} ${dir}/${sample}1.${suffix},${dir}/${sample}2.${suffix},${dir}/${sample}3.${suffix} --Output ${treat}_BS > ${treat}_BS.log &

java -jar -Xmx50g /home/yichun/tools/cash_v2.2.1/cash.jar --MergePval G --Combine True --DisplayAllEvent True --StrandSpecific NONE --SpliceCons True --JuncAllSample 25 --JuncOneGroup 10 --minAnchorLen 5 --minIntronLen 25 --minJuncReadsForNewIso 10 --runSepChr True --GTF /home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gff3 --Case:BS12h ${dir}/BS12h1.${suffix},${dir}/BS12h2.${suffix},${dir}/BS12h3.${suffix} --Control:${treat} ${dir}/${sample}1.${suffix},${dir}/${sample}2.${suffix},${dir}/${sample}3.${suffix} --Output ${treat}_BS12h > ${treat}_BS12h.log &

java -jar -Xmx50g /home/yichun/tools/cash_v2.2.1/cash.jar --MergePval G --Combine True --DisplayAllEvent True --StrandSpecific NONE --SpliceCons True --JuncAllSample 25 --JuncOneGroup 10 --minAnchorLen 5 --minIntronLen 25 --minJuncReadsForNewIso 10 --runSepChr True --GTF /home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gff3 --Case:BS24h ${dir}/BS24h1.${suffix},${dir}/BS24h2.${suffix},${dir}/BS24h3.${suffix} --Control:${treat} ${dir}/${sample}1.${suffix},${dir}/${sample}2.${suffix},${dir}/${sample}3.${suffix} --Output ${treat}_BS24h > ${treat}_BS24h.log &

wait
"
done

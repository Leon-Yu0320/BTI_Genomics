#!/bin/bash

genome_DIR="/home/liangyu/1_Project/8_Sorghum/6_Annotation"
work_dir="/home/liangyu/1_Project/8_Sorghum/7_Motif"

#grep 'mRNA' $genome_DIR/Sorghum.gff_primary_transcripts.gff3 > $genome_DIR/Sorghum_PT-mRNA.gff3

#for ID in $(cat $work_dir/Gene_model.list);
#do
        #grep "$ID" $genome_DIR/Sorghum_PT-mRNA.gff3 | cut -f 1,4,5,9 >> $work_dir/coordindate.list  

#done

#cut -d ':' -f1 coordindate.list > coordindate_clean.list 
#sed -i 's/ID=//g' coordindate_clean.list 

#cat coordindate_clean.list |awk -F ' ' '{print($2-1000)}' > temp
#paste coordindate_clean.list temp > coordinate.final

#awk '{ print $1 ":" $5 "-" $2}' coordinate.final >  coordinate.txt

for LINE in $(cat $work_dir/coordinate.txt);
do 
        samtools faidx $genome_DIR/Sbicolor_454.fa $LINE >> Promoter_seq.fasta
done

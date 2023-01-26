#!/bin/bash

Ref_dir="/home/liangyu/0_data/1_Genome"
workdir="/home/liangyu/1_Project/5_Cotton/3_Varaiants_calling/0_trail" 
BAM_dir="/home/liangyu/1_Project/5_Cotton/3_Varaiants_calling/0_BAM" 

# Other settings
nt=32 #number of threads to use in computation



#Build index and dict
cd $Ref_dir
samtools faidx $Ref_dir/Ghirsutum_527_v2.0.fa
#bwa index $Ref_genome 
picard CreateSequenceDictionary -R $Ref_dir/Ghirsutum_527_v2.0.fa -O $Ref_dir/Ghirsutum_527_v2.0.dict

mkdir $workdir
cd $workdir

#loop settings
for file in 203.R1 204.R1 205.R1 206.R1 207.R1 208.R1 209.R1 210.R1 211.R1 212.R1 213.R1 214.R1 215.R1 216.R1 217.R1 218.R1 219.R1 220.R1 221.R1 222.R1 223.R1 224.R1 225.R1 226.R1 228.R1 229.R1 230.R1 231.R1 232.R1 233.R1 234.R1 235.R1 236.R1 237.R1 238.R1 239.R1 240.R1 241.R1 242.R1 244.R1 245.R1;
do
 	mkdir $file
  	cd $file

	#Mark and remove duplictaes from Bam
	picard MarkDuplicates -REMOVE_DUPLICATES true -I "/home/liangyu/1_Project/5_Cotton/3_Varaiants_calling/0_BAM/${file}.sorted.bam" -O ${file}.sort.dup.bam -M ${file}.metrics.txt 

	picard AddOrReplaceReadGroups -I ${file}.sort.dup.bam -O ${file}.sort.dup.Gr.bam -LB cotton -PL illumina -PU WL -SM ${file}

	#index bam file
	samtools index ${file}.sort.dup.Gr.bam

	#Transform the NCigarreads
	java -Xmx128G -jar /home/liangyu/anaconda3/opt/gatk-3.8/GenomeAnalysisTK.jar -T SplitNCigarReads -R $Ref_dir/Ghirsutum_527_v2.0.fa -I ${file}.sort.dup.Gr.bam -o ${file}.sort.dup.Gr_NC.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

	java -Xmx128G -jar /home/liangyu/anaconda3/opt/gatk-3.8/GenomeAnalysisTK.jar -T HaplotypeCaller -R $Ref_dir/Ghirsutum_527_v2.0.fa -I ${file}.sort.dup.Gr_NC.bam -nct 3 --emitRefConfidence GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -o ${file}.gvcf
	mv ${file}.gvcf $workdir
	cd ..
done

cd $workdir
ls | grep ".gvcf" | sort > gvcf.list

gatk3 -T CombineGVCFs -R $Ref_dir/Ghirsutum_527_v2.0.fa --variant $workdir/gvcf.list -o cotton.cohort.g.vcf
gatk3 -T GenotypeGVCFs -R $Ref_dir/Ghirsutum_527_v2.0.fa --variant $workdir/cohort.g.vcf -o cotton_final.vcf


#Filter VCF based on quality, depth, minor allele frequecny
vcftools --vcf $workdir/cotton_final.vcf --max-missing 0.90 --maf 0.05 --min-alleles 2 --max-alleles 2 --minQ 20 --minDP 10 --recode --recode-INFO-all --out $workdir/cotton_final_filter.vcf
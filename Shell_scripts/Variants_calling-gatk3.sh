#! bin/bash


#Build index and dict
samtools faidx $Ref_dir/Ghirsutum_527_v2.0.fa
bwa index $Ref_dir/Ghirsutum_527_v2.0.fa 
picard CreateSequenceDictionary -R $Ref_dir/Ghirsutum_527_v2.0.fa -O $Ref_dir/Ghirsutum_527_v2.0.dict

### set the thread numbers
nt=14

#loop settings
for file in $(cat $workdir/sample.list);
do
  mkdir ${workdir}/${file}

  ## Map reads by bwa piepline
  bwa mem -M -R "@RG\tID:$file\tSM:$file\tPL:ILLUMINA" \
    -t $nt  \
    $Ref_dir/Ghirsutum_527_v2.0.fa \
    ${BAM_dir}/${file}_forward_paired.fq.gz \
    ${BAM_dir}/${file}_reverse_paired.fq.gz | samtools view -@10 -bS - -o ${workdir}/${file}/${file}.reorder.bam 

  samtools view -b -F 4 ${workdir}/${file}/${file}.reorder.bam > ${workdir}/${file}/${file}.cut.bam

  ## Sort bam files
  picard SortSam \
   INPUT=${workdir}/${file}/${file}.cut.bam \
   OUTPUT=${workdir}/${file}/${file}.sort.bam SORT_ORDER=coordinate

    #Mark and remove duplictaes from Bam
   picard MarkDuplicates -REMOVE_DUPLICATES true \
    -I ${workdir}/${file}/${file}.sort.bam \
    -O ${workdir}/${file}/${file}.sort.dup.bam \
    -M ${workdir}/${file}/${file}.metrics.txt 

 #add group numbers and infotrmation
  picard AddOrReplaceReadGroups \
    -I ${workdir}/${file}/${file}.sort.dup.bam \
    -O ${workdir}/${file}/${file}.sort.dup.Gr.bam -LB cotton -PL illumina -PU Diversity -SM ${file} &

   #index bam file
   samtools index ${workdir}/${file}/${file}.sort.dup.Gr.bam

  java -Xmx128G -jar $GATK_dir/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R $Ref_dir/Ghirsutum_527_v2.0.fa \
    -I ${workdir}/${file}/${file}.sort.dup.Gr.bam \
    -nct $nt --emitRefConfidence GVCF \
    -variant_index_type LINEAR \
    -variant_index_parameter 128000 \
    -o ${workdir}/${file}/${file}.gvcf

  mv ${workdir}/${file}/${file}.gvcf ${workdir}/

done

ls ${workdir} | grep ".gvcf" | sort > gvcf.list

### Generate the cohort VCFfiles
java -Xmx128G -jar $GATK_dir/GenomeAnalysisTK.jar \
 -T CombineGVCFs -nct $nt \
 -R $Ref_dir/Ghirsutum_527_v2.0.fa \
 --variant ${workdir}/gvcf.list 
 -o cotton.cohort.gvcf \

### Transform the VCFs format
java -Xmx128G -jar $GATK_dir/GenomeAnalysisTK.jar \
 -T GenotypeGVCFs -nct $nt \
 -R $Ref_dir/Ghirsutum_527_v2.0.fa \
 --variant $workdir/cohort.g.vcf -o cotton_DNA_genotype.vcf

### Filter the variants of multiple samples


workdir="/home/liangyu/1_Project/4_spinach/5_VCFs" 

### STEP.1 remove loci with missing data based on cut-off to reduce bias of total mapping depth
#keep loci with up to 10% missing data for plot
vcftools --vcf $workdir/F1_map_variants.vcf \
	--max-missing 0.90 --min-alleles 2 \
	--recode --recode-INFO-all \
	--out $workdir/F1_plot.vcf


### STEP.2 PLOT DP, QD, and QUAL for each loci 
###Extrat VCF depth information (This is the combined depth among samples)
egrep -v "^#" $workdir/F1_plot.vcf.recode.vcf | \
cut -f 8 | \
sed 's/^.*;DP=\([0-9]*\);.*$/\1/' > $workdir/F1_map_DP.txt


###Extrat VCF QualitybyDepth information 
egrep -v "^#" $workdir/F1_plot.vcf.recode.vcf | \
cut -f 8 | \
sed 's/^.*;QD=\([0-9]*.[0-9]*\);.*$/\1/' > $workdir/F1_map_QD.txt

###Extrat VCF Quality information
egrep -v "^#" $workdir/F1_plot.vcf.recode.vcf | \
cut -f 6 > $workdir/F1_map_QUAL.txt

#Plot three component
Rscript 1_Plot_depth_quality.R --DP $workdir/F1_map_DP.txt --QUAL $workdir/F1_map_QUAL.txt --QD $workdir/F1_map_QD.txt


### STEP.3 Remove loci based on cutoff derived from plot information 
## MinDP and MaxDP is determined by average total depth (total depth/population size)
## minGQ=20: 1% chance that the call is incorrect
## minQ=20: 1 % chance that there is no variant at the site

vcftools --vcf $workdir/F1_plot.vcf.recode.vcf  \
	--minDP 5 --maxDP 25 \
	--minGQ 20 --minQ 20 \
	--recode --recode-INFO-all \
	--out $workdir/F1.vcf

#filter based qualitybydepth (QD) (QD is the normalized quality based on mapping depth usually set as 2)
vcffilter -f "QD > 5" $workdir/F1.vcf.recode.vcf > $workdir/F1_final_filter.vcf
## Line NO. 1079049


### STEP 4. extract genotype information only 
cat $workdir/F1_final_filter.vcf | perl while (<>) { s!:\S+!!g; print; } > F1_final_filter.clean 

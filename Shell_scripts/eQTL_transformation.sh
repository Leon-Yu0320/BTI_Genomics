#! bin/bash

work_dir="/home/liangyu/1_Project/5_Cotton/10_eQTL"

### STEP.1 PLOT DP, QD, and QUAL for each loci 
###Extrat VCF depth information (This is the combined depth among samples)
egrep -v "^#" $workdir/cotton_final_2021_WL.vcf  | \
cut -f 8 | \
sed 's/^.*;DP=\([0-9]*\);.*$/\1/' > $workdir/Gh_WL_DP.txt

###Extrat VCF QualitybyDepth information 
egrep -v "^#" $workdir/cotton_final_2021_WL.vcf  | \
cut -f 8 | \
sed 's/^.*;QD=\([0-9]*.[0-9]*\);.*$/\1/' > $workdir/Gh_WL_QD.txt

###Extrat VCF Quality information
egrep -v "^#" $workdir/cotton_final_2021_WL.vcf  | \
cut -f 6 > $workdir/Gh_WL_QUAL.txt

#Plot three components
Rscript $workdir/1_Plot_depth_quality.R \
    --DP $workdir/Gh_WL_DP.txt \
    --QUAL $workdir/Gh_WL_QD.txt \
    --QD $workdir/F1_map_QD.txt


#filter by missing rate and mapping depth
vcftools --vcf $work_dir/cotton_final_2021_WL.vcf 
    --minDP 5 
    --max-missing 0.9 
    --recode --recode-INFO-all 
    --out $work_dir/cotton_final_2021_WL.filter.vcf


#### Transform into 0,1,2
vcftools --vcf $work_dir/cotton_final_2021_WL.filter.vcf.recode.vcf --012 --out $work_dir/snp_matrix

#Add sample ID and SNPs ID for each matrix
awk '{$1=null;print $0}' $work_dir/snp_matrix.012 > $work_dir/tmp.txt
paste -d" " $work_dir/snp_matrix.012.indv $work_dir/tmp.txt | sed 's/  / /g' > $work_dir/snp_matrix_indv.txt

#Extract SNP IDs from vcf 
echo "ID" > $work_dir/snpid.txt
grep -v "^#" $work_dir/cotton_final_2021_WL.filter.vcf.recode.vcf  | awk '{print $1"_"$2}' >> $work_dir/snpid.txt
cat $work_dir/snpid.txt | tr "\n" " " > $work_dir/snpid2.txt
echo "" >> snpid2.txt

#combine files
cat snpid2.txt snp_matrix_indv.txt> genotype.txt

#transpose the matrix
awk '{i=1;while(i <= NF){col[i]=col[i] $i " ";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' genotype.txt | sed 's/[ \t]*$//g' > cotton_final_2021_WL_transpose.txt

!# /bin/bash/

IFS=$'\n';
for LINE in $(cat $work_dir/Module.list);
do 
  NAME=$(echo $LINE | cut -d "." -f1 | cut -d "-" -f3)
	cat ${seq_dir}/$LINE | cut -d "," -f1 | grep -v 'node' > ${work_dir}/$NAME.list 

  ###extract promoter sequences to differnet moudles
  faSomeRecords ${work_dir}/Promoter_seq.fasta ${work_dir}/$NAME.list ${work_dir}/Sequences/$NAME.promoter.fasta
	
  ###Perform motif enrichment analysis 
	ame --verbose 1 --oc ${work_dir}/Results/$NAME --o ${work_dir}/Results/$NAME \
	--scoring avg --method fisher --hit-lo-fraction 0.25 --evalue-report-threshold 10.0 \
	--control --shuffle-- --kmer 2 ${work_dir}/Sequences/$NAME.promoter.fasta $work_dir/ArabidopsisDAPv1.meme

  ###revise name of files
	mv ${work_dir}/Results/$NAME/ame.tsv ${work_dir}/Results/${NAME}/$NAME.tsv
	mv ${work_dir}/Results/$NAME/ame.html ${work_dir}/Results/${NAME}/$NAME.html
	mv ${work_dir}/Results/$NAME/sequences.tsv ${work_dir}/Results/${NAME}/$NAME.sequence

		 
  ###Extract the gene ID with DAPseq signal from database 
	awk -v $NAME '{print $0}' $work_dir/AME-MAD85_TF2DAPseq > ${work_dir}/Results/${NAME}/$NAME.TF2DAPseq

	for i in $(cat ${work_dir}/Results/${NAME}/$NAME.TF2DAPseq | cut -f4);
	do
		grep $i ${work_dir}/Results/${NAME}/$NAME.tsv >> ${work_dir}/Results/${NAME}/$NAME.select
	done

  ###extract target sequences based on enriched motifs and only targets marked with "true positive will be kept"
	for i in $(cat ${work_dir}/Results/${NAME}/$NAME.select | cut -f3);
	do 
		grep $i ${work_dir}/Results/${NAME}/$NAME.sequence | awk '$6 == "tp" {print $0}'  >> ${work_dir}/Results/${NAME}/$NAME.peak
	done

  ###Create the clean file with information for cytoscape downstream analysis
	cut -f2,3 ${work_dir}/Results/${NAME}/$NAME.peak | cut -d "." -f2,3 | awk '$3 = "DNA-binding" {print $0}' > ${work_dir}/Results/${NAME}/$NAME.clean

done

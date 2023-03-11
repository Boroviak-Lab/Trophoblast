#~/./trim_galore --paired --fastqc cmTSC1.1.fq.gz  cmTSC1.2.fq.gz --output_dir /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed
#~/./trim_galore --paired --fastqc cmTSC2.1.fq.gz  cmTSC2.2.fq.gz --output_dir /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed
#~/./trim_galore --paired --fastqc cmTSC3.1.fq.gz  cmTSC3.2.fq.gz --output_dir /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed
#~/./trim_galore --paired --fastqc Naive1.1.fq.gz  Naive1.2.fq.gz --output_dir /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed
#~/./trim_galore --paired --fastqc Naive2.1.fq.gz  Naive2.2.fq.gz --output_dir /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed
#~/./trim_galore --paired --fastqc Naive3.1.fq.gz  Naive3.2.fq.gz --output_dir /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed
#~/./trim_galore --paired --fastqc Primed1.1.fq.gz  Primed1.2.fq.gz --output_dir /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed
#~/./trim_galore --paired --fastqc Primed2.1.fq.gz  Primed2.2.fq.gz --output_dir /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed
#~/./trim_galore --paired --fastqc Primed3.1.fq.gz  Primed3.2.fq.gz --output_dir /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed

#cd ~/bismark_v0.16.3 
#bismark_genome_preparation /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Sequence/WholeGenomeFasta/

#sdas

#cd ~/bismark_v0.16.3
#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Sequence/WholeGenomeFasta/ -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC1.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC1.2_val_2.fq.gz
#sadsadsa

#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC1.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC1.2_val_2.fq.gz
#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive1.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive1.2_val_2.fq.gz
#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Primed1.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Primed1.2_val_2.fq.gz


#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC2.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC2.2_val_2.fq.gz
#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive2.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive2.2_val_2.fq.gz
#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Primed2.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Primed2.2_val_2.fq.gz

#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Bisulfite_Genome -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC3.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC3.2_val_2.fq.gz
#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive3.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive3.2_val_2.fq.gz
#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Primed3.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Primed3.2_val_2.fq.gz

#~/bismark_v0.16.3/./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/allMeth/Processed/ 
#--genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -1 /mnt/scratch/gurdon/cap76/Thorsten/
#allMeth/fastqc_trimmed_results/${FILES[$SLURM_ARRAY_TASK_ID-1]}_1_val_1.fq.gz -2 /mnt/scratch/
#gurdon/cap76/Thorsten/allMeth/fastqc_trimmed_results/${FILES[$SLURM_ARRAY_TASK_ID-1]}_2_val_2.fq.gz

#~/bismark_v0.16.3/./deduplicate_bismark --bam -s /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC1.1_val_1_bismark_bt2_pe.bam
#~/bismark_v0.16.3/./deduplicate_bismark --bam -s /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC2.1_val_1_bismark_bt2_pe.bam
#~/bismark_v0.16.3/./deduplicate_bismark --bam -p /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC3.1_val_1_bismark_bt2_pe.bam

#~/bismark_v0.16.3/./deduplicate_bismark --bam -s /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive1.1_val_1_bismark_bt2_pe.bam
#~/bismark_v0.16.3/./deduplicate_bismark --bam -s /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive2.1_val_1_bismark_bt2_pe.bam
#~/bismark_v0.16.3/./deduplicate_bismark --bam -s /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive3.1_val_1_bismark_bt2_pe.bam

#~/bismark_v0.16.3/./deduplicate_bismark --bam -s /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Primed1.1_val_1_bismark_bt2_pe.bam
#~/bismark_v0.16.3/./deduplicate_bismark --bam -s /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Primed2.1_val_1_bismark_bt2_pe.bam
#~/bismark_v0.16.3/./deduplicate_bismark --bam -s /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Primed3.1_val_1_bismark_bt2_pe.bam



#/mnt/scratch/gurdon/cap76/Thorsten/allMeth/Processed/${FILES[$SLURM_ARRAY_TASK_ID-1]}_1_val_1_bismark_bt2_pe.bam ; 
#~/bismark_v0.16.3/./bismark_methylation_extractor -o /mnt/scratch/gurdon/cap76/Thorsten/allMeth/Processed/ --scaffolds --comprehensive --bedGraph --gzip /mnt/scratch/gurdon/cap76/Thorsten/allMeth/Processed/${FILES[$SLURM_ARRAY_TASK_ID-1]}_1_val_1_bismark_bt2_pe.deduplicated.bam"

~/bismark_v0.16.3/./bismark_methylation_extractor  --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/BEDGRAPHS --scaffolds --cytosine_report --CX --bedGraph --gzip /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC1.1_val_1_bismark_bt2_pe.deduplicated.bam
~/bismark_v0.16.3/./bismark_methylation_extractor  --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/BEDGRAPHS --scaffolds --cytosine_report --CX --bedGraph --gzip /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive1.1_val_1_bismark_bt2_pe.deduplicated.bam
~/bismark_v0.16.3/./bismark_methylation_extractor   --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/BEDGRAPHS --scaffolds --cytosine_report --CX --bedGraph --gzip /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Primed1.1_val_1_bismark_bt2_pe.deduplicated.bam

~/bismark_v0.16.3/./bismark_methylation_extractor  --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/BEDGRAPHS --scaffolds --cytosine_report --CX --bedGraph --gzip /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC2.1_val_1_bismark_bt2_pe.deduplicated.bam
~/bismark_v0.16.3/./bismark_methylation_extractor  --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/BEDGRAPHS --scaffolds --cytosine_report --CX --bedGraph --gzip /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive2.1_val_1_bismark_bt2_pe.deduplicated.bam
~/bismark_v0.16.3/./bismark_methylation_extractor  --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/BEDGRAPHS --scaffolds --cytosine_report --CX --bedGraph --gzip /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Primed2.1_val_1_bismark_bt2_pe.deduplicated.bam

~/bismark_v0.16.3/./bismark_methylation_extractor  --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/BEDGRAPHS --scaffolds --cytosine_report --CX --bedGraph --gzip /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC3.1_val_1_bismark_bt2_pe.deduplicated.bam
~/bismark_v0.16.3/./bismark_methylation_extractor  --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/BEDGRAPHS --scaffolds --cytosine_report --CX --bedGraph --gzip /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive3.1_val_1_bismark_bt2_pe.deduplicated.bam
~/bismark_v0.16.3/./bismark_methylation_extractor  --genome /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/ -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/BEDGRAPHS --scaffolds --cytosine_report --CX --bedGraph --gzip /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Primed3.1_val_1_bismark_bt2_pe.deduplicated.bam

#Primed3.1_val_1_bismark_bt2_pe.deduplicated.bam
#~/bismark_v0.16.3/./bismark_methylation_extractor


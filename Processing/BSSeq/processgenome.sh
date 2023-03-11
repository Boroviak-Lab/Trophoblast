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

#./bismark_genome_preparation 
~/bismark_v0.16.3/./bismark_genome_preparation /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/


# ~/bismark_v0.16.3/bismark_genome_preparation
# -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Sequence/WholeGenomeFasta/ -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC2.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC2.2_val_2.fq.gz


#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/surani/cap76/Genomes/UCSC/hg19 -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC3.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC3.2_val_2.fq.gz

#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Sequence/WholeGenomeFasta/ -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive1.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive1.2_val_2.fq.gz
#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Sequence/WholeGenomeFasta/ -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive2.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive2.2_val_2.fq.gz
#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Sequence/WholeGenomeFasta/ -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive3.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/Naive3.2_val_2.fq.gz

#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Sequence/WholeGenomeFasta/ -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC1.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC1.2_val_2.fq.gz
#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Sequence/WholeGenomeFasta/ -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC2.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC2.2_val_2.fq.gz
#./bismark -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --genome /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Sequence/WholeGenomeFasta/ -1 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC3.1_val_1.fq.gz -2 /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC3.2_val_2.fq.gz



#./deduplicate_bismark --bam /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC1.1_val_1_bismark_bt2_pe.bam
#./bismark_methylation_extractor -o /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed --comprehensive --bedGraph --gzip /mnt/scratch/gurdon/cap76/Thorsten/Dylan/BSSeq/Processed/cmTSC1.1_val_1_bismark_bt2.deduplicated.bam




#!./bismark_methylation_extractor -o /mnt/scratch/surani/cap76/BEDGRAPHS --comprehensive --bedGraph --gzip /mnt/scratch/surani/cap76/Processed/SRR17773$i*_trimmed_bismark_bt2.deduplicated.bam"



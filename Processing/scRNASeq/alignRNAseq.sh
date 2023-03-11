#!/bin/bash
#!
#! Example SLURM job script for Gurdon Institute Cluster
#! Last updated: Sat Apr 18 13:05:53 BST 2015
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! sbatch directives begin here ###############################
#! Name of the job:
#SBATCH -J PGCmeth
#! How many (MPI) tasks will there be in total? (<= nodes*16)
#SBATCH -n 6
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
#SBATCH --no-requeue
#! Job farming
#SBATCH --array=1-96%20

#769
#mkdir /mnt/scratch/surani/cap76/TEST3
#STAR --runMode genomeGenerate --limitGenomeGenerateRAM 48000000000  --genomeDir /mnt/scratch/surani/cap76/TEST3 --genomeFastaFiles cjacchus.fa --sjdbGTFfile cjacchus.gtf --sjdbOverhang 149

cd /mnt/scratch/gurdon/cap76/Thorsten/CRUK_210218_A00489_0793_AH3FMJDRXY/

#FILES=( *.bw )
#SLX-19811.i702_i503.HKJ5LBBXY.s_1.r_1.fq.gz
#SLX-19811.i702_i503.HKJ5LBBXY.s_1.r_1.fq.gz

#m=( $(ls *.cram | sed -e 's/\.cram$//' ) )
#FILES=( $(ls *.cram | sed -e 's/\.cram$//' ) )
#a=($(echo $t | tr ',' "\n"))

#ERR2851911_1.fastq.gz

#m=( $(ls *_1.fq.gz | sed -e 's/\_1.fq.gz$//' ) )
#FILES=( $(ls *_1.fq.gz | sed -e 's/\_1.fq.gz$//' ) )
#a=($(echo $t | tr ',' "\n"))



FILES=( $(ls *_1.fq.gz | sed -e 's/\_1.fq.gz$//' ) )


#ERR2851716_1.fastq.gz
mkdir fastqc_results2 
mkdir fastqc_trimmed_results2 
mkdir STAR_results2 

cp /mnt/scratch/gurdon/cap76/Thorsten/MissingTranscriptomics/wigToBigWig STAR_results2
cp /mnt/scratch/gurdon/cap76/Thorsten/MissingTranscriptomics/wigToBigWig fastqc_trimmed_results2

T=8 
wdir='/mnt/scratch/gurdon/cap76/Thorsten/CRUK_210218_A00489_0793_AH3FMJDRXY/'
FASTA='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.fa'
GTF='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/newGTF_modified_extended5kb2.gtf'
INDICES='/mnt/scratch/surani/cap76/TEST2'
CHRSIZES='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.sizes'
GENEMODEL='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.bed'

mystring="
fastqc -o ${wdir}fastqc_results2 -t $T ${FILES[$SLURM_ARRAY_TASK_ID-1]}_1.fq.gz ${FILES[$SLURM_ARRAY_TASK_ID-1]}_2.fq.gz ; 
echo 'FASTQC has finished'  >> ${wdir}o.smartseq2 ; ~/trim_galore --paired  -o ${wdir}fastqc_trimmed_results2 ${FILES[$SLURM_ARRAY_TASK_ID-1]}_1.fq.gz ${FILES[$SLURM_ARRAY_TASK_ID-1]}_2.fq.gz ; echo 'trim_galore has finished' >> ${wdir}o.smartseq2 ; cd ${wdir}fastqc_trimmed_results2 ; cat ${FILES[$SLURM_ARRAY_TASK_ID-1]}_*.fastq.gz_trimming_report.txt | grep -A2 'processed' > QC_fastq.${FILES[$SLURM_ARRAY_TASK_ID-1]}.txt ; STAR --runThreadN $T --runMode alignReads --genomeDir $INDICES --readFilesIn ${FILES[$SLURM_ARRAY_TASK_ID-1]}_1_val_1.fq.gz ${FILES[$SLURM_ARRAY_TASK_ID-1]}_2_val_2.fq.gz --readFilesCommand gunzip -c --outFileNamePrefix ${wdir}STAR_results2/${FILES[$SLURM_ARRAY_TASK_ID-1]} --sjdbGTFfile $GTF --sjdbOverhang 149 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outWigType wiggle --outWigStrand Unstranded --outWigNorm RPM ; echo 'STAR has finished' >> ${wdir}o.smartseq2 ; cd .. ; cd STAR_results2 ; cat  ${FILES[$SLURM_ARRAY_TASK_ID-1]}Log.final.out | grep 'Uniquely mapped reads' > QC_uniquelyMappedReads.${FILES[$SLURM_ARRAY_TASK_ID-1]}.txt ; samtools index ${FILES[$SLURM_ARRAY_TASK_ID-1]}Aligned.sortedByCoord.out.bam ; cd .. ; echo 'Smart-SEQ2 analysis completed...' >> ${wdir}o.smartseq2 ; echo ${FILES[$SLURM_ARRAY_TASK_ID-1]} >> ${wdir}o.smartseq2 ; echo '--------------------------------' >> ${wdir}o.smartseq2"

echo $mystring

eval $mystring

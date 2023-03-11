GTF='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/newGTF_modified_extended5kb2.gtf'
#GTF2='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.gtf'
#GTF3='/mnt/scratch/gurdon/cap76/Thorsten/newGTF_modified_extended1kb2rr.gtf'

A=$(cat 03_selectSamples_and_featurecounts1.txt)

#featureCounts 
../../CRUK/subread-2.0.1-Linux-x86_64/bin/featureCounts -p -T 16 -t exon -g gene_id -a $GTF -o featurecountsAll_extended_1B.txt $(echo $A)
#featureCounts -p -T 16 -t exon -g gene_id -a $GTF2 -o featurecountsAll_base_1.txt $(echo $A)
#featureCounts -p -T 16 -t exon -g gene_id -a $GTF3 -o featurecountsAll_extendedB_1.txt $(echo $A)



#A=$(cat bamlist_2.txt)
#
#featureCounts -p -T 16 -t exon -g gene_id -a $GTF -o featurecountsAll_extended_2.txt $(echo $A)
#featureCounts -p -T 16 -t exon -g gene_id -a $GTF2 -o featurecountsAll_base_2.txt $(echo $A)
#featureCounts -p -T 16 -t exon -g gene_id -a $GTF3 -o featurecountsAll_extendedB_2.txt $(echo $A)


#featureCounts -p -T 16 -t exon -g gene_id -a /mnt/scratch/gurdon/cap76/Thorsten/newGTF_modified_extended1kb.gtf -o featurecountsAll_extendedB_2.txt $(echo $(cat bamlist_2.txt))



#featureCounts -p -T 16 -t exon -g gene_id -a $GTF -o featurecountsAll_extended_miss.txt 3536STDY8289482Aligned.sortedByCoord.out.bam
#featureCounts -p -T 16 -t exon -g gene_id -a $GTF2 -o featurecountsAll_base_miss.txt 3536STDY8289482Aligned.sortedByCoord.out.bam

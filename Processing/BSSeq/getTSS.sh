


grep -v "^#" newGTF_modified_extended5kb2.gtf | cut -f3 | sort | uniq -c | sort -k1rn
awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' newGTF_modified_extended5kb2.gtf | sortBed | mergeBed -i - > exon_merged.bed
awk 'BEGIN{OFS="\t";} $3=="transcript" {print $1,$4-1,$5}' newGTF_modified_extended5kb2.gtf | sortBed | mergeBed -i - > transcript_merged.bed
bedtools subtract -a transcript_merged.bed -b exon_merged.bed >  introns.bed



#Intronic
awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' newGTF_modified_extended5kb2.gtf | sortBed | subtractBed -a stdin -b exon_merged.bed > intron.bed
#intersectBed -a gencode_v${v}_exon_merged.bed.gz -b gencode_v${v}_intron.bed.gz

#Intergenic
awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' newGTF_modified_extended5kb2.gtf > genes.bed
awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' genes.bed | sortBed | complementBed -i stdin -g cjacchus.sizes > intergenic.bed

bedtools subtract -a intergenic.bed -b genes.bed > intergenic.us.bed

#FASTA='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.fa'
#GTF='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/newGTF_modified_extended5kb2.gtf'


#awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' |

#complementBed -i stdin -g hg19.genome |
#gzip > gencode_v${v}_intergenic.bed.gz


#awk 'BEGIN{OFS="\t";} $3=="transcript" {print $1,$4-1,$5}' newGTF_modified_extended5kb2.gtf | sortBed | mergeBed -i - > transcript_merged.bed
#awk 'BEGIN{OFS="\t";} {print $1,$4-1,$5}' newGTF_modified_extended5kb2.gtf > sortgtf.bed 
#bedtools sort -i /mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.sizes > chr.sorted
#sortBed /mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.sizes > chr.sorted
#bedtools complement -i intergenic.bed -g sortgtf.bed > intergenic_sorted.bed
#

promoter.pl newGTF_modified_extended5kb2.gtf 200 | head




awk 'BEGIN {FS=OFS="\t"}; { print $2,$3,$4+1,$6,$1 }' TSS.bed > TSS_s.bed


bedtools slop -i TSS_s.bed  -g /mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.sizes -b 1000 > TSS_1kb.sh
bedtools slop -i TSS_s.bed  -g /mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.sizes -b 5000 > TSS_5kb.sh
#bedtools slop -i TSS.bed  -g /mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.sizes -b 1000 > TSS_1kb.sh


#awk 'BEGIN{OFS="\t";FS="\t"} {if ($7 == "+") print $1,$4-2000,$4+50,$9,".",$7; else print $1,$5-2000,$5+50,$9,".",$7}' TSS.bed > TSS.bed




awk 'BEGIN {FS=OFS="\t"}; { print $2,$3,$4+1,$6,$1 }' gene_body.bed > gene_body_s.bed


#bedtools slop -i TSS_s.bed  -g /mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.sizes -b 1000 > TSS_1kb.sh
#bedtools slop -i TSS_s.bed  -g /mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.sizes -b 5000 > TSS_5kb.sh


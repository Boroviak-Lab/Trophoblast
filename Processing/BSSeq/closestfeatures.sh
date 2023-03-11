#Trim header first - add intoo scrip
awk -v FS='\t' -v OFS='\t' '{ { print $1, $2, $3, $5 } }' AllTB_Primed_Gain.txt > AllTB_Primed_Gains.txt
awk -v FS='\t' -v OFS='\t' '{ { print $1, $2, $3, $5 } }' AllTB_Primed_Loss.txt > AllTB_Primed_Losss.txt
awk -v FS='\t' -v OFS='\t' '{ { print $1, $2, $3, $5 } }' AllTB_Niave_Gain.txt > AllTB_Niave_Gains.txt
awk -v FS='\t' -v OFS='\t' '{ { print $1, $2, $3, $5 } }' AllTB_Niave_Loss.txt > AllTB_Niave_Losss.txt
awk -v FS='\t' -v OFS='\t' '{ { print $1, $2, $3, $5 } }' AllPrimed_Niave_Gain.txt > AllPrimed_Niave_Gains.txt
awk -v FS='\t' -v OFS='\t' '{ { print $1, $2, $3, $5 } }' AllPrimed_Niave_Loss.txt > AllPrimed_Niave_Losss.txt

#AllPrimed_Niave_Loss.txt


#bedtools sort -i ../TSS_1kb.sh > Promoter_5kb.s.bed
#bedtools sort -i ../CpG.bed > CpG.s.bed
#bedtools sort -i ../intergenic.bed > intergenic.s.bed

#bedtools sort -i gene_body_s.bed



awk '{ if ($1 == 1 || $1 == 2 || $1 == 3 || $1 == 4 || $1 == 5 || $1 == 6 || $1 == 7 || $1 == 8 || $1 == 9 || $1 == 10 || $1 == 11 || $1 == 12 || $1 == 13 || $1 == 14 || $1 == 15 || $1 == 16 || $1 == 17 || $1 == 18 || $1 == 19 || $1 == 20 || $1 == 21 || $1 == 22) { print } }' ../intergenic.us.bed > intergenic.s.1-22.txt
awk '{ if ($1 == 1 || $1 == 2 || $1 == 3 || $1 == 4 || $1 == 5 || $1 == 6 || $1 == 7 || $1 == 8 || $1 == 9 || $1 == 10 || $1 == 11 || $1 == 12 || $1 == 13 || $1 == 14 || $1 == 15 || $1 == 16 || $1 == 17 || $1 == 18 || $1 == 19 || $1 == 20 || $1 == 21 || $1 == 22) { print } }' ../introns.bed > intron.s.1-22.txt
awk '{ if ($1 == 1 || $1 == 2 || $1 == 3 || $1 == 4 || $1 == 5 || $1 == 6 || $1 == 7 || $1 == 8 || $1 == 9 || $1 == 10 || $1 == 11 || $1 == 12 || $1 == 13 || $1 == 14 || $1 == 15 || $1 == 16 || $1 == 17 || $1 == 18 || $1 == 19 || $1 == 20 || $1 == 21 || $1 == 22) { print } }' ../exon_merged.bed > exon.s.1-22.txt

#../intergenic.us.bed
awk '{ if ($1 == 1 || $1 == 2 || $1 == 3 || $1 == 4 || $1 == 5 || $1 == 6 || $1 == 7 || $1 == 8 || $1 == 9 || $1 == 10 || $1 == 11 || $1 == 12 || $1 == 13 || $1 == 14 || $1 == 15 || $1 == 16 || $1 == 17 || $1 == 18 || $1 == 19 || $1 == 20 || $1 == 21 || $1 == 22) { print } }' Promoter_5kb.s.bed > Prom_1kb.s.1-22.txt
awk '{ if ($1 == 1 || $1 == 2 || $1 == 3 || $1 == 4 || $1 == 5 || $1 == 6 || $1 == 7 || $1 == 8 || $1 == 9 || $1 == 10 || $1 == 11 || $1 == 12 || $1 == 13 || $1 == 14 || $1 == 15 || $1 == 16 || $1 == 17 || $1 == 18 || $1 == 19 || $1 == 20 || $1 == 21 || $1 == 22) { print } }' ../gene_body_s.bed > GB_1kb.s.1-22.txt
awk '{ if ($1 == 1 || $1 == 2 || $1 == 3 || $1 == 4 || $1 == 5 || $1 == 6 || $1 == 7 || $1 == 8 || $1 == 9 || $1 == 10 || $1 == 11 || $1 == 12 || $1 == 13 || $1 == 14 || $1 == 15 || $1 == 16 || $1 == 17 || $1 == 18 || $1 == 19 || $1 == 20 || $1 == 21 || $1 == 22) { print } }' CpG.s.bed > CpG.s.1-22.txt
awk '{ if ($1 == 1 || $1 == 2 || $1 == 3 || $1 == 4 || $1 == 5 || $1 == 6 || $1 == 7 || $1 == 8 || $1 == 9 || $1 == 10 || $1 == 11 || $1 == 12 || $1 == 13 || $1 == 14 || $1 == 15 || $1 == 16 || $1 == 17 || $1 == 18 || $1 == 19 || $1 == 20 || $1 == 21 || $1 == 22) { print } }' intergenic.s.bed > intergenic.s.1-22.txt

bedtools sort -i AllTB_Primed_Gains.txt > AllTB_Primed_Gains.s.txt
bedtools sort -i AllTB_Primed_Losss.txt >  AllTB_Primed_Loss.s.txt
bedtools sort -i AllTB_Niave_Gains.txt > AllTB_Niave_Gains.s.txt
bedtools sort -i AllTB_Niave_Losss.txt > AllTB_Niave_Loss.s.txt

bedtools sort -i AllPrimed_Niave_Gains.txt > AllPrimed_Niave_Gains.s.txt
bedtools sort -i AllPrimed_Niave_Losss.txt > AllPrimed_Niave_Loss.s.txt



bedtools closest -a AllTB_Primed_Gains.s.txt -b Prom_1kb.s.1-22.txt -d > AllTB_Primed_Gains_anno.txt
bedtools closest -a AllTB_Primed_Loss.s.txt -b Prom_1kb.s.1-22.txt -d > AllTB_Primed_Loss_anno.txt
bedtools closest -a AllTB_Niave_Gains.s.txt -b Prom_1kb.s.1-22.txt -d  > AllTB_Niave_Gains_anno.txt
bedtools closest -a AllTB_Niave_Loss.s.txt -b Prom_1kb.s.1-22.txt -d > AllTB_Niave_Loss_anno.txt

bedtools closest -a AllPrimed_Niave_Gains.s.txt -b Prom_1kb.s.1-22.txt -d  > AllPrimed_Niave_Gains_anno.txt
bedtools closest -a AllPrimed_Niave_Loss.s.txt -b Prom_1kb.s.1-22.txt -d > AllPrimed_Niave_Loss_anno.txt


#bedtools closest -a AllTB_Primed_Gain.s.txt -b PromConv_1kb.s.1-22.txt ConvGB_1kb.s.1-22.txt CpG.s.1-22.txt intergenic.s.1-22.txt > AllTB_Primed_Gain_closest.txt


bedtools intersect -a AllTB_Primed_Gains.txt -b intergenic.snp.1-22.txt > AllTB_Primed_Gains_overlapIG.bed
bedtools intersect -a AllTB_Primed_Losss.txt -b intergenic.snp.1-22.txt > AllTB_Primed_Losss_overlapIG.bed
bedtools intersect -a AllTB_Niave_Gains.txt -b intergenic.snp.1-22.txt > AllTB_Niave_Gains_overlapIG.bed
bedtools intersect -a AllTB_Niave_Losss.txt -b intergenic.snp.1-22.txt > AllTB_Niave_Losss_overlapIG.bed

bedtools intersect -a AllTB_Primed_Gains.txt -b Prom_1kb.s.1-22.txt > AllTB_Primed_Gains_overlapPr.bed
bedtools intersect -a AllTB_Primed_Losss.txt -b Prom_1kb.s.1-22.txt > AllTB_Primed_Losss_overlapPr.bed
bedtools intersect -a AllTB_Niave_Gains.txt -b Prom_1kb.s.1-22.txt > AllTB_Niave_Gains_overlapPr.bed
bedtools intersect -a AllTB_Niave_Losss.txt -b Prom_1kb.s.1-22.txt > AllTB_Niave_Losss_overlapPr.bed

bedtools intersect -a AllTB_Primed_Gains.txt -b CpG.s.1-22.txt > AllTB_Primed_Gains_overlapCPG.bed
bedtools intersect -a AllTB_Primed_Losss.txt -b CpG.s.1-22.txt > AllTB_Primed_Losss_overlapCPG.bed
bedtools intersect -a AllTB_Niave_Gains.txt -b CpG.s.1-22.txt > AllTB_Niave_Gains_overlapCPG.bed
bedtools intersect -a AllTB_Niave_Losss.txt -b CpG.s.1-22.txt > AllTB_Niave_Losss_overlapCPG.bed

bedtools intersect -a AllTB_Primed_Gains.txt -b GB_1kb.s.1-22.txt > AllTB_Primed_Gains_overlapGB.bed
bedtools intersect -a AllTB_Primed_Losss.txt -b GB_1kb.s.1-22.txt > AllTB_Primed_Losss_overlapGB.bed
bedtools intersect -a AllTB_Niave_Gains.txt -b GB_1kb.s.1-22.txt > AllTB_Niave_Gains_overlapGB.bed
bedtools intersect -a AllTB_Niave_Losss.txt -b GB_1kb.s.1-22.txt > AllTB_Niave_Losss_overlapGB.bed


bedtools intersect -a AllPrimed_Niave_Gains.txt -b intergenic.snp.1-22.txt > AllPrimed_Niave_Gains_overlapIG.bed
bedtools intersect -a AllPrimed_Niave_Losss.txt -b intergenic.snp.1-22.txt > AllPrimed_Niave_Losss_overlapIG.bed
bedtools intersect -a AllPrimed_Niave_Gains.txt -b Prom_1kb.s.1-22.txt > AllPrimed_Niave_Gains_overlapPr.bed
bedtools intersect -a AllPrimed_Niave_Losss.txt -b Prom_1kb.s.1-22.txt > AllPrimed_Niave_Losss_overlapPr.bed
bedtools intersect -a AllPrimed_Niave_Gains.txt -b CpG.s.1-22.txt > AllPrimed_Niave_Gains_overlapCPG.bed
bedtools intersect -a AllPrimed_Niave_Losss.txt -b CpG.s.1-22.txt > AllPrimed_Niave_Losss_overlapCPG.bed
bedtools intersect -a AllPrimed_Niave_Gains.txt -b GB_1kb.s.1-22.txt > AllPrimed_Niave_Gains_overlapGB.bed
bedtools intersect -a AllPrimed_Niave_Losss.txt -b GB_1kb.s.1-22.txt > AllPrimed_Niave_Losss_overlapGB.bed

bedtools intersect -a AllPrimed_Niave_Gains.txt -b intron.s.1-22.txt > AllPrimed_Niave_Gains_overlapINT.bed
bedtools intersect -a AllPrimed_Niave_Losss.txt -b intron.s.1-22.txt > AllPrimed_Niave_Losss_overlapINT.bed
bedtools intersect -a AllPrimed_Niave_Gains.txt -b exon.s.1-22.txt > AllPrimed_Niave_Gains_overlapEX.bed
bedtools intersect -a AllPrimed_Niave_Losss.txt -b exon.s.1-22.txt > AllPrimed_Niave_Losss_overlapEX.bed

bedtools intersect -a AllTB_Primed_Gains.txt -b intron.s.1-22.txt > AllTB_Primed_Gains_overlapINT.bed
bedtools intersect -a AllTB_Primed_Losss.txt -b intron.s.1-22.txt > AllTB_Primed_Losss_overlapINT.bed
bedtools intersect -a AllTB_Primed_Gains.txt -b exon.s.1-22.txt > AllTB_Primed_Gains_overlapEX.bed
bedtools intersect -a AllTB_Primed_Losss.txt -b exon.s.1-22.txt > AllTB_Primed_Losss_overlapEX.bed

bedtools intersect -a AllTB_Niave_Gains.txt -b intron.s.1-22.txt > AllTB_Niave_Gains_overlapINT.bed
bedtools intersect -a AllTB_Niave_Losss.txt -b intron.s.1-22.txt > AllTB_Niave_Losss_overlapINT.bed
bedtools intersect -a AllTB_Niave_Gains.txt -b exon.s.1-22.txt > AllTB_Niave_Gains_overlapEX.bed
bedtools intersect -a AllTB_Niave_Losss.txt -b exon.s.1-22.txt > AllTB_Niave_Losss_overlapEX.bed


awk '{ if ($1 == 1 || $1 == 2 || $1 == 3 || $1 == 4 || $1 == 5 || $1 == 6 || $1 == 7 || $1 == 8 || $1 == 9 || $1 == 10 || $1 == 11 || $1 == 12 || $1 == 13 || $1 == 14 || $1 == 15 || $1 == 16 || $1 == 17 || $1 == 18 || $1 == 19 || $1 == 20 || $1 == 21 || $1 == 22) { print } }' ../intergenic.us.bed > intergenic.s.1-22.txt
awk '{ if ($1 == 1 || $1 == 2 || $1 == 3 || $1 == 4 || $1 == 5 || $1 == 6 || $1 == 7 || $1 == 8 || $1 == 9 || $1 == 10 || $1 == 11 || $1 == 12 || $1 == 13 || $1 == 14 || $1 == 15 || $1 == 16 || $1 == 17 || $1 == 18 || $1 == 19 || $1 == 20 || $1 == 21 || $1 == 22) { print } }' ../introns.bed > intron.s.1-22.txt
awk '{ if ($1 == 1 || $1 == 2 || $1 == 3 || $1 == 4 || $1 == 5 || $1 == 6 || $1 == 7 || $1 == 8 || $1 == 9 || $1 == 10 || $1 == 11 || $1 == 12 || $1 == 13 || $1 == 14 || $1 == 15 || $1 == 16 || $1 == 17 || $1 == 18 || $1 == 19 || $1 == 20 || $1 == 21 || $1 == 22) { print } }' ../exon_merged.bed > exon.s.1-22.txt


#GB_1kb.s.1-22.txt

#bedtools intersect -a AllTB_Primed_Gain.s.txt -b PromConv_1kb.s.1-22.txt > AllTB_Primed_Gain_PromConv.bed
#bedtools intersect -a AllTB_Primed_Gain.s.txt -b ConvGB_1kb.s.1-22.txt  > AllTB_Primed_Gain_ConvGB.bed
#bedtools intersect -a AllTB_Primed_Gain.s.txt -b CpG.s.1-22.txt > AllTB_Primed_Gain_CpG.bed
#bedtools intersect -a AllTB_Primed_Gain.s.txt -b intergenic.s.1-22.txt > AllTB_Primed_Gain_intergenic.bed

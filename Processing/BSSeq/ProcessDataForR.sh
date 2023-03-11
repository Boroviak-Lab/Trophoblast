#Script to generate the profile plots for DE regions

##Put into a bed format
#awk -v FS='\t' -v OFS='\t' '{ { print $1, $2, $3, $5 } }' AllTB_Primed_Gain.txt > AllTB_Primed_Gains.bed
#awk -v FS='\t' -v OFS='\t' '{ { print $1, $2, $3, $5 } }' AllTB_Primed_Loss.txt > AllTB_Primed_Losss.bed
#awk -v FS='\t' -v OFS='\t' '{ { print $1, $2, $3, $5 } }' AllTB_Niave_Gain.txt > AllTB_Niave_Gains.bed
#awk -v FS='\t' -v OFS='\t' '{ { print $1, $2, $3, $5 } }' AllTB_Niave_Loss.txt > AllTB_Niave_Losss.bed
#awk -v FS='\t' -v OFS='\t' '{ { print $1, $2, $3, $5 } }' AllPrimed_Niave_Gain.txt > AllPrimed_Niave_Gains.bed
#awk -v FS='\t' -v OFS='\t' '{ { print $1, $2, $3, $5 } }' AllPrimed_Niave_Loss.txt > AllPrimed_Niave_Losss.bed

##Merge, sort, merge peaks to create a combined file
#cat AllTB_Primed_Gains.bed AllTB_Primed_Losss.bed AllTB_Niave_Gains.bed AllTB_Niave_Losss.bed AllPrimed_Niave_Gains.bed AllPrimed_Niave_Losss.bed > allmerge.bed
#sort -k1,1 -k2,2n allmerge.bed > allmerge.s.bed
#bedtools merge -d 1000 -i allmerge.s.bed > allMerge.bed

##Anotate closes promoter
#bedtools closest -a allMerge.bed -b Prom_1kb.s.1-22.txt -d > allMerge_anno.bed

##Get signal
computeMatrix scale-regions -S Primed1.bw Primed2.bw Primed3.bw Niave1.bw Niave2.bw Niave3.bw TSC1.bw TSC2.bw TSC3.bw -R allMerge.bed --binSize 1 --beforeRegionStartLength 0 --regionBodyLength 1 --afterRegionStartLength 0 --skipZeros -o AllDEInVit.mat.gz
##computeMatrix scale-regions -S Primed1.bw Primed2.bw Primed3.bw Niave1.bw Niave2.bw Niave3.bw TSC1.bw TSC2.bw TSC3.bw -R Prom_1kb.s.1-22.txt --binSize 1 --beforeRegionStartLength 0 --regionBodyLength 1 --afterRegionStartLength 0 --skipZeros -o AllInVit.mat.gz
##computeMatrix scale-regions -S EmDisc.bw Am.bw VE.bw ExMes.bw Tb.bw Gland.bw -R allMerge.bed --binSize 1 --beforeRegionStartLength 0 --regionBodyLength 1 --afterRegionStartLength 0 --skipZeros -o AllDEInViv.mat.gz


CHRSIZES='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.sizes'

#sort -k1,1 -k2,2n -k3,3n -k6,6 -u AllTB_Niave_Gains_anno.txt > AllTB_Niave_Gains_anno.u.bed
#sort -k1,1 -k2,2n -k3,3n -k6,6 -u AllTB_Primed_Loss.txt > AllTB_Primed_Loss.u.bed
#sort -k1,1 -k2,2n -k3,3n -k6,6 -u AllTB_Niave_Gains_anno.txt > AllTB_Niave_Gains_anno.u.bed
#sort -k1,1 -k2,2n -k3,3n -k6,6 -u AllTB_Niave_Loss_anno.txt > AllTB_Niave_Loss_anno.u.bed
#sort -k1,1 -k2,2n -k3,3n -k6,6 -u AllPrimed_Niave_Gains_anno.txt > AllPrimed_Niave_Gains_anno.u.bed
#sort -k1,1 -k2,2n -k3,3n -k6,6 -u  AllPrimed_Niave_Loss_anno.txt > AllPrimed_Niave_Loss_anno.u.bed

computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R AllTB_Primed_Gain.u.bed --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o tbpmatrixConvG.mat.gz
plotProfile -m tbpmatrixConvG.mat.gz -out PT_lineG1.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup

computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R AllTB_Primed_Loss.u.bed --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o tbpmatrixConvL.mat.gz
plotProfile -m tbpmatrixConvL.mat.gz -out PT_lineL1.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup

computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R AllTB_Niave_Gains_anno.u.bed --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o tbpmatrixConvG_2.mat.gz
plotProfile -m tbpmatrixConvG_2.mat.gz -out PT_lineG1_2.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup

computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R AllTB_Niave_Loss_anno.u.bed --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o tbpmatrixConvL_2.mat.gz
plotProfile -m tbpmatrixConvL_2.mat.gz -out PT_lineL1_2.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup

computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R AllPrimed_Niave_Gains_anno.u.bed --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o tbpmatrixConvG_3.mat.gz
plotProfile -m tbpmatrixConvG_3.mat.gz -out PT_lineG1_3.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup

computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R AllPrimed_Niave_Loss_anno.u.bed --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o tbpmatrixConvL_3.mat.gz
plotProfile -m tbpmatrixConvL_3.mat.gz -out PT_lineL1_3.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup



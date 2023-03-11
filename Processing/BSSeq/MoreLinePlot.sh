CHRSIZES='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.sizes'


#computeMatrix scale-regions -S Primed1.bw Primed2.bw Primed3.bw Niave1.bw Niave2.bw Niave3.bw TSC1.bw TSC2.bw TSC3.bw -R ../TSS_5kb.sh --binSize 1 --beforeRegionStartLength 0 --regionBodyLength 1 --afterRegionStartLength 0 --skipZeros -o AllInVit.mat.gz
#computeMatrix scale-regions -S EmDisc.bw Am.bw VE.bw ExMes.bw Tb.bw Gland.bw -R ../TSS_5kb.sh --binSize 1 --beforeRegionStartLength 0 --regionBodyLength 1 --afterRegionStartLength 0 --skipZeros -o AllInViv.mat.gz

#bedtools merge -d 1000 -i allPeaks.bed > allMerge.bed


#computeMatrix scale-regions -S Primed1.bw Primed2.bw Primed3.bw Niave1.bw Niave2.bw Niave3.bw TSC1.bw TSC2.bw TSC3.bw -R allMerge.bed --binSize 1 --beforeRegionStartLength 0 --regionBodyLength 1 --afterRegionStartLength 0 --skipZeros -o AllDEInVit.mat.gz
#computeMatrix scale-regions -S EmDisc.bw Am.bw VE.bw ExMes.bw Tb.bw Gland.bw -R allMerge.bed --binSize 1 --beforeRegionStartLength 0 --regionBodyLength 1 --afterRegionStartLength 0 --skipZeros -o AllDEInViv.mat.gz


#computeMatrix scale-regions -S Primed1.bw Primed2.bw Primed3.bw Niave1.bw Niave2.bw Niave3.bw TSC1.bw TSC2.bw TSC3.bw EmDisc.bw Tb.bw cMor.bw -R ../TSS_5kb.sh --binSize 1 --beforeRegionStartLength 0 --regionBodyLength 1 --afterRegionStartLength 0 --skipZeros -o AllIn.mat.gz


computeMatrix scale-regions -S Primed1.bw Primed2.bw Primed3.bw Niave1.bw Niave2.bw Niave3.bw TSC1.bw TSC2.bw TSC3.bw EmDisc.bw Tb.bw cMor.bw -R allMerge.bed --binSize 100 --beforeRegionStartLength 1000 --regionBodyLength 100 --afterRegionStartLength 1000 --skipZeros -o AllDEInExt.mat.gz



#m.bw    EightCell.bw    ElevenCell.bw  ExMes.bw     Gland.bw     Niave1.bw  Niave3.bw   Primed2.bw  RemodCS7.bw   SYS.bw     Tb.bw             TrueEightCellIV.bw  TSC1.bw  TSC3.bw
#cMor.bw  EightCellIV.bw  EmDisc.bw      FourCell.bw  GlandCS7.bw  Niave2.bw  Primed1.bw  Primed3.bw  StromaCS7.bw  SYSCS7.bw  TrueEightCell.bw  TrueTenCellIV.bw    TSC2.bw  VE.bw
#c

#bedtools intersect -a AllPrimed_Niave_Gains.txt -b intergenic.snp.1-22.txt > AllPrimed_Niave_Gains_overlapIG.bed
#bedtools intersect -a AllPrimed_Niave_Losss.txt -b intergenic.snp.1-22.txt > AllPrimed_Niave_Losss_overlapIG.bed

#AllPrimed_Niave_Gains_overlapIG.bed
#AllPrimed_Niave_Gains_overlapIG.bed

#bedtools annotate -i /mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.gtf -files AllDEIn.mat > AllDEInAno.mat
#bedtools annotate -i /mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.gtf -files AllIn.mat > AllInAno.mat

#bedtools sort -i ../TSS_1kb.sh > TSS_1kb.s.sh

#bedtools sort -i AllDEIn.mat > AllDEInAno.s.mat
#bedtools sort -i  AllIn.mat >  AllIn.s.mat

#bedtools closest -a AllDEInAno.s.mat -b TSS_1kb.s.sh > AllDEInAno.mat
#bedtools closest -a AllIn.s.mat -b TSS_1kb.s.sh > AllInAno.mat

#<FILE> \
#                           -b <FILE1, FILE2, ..., FILEN>


#AllDEIn.mat
#AllIn.mat

#bedtools closest -a AllTB_Primed_Gains.s.txt -b Prom_1kb.s.1-22.txt -d > AllTB_Primed_Gains_anno.txt
#bedtools closest -a AllTB_Primed_Loss.s.txt -b Prom_1kb.s.1-22.txt -d > AllTB_Primed_Loss_anno.txt
#bedtools closest -a AllTB_Niave_Gains.s.txt -b Prom_1kb.s.1-22.txt -d  > AllTB_Niave_Gains_anno.txt
#bedtools closest -a AllTB_Niave_Loss.s.txt -b Prom_1kb.s.1-22.txt -d > AllTB_Niave_Loss_anno.txt

#bedtools closest -a AllPrimed_Niave_Gains.s.txt -b Prom_1kb.s.1-22.txt -d  > AllPrimed_Niave_Gains_anno.txt
#bedtools closest -a AllPrimed_Niave_Loss.s.txt -b Prom_1kb.s.1-22.txt -d > AllPrimed_Niave_Loss_anno.txt

#cat AllTB_Primed_Gains.s.txt AllTB_Primed_Loss.s.txt AllTB_Niave_Gains.s.txt AllTB_Niave_Loss.s.txt AllPrimed_Niave_Gains.s.txt AllPrimed_Niave_Loss.s.txt > allPeaks.bed

#bedtools merge -i allPeaks.bed > allMerge.bed
#cat All8C_Mor_Gain.txt All8C_Mor_Loss.txt AllEmD_Tb_Gain.txt AllEmD_Tb_Loss.txt AllcMor_EmDisc_Gain.txt AllcMor_EmDisc_Loss.txt AllIVF_Gain.txt AllIVF_Loss.txt > allPeaksInVivo.bed
#bedtools merge -d 1000 -i allPeaksInVivo.bed > allPeaksInVivo.merge.bed
#computeMatrix scale-regions -S Primed1.bw Primed2.bw Primed3.bw Niave1.bw Niave2.bw Niave3.bw TSC1.bw TSC2.bw TSC3.bw EmDisc.bw Tb.bw cMor.bw -R allPeaksInVivo.bed --binSize 1 --beforeRegionStartLength 0 --regionBodyLength 1 --afterRegionStartLength 0 --skipZeros -o AllDEInVivoPeaks.mat


#<BED/GFF/VCF/BAM>
#bedtools closest -a AllIn.s.mat -b TSS_1kb.s.sh > AllInAno.mat


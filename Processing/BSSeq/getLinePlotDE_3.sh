CHRSIZES='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.sizes'

#sort -k1,1 -k2,2n -k3,3n -k6,6 -u AllTB_Primed_Gain.txt > AllTB_Primed_Gain.u.bed
#sort -k1,1 -k2,2n -k3,3n -k6,6 -u AllTB_Primed_Loss.txt > AllTB_Primed_Loss.u.bed

#sort -k1,1 -k2,2n -k3,3n -k6,6 -u AllTB_Niave_Gains_anno.txt > AllTB_Niave_Gains_anno.u.bed
#sort -k1,1 -k2,2n -k3,3n -k6,6 -u AllTB_Niave_Loss_anno.txt > AllTB_Niave_Loss_anno.u.bed
#sort -k1,1 -k2,2n -k3,3n -k6,6 -u AllPrimed_Niave_Gains_anno.txt > AllPrimed_Niave_Gains_anno.u.bed
#sort -k1,1 -k2,2n -k3,3n -k6,6 -u  AllPrimed_Niave_Loss_anno.txt > AllPrimed_Niave_Loss_anno.u.bed


computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R AllTB_Primed_Gain.u.bed --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o tbpmatrixConvG.mat.gz
plotProfile -m tbpmatrixConvG.mat.gz -out PT_lineG1.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup

computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R AllTB_Primed_Loss.u.bed --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o tbpmatrixConvL.mat.gz
plotProfile -m tbpmatrixConvL.mat.gz -out PT_lineL1.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup

#computeMatrix scale-regions -S EmDisc.bw Tb.bw Gland.bw -R AllTB_Primed_Gain.txt --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o matrixConviv1.mat.gz
#plotProfile -m matrixConviv1.mat.gz -out PT_lineIVG1.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup

#computeMatrix scale-regions -S EmDisc.bw Tb.bw Gland.bw -R AllTB_Primed_Loss.txt --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o matrixConviv.mat.gz
#plotProfile -m matrixConviv.mat.gz -out PT_lineIVL1.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup

#AllTB_Niave_Gain_anno.txt
#AllTB_Niave_Loss_anno.txt

#AllPrimed_Niave_Gain_anno.txt
#AllPrimed_Niave_Loss_anno.txt

computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R AllTB_Niave_Gains_anno.u.bed --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o tbpmatrixConvG_2.mat.gz
plotProfile -m tbpmatrixConvG_2.mat.gz -out PT_lineG1_2.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup

computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R AllTB_Niave_Loss_anno.u.bed --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o tbpmatrixConvL_2.mat.gz
plotProfile -m tbpmatrixConvL_2.mat.gz -out PT_lineL1_2.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup

#computeMatrix scale-regions -S EmDisc.bw Tb.bw Gland.bw -R Output/AllTB_Niave_Gains_anno.txt --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o matrixConviv1_2.mat.gz
#plotProfile -m matrixConviv1_2.mat.gz -out PT_lineIVG1_2.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup

#computeMatrix scale-regions -S EmDisc.bw Tb.bw Gland.bw -R Output/AllTB_Niave_Loss_anno.txt --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o matrixConviv_2.mat.gz
#plotProfile -m matrixConviv_2.mat.gz -out PT_lineIVL1_2.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup

#AllPrimed_Niave_Gain_anno.txt
#AllPrimed_Niave_Loss_anno.txt

computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R AllPrimed_Niave_Gains_anno.u.bed --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o tbpmatrixConvG_3.mat.gz
plotProfile -m tbpmatrixConvG_3.mat.gz -out PT_lineG1_3.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup

computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R AllPrimed_Niave_Loss_anno.u.bed --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o tbpmatrixConvL_3.mat.gz
plotProfile -m tbpmatrixConvL_3.mat.gz -out PT_lineL1_3.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup

#computeMatrix scale-regions -S EmDisc.bw Tb.bw Gland.bw -R Output/AllPrimed_Niave_Gains_anno.txt --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o matrixConviv1_3.mat.gz
#plotProfile -m matrixConviv1_3.mat.gz -out PT_lineIVG1_3.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup

#computeMatrix scale-regions -S EmDisc.bw Tb.bw Gland.bw -R Output/AllPrimed_Niave_Loss_anno.txt --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o matrixConviv_3.mat.gz
#plotProfile -m matrixConviv_3.mat.gz -out PT_lineIVL1_3.pdf --numPlotsPerRow 1 --startLabel "DMR-1kbp" --endLabel "DMR + 1kbp" --perGroup



#bedtools annotate -i .bed -files FILE1 FILE2 FILE3 ... FILEn


#computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R Output/AllTB_Niave_Gain_anno.txt --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o tbpmatrixConvG_2.mat.gz


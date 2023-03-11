CHRSIZES='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.sizes'
#computeMatrix scale-regions -S EmDisc.bw Am.bw VE.bw ExMes.bw Tb.bw Gland.bw -R ../CpG.bed --beforeRegionStartLength 1000 --regionBodyLength 5000 --afterRegionStartLength 1000 --skipZeros -o InVivo3.mat.gz
EmDisc.bw
computeMatrix scale-regions -S EmDisc.bw Tb.bw Gland.bw -R AllTB_Primed_Gain.txt --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o matrixConviv1.mat.gz
plotProfile -m matrixConviv1.mat.gz -out PT_lineIVG1.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

computeMatrix scale-regions -S EmDisc.bw Tb.bw Gland.bw -R AllTB_Primed_Gain.txt --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o matrixConviv.mat.gz
plotProfile -m matrixConviv.mat.gz -out PT_lineIVL1.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup



#bedtools annotate -i .bed -files FILE1 FILE2 FILE3 ... FILEn

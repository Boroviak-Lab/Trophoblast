CHRSIZES='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.sizes'

#computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R Prom_1kb.s.1-22.txt --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o Prom.mat.gz
#plotProfile -m Prom.mat.gz -out Prom_line.pdf --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R CpG.s.1-22.txt --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o CPG.mat.gz
#plotProfile -m CPG.mat.gz -out CPG_line.pdf --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup
#plotProfile -m IEX.mat.gz -out EX_line.pdf --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup


computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R intergenic.snp.1-22.txt --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o IG.mat.gz
plotProfile -m IG.mat.gz -out IG_line.pdf --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R GB_1kb.s.1-22.txt --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o  GB.mat.gz
#plotProfile -m GB.mat.gz -out GB_line.pdf --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup


#computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R intron.s.1-22.txt --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o INT.mat.gz
#plotProfile -m INT.mat.gz -out INT_line.pdf --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup


#computeMatrix scale-regions -S Niave1.bw Niave2.bw Niave3.bw Primed1.bw Primed2.bw Primed3.bw TSC1.bw TSC2.bw TSC3.bw -R exon.s.1-22.txt --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros -o EX.mat.gz
#plotProfile -m EX.mat.gz -out IEX_line.pdf --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup



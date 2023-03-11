CHRSIZES='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.sizes'
#cp /mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/wigToBigWig ./

#./bedGraphToBigWig Primed2.1_val_1_bismark_bt2_pe.deduplicated.bedGraph $CHRSIZES  Primed2.bw
#./bedGraphToBigWig Naive2.1_val_1_bismark_bt2_pe.deduplicated.bedGraph $CHRSIZES  Niave2.bw
#./bedGraphToBigWig cmTSC2.1_val_1_bismark_bt2_pe.deduplicated.bedGraph $CHRSIZES  TSC2.bw
#./bedGraphToBigWig Primed1.1_val_1_bismark_bt2_pe.deduplicated.bedGraph $CHRSIZES  Primed1.bw
#./bedGraphToBigWig Naive1.1_val_1_bismark_bt2_pe.deduplicated.bedGraph $CHRSIZES  Niave1.bw
#./bedGraphToBigWig cmTSC1.1_val_1_bismark_bt2_pe.deduplicated.bedGraph $CHRSIZES  TSC1.bw

#./bedGraphToBigWig Primed3.1_val_1_bismark_bt2_pe.deduplicated.bedGraph $CHRSIZES  Primed3.bw
#./bedGraphToBigWig Naive3.1_val_1_bismark_bt2_pe.deduplicated.bedGraph $CHRSIZES  Niave3.bw
#./bedGraphToBigWig cmTSC3.1_val_1_bismark_bt2_pe.deduplicated.bedGraph $CHRSIZES  TSC3.bw

#./bedGraphToBigWig ../../../MethylationAll/BAM/Am_CS5.sorted.bedGraph $CHRSIZES Am.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/ExMes_CS5.sorted.bedGraph $CHRSIZES ExMes.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/SYS_CS5.sorted.bedGraph $CHRSIZES SYS.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/VE_CS5.sorted.bedGraph $CHRSIZES VE.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/EmDisc_CS5.sorted.bedGraph $CHRSIZES EmDisc.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/Gland_CS5.sorted.bedGraph $CHRSIZES Gland.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/Tb_CS5.sorted.bedGraph $CHRSIZES Tb.bw



#./bedGraphToBigWig ../../../MethylationAll/BAM/Am_CS5.sorted.bedGraph $CHRSIZES 
#./bedGraphToBigWig ../../../MethylationAll/BAM/FourCell.sorted.bedGraph $CHRSIZES FourCell.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/Tb_CS5.sorted.bedGraph $CHRSIZES
#./bedGraphToBigWig ../../../MethylationAll/BAM/cMor.sorted.bedGraph $CHRSIZES cMor.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/Gland_CS5.sorted.bedGraph $CHRSIZES
#./bedGraphToBigWig ../../../MethylationAll/BAM/TrueEightCellIV.sorted.bedGraph $CHRSIZES TrueEightCellIV.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/EightCellIV.sorted.bedGraph $CHRSIZES EightCellIV.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/GlandCS7.sorted.bedGraph $CHRSIZES GlandCS7.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/TrueEightCell.sorted.bedGraph $CHRSIZES TrueEightCell.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/EightCell.sorted.bedGraph $CHRSIZES EightCell.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/RemodCS7.sorted.bedGraph $CHRSIZES RemodCS7.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/TrueTenCellIV.sorted.bedGraph $CHRSIZES TrueTenCellIV.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/ElevenCell.sorted.bedGraph $CHRSIZES ElevenCell.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/StromaCS7.sorted.bedGraph $CHRSIZES StromaCS7.bw
#./bedGraphToBigWig ../../../MethylationAll/BAM/VE_CS5.sorted.bedGraph $CHRSIZES
#./bedGraphToBigWig ../../../MethylationAll/BAM/EmDisc_CS5.sorted.bedGraph $CHRSIZES
#./bedGraphToBigWig ../../../MethylationAll/BAM/SYS_CS5.sorted.bedGraph $CHRSIZES
#./bedGraphToBigWig ../../../MethylationAll/BAM/ExMes_CS5.sorted.bedGraph $CHRSIZES
#./bedGraphToBigWig ../../../MethylationAll/BAM/SYSCS7.sorted.bedGraph $CHRSIZES SYSCS7.bw


#multiBigwigSummary bins -b FourCell.bw TrueEightCellIV.bw TrueEightCell.bw TrueTenCellIV.bw ElevenCell.bw cMor.bw -out scores_per_binv3.npz --outRawCounts scores_per_binv3.tab
#plotCorrelation -in scores_per_binv3.npz --corMethod spearman --skipZeros  --plotTitle "Spearman Correlation of Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o heatmap_SpearmanCorr_readCountsv3.png --outFileCorMatrix SpearmanCorr_readCountsv3.tab
#plotPCA -in scores_per_binv3.npz -o PCA_readCounts3.png -T "PCA of read counts"



multiBigwigSummary bins -b FourCell.bw TrueEightCell.bw cMor.bw  EmDisc.bw Am.bw VE.bw ExMes.bw Tb.bw Gland.bw  GlandCS7.bw RemodCS7.bw StromaCS7.bw -out scores_per_binv4.npz --outRawCounts scores_per_binv4.tab
plotCorrelation -in scores_per_binv4.npz --corMethod spearman --skipZeros  --plotTitle "Spearman Correlation of Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o heatmap_SpearmanCorr_readCountsv4.png --outFileCorMatrix SpearmanCorr_readCountsv4.tab
plotPCA -in scores_per_binv4.npz -o PCA_readCounts4.png -T "PCA of read counts"



#multiBigwigSummary bins -b EmDisc.bw Am.bw VE.bw ExMes.bw Tb.bw Gland.bw FourCell.bw cMor.bw TrueEightCellIV.bw TrueEightCell.bw TrueTenCellIV.bw ElevenCell.bw GlandCS7.bw RemodCS7.bw StromaCS7.bw SYSCS7.bw -out scores_per_binv2.npz --outRawCounts scores_per_binv2.tab
#plotCorrelation -in scores_per_binv2.npz --corMethod spearman --skipZeros  --plotTitle "Spearman Correlation of Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o heatmap_SpearmanCorr_readCountsv2.png --outFileCorMatrix SpearmanCorr_readCountsv2.tab



#plotPCA -in scores_per_binv2.npz -o PCA_readCounts.png -T "PCA of read counts"

#multiBigwigSummary bins -b Primed1.bw Primed2.bw Primed3.bw Niave1.bw Niave2.bw Niave3.bw TSC1.bw TSC2.bw TSC3.bw --labels P1 P2 P3 N1 N2 N3 T1 T2 T3 -out scores_per_bin2.npz --outRawCounts scores_per_bin2.tab
#plotCorrelation -in scores_per_bin2.npz --corMethod spearman --skipZeros  --plotTitle "Spearman Correlation of Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o heatmap_SpearmanCorr_readCounts2.png --outFileCorMatrix SpearmanCorr_readCounts2.tab

#plotPCA -in scores_per_bin2.npz -o PCA_readCounts1.png -T "PCA of read counts"

#computeMatrix scale-regions -S EmDisc.bw Am.bw VE.bw ExMes.bw Tb.bw Gland.bw -R PromConv_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o InVivo1.mat.gz
#plotProfile -m InVivo1.mat.gz -out InVivo1_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S EmDisc.bw Am.bw VE.bw ExMes.bw Tb.bw Gland.bw -R ConvGB_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o InVivo2.mat.gz
#plotProfile -m InVivo2.mat.gz -out InVivo2_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S EmDisc.bw Am.bw VE.bw ExMes.bw Tb.bw Gland.bw -R ../CpG.bed --beforeRegionStartLength 1000 --regionBodyLength 5000 --afterRegionStartLength 1000 --skipZeros -o InVivo3.mat.gz
#plotProfile -m InVivo3.mat.gz -out InVivo3_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S EmDisc.bw Am.bw VE.bw ExMes.bw Tb.bw Gland.bw -R intergenic_line.png --beforeRegionStartLength 1000 --regionBodyLength 5000 --afterRegionStartLength 1000 --skipZeros -o InVivo4.mat.gz
#plotProfile -m InVivo4.mat.gz -out InVivo4_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#../CpG.bed  ../exon_merged.bed  ../gene_body.bed  ../intergenic.bed  ../intron.bed  ../Promoter_200.bed  ../Promoter_500.bed  ../Promoter_5kb.bed  ../TSS.bed
#Niave1.bw  Niave2.bw  Primed1.bw  Primed2.bw  TSC1.bw  TSC2.bw
#--startLabel "TSS-1kbp" --endLabel "TSS + 1kbp"

#computeMatrix scale-regions -S Niave1.bw Primed1.bw TSC1.bw -R PromConv_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o matrixConv.mat.gz
#plotProfile -m matrixConv.mat.gz -out Conv_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S Niave1.bw Primed1.bw TSC1.bw -R NiaveConv_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o matrixNiave.mat.gz
#plotProfile -m matrixNiave.mat.gz -out Niave_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S Niave1.bw Primed1.bw TSC1.bw -R TSCConv_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o matrixTSC.mat.gz
#plotProfile -m matrixTSC.mat.gz -out TSC_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup



#computeMatrix scale-regions -S Niave1.bw Primed1.bw TSC1.bw -R ConvGB_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o matrixConvGB.mat.gz
#plotProfile -m matrixConvGB.mat.gz -out ConvGB_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S Niave1.bw Primed1.bw TSC1.bw -R TSCGB_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o matrixTSCGB.mat.gz
#plotProfile -m matrixTSCGB.mat.gz -out TSCGB_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S Niave1.bw Primed1.bw TSC1.bw -R NiaveGB_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o matrixNiaveGB.mat.gz
#plotProfile -m matrixNiaveGB.mat.gz -out NiaveGB_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S Primed1.bw -R PromConv_1kb.txt NiaveConv_1kb.txt TSCConv_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o matrixConv1.mat.gz
#plotProfile -m matrixConv1.mat.gz -out Conv_line2.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" 

#computeMatrix scale-regions -S Niave1.bw -R PromConv_1kb.txt NiaveConv_1kb.txt TSCConv_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o matrixNiave1.mat.gz
#plotProfile -m matrixNiave1.mat.gz -out Niave_line2.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" 

#computeMatrix scale-regions -S TSC1.bw -R PromConv_1kb.txt NiaveConv_1kb.txt TSCConv_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o matrixTSC1.mat.gz
#plotProfile -m matrixTSC1.mat.gz -out TSC_line2.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" 


#computeMatrix scale-regions -S Primed1.bw -R ConvGB_1kb.txt NiaveGB_1kb.txt TSCGB_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o matrixConvGB2.mat.gz
#plotProfile -m matrixConvGB2.mat.gz -out ConvGB_line2.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" 

#computeMatrix scale-regions -S TSC1.bw -R ConvGB_1kb.txt NiaveGB_1kb.txt TSCGB_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o matrixTSCGB2.mat.gz
#plotProfile -m matrixTSCGB2.mat.gz -out TSCGB_line2.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" 

#computeMatrix scale-regions -S Niave1.bw -R ConvGB_1kb.txt NiaveGB_1kb.txt TSCGB_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o matrixNiaveGB2.mat.gz
#plotProfile -m matrixNiaveGB2.mat.gz -out NiaveGB_line2.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" 


#computeMatrix scale-regions -S Niave1.bw Primed1.bw TSC1.bw -R ../CpG.bed --beforeRegionStartLength 1000 --regionBodyLength 5000 --afterRegionStartLength 1000 --skipZeros -o matrixGpG.mat.gz
#plotProfile -m matrixGpG.mat.gz -out GpG_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S Niave1.bw Primed1.bw TSC1.bw -R ../intergenic.bed --beforeRegionStartLength 1000 --regionBodyLength 5000 --afterRegionStartLength 1000 --skipZeros -o matrixintergenic.mat.gz
#plotProfile -m matrixintergenic.mat.gz -out intergenic_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup


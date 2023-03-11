CHRSIZES='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.sizes'
#cp /mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/wigToBigWig ./


computeMatrix scale-regions -S cMor.bw EmDisc.bw Tb.bw -R ../TSS_5kb.sh --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o InVivo1.mat.gz
plotProfile -m InVivo1.mat.gz -out InViv_Prom_linefin.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup
computeMatrix scale-regions -S cMor.bw EmDisc.bw Tb.bw -R ../intergenic.us.bed --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o InVivo2.mat.gz
plotProfile -m InVivo2.mat.gz -out InViv_Iter_linefin.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup
computeMatrix scale-regions -S cMor.bw EmDisc.bw Tb.bw -R ../introns.bed --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o InVivo3.mat.gz
plotProfile -m InVivo3.mat.gz -out InViv_Inton_linefin.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup
computeMatrix scale-regions -S cMor.bw EmDisc.bw Tb.bw -R ../gene_body_s.bed --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o InVivo4.mat.gz
plotProfile -m InVivo4.mat.gz -out InViv_GB_linefin.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup
computeMatrix scale-regions -S cMor.bw EmDisc.bw Tb.bw -R ../exon_merged.bed --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o InVivo5.mat.gz
plotProfile -m InVivo5.mat.gz -out InViv_Exon_linefin.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup
computeMatrix scale-regions -S cMor.bw EmDisc.bw Tb.bw -R ../CpG.bed --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o InVivo6.mat.gz
plotProfile -m InVivo6.mat.gz -out InViv_CpG_linefin.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup



#../TSS_5kb.sh
#../intergenic.us.bed
#../introns.bed
#../gene_body_s.bed
#../exon_merged.bed
#../CpG.bed

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



#FourCell.sorted.bam TrueTenCellIV.sorted.bam TrueEightCellIV.sorted.bam ElevenCell.sorted.bam EightCellIV.sorted.bam cMor.sorted.bam SYSCS7.sorted.bam StromaCS7.sorted.bam RemodCS7.sorted.bam GlandCS7.sorted.bam

#computeMatrix scale-regions -S ../../../MethylationAll/BAM/FourCell.sorted.bam ../../../MethylationAll/BAM/TrueTenCellIV.sorted.bam ../../../MethylationAll/BAM/TrueEightCellIV.sorted.bam ../../../MethylationAll/BAM/ElevenCell.sorted.bam ../../../MethylationAll/BAM/EightCellIV.sorted.bam ../../../MethylationAll/BAM/cMor.sorted.bam ../../../MethylationAll/BAM/SYSCS7.sorted.bam ../../../MethylationAll/BAM/StromaCS7.sorted.bam ../../../MethylationAll/BAM/RemodCS7.sorted.bam ../../../MethylationAll/BAM/GlandCS7.sorted.bam -R PromConv_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o InVivo5.mat.gz
#plotProfile -m InVivo5.mat.gz -out InViv5_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S ../../../MethylationAll/BAM/FourCell.sorted.bam ../../../MethylationAll/BAM/TrueTenCellIV.sorted.bam ../../../MethylationAll/BAM/TrueEightCellIV.sorted.bam ../../../MethylationAll/BAM/ElevenCell.sorted.bam ../../../MethylationAll/BAM/EightCellIV.sorted.bam ../../../MethylationAll/BAM/cMor.sorted.bam ../../../MethylationAll/BAM/SYSCS7.sorted.bam ../../../MethylationAll/BAM/StromaCS7.sorted.bam ../../../MethylationAll/BAM/RemodCS7.sorted.bam ../../../MethylationAll/BAM/GlandCS7.sorted.bam -R ConvGB_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o InVivo6.mat.gz
#plotProfile -m InVivo6.mat.gz -out InVivo6_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S ../../../MethylationAll/BAM/FourCell.sorted.bam ../../../MethylationAll/BAM/TrueTenCellIV.sorted.bam ../../../MethylationAll/BAM/TrueEightCellIV.sorted.bam ../../../MethylationAll/BAM/ElevenCell.sorted.bam ../../../MethylationAll/BAM/EightCellIV.sorted.bam ../../../MethylationAll/BAM/cMor.sorted.bam ../../../MethylationAll/BAM/SYSCS7.sorted.bam ../../../MethylationAll/BAM/StromaCS7.sorted.bam ../../../MethylationAll/BAM/RemodCS7.sorted.bam ../../../MethylationAll/BAM/GlandCS7.sorted.bam -R ../CpG.bed --beforeRegionStartLength 1000 --regionBodyLength 5000 --afterRegionStartLength 1000 --skipZeros -o InVivo7.mat.gz
#plotProfile -m InVivo7.mat.gz -out InVivo7_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S ../../../MethylationAll/BAM/FourCell.sorted.bam ../../../MethylationAll/BAM/TrueTenCellIV.sorted.bam ../../../MethylationAll/BAM/TrueEightCellIV.sorted.bam ../../../MethylationAll/BAM/ElevenCell.sorted.bam ../../../MethylationAll/BAM/EightCellIV.sorted.bam ../../../MethylationAll/BAM/cMor.sorted.bam ../../../MethylationAll/BAM/SYSCS7.sorted.bam ../../../MethylationAll/BAM/StromaCS7.sorted.bam ../../../MethylationAll/BAM/RemodCS7.sorted.bam ../../../MethylationAll/BAM/GlandCS7.sorted.bam -R intergenic_line.png --beforeRegionStartLength 1000 --regionBodyLength 5000 --afterRegionStartLength 1000 --skipZeros -o InVivo8.mat.gz
#plotProfile -m InVivo8.mat.gz -out InVivo8_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup


#computeMatrix scale-regions -S EmDisc.bw Am.bw VE.bw ExMes.bw Tb.bw Gland.bw -R PromConv_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o InVivo1.mat.gz
#plotProfile -m InVivo1.mat.gz -out InVivo1_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S EmDisc.bw Am.bw VE.bw ExMes.bw Tb.bw Gland.bw -R ConvGB_1kb.txt --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o InVivo2.mat.gz
#plotProfile -m InVivo2.mat.gz -out InVivo2_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S EmDisc.bw Am.bw VE.bw ExMes.bw Tb.bw Gland.bw -R ../CpG.bed --beforeRegionStartLength 1000 --regionBodyLength 5000 --afterRegionStartLength 1000 --skipZeros -o InVivo3.mat.gz
#plotProfile -m InVivo3.mat.gz -out InVivo3_line.png --numPlotsPerRow 1 --startLabel "TSS-1kbp" --endLabel "TSS + 1kbp" --perGroup

#computeMatrix scale-regions -S EmDisc.bw Am.bw VE.bw ExMes.bw Tb.bw Gland.bw -R ../intergenic.bed --beforeRegionStartLength 1000 --regionBodyLength 5000 --afterRegionStartLength 1000 --skipZeros -o InVivo4.mat.gz
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


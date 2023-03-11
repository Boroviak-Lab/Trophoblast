library(DMRcaller)

CHR <- 2
methylationData_1 <- readBismark(paste('cmTSC1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('cmTSC2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_3 <- readBismark(paste('cmTSC3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_4 <- readBismark(paste('Naive1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_5 <- readBismark(paste('Naive2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_6 <- readBismark(paste('Naive3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationDataList <- joinReplicates(methylationData_1, methylationData_2)
methylationDataList <- joinReplicates(methylationDataList, methylationData_3)
methylationDataList <- joinReplicates(methylationDataList, methylationData_4)
methylationDataList <- joinReplicates(methylationDataList, methylationData_5)
methylationDataList <- joinReplicates(methylationDataList, methylationData_6)
saveRDS(methylationDataList,file=paste("methylationDataList_TbvNiave_",CHR,".rds",sep=""))
condition <- c("a", "a", "a", "b", "b", "b")
DMRsReplicatesBinsCG <- computeDMRsReplicates(methylationData = methylationDataList,
condition = condition, regions = NULL, context = "CG", method = "bins",
binSize = 1000, test = "betareg", pseudocountM = 1, pseudocountN = 2, pValueThreshold = 0.01,
minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_TbvNiave_", CHR, ".rds",sep="") )


CHR <- 22
methylationData_1 <- readBismark(paste('cmTSC1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('cmTSC2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_3 <- readBismark(paste('cmTSC3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_4 <- readBismark(paste('Naive1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_5 <- readBismark(paste('Naive2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_6 <- readBismark(paste('Naive3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationDataList <- joinReplicates(methylationData_1, methylationData_2)
methylationDataList <- joinReplicates(methylationDataList, methylationData_3)
methylationDataList <- joinReplicates(methylationDataList, methylationData_4)
methylationDataList <- joinReplicates(methylationDataList, methylationData_5)
methylationDataList <- joinReplicates(methylationDataList, methylationData_6)
saveRDS(methylationDataList,file=paste("methylationDataList_TbvNiave_",CHR,".rds",sep=""))
condition <- c("a", "a", "a", "b", "b", "b")
DMRsReplicatesBinsCG <- computeDMRsReplicates(methylationData = methylationDataList,
condition = condition, regions = NULL, context = "CG", method = "bins",
binSize = 1000, test = "betareg", pseudocountM = 1, pseudocountN = 2, pValueThreshold = 0.01,
minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_TbvNiave_", CHR, ".rds",sep="") )


library(DMRcaller)

CHR <- 1
methylationData_1 <- readBismark(paste('cmTSC1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('cmTSC2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_3 <- readBismark(paste('cmTSC3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_4 <- readBismark(paste('Primed1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_5 <- readBismark(paste('Primed2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_6 <- readBismark(paste('Primed3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationDataList <- joinReplicates(methylationData_1, methylationData_2)
methylationDataList <- joinReplicates(methylationDataList, methylationData_3)
methylationDataList <- joinReplicates(methylationDataList, methylationData_4)
methylationDataList <- joinReplicates(methylationDataList, methylationData_5)
methylationDataList <- joinReplicates(methylationDataList, methylationData_6)
saveRDS(methylationDataList,file=paste("methylationDataList_TbvPrimed_",CHR,".rds",sep=""))
condition <- c("a", "a", "a", "b", "b", "b")
DMRsReplicatesBinsCG <- computeDMRsReplicates(methylationData = methylationDataList,
condition = condition, regions = NULL, context = "CG", method = "bins",
binSize = 1000, test = "betareg", pseudocountM = 1, pseudocountN = 2, pValueThreshold = 0.01,
minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_TbvPrimed_", CHR, ".rds",sep="") )
#pdf(paste("Profiles_TbvPrimed_", CHR, ".pdf",sep=""))
#par(mar=c(4, 4, 3, 1)+0.1)
#plotMethylationProfileFromData(methylationDataList,conditionsNames = condition, windowSize = 10000, autoscale = FALSE, context = c("CG"))
#dev.off()

CHR <- 2
methylationData_1 <- readBismark(paste('cmTSC1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('cmTSC2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_3 <- readBismark(paste('cmTSC3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_4 <- readBismark(paste('Primed1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_5 <- readBismark(paste('Primed2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_6 <- readBismark(paste('Primed3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationDataList <- joinReplicates(methylationData_1, methylationData_2)
methylationDataList <- joinReplicates(methylationDataList, methylationData_3)
methylationDataList <- joinReplicates(methylationDataList, methylationData_4)
methylationDataList <- joinReplicates(methylationDataList, methylationData_5)
methylationDataList <- joinReplicates(methylationDataList, methylationData_6)
saveRDS(methylationDataList,file=paste("methylationDataList_TbvPrimed_",CHR,".rds",sep=""))
condition <- c("a", "a", "a", "b", "b", "b")
DMRsReplicatesBinsCG <- computeDMRsReplicates(methylationData = methylationDataList,
condition = condition, regions = NULL, context = "CG", method = "bins",
binSize = 1000, test = "betareg", pseudocountM = 1, pseudocountN = 2, pValueThreshold = 0.01,
minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_TbvPrimed_", CHR, ".rds",sep="") )
#pdf(paste("Profiles_TbvPrimed_", CHR, ".pdf",sep=""))
#par(mar=c(4, 4, 3, 1)+0.1)
#plotMethylationProfileFromData(methylationDataList,conditionsNames = condition, windowSize = 10000, autoscale = FALSE, context = c("CG"))
#dev.off()

CHR <- 3
methylationData_1 <- readBismark(paste('cmTSC1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('cmTSC2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_3 <- readBismark(paste('cmTSC3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_4 <- readBismark(paste('Primed1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_5 <- readBismark(paste('Primed2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_6 <- readBismark(paste('Primed3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationDataList <- joinReplicates(methylationData_1, methylationData_2)
methylationDataList <- joinReplicates(methylationDataList, methylationData_3)
methylationDataList <- joinReplicates(methylationDataList, methylationData_4)
methylationDataList <- joinReplicates(methylationDataList, methylationData_5)
methylationDataList <- joinReplicates(methylationDataList, methylationData_6)
saveRDS(methylationDataList,file=paste("methylationDataList_TbvPrimed_",CHR,".rds",sep=""))
condition <- c("a", "a", "a", "b", "b", "b")
DMRsReplicatesBinsCG <- computeDMRsReplicates(methylationData = methylationDataList,
condition = condition, regions = NULL, context = "CG", method = "bins",
binSize = 1000, test = "betareg", pseudocountM = 1, pseudocountN = 2, pValueThreshold = 0.01,
minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_TbvPrimed_", CHR, ".rds",sep="") )
#pdf(paste("Profiles_TbvPrimed_", CHR, ".pdf",sep=""))
#par(mar=c(4, 4, 3, 1)+0.1)
#plotMethylationProfileFromData(methylationDataList,conditionsNames = condition, windowSize = 10000, autoscale = FALSE, context = c("CG"))
#dev.off()

CHR <- 4
methylationData_1 <- readBismark(paste('cmTSC1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('cmTSC2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_3 <- readBismark(paste('cmTSC3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_4 <- readBismark(paste('Primed1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_5 <- readBismark(paste('Primed2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_6 <- readBismark(paste('Primed3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationDataList <- joinReplicates(methylationData_1, methylationData_2)
methylationDataList <- joinReplicates(methylationDataList, methylationData_3)
methylationDataList <- joinReplicates(methylationDataList, methylationData_4)
methylationDataList <- joinReplicates(methylationDataList, methylationData_5)
methylationDataList <- joinReplicates(methylationDataList, methylationData_6)
saveRDS(methylationDataList,file=paste("methylationDataList_TbvPrimed_",CHR,".rds",sep=""))
condition <- c("a", "a", "a", "b", "b", "b")
DMRsReplicatesBinsCG <- computeDMRsReplicates(methylationData = methylationDataList,
condition = condition, regions = NULL, context = "CG", method = "bins",
binSize = 1000, test = "betareg", pseudocountM = 1, pseudocountN = 2, pValueThreshold = 0.01,
minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_TbvPrimed_", CHR, ".rds",sep="") )
#pdf(paste("Profiles_TbvPrimed_", CHR, ".pdf",sep=""))
#par(mar=c(4, 4, 3, 1)+0.1)
#plotMethylationProfileFromData(methylationDataList,conditionsNames = condition, windowSize = 10000, autoscale = FALSE, context = c("CG"))
#dev.off()

CHR <- 5
methylationData_1 <- readBismark(paste('cmTSC1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('cmTSC2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_3 <- readBismark(paste('cmTSC3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_4 <- readBismark(paste('Primed1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_5 <- readBismark(paste('Primed2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_6 <- readBismark(paste('Primed3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationDataList <- joinReplicates(methylationData_1, methylationData_2)
methylationDataList <- joinReplicates(methylationDataList, methylationData_3)
methylationDataList <- joinReplicates(methylationDataList, methylationData_4)
methylationDataList <- joinReplicates(methylationDataList, methylationData_5)
methylationDataList <- joinReplicates(methylationDataList, methylationData_6)
saveRDS(methylationDataList,file=paste("methylationDataList_TbvPrimed_",CHR,".rds",sep=""))
condition <- c("a", "a", "a", "b", "b", "b")
DMRsReplicatesBinsCG <- computeDMRsReplicates(methylationData = methylationDataList,
condition = condition, regions = NULL, context = "CG", method = "bins",
binSize = 1000, test = "betareg", pseudocountM = 1, pseudocountN = 2, pValueThreshold = 0.01,
minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_TbvPrimed_", CHR, ".rds",sep="") )
#pdf(paste("Profiles_TbvPrimed_", CHR, ".pdf",sep=""))
#par(mar=c(4, 4, 3, 1)+0.1)
#plotMethylationProfileFromData(methylationDataList,conditionsNames = condition, windowSize = 10000, autoscale = FALSE, context = c("CG"))
#dev.off()

CHR <- 6
methylationData_1 <- readBismark(paste('cmTSC1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('cmTSC2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_3 <- readBismark(paste('cmTSC3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_4 <- readBismark(paste('Primed1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_5 <- readBismark(paste('Primed2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_6 <- readBismark(paste('Primed3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationDataList <- joinReplicates(methylationData_1, methylationData_2)
methylationDataList <- joinReplicates(methylationDataList, methylationData_3)
methylationDataList <- joinReplicates(methylationDataList, methylationData_4)
methylationDataList <- joinReplicates(methylationDataList, methylationData_5)
methylationDataList <- joinReplicates(methylationDataList, methylationData_6)
saveRDS(methylationDataList,file=paste("methylationDataList_TbvPrimed_",CHR,".rds",sep=""))
condition <- c("a", "a", "a", "b", "b", "b")
DMRsReplicatesBinsCG <- computeDMRsReplicates(methylationData = methylationDataList,
condition = condition, regions = NULL, context = "CG", method = "bins",
binSize = 1000, test = "betareg", pseudocountM = 1, pseudocountN = 2, pValueThreshold = 0.01,
minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_TbvPrimed_", CHR, ".rds",sep="") )
#pdf(paste("Profiles_TbvPrimed_", CHR, ".pdf",sep=""))
#par(mar=c(4, 4, 3, 1)+0.1)
#plotMethylationProfileFromData(methylationDataList,conditionsNames = condition, windowSize = 10000, autoscale = FALSE, context = c("CG"))
#dev.off()

CHR <- 7
methylationData_1 <- readBismark(paste('cmTSC1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('cmTSC2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_3 <- readBismark(paste('cmTSC3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_4 <- readBismark(paste('Primed1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_5 <- readBismark(paste('Primed2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_6 <- readBismark(paste('Primed3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationDataList <- joinReplicates(methylationData_1, methylationData_2)
methylationDataList <- joinReplicates(methylationDataList, methylationData_3)
methylationDataList <- joinReplicates(methylationDataList, methylationData_4)
methylationDataList <- joinReplicates(methylationDataList, methylationData_5)
methylationDataList <- joinReplicates(methylationDataList, methylationData_6)
saveRDS(methylationDataList,file=paste("methylationDataList_TbvPrimed_",CHR,".rds",sep=""))
condition <- c("a", "a", "a", "b", "b", "b")
DMRsReplicatesBinsCG <- computeDMRsReplicates(methylationData = methylationDataList,
condition = condition, regions = NULL, context = "CG", method = "bins",
binSize = 1000, test = "betareg", pseudocountM = 1, pseudocountN = 2, pValueThreshold = 0.01,
minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_TbvPrimed_", CHR, ".rds",sep="") )
#pdf(paste("Profiles_TbvPrimed_", CHR, ".pdf",sep=""))
#par(mar=c(4, 4, 3, 1)+0.1)
#plotMethylationProfileFromData(methylationDataList,conditionsNames = condition, windowSize = 10000, autoscale = FALSE, context = c("CG"))
#dev.off()

CHR <- 8
methylationData_1 <- readBismark(paste('cmTSC1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('cmTSC2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_3 <- readBismark(paste('cmTSC3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_4 <- readBismark(paste('Primed1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_5 <- readBismark(paste('Primed2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_6 <- readBismark(paste('Primed3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationDataList <- joinReplicates(methylationData_1, methylationData_2)
methylationDataList <- joinReplicates(methylationDataList, methylationData_3)
methylationDataList <- joinReplicates(methylationDataList, methylationData_4)
methylationDataList <- joinReplicates(methylationDataList, methylationData_5)
methylationDataList <- joinReplicates(methylationDataList, methylationData_6)
saveRDS(methylationDataList,file=paste("methylationDataList_TbvPrimed_",CHR,".rds",sep=""))
condition <- c("a", "a", "a", "b", "b", "b")
DMRsReplicatesBinsCG <- computeDMRsReplicates(methylationData = methylationDataList,
condition = condition, regions = NULL, context = "CG", method = "bins",
binSize = 1000, test = "betareg", pseudocountM = 1, pseudocountN = 2, pValueThreshold = 0.01,
minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_TbvPrimed_", CHR, ".rds",sep="") )
#pdf(paste("Profiles_TbvPrimed_", CHR, ".pdf",sep=""))
#par(mar=c(4, 4, 3, 1)+0.1)
#plotMethylationProfileFromData(methylationDataList,conditionsNames = condition, windowSize = 10000, autoscale = FALSE, context = c("CG"))
#dev.off()

CHR <- 9
methylationData_1 <- readBismark(paste('cmTSC1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('cmTSC2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_3 <- readBismark(paste('cmTSC3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_4 <- readBismark(paste('Primed1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_5 <- readBismark(paste('Primed2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_6 <- readBismark(paste('Primed3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationDataList <- joinReplicates(methylationData_1, methylationData_2)
methylationDataList <- joinReplicates(methylationDataList, methylationData_3)
methylationDataList <- joinReplicates(methylationDataList, methylationData_4)
methylationDataList <- joinReplicates(methylationDataList, methylationData_5)
methylationDataList <- joinReplicates(methylationDataList, methylationData_6)
saveRDS(methylationDataList,file=paste("methylationDataList_TbvPrimed_",CHR,".rds",sep=""))
condition <- c("a", "a", "a", "b", "b", "b")
DMRsReplicatesBinsCG <- computeDMRsReplicates(methylationData = methylationDataList,
condition = condition, regions = NULL, context = "CG", method = "bins",
binSize = 1000, test = "betareg", pseudocountM = 1, pseudocountN = 2, pValueThreshold = 0.01,
minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_TbvPrimed_", CHR, ".rds",sep="") )
#pdf(paste("Profiles_TbvPrimed_", CHR, ".pdf",sep=""))
#par(mar=c(4, 4, 3, 1)+0.1)
#plotMethylationProfileFromData(methylationDataList,conditionsNames = condition, windowSize = 10000, autoscale = FALSE, context = c("CG"))
#dev.off()

CHR <- 10
methylationData_1 <- readBismark(paste('cmTSC1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('cmTSC2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_3 <- readBismark(paste('cmTSC3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_4 <- readBismark(paste('Primed1.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_5 <- readBismark(paste('Primed2.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_6 <- readBismark(paste('Primed3.1_val_1_bismark_bt2_pe.deduplicated.CX_report_CHR', CHR, '.txt', sep="") )
methylationDataList <- joinReplicates(methylationData_1, methylationData_2)
methylationDataList <- joinReplicates(methylationDataList, methylationData_3)
methylationDataList <- joinReplicates(methylationDataList, methylationData_4)
methylationDataList <- joinReplicates(methylationDataList, methylationData_5)
methylationDataList <- joinReplicates(methylationDataList, methylationData_6)
saveRDS(methylationDataList,file=paste("methylationDataList_TbvPrimed_",CHR,".rds",sep=""))
condition <- c("a", "a", "a", "b", "b", "b")
DMRsReplicatesBinsCG <- computeDMRsReplicates(methylationData = methylationDataList,
condition = condition, regions = NULL, context = "CG", method = "bins",
binSize = 1000, test = "betareg", pseudocountM = 1, pseudocountN = 2, pValueThreshold = 0.01,
minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_TbvPrimed_", CHR, ".rds",sep="") )
#pdf(paste("Profiles_TbvPrimed_", CHR, ".pdf",sep=""))
#par(mar=c(4, 4, 3, 1)+0.1)
#plotMethylationProfileFromData(methylationDataList,conditionsNames = condition, windowSize = 10000, autoscale = FALSE, context = c("CG"))
#dev.off()

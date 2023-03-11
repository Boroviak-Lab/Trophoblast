library(DMRcaller)


CHR <- 20

#methylationData_1 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
#methylationData_2 <- readBismark(paste('cMor.sorted.CX_report_CHR', CHR, '.txt', sep="") )
#DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score",
#pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
#saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_8C_cMor_", CHR, ".rds",sep="") )
methylationData_1 <- readBismark(paste('cMor.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('EmDisc_CS5.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score",
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_cMor_EmDisc", CHR, ".rds",sep="") )
methylationData_1 <- readBismark(paste('EmDisc_CS5.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('Tb_CS5.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score",
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_EmDisc_Tb", CHR, ".rds",sep="") )

CHR <- 21
methylationData_1 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('cMor.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score",
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_8C_cMor_", CHR, ".rds",sep="") )
methylationData_1 <- readBismark(paste('cMor.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('EmDisc_CS5.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score",
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_cMor_EmDisc", CHR, ".rds",sep="") )
methylationData_1 <- readBismark(paste('EmDisc_CS5.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('Tb_CS5.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score",
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_EmDisc_Tb", CHR, ".rds",sep="") )

CHR <- 22
methylationData_1 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('cMor.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score",
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_8C_cMor_", CHR, ".rds",sep="") )
methylationData_1 <- readBismark(paste('cMor.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('EmDisc_CS5.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score",
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_cMor_EmDisc", CHR, ".rds",sep="") )
methylationData_1 <- readBismark(paste('EmDisc_CS5.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('Tb_CS5.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score",
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_EmDisc_Tb", CHR, ".rds",sep="") )


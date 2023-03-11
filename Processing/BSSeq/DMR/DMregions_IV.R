library(DMRcaller)

#CHR <- 1

#methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
#methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
#DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
#pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
#saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

#CHR <- 2
#methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
#methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
#DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
#pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
#saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 3
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 4
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 5
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 6
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 7
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 8
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 9
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 10
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 11
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 12
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 13
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 14
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 15
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 16
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 17
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 18
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 19
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 20
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 21
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score", 
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

CHR <- 22
methylationData_1 <- readBismark(paste('TrueEightCellIV.sorted.CX_report_CHR', CHR, '.txt', sep="") )
methylationData_2 <- readBismark(paste('TrueEightCell.sorted.CX_report_CHR', CHR, '.txt', sep="") )
DMRsReplicatesBinsCG <- computeDMRs(methylationData_1,methylationData_2, regions = NULL, context = "CG", method = "bins",binSize = 100, test = "score",
pValueThreshold = 0.01, minCytosinesCount = 4, minProportionDifference = 0.1, minGap = 0, minSize = 50, minReadsPerCytosine = 4, cores = 8)
saveRDS(DMRsReplicatesBinsCG,file=paste("DMR_IVF_", CHR, ".rds",sep="") )

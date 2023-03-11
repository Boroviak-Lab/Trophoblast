
D1 <- readRDS('DMR_8C_cMor_1.rds')
D2 <- readRDS('DMR_8C_cMor_2.rds')
D3 <- readRDS('DMR_8C_cMor_3.rds')
D4 <- readRDS('DMR_8C_cMor_4.rds')
D5 <- readRDS('DMR_8C_cMor_5.rds')
D6 <- readRDS('DMR_8C_cMor_6.rds')
D7 <- readRDS('DMR_8C_cMor_7.rds')
D8 <- readRDS('DMR_8C_cMor_8.rds')
D9 <- readRDS('DMR_8C_cMor_9.rds')
D10 <- readRDS('DMR_8C_cMor_10.rds')
D11 <- readRDS('DMR_8C_cMor_11.rds')
D12 <- readRDS('DMR_8C_cMor_12.rds')
D13 <- readRDS('DMR_8C_cMor_13.rds')
D14 <- readRDS('DMR_8C_cMor_14.rds')
D15 <- readRDS('DMR_8C_cMor_15.rds')
D16 <- readRDS('DMR_8C_cMor_16.rds')
D17 <- readRDS('DMR_8C_cMor_17.rds')
D18 <- readRDS('DMR_8C_cMor_18.rds')
D19 <- readRDS('DMR_8C_cMor_19.rds')
D20 <- readRDS('DMR_8C_cMor_20.rds')
D21 <- readRDS('DMR_8C_cMor_21.rds')
D22 <- readRDS('DMR_8C_cMor_22.rds')

MD <- rbind(as.data.frame(D1),as.data.frame(D1),as.data.frame(D3),as.data.frame(D4),as.data.frame(D5),as.data.frame(D6),
as.data.frame(D7),as.data.frame(D8),as.data.frame(D9),as.data.frame(D10),as.data.frame(D11),as.data.frame(D12),
as.data.frame(D13),as.data.frame(D14),as.data.frame(D15),as.data.frame(D16),as.data.frame(D17),as.data.frame(D18),
as.data.frame(D19),as.data.frame(D20),as.data.frame(D21),as.data.frame(D22))

I1 <- which(MD$regionType=="gain")
I2 <- which(MD$regionType=="loss")

write.table(MD[I1,c("seqnames","start","end","direction","pValue") ], file = "All8C_Mor_Gain.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
write.table(MD[I2,c("seqnames","start","end","direction","pValue") ], file = "All8C_Mor_Loss.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)


D1 <- readRDS('DMR_EmDisc_Tb1.rds')
D2 <- readRDS('DMR_EmDisc_Tb2.rds')
D3 <- readRDS('DMR_EmDisc_Tb3.rds')
D4 <- readRDS('DMR_EmDisc_Tb4.rds')
D5 <- readRDS('DMR_EmDisc_Tb5.rds')
D6 <- readRDS('DMR_EmDisc_Tb6.rds')
D7 <- readRDS('DMR_EmDisc_Tb7.rds')
D8 <- readRDS('DMR_EmDisc_Tb8.rds')
D9 <- readRDS('DMR_EmDisc_Tb9.rds')
D10 <- readRDS('DMR_EmDisc_Tb10.rds')
D11 <- readRDS('DMR_EmDisc_Tb11.rds')
D12 <- readRDS('DMR_EmDisc_Tb12.rds')
D13 <- readRDS('DMR_EmDisc_Tb13.rds')
D14 <- readRDS('DMR_EmDisc_Tb14.rds')
D15 <- readRDS('DMR_EmDisc_Tb15.rds')
D16 <- readRDS('DMR_EmDisc_Tb16.rds')
D17 <- readRDS('DMR_EmDisc_Tb17.rds')
D18 <- readRDS('DMR_EmDisc_Tb18.rds')
D19 <- readRDS('DMR_EmDisc_Tb19.rds')
D20 <- readRDS('DMR_EmDisc_Tb20.rds')
D21 <- readRDS('DMR_EmDisc_Tb21.rds')
D22 <- readRDS('DMR_EmDisc_Tb22.rds')

MD <- rbind(as.data.frame(D1),as.data.frame(D1),as.data.frame(D3),as.data.frame(D4),as.data.frame(D5),as.data.frame(D6),
as.data.frame(D7),as.data.frame(D8),as.data.frame(D9),as.data.frame(D10),as.data.frame(D11),as.data.frame(D12),
as.data.frame(D13),as.data.frame(D14),as.data.frame(D15),as.data.frame(D16),as.data.frame(D17),as.data.frame(D18),
as.data.frame(D19),as.data.frame(D20),as.data.frame(D21),as.data.frame(D22))

I1 <- which(MD$regionType=="gain")
I2 <- which(MD$regionType=="loss")

write.table(MD[I1,c("seqnames","start","end","direction","pValue") ], file = "AllEmD_Tb_Gain.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
write.table(MD[I2,c("seqnames","start","end","direction","pValue") ], file = "AllEmD_Tb_Loss.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)



D1 <- readRDS('DMR_cMor_EmDisc1.rds')
D2 <- readRDS('DMR_cMor_EmDisc2.rds')
D3 <- readRDS('DMR_cMor_EmDisc3.rds')
D4 <- readRDS('DMR_cMor_EmDisc4.rds')
D5 <- readRDS('DMR_cMor_EmDisc5.rds')
D6 <- readRDS('DMR_cMor_EmDisc6.rds')
D7 <- readRDS('DMR_cMor_EmDisc7.rds')
D8 <- readRDS('DMR_cMor_EmDisc8.rds')
D9 <- readRDS('DMR_cMor_EmDisc9.rds')
D10 <- readRDS('DMR_cMor_EmDisc10.rds')
D11 <- readRDS('DMR_cMor_EmDisc11.rds')
D12 <- readRDS('DMR_cMor_EmDisc12.rds')
D13 <- readRDS('DMR_cMor_EmDisc13.rds')
D14 <- readRDS('DMR_cMor_EmDisc14.rds')
D15 <- readRDS('DMR_cMor_EmDisc15.rds')
D16 <- readRDS('DMR_cMor_EmDisc16.rds')
D17 <- readRDS('DMR_cMor_EmDisc17.rds')
D18 <- readRDS('DMR_cMor_EmDisc18.rds')
D19 <- readRDS('DMR_cMor_EmDisc19.rds')
D20 <- readRDS('DMR_cMor_EmDisc20.rds')
D21 <- readRDS('DMR_cMor_EmDisc21.rds')
D22 <- readRDS('DMR_cMor_EmDisc22.rds')

MD <- rbind(as.data.frame(D1),as.data.frame(D1),as.data.frame(D3),as.data.frame(D4),as.data.frame(D5),as.data.frame(D6),
as.data.frame(D7),as.data.frame(D8),as.data.frame(D9),as.data.frame(D10),as.data.frame(D11),as.data.frame(D12),
as.data.frame(D13),as.data.frame(D14),as.data.frame(D15),as.data.frame(D16),as.data.frame(D17),as.data.frame(D18),
as.data.frame(D19),as.data.frame(D20),as.data.frame(D21),as.data.frame(D22))

I1 <- which(MD$regionType=="gain")
I2 <- which(MD$regionType=="loss")

write.table(MD[I1,c("seqnames","start","end","direction","pValue") ], file = "AllcMor_EmDisc_Gain.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
write.table(MD[I2,c("seqnames","start","end","direction","pValue") ], file = "AllcMor_EmDisc_Loss.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)



D1 <- readRDS('DMR_IVF_1.rds')
D2 <- readRDS('DMR_IVF_2.rds')
D3 <- readRDS('DMR_IVF_3.rds')
D4 <- readRDS('DMR_IVF_4.rds')
D5 <- readRDS('DMR_IVF_5.rds')
D6 <- readRDS('DMR_IVF_6.rds')
D7 <- readRDS('DMR_IVF_7.rds')
D8 <- readRDS('DMR_IVF_8.rds')
D9 <- readRDS('DMR_IVF_9.rds')
D10 <- readRDS('DMR_IVF_10.rds')
D11 <- readRDS('DMR_IVF_11.rds')
D12 <- readRDS('DMR_IVF_12.rds')
D13 <- readRDS('DMR_IVF_13.rds')
D14 <- readRDS('DMR_IVF_14.rds')
D15 <- readRDS('DMR_IVF_15.rds')
D16 <- readRDS('DMR_IVF_16.rds')
D17 <- readRDS('DMR_IVF_17.rds')
D18 <- readRDS('DMR_IVF_18.rds')
D19 <- readRDS('DMR_IVF_19.rds')
D20 <- readRDS('DMR_IVF_20.rds')
D21 <- readRDS('DMR_IVF_21.rds')
D22 <- readRDS('DMR_IVF_22.rds')

MD <- rbind(as.data.frame(D1),as.data.frame(D1),as.data.frame(D3),as.data.frame(D4),as.data.frame(D5),as.data.frame(D6),
as.data.frame(D7),as.data.frame(D8),as.data.frame(D9),as.data.frame(D10),as.data.frame(D11),as.data.frame(D12),
as.data.frame(D13),as.data.frame(D14),as.data.frame(D15),as.data.frame(D16),as.data.frame(D17),as.data.frame(D18),
as.data.frame(D19),as.data.frame(D20),as.data.frame(D21),as.data.frame(D22))

I1 <- which(MD$regionType=="gain")
I2 <- which(MD$regionType=="loss")

write.table(MD[I1,c("seqnames","start","end","direction","pValue") ], file = "AllIVF_Gain.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
write.table(MD[I2,c("seqnames","start","end","direction","pValue") ], file = "AllVF_Loss.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)


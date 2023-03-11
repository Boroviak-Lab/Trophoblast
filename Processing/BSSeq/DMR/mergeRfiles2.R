D1 <- readRDS('DMR_TbvNiave_1.rds')
D2 <- readRDS('DMR_TbvNiave_2.rds')
D3 <- readRDS('DMR_TbvNiave_3.rds')
D4 <- readRDS('DMR_TbvNiave_4.rds')
D5 <- readRDS('DMR_TbvNiave_5.rds')
D6 <- readRDS('DMR_TbvNiave_6.rds')
D7 <- readRDS('DMR_TbvNiave_7.rds')
D8 <- readRDS('DMR_TbvNiave_8.rds')
D9 <- readRDS('DMR_TbvNiave_9.rds')
D10 <- readRDS('DMR_TbvNiave_10.rds')
D11 <- readRDS('DMR_TbvNiave_11.rds')
D12 <- readRDS('DMR_TbvNiave_12.rds')
D13 <- readRDS('DMR_TbvNiave_13.rds')
D14 <- readRDS('DMR_TbvNiave_14.rds')
D15 <- readRDS('DMR_TbvNiave_15.rds')
D16 <- readRDS('DMR_TbvNiave_16.rds')
D17 <- readRDS('DMR_TbvNiave_17.rds')
D18 <- readRDS('DMR_TbvNiave_18.rds')
D19 <- readRDS('DMR_TbvNiave_19.rds')
D20 <- readRDS('DMR_TbvNiave_20.rds')
D21 <- readRDS('DMR_TbvNiave_21.rds')
D22 <- readRDS('DMR_TbvNiave_22.rds')

MD <- rbind(as.data.frame(D1),as.data.frame(D1),as.data.frame(D3),as.data.frame(D4),as.data.frame(D5),as.data.frame(D6),
as.data.frame(D7),as.data.frame(D8),as.data.frame(D9),as.data.frame(D10),as.data.frame(D11),as.data.frame(D12),
as.data.frame(D13),as.data.frame(D14),as.data.frame(D15),as.data.frame(D16),as.data.frame(D17),as.data.frame(D18),
as.data.frame(D19),as.data.frame(D20),as.data.frame(D21),as.data.frame(D22))

I1 <- which(MD$regionType=="gain")
I2 <- which(MD$regionType=="loss")

write.table(MD[I1,c("seqnames","start","end","direction","pValue") ], file = "AllTB_Niave_Gain.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
write.table(MD[I2,c("seqnames","start","end","direction","pValue") ], file = "AllTB_Niave_Loss.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)


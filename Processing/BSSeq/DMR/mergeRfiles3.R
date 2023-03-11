D1 <- readRDS('DMR_PrimedvNiave_1.rds')
D2 <- readRDS('DMR_PrimedvNiave_2.rds')
D3 <- readRDS('DMR_PrimedvNiave_3.rds')
D4 <- readRDS('DMR_PrimedvNiave_4.rds')
D5 <- readRDS('DMR_PrimedvNiave_5.rds')
D6 <- readRDS('DMR_PrimedvNiave_6.rds')
D7 <- readRDS('DMR_PrimedvNiave_7.rds')
D8 <- readRDS('DMR_PrimedvNiave_8.rds')
D9 <- readRDS('DMR_PrimedvNiave_9.rds')
D10 <- readRDS('DMR_PrimedvNiave_10.rds')
D11 <- readRDS('DMR_PrimedvNiave_11.rds')
D12 <- readRDS('DMR_PrimedvNiave_12.rds')
D13 <- readRDS('DMR_PrimedvNiave_13.rds')
D14 <- readRDS('DMR_PrimedvNiave_14.rds')
D15 <- readRDS('DMR_PrimedvNiave_15.rds')
D16 <- readRDS('DMR_PrimedvNiave_16.rds')
D17 <- readRDS('DMR_PrimedvNiave_17.rds')
D18 <- readRDS('DMR_PrimedvNiave_18.rds')
D19 <- readRDS('DMR_PrimedvNiave_19.rds')
D20 <- readRDS('DMR_PrimedvNiave_20.rds')
D21 <- readRDS('DMR_PrimedvNiave_21.rds')
D22 <- readRDS('DMR_PrimedvNiave_22.rds')

MD <- rbind(as.data.frame(D1),as.data.frame(D1),as.data.frame(D3),as.data.frame(D4),as.data.frame(D5),as.data.frame(D6),
as.data.frame(D7),as.data.frame(D8),as.data.frame(D9),as.data.frame(D10),as.data.frame(D11),as.data.frame(D12),
as.data.frame(D13),as.data.frame(D14),as.data.frame(D15),as.data.frame(D16),as.data.frame(D17),as.data.frame(D18),
as.data.frame(D19),as.data.frame(D20),as.data.frame(D21),as.data.frame(D22))

I1 <- which(MD$regionType=="gain")
I2 <- which(MD$regionType=="loss")

write.table(MD[I1,c("seqnames","start","end","direction","pValue") ], file = "AllPrimed_Niave_Gain.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
write.table(MD[I2,c("seqnames","start","end","direction","pValue") ], file = "AllPrimed_Niave_Loss.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)


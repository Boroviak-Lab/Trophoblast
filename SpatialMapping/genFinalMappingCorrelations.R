library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)
library(pheatmap)
library(reshape2)

set.seed(1)

saveext = "Output"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

marmoset_dataInVivo <- readRDS(file="Data/InVivoRef.rds")
marmoset_data_InVitro <- readRDS(file="Data/InVitro.rds")

mammal.anchors <- FindIntegrationAnchors(object.list = list(marmoset_dataInVivo, marmoset_data_InVitro), dims = 1:20, anchor.features = 5000, k.filter = 100)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

newID <- as.character(Idents(mammal.combined))
newID[which(newID=="EmDisc_CS5_Am")] <- "EmDisc_CS5"
newID[which(newID=="EmDisc_gast_CS6 ")] <- "EmDisc_CS6"
newID[which(newID=="EmDisc_Gast_CS6")] <- "EmDisc_CS6"
newID[which(newID=="EmDisc_stalk_CS6")] <- "Stalk_CS6"
newID[which(newID=="EmDisc_stalk_CS7")] <- "Stalk_CS7"
newID[which(newID=="Stalk_CS7_PGC")] <- "PGC_CS7"
newID[which(newID=="EmDisc_CS7_PGC")] <- "PGC_CS7"
newID[which(newID=="ExMes_stalk_CS7")] <- "Stalk_CS7"
newID[which(newID=="EmDisc_CS7_Am")] <- "EmDisc_CS7"
newID[which(newID=="EmDisc_CS7_Am")] <- "EmDisc_CS7"
newID[which(newID=="Am_CS6_EmDisc")] <- "Am_CS6"
newID[which(newID=="EmDisc_CS6_Gast")] <- "EmDisc_CS6"
newID[which(newID=="EmDisc_gast_CS6")] <- "EmDisc_CS6"
newID[which(newID=="EmDisc_CS6_Am")] <- "EmDisc_CS6"
newID[which(newID=="EmDisc_Stalk_CS6")] <- "Stalk_CS6"
newID[which(newID=="EmDisc_Gast_CS7")] <- "EmDisc_CS7"
newID[which(newID=="EmDisc_CS7_Gast")] <- "EmDisc_CS7"
newID[which(newID=="EmDisc_CS6_Gast")] <- "EmDisc_CS6"
newID[which(newID=="Am_CS7_EmDisc")] <- "Am_CS7"
newID[which(newID=="cMor_CS3")] <- "cMor_CS2"

newID[which(newID=="ESC_CMESC")] <- "Conventional1"
newID[which(newID=="ESC_conv2")] <- "Conventional2"
newID[which(newID=="EMS3_PLAXA")] <- "PLAXA1"


Idents(mammal.combined) <- newID

cType <- c("OKAEP5esc","CM_TSF2","CM_TSC","cmTSCPAVSbub","TbLC","PAVSCHIR","cmTSCPAVSiwp2","cmTSCPAVS","PLAXA","PLAXA1","HYPO","newTSP3","Conventional","Conventional1","Conventional2","OKAEP5esc","Am","BMP_noMEF","2307","2308","Other","Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS2","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","Stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","Am","Amnoid_bead","BMP_MEF","BMP_noMEF","EmD","EmDisc","ActA_MEF","ActA_noMEF","SB43_MEF","CHIR_MEF","FGF_noMEF","Am_CS5_PGC","Am_CS5_ExMesPGC","EmDisc_CS5_Am","PGC_CS6","Stalk_CS7_PGC","EmDisc_CS7_PGC","EmDisc_CS7_Am","Am_CS6_EmDisc","EmDisc_CS5_Gast","EmDisc_CS6_Am","EmDisc_CS6_Gast","EmDisc_CS7_Gast","Stalk_CS7_PGC","Gland_CS6_","ReStroma_CS5","ReStroma_CS6","ReStroma_CS7","Stroma_CS5","Stroma_CS6","Stroma_CS7","VE_CS6","VE_CS7","Am_CS7_EmDisc","EmDisc_Gast_CS6","EmDisc_Gast_CS7","Stalk_CS6")
BaseCol <- c("#a8870f","black","green","#e607a5","#a8870f","#8a2f6f","#8017c2","#BF0489","#00BFBF","#00BFBF","#E6B500","#7108a6","#0233BF","#0233BF","#0233BF","#c49a00","#877bd6","#1A0873","lightgrey","lightgrey","lightgrey","#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#7b3294","#c2a5cf","#a6dba0","#008837","#ca0020","#f4a582","#92c5de","#0571b0","#d8b365","#5ab4ac","#4d4d4d","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#0233BF","#5F54C7","#0c9cf5","#0767DA","#0767DA","#0233BF","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0",'#D74404',"black","#5F54C7","#0767DA","#0233BF","#c49a00")


colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}

coluse <- BaseCol[colind]
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_All_inteset_5kfilt50-allESC-wall_rr2.pdf",sep=""),width = 20, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_All_inteset_5kfilt50-allESC-wall_rr2.pdf",sep=""),width = 20, height = 8, useDingbats = FALSE)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "tsne", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/TSNE_All_inteset_5kfilt50-allESC-wall_rr2.pdf",sep=""),width = 20, height = 8, useDingbats = FALSE)


DefaultAssay(mammal.combined) <- "RNA"

FeaturePlot(mammal.combined,  reduction = "umap", features = "JUP", split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/UMAP_JUP_alignstages_rr.pdf",sep=""),width = 20, height = 8)

FeaturePlot(mammal.combined,  reduction = "pca", features = "JUP", split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/PCA_JUP_alignstages_rr.pdf",sep=""),width = 20, height = 8)


DefaultAssay(mammal.combined) <- "integrated"
Data2 <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20, k.param = 100)

KNN_K100 <- Data2@graphs$integrated_nn
SNN_K100 <- Data2@graphs$integrated_snn

Data2 <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20, k.param = 50)

KNN_K50 <- Data2@graphs$integrated_nn
SNN_K50 <- Data2@graphs$integrated_snn

Data2 <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20, k.param = 150)

KNN_K150 <- Data2@graphs$integrated_nn
SNN_K150 <- Data2@graphs$integrated_snn

saveRDS(KNN_K100,file=paste(saveext,"KNN_K100_set1_wall_rr.rds",sep=""))
saveRDS(SNN_K100,file=paste(saveext,"SNN_K100_set1_wall_rr.rds",sep=""))
saveRDS(KNN_K150,file=paste(saveext,"KNN_K150_set1_wall_rr.rds",sep=""))
saveRDS(SNN_K150,file=paste(saveext,"SNN_K150_set1_wall_rr.rds",sep=""))
saveRDS(KNN_K50,file=paste(saveext,"KNN_K50_set1_wall_rr.rds",sep=""))
saveRDS(SNN_K50,file=paste(saveext,"SNN_K50_set1_wall_rr.rds",sep=""))

library(gtools)
k = 10
inds <- which(mammal.combined$Dataset=="2) Marmoset in vivo") #seq(from=1,to=length(colnames(marmoset_dataInVivoAlt1)),by=1)
inds1 <- which(mammal.combined$Dataset=="InVitro") #seq(from=1,to=length(colnames(marmoset_dataInVivoAlt1)),by=1)

C1 <- matrix(, ncol = k, nrow = length(Idents(mammal.combined)))
S1 <- matrix(, ncol = k, nrow = length(Idents(mammal.combined)))
C1p <- matrix(, ncol = k, nrow = length(Idents(mammal.combined)))

lC1 <- matrix(, ncol = k, nrow = length(Idents(mammal.combined)))
lC1p <- matrix(, ncol = k, nrow = length(Idents(mammal.combined)))


Locations <- as.data.frame(mammal.combined$LOC)
for (i in 1:length(Idents(mammal.combined)) ) {
  a1 <- sort(SNN_K50[i,inds],decreasing = TRUE)[1:k]
  C1[i,] <- names(a1)
  S1[i,] <- a1
  lC1[i,] <- as.character(Locations[names(a1),])
  aperm <- SNN_K50[i,]
  names(aperm) <- permute(names(aperm))
  a2 <- sort(aperm[inds],decreasing = TRUE)[1:k]
  C1p[i,] <- names(a2)
  
  lC1p[i,] <- as.character(Locations[names(a2),])
}
write.table(Idents(mammal.combined)[inds1],file=paste(saveext,'AllMapIdents_set1_all_rr.csv',sep=""),sep=",")
write.table(mammal.combined$LOC[inds1],file=paste(saveext,'AllMapIdents_set1_all_rr.csv',sep=""),sep=",")


write.table(data.frame(as.data.frame(colnames(mammal.combined)[inds1]),as.data.frame(mammal.combined$LOC[inds1]),as.data.frame(Idents(mammal.combined)[inds1]), as.data.frame(C1[inds1,]), as.data.frame(lC1[inds1,]), 1 ),file=paste(saveext,'C1_set1_rr.csv',sep=""),sep=",",row.names = F)
write.table(data.frame(as.data.frame(colnames(mammal.combined)[inds1]),as.data.frame(mammal.combined$LOC[inds1]),as.data.frame(Idents(mammal.combined)[inds1]), as.data.frame(C1p[inds1,]),as.data.frame(lC1p[inds1,]),1 ),file=paste(saveext,'C1perm_set1_rr.csv',sep=""),sep=",",row.names = F)
write.table(S1[inds1,],file=paste(saveext,'S1_set1_all_rr.csv',sep=""),sep=",")

write.table(mammal.combined$LOC,file=paste(saveext,'AllLocKey_all_rr.csv',sep=""),sep=",")
write.table(Idents(mammal.combined),file=paste(saveext,'AllLocIDKey_all_rr.csv',sep=""),sep=",")

#Theses are files that work in the MATLAB counterparat
DefaultAssay(mammal.combined) <- "integrated"

d0 <- GetAssayData(mammal.combined, assay = "integrated")
d1 <- d0[,which(mammal.combined$Dataset=="2) Marmoset in vivo")]
d2 <- d0[,which(mammal.combined$Dataset=="InVitro")]

Locs1 <- mammal.combined$LOC[which(mammal.combined$Dataset=="2) Marmoset in vivo")] #subset(mammal.combined, cells=colnames(marmoset_dataInVivo))$LOC
Locs2 <- mammal.combined$LOC[which(mammal.combined$Dataset=="InVitro")] #subset(mammal.combined, cells=colnames(marmoset_data_InVitroAlt2))$LOC

ID1 <- Idents(mammal.combined)[which(mammal.combined$Dataset=="2) Marmoset in vivo")] #Idents(subset(mammal.combined, cells=colnames(marmoset_dataInVivo)))
ID2 <- Idents(mammal.combined)[which(mammal.combined$Dataset=="InVitro")] #Idents(subset(mammal.combined, cells=colnames(marmoset_data_InVitroAlt2)))

write.csv(data.frame(Loc=Locs1,Type=ID1,z=1), file=paste(saveext,"/Correlation_int_locations1_wall_rr2.csv",sep=""))
write.csv(data.frame(Loc=Locs2,Type=ID2,z=1), file=paste(saveext,"/Correlation_int_locations2_wall_rr2.csv",sep=""))

C <- cor(as.data.frame(d1),as.data.frame(d2),method="pearson")
write.csv(as.data.frame(C), file=paste(saveext,"/Correlation_int_set1_wall_rr2.csv",sep=""))

path_data=#working directory
setwd(path_data)

set.seed(1234567)  
library(RColorBrewer)
library(ggplot2)
library(mclust)
library(pcaMethods)
library(umap)
library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)
library(ComplexHeatmap)
library('ggrepel')

cols <- c(
  ## preimplantation
  "Zy_CS1" = "#56E600",
  "4-cell_CS2" = "#48BF00",
  "8-cell_CS2" = "#00BF30",
  "cMor_CS3" = "#009926",
  "ICM_CS3" = "#00E6E6",
  "Epi_CS3" = "#00BFBF",
  "Hyp_CS3" = "#e6cb00",
  "Tb_CS3" = "#BF0489", 
  "TE" = "#BF0489", 
  
  ## Postimplantation 
  #Embryonic
  "EmDisc_CS5" =	"#0767DA",
  "EmDisc_CS6" = "#0767DA",
  "EmDisc_CS7" = "#0767DA",
  "PGC_CS5" = "#BCBF02",
  "PGC_CS6" = "#BCBF02",
  "PGC_CS7" = "#BCBF02",
  
  # Extraembryonic
  "Am_CS5" = "#5F54C7",
  "Am_CS6" = "#5F54C7",
  "Am_CS7" = "#5F54C7",
  "Am" = "#5F54C7",
  
  "VE_CS5" = "#F04C04",
  "VE_CS6" = "#F04C04",
  "VE_CS7" = "#F04C04",
  
  "SYS_CS5" = "#d17600",
  "SYS_CS6" = "#d17600",
  "SYS_CS7" = "#d17600",
  "ExMes_CS5" = "#c49a00",
  "ExMes_CS6" = "#c49a00",
  "ExMes_CS7" = "#c49a00",
  "ExMes_stalk_CS7" = "#967700",
  "Stalk_CS6" = "#967700",
  #Trophoblast
  "CTB"="#b832c2",
  "Tb_CS5"="#b832c2",
  "Tb_CS6"="#b832c2",
  "Tb_CS7"="#b832c2",
  "Tb_embryonal"="#b832c2",
  "EVT" = "#015244",
  "STB"="red",
  "CTB_twin"="#824087",
  "Tb_luminal"="#824087",
  "SYNC"="#921FE6",
  "Tb_abembryonal_CS5" = "#331157",
  "Tb_abembryonal_CS7" = "#331157",
  
  ## Maternal
  "Gland.E15" = "grey",
  "Stroma.E15"="grey",
  "ReGland.E15"="grey",
  "ReStroma.E15"="grey",
  "Myo.E15"="grey",
  "Gland.CS7." ="grey",
  "Stroma.CS7."="grey",
  "ReGland.CS7."="grey",
  "ReStroma.CS7."="grey",
  "Myo.CS7."="grey",
  
  ## in vitro culture
  "Conventional" = "#0233BF",
  "ESC_conv2" = "#0233BF",
  "ESC_CMESC" = "#0233BF",
  "EmDisc_#1" = "#029cbf",
  "EmDisc_#2" = "#029cbf",
  "PLAXA" = "#057699",
  "EMS3_PLAXA" = "#00BFBF",
  "PLAXC" = rgb(242/255, 184/255, 40/255, 1),
  "PLA" = "grey",
  
  ## in vitro culture
  "Amnioid_esc" = "grey",
  "Amnioid_#1" = "#0233BF",
  "Amnioid_#2" = "#0233BF",
  "Amnoid_bead" = "pink",
  "BMP_MEF" = "#3b02bf",
  "BMP_noMEF" = "#3b02bf",
  "HYPO" = "#8B8000",
  "NOG" = "magenta",
  "OKAEP5esc" = "pink",
  "newTSP3" = "#FF00FF",
  "cmTSCPAVS" = "green",
  "cmTSCPAVSiwp2" = "yellow",
  "TbCL" = "red"
    
  
)

#load in vivo dataset
invivo2_DATA<- readRDS('Amnioids_aligned.rds')
invivo_DATAraw<- readRDS('InVivo2.rds')
invivo_DATA<- subset(x=invivo_DATAraw,idents=c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Hyp_CS3","Tb_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","Stalk_CS6","SYS_CS5","SYS_CS6","SYS_CS7","PGC_CS7","Tb_CS5","Tb_CS6","Tb_CS7","Tb_abembryonal_CS7"))
invivoAMinVIVO_DATA<- subset(x=invivo2_DATA,idents=unique(Idents(invivo2_DATA))[grep("EmDisc_CS5|EmDisc_CS6|EmDisc_CS7|Am_CS5|Am_CS6|Am_CS7",unique(Idents(invivo2_DATA)))])
invivo_DATA <- subset(invivo_DATA, subset = nFeature_RNA > 0)

#manual reannotate abembryonal
toadd_abemb<-Idents(invivo_DATA)
toadd_abemb[grep("E25A_195|E25A_190_N",invivo_DATA$LOC)]<-"Tb_abembryonal_CS7"
Idents(invivo_DATA)<-toadd_abemb

#normalization and basic dimensonality reduction
invivo_DATA <- NormalizeData(invivo_DATA, verbose = FALSE)
invivo_DATA <- FindVariableFeatures(invivo_DATA, selection.method = "vst", nfeatures = 10000)
invivo_DATA <- ScaleData(invivo_DATA, verbose = FALSE)
invivo_DATA <- RunPCA(invivo_DATA, npcs = 20, verbose = FALSE)
invivo_DATA <- RunUMAP(invivo_DATA, reduction = "pca", dims = 1:20)
invivo_DATA <- RunTSNE(invivo_DATA, reduction = "pca", dims = 1:20)

#Visualize dataset and check clustering visually with lineage specific genes
DimPlot(invivo_DATA, reduction = "pca",dims = c(1, 2),pt.size=10,cols=cols) 
#FeaturePlot(invivo_DATA,features = c("GATA6","PDGFRA","GATA3","POU5F1"),pt.size=5)
ggsave("Pre_Post_alllineages.pdf",width =14, height = 16)

#preimplantation subset and visualize
unique(Idents(invivo_DATA))[grep("CS3|cs2",unique(Idents(invivo_DATA)))]
preimp<-subset(x=invivo_DATA,idents=unique(Idents(invivo_DATA))[grep("CS3|CS2",unique(Idents(invivo_DATA)))])
preimp <- subset(preimp, subset = nFeature_RNA > 0)
preimp <- NormalizeData(preimp, verbose = FALSE)
preimp <- FindVariableFeatures(preimp, selection.method = "vst", nfeatures = 10000)
preimp <- ScaleData(preimp, verbose = FALSE)
preimp <- RunPCA(preimp, npcs = 20, verbose = FALSE)
preimp <- RunUMAP(preimp, reduction = "pca", dims = 1:20)
preimp <- RunTSNE(preimp, reduction = "pca", dims = 1:20)
DimPlot(preimp, reduction = "pca",dims = c(1, 2),pt.size=4,cols=cols) + NoLegend()
DimPlot(preimp, reduction = "pca",dims = c(1, 2),pt.size=4,cols=cols) + NoLegend()
ggsave("Prealllineages_PCA.pdf",width =8, height = 8)


#Without cleavage stages
preimp_nocleave<-subset(x=invivo_DATA,idents=unique(Idents(invivo_DATA))[grep("CS3",unique(Idents(invivo_DATA)))])
preimp_nocleave <- subset(preimp_nocleave, subset = nFeature_RNA > 0)
preimp_nocleave <- NormalizeData(preimp_nocleave, verbose = FALSE)
preimp_nocleave <- FindVariableFeatures(preimp_nocleave, selection.method = "vst", nfeatures = 10000)
preimp_nocleave <- ScaleData(preimp_nocleave, verbose = FALSE)
preimp_nocleave <- RunPCA(preimp_nocleave, npcs = 20, verbose = FALSE)
preimp_nocleave <- RunUMAP(preimp_nocleave, reduction = "pca", dims = 1:20)
preimp_nocleave <- RunTSNE(preimp_nocleave, reduction = "pca", dims = 1:20)
DimPlot(preimp_nocleave, reduction = "pca",dims = c(1, 2),pt.size=4,cols=cols) + NoLegend()
ggsave("Prenocleave_PCA.pdf",width =8, height = 8)

#order lineages
new_order<-levels(preimp_nocleave)[c(1,4,2,3,5)]
preimp_nocleave@active.ident <- factor(x = preimp_nocleave@active.ident, levels = new_order)

#Find lineage specific genes
pre_ALL<-FindAllMarkers(preimp_nocleave)
write.csv(pre_ALL, file="C:/Users/Dylan/Documents/cambridge/THESIS/rnaseq/DE_PRE.csv")
top50<-c()
for (i in 1:length(unique(pre_ALL$cluster))){
  topgenes<-pre_ALL[which(pre_ALL$cluster == unique(pre_ALL$cluster)[i]),]
  write.csv(topgenes, file=paste("C:/Users/Dylan/Documents/cambridge/THESIS/rnaseq/top_",unique(pre_ALL$cluster)[i],".csv"),sep = "")
  top50<-c(top50,topgenes[c(1:100),7])
}

#visualize top 50 unbiased lineage markers for each lineage
DoHeatmap(object = preimp_nocleave,features=top50,group.by='ident') +
  scale_fill_gradientn(colors = c("#00B0F0", "white", "#FF0000"))
ggsave("preimp_heatmap.pdf",width =10, height = 12)

#example code for gene specific visualization
FeaturePlot(preimp_nocleave,features = c("JAM2","FABP3","JAM2","MUC16","NAV2","FXYD3"),pt.size=2)
ggsave("C:/Users/Dylan/Documents/cambridge/THESIS/rnaseq/GOI_pre_troph_umap.pdf",width =7, height = 7)
VlnPlot(object=preimp_nocleave,assay="RNA",features=c("JAM2","FABP3","MIOX","GATA2","GATA3","FXYD3"),pt.size=1.5,cols=cols)
ggsave("GOI_pre_troph_VIOLIN.pdf",width =13, height = 14)

#POSTIMPLANTATION

#remove maternal and cleavage
prepost_withmat<-subset(x=invivo_DATA,idents=unique(Idents(invivo_DATA))[grep("CS3|CS5|CS6|CS7",unique(Idents(invivo_DATA)))])
prepost_without<-subset(x=prepost_withmat,idents=unique(Idents(prepost_withmat))[grep("Gland|Stroma|Myo|Zy|cMor",unique(Idents(prepost_withmat)),invert=TRUE)])
prepost_without <- NormalizeData(prepost_without, verbose = FALSE)
prepost_without <- FindVariableFeatures(prepost_without, selection.method = "vst", nfeatures = 10000)
prepost_without <- ScaleData(prepost_without, verbose = FALSE)
prepost_without <- RunPCA(prepost_without, npcs = 20, verbose = FALSE)
prepost_without <- RunUMAP(prepost_without, reduction = "pca", dims = 1:20)
prepost_without <- RunTSNE(prepost_without, reduction = "pca", dims = 1:20)
DimPlot(prepost_without, reduction = "pca",dims = c(1, 2),pt.size=4,cols=cols) + NoLegend()
ggsave(paste(path_data,"Prepost_PCA.pdf",sep=""),width =8, height = 8)

#GENERAL PRE-POST
#Reaanotate pre and postimplantation trophoblast lineages
groupedprepost<-prepost_without
groupedprepost_idents<-as.character(Idents(groupedprepost))
groupedprepost_idents[-grep("CS3",Idents(prepost_withmat))]<-as.character("POST")
groupedprepost_idents[grep("CS3",Idents(prepost_withmat))]<-as.character("PRE")
Idents(groupedprepost)<-groupedprepost_idents

#find markers between pre and postimplantation trophoblast
PRE_POST<-FindMarkers(groupedprepost,"POST","PRE",logfc.threshold = 1, test.use = "wilcox", min.pct = 0.90)
PRE_POST_POST<-PRE_POST[order(PRE_POST$avg_log2FC, decreasing=TRUE),]
PRE_POST_PRE<-PRE_POST[order(PRE_POST$avg_log2FC, decreasing=FALSE),]
top_markers<-c(rownames(PRE_POST_PRE[c(1:50),]),rownames(PRE_POST_POST[c(1:50),]))

#set lineage order
lineage_orderall<-c("cMor_CS3",
                    "ICM_CS3",
                    "Epi_CS3",
                    "Tb_CS3",
                    "Hyp_CS3",
                    "EmDisc_CS5",
                    "EmDisc_CS6",
                    "EmDisc_CS7",
                    "PGC_CS7",
                    "Am_CS5",
                    "Am_CS6",
                    "Am_CS7",
                    "VE_CS5",
                    "VE_CS6",
                    "VE_CS7",
                    "SYS_CS5",
                    "SYS_CS6",
                    "SYS_CS7",
                    "ExMes_CS5",
                    "ExMes_CS6",
                    "ExMes_CS7",
                    "Stalk_CS6",
                    "Tb_CS5",
                    "Tb_CS6",
                    "Tb_CS7",
                    "Tb_abembryonal_CS7")

new_order<-levels(prepost_withmat)[match(lineage_orderall,(levels(prepost_withmat)))]
prepost_withmat@active.ident <- factor(x = prepost_withmat@active.ident, levels = new_order)
DoHeatmap(object = prepost_withmat ,features=top_markers,group.by='ident',group.colors=cols) +
  scale_fill_gradientn(colors = c("#00B0F0", "white", "#FF2421"))
ggsave("PRE_POST_heatmap.pdf",width =12, height = 8)

#subset trophoblast lineages and analysis
justTS=subset(x=invivo_DATA,idents=c("Tb_CS3","Tb_CS5","Tb_CS6","Tb_CS7","Tb_abembryonal_CS7"))
justTS <- NormalizeData(justTS, verbose = FALSE)
justTS <- FindVariableFeatures(justTS, selection.method = "vst", nfeatures = 10000)
justTS <- ScaleData(justTS, verbose = FALSE)
justTS <- RunPCA(justTS, npcs = 20, verbose = FALSE)
justTS <- RunUMAP(justTS, reduction = "pca", dims = 1:20)
justTS <- RunTSNE(justTS, reduction = "pca", dims = 1:20)

cols_tb <-  c(
  ## preimplantation
  "Tb_CS3" = "#BF0489", 
  "Tb_CS5"="#b832c2",
  "Tb_CS6"="#802487",
  "Tb_CS7"="#49104d",
  "Tb_abembryonal_CS7" = "#150b33"
)
DimPlot(justTS, reduction = "pca",dims = c(1, 2),pt.size=5,cols=cols_tb)

#label all embryonal torphoblast from different stages
justTS_postmerged<-justTS
justts_idents<-as.character(Idents(justTS_postmerged))
justts_idents[which(justts_idents=="Tb_CS5"|justts_idents=="Tb_CS6"|justts_idents=="Tb_CS7")]<-as.character("Tb_embryonal")
Idents(justTS_postmerged)<-justts_idents


#DEGS and visualize
new_order<-levels(justTS)[c(2,3,4,5,1)]
justTS@active.ident <- factor(x = justTS@active.ident, levels = new_order)
FeaturePlot(justTS,features = c("MIOX","FABP3","JAM2","MLLT1","CGB3","PRR9"),pt.size=2)
VlnPlot(object=justTS,assay="RNA",features=c("POU5F1","CDX2","DAB2","TFAP2C","GATA3","CGA","CGB3","NR2F2","GCM1"),pt.size=1.5,cols=cols)

## Annotation for twin luminal and embryonal trophoblast
crosspsecies_annotations <- read.csv('crossspecies_annotation.csv', head=TRUE,stringsAsFactors = FALSE)
crosspsecies_annotations_marmoset <- crosspsecies_annotations[which(crosspsecies_annotations$Dataset=='Marmoset'),]
justTS_postmerged<-justTS
justts_idents<-as.character(Idents(justTS_postmerged))
justts_idents[which(justts_idents=="Tb_CS5"|justts_idents=="Tb_CS6"|justts_idents=="Tb_CS7")]<-as.character("Tb_embryonal")
Idents(justTS_postmerged)<-justts_idents
newts_idents<-as.character(Idents(justTS_postmerged))
newts_idents[match(crosspsecies_annotations_marmoset$Cell[which(crosspsecies_annotations_marmoset$Primary.cluster.anotation=="CTB_twin")],colnames(justTS_postmerged))]<-as.character(crosspsecies_annotations_marmoset$Primary.cluster.anotation)[which(crosspsecies_annotations_marmoset$Primary.cluster.anotation=="CTB_twin")]
Idents(prepost_without)<-justts_idents
Idents(justTS_postmerged)<-newts_idents

#find deferentially expressed genes
TE_CTB<-FindMarkers(prepost_without,"Tb_embryonal","Tb_CS3",logfc.threshold = 0, test.use = "wilcox", min.pct = 0.50)
TE_CTBtwin<-FindMarkers(prepost_without,"CTB_twin","Tb_CS3",logfc.threshold = 0, test.use = "wilcox", min.pct = 0.50)
CTB_SYNtwin<-FindMarkers(prepost_without,"Tb_embryonal","CTB_twin",logfc.threshold = 0, test.use = "wilcox", min.pct = 0.50)
Am_CTB<-FindMarkers(prepost_without,"Am_CS6","Tb_embryonal",logfc.threshold = 0, test.use = "wilcox", min.pct = 0.50)
Am_CTBtwin<-FindMarkers(prepost_without,"Am_CS6","CTB_twin",logfc.threshold = 0, test.use = "wilcox", min.pct = 0.50)

#AMNION analysis

#colours
cols_tbam <-  c(
  ## preimplantation
  "Tb_CS3" = "#BF0489", 
  "Tb_luminal"="#824087",
  "Tb_embryonal"="#b832c2",
  "Tb_abembryonal_CS7" = "#150b33",
  "Am_CS5" = "#5F54C7",
  "Am_CS6" = "#292175",
  "Am_CS7" = "#0b082e"
)


##amnion analysis and comparison to trophoblast
justam=subset(x=invivo_DATA,idents=c("Am_CS5","Am_CS6","Am_CS7"))
amnion_tb<-merge(justam,y=c(justTS_postmerged))
amnion_tb <- subset(amnion_tb, subset = nFeature_RNA > 0)
amnion_tb <- NormalizeData(amnion_tb, verbose = FALSE)
amnion_tb <- FindVariableFeatures(amnion_tb, selection.method = "vst", nfeatures = 10000)
amnion_tb <- ScaleData(amnion_tb, verbose = FALSE)
amnion_tb <- RunPCA(amnion_tb, npcs = 20, verbose = FALSE)
amnion_tb <- RunUMAP(amnion_tb, reduction = "pca", dims = 1:20)
amnion_tb <- RunTSNE(amnion_tb, reduction = "pca", dims = 1:20)
DimPlot(amnion_tb , reduction = "umap",dims = c(1, 2),pt.size=10, cols_tbam)  

amnion_tb_grouped<-amnion_tb
amnion_tb_grouped_idents<-as.character(Idents(amnion_tb_grouped))
amnion_tb_grouped_idents[which(amnion_tb_grouped_idents=="Am_CS5"|amnion_tb_grouped_idents=="Am_CS6"|amnion_tb_grouped_idents=="Am_CS7")]<-as.character("Am_post")
Idents(amnion_tb_grouped)<-amnion_tb_grouped_idents

#find differentially expressed genes between amnion and trophoblast lineages
Am_ALL<-FindAllMarkers(amnion_tb_grouped)
celltypes<-c("Tb_CS3","Tb_embryonal","Tb_luminal","Tb_abembryonal_CS7","Am_post")
heatmap_invivo_genes<-c()
for (i in 1:length(celltypes)){
  temp_cell<-celltypes[i]
  temp_data<-Am_ALL[which(Am_ALL$cluster==temp_cell),]
  temp_data<-temp_data[temp_data$p_val_adj<0.05,]
  temp_data<-temp_data[order(-temp_data$avg_log2FC),]
  heatmap_invivo_genes<-c(heatmap_invivo_genes,temp_data$gene[1:50])
  
}
#reorder
new_order<-levels(amnion_tb)[c(5,6,7,4,1,2,3)]
amnion_tb@active.ident <- factor(x = amnion_tb@active.ident, levels = new_order)

#visualize top regulated genes
DoHeatmap(object = amnion_tb,features=heatmap_invivo_genes,group.by='ident',group.colors=cols_tbam) +
  scale_fill_gradientn(colors = c("#00B0F0", "white", "#FF0000"))
ggsave("Tballam_embryonal_heatmap.pdf",width =10, height = 10)

#load in vitro trophoblast data and normalize

invivoAMinVITRO_DATA<- subset(x=invivo2_DATA,idents=unique(Idents(invivo2_DATA))[grep("PAVS|newTSp3|Okaeescp5",unique(Idents(invivo2_DATA)))])
invitro_DATA<-invivoAMinVITRO_DATA
invitro_DATA$species <- "In Vitro"
invitro_DATA <- subset(invitro_DATA, subset = nFeature_RNA > 0)
invitro_DATA <- NormalizeData(invitro_DATA, verbose = FALSE)
invitro_DATA <- FindVariableFeatures(invitro_DATA, selection.method = "vst", nfeatures = 10000)
invitro_DATA <- ScaleData(invitro_DATA, verbose = FALSE)
invitro_DATA <- RunPCA(invitro_DATA, npcs = 20, verbose = FALSE)
invitro_DATA <- RunUMAP(invitro_DATA, reduction = "pca", dims = 1:20)
invitro_DATA <- RunTSNE(invitro_DATA, reduction = "pca", dims = 1:20)
DimPlot(invitro_DATA, reduction = "pca",dims = c(1, 2),pt.size=10) 

#merge in vitro and in vivo datasets using canonical correlation analysis
mammal.anchors <- FindIntegrationAnchors(object.list = list(invitro_DATA,prepost_without), dims = 1:20, anchor.features = 10000,k.filter = 100)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined<-ScaleData(mammal.combined, assay = 'RNA', verbose = FALSE)
DimPlot(mammal.combined, reduction = "pca",dims = c(1, 2), cols=cols ,split.by="species",  label = FALSE,label.size = 4,pt.size =3)# + NoLegend()

#example of single marker visualization
FeaturePlot(object=mammal.combined, reduction = "pca",features=c("GATA3","WNT3A","BMP4"),pt.size=4,split.by="species")
VlnPlot(object=mammal.combined,assay="RNA",features=c("GATA3","WNT3A","BMP4"),pt.size=1.5)
ggsave("invitro_invivo_PCA_merge.pdf",width =24, height = 16)

#Merge in vitro with amnion data 

invitroAMNION<-subset(x=mammal.combined,idents=c("newTSP3",
                                                      "cmTSCPAVS",
                                                      "Tb_CS3", 
                                                      "CTB_twin",
                                                      "Tb_embryonal",
                                                      "Tb_abembryonal_CS7",
                                                      "Am_CS5",
                                                      "Am_CS6",
                                                      "Am_CS7"))
#order
new_order<-levels(invitroAMNION)[c(8,5,7,9,1,2,3,4,6)]
invitroAMNION@active.ident <- factor(x = invitroAMNION@active.ident, levels = new_order)

DoHeatmap(object = invitroAMNION,features=heatmap_invivo_genes,group.by='ident',group.colors=cols) +
  scale_fill_gradientn(colors = c("#00B0F0", "white", "#FF0000"))
ggsave("Tballam_embryonal_heatmap.pdf",width =10, height = 10)

#correlation analysis
marmoset.combined.av<-matrix(,nrow=nrow(fig2_invitroAMNION[["integrated"]]@scale.data),ncol=length(unique(Idents(fig2_invitroAMNION))))
for (i in 1:length(unique(Idents(fig2_invitroAMNION)))){
  tis_temp<-unique(Idents(fig2_invitroAMNION))[i]
  if(length(which(Idents(fig2_invitroAMNION)==tis_temp))>1){
    marmoset.combined.av[,i]<-rowMeans(fig2_invitroAMNION[["integrated"]]@scale.data[,which(Idents(fig2_invitroAMNION)==tis_temp)],na.rm=TRUE)
    
  }
  else {
    marmoset.combined.av[,i]<-fig2_invitroAMNION[["integrated"]]@scale.data[,which(Idents(fig2_invitroAMNION)==tis_temp)]
  }
  
}
colnames(marmoset.combined.av)<-unique(Idents(fig2_okae))
rownames(marmoset.combined.av)<-rownames(fig2_okae[["integrated"]]@scale.data)
M<-cor(as.matrix(marmoset.combined.av))
M_2<-(M+1)/2
corrplot(M, method="color", order="hclust",outline='white',
         col=rev(c(brewer.pal(n=8, name="RdYlBu"),rep("white",2))))

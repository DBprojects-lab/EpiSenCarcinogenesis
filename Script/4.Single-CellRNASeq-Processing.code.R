library(Seurat)
library(dplyr)
library(hdf5r)
library(data.table)
#############################################################
######Tissue single-cell data processing--Adrenal_gland######
#############################################################
###https://figshare.com/articles/dataset/HCL_DGE_Data/7235471#Download data
setwd("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Adrenal")
count1 <- read.table("AdultAdrenalGland2.rmbatchdge.txt.gz", header = TRUE)
anndata1 <- read.csv("Adult-Adrenal-Gland2_rmbatchAnno.csv",header = T)
count2 <- read.table("AdultAdrenalGland3.rmbatchdge.txt.gz", header = TRUE)
anndata2 <- read.csv("Adult-Adrenal-Gland3_rmbatchAnno.csv",header = T)
#Combined sample
inter_gene<- intersect(rownames(count1),rownames(count2))
count1 <- count1[inter_gene,]
count2 <- count2[inter_gene,]
allcount <- cbind(count1,count2)
allanndata <- rbind(anndata1,anndata2)
allanndata <- allanndata[,c(4,14)]
head(colnames(allcount))
table(allanndata$Celltype)
allanndata$Celltype <- apply(as.matrix(allanndata$Celltype),1,function(x){
  x <- strsplit(x,"_")[[1]][1]
  return(x)})
table(allanndata$Celltype)
allanndata$Celltype[which(allanndata$Celltype == "lymphatic endothelial cell")] <- "Endothelial cell"
allanndata$Celltype[which(allanndata$Celltype == "CD8 T cell")] <- "T cell"
allanndata$Celltype[which(allanndata$Celltype == "B cell (Plasmocyte)")] <- "Plasmocyte"
allcount <- allcount[,allanndata$Cell_id]
save(allcount,file = "Adrenal_count.rda")
save(allanndata,file = "Adrenal_allanndata.rda")

###seruat-
sce.all <- CreateSeuratObject(counts = allcount)
sce.all@meta.data$Celltype  <- allanndata$Celltype
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sce.all=subset(sce.all,subset = nFeature_RNA > 200 & nFeature_RNA <2000  & percent.mt<20)
sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize") 
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst")
all.genes <- rownames(sce.all)
sce.all <- ScaleData(sce.all, features = all.genes)
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
ElbowPlot(sce.all,ndims = 50)
sce.all<- RunUMAP(sce.all, reduction = "pca", dims = 1:20)
sce.all<- FindNeighbors(sce.all, reduction = "pca", dims = 1:20)
sce.all<- FindClusters(sce.all, resolution = 0.7)
pdf("FirstClustering.pdf",width=9)
DimPlot(sce.all,reduction="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
diff.wilcox = FindAllMarkers(sce.all, logfc.threshold = 0.25, min.pct = 0.1, 
                             only.pos = TRUE, test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(sce.all))
##Cell type annotations
pdf("topMarkerGeneFeature.pdf",width=12,height=20)
DoHeatmap(sce.all, features = top10, group.by = "seurat_clusters", group.bar = T, size = 4)
dev.off()
#B cell 14
#Endothelial cell 5,13
#Fibroblast 15
#Smooth muscle cell 16
#T cell 6
#NK cell 8，12
#Macrophage 9,10,11
#Zona fasciculata cell  0,1,2,3,7
sce.all <- subset(sce.all,subset = seurat_clusters!=1)
sce.all <- subset(sce.all,subset = seurat_clusters!=4)
sce.all <- RenameIdents(sce.all,
                        `0` = "Zona fasciculata cell",`2` = "Zona fasciculata cell",`3` = "Zona fasciculata cell",`5` = "Endothelial cell",`6` = "T cell",`7` = "Zona fasciculata cell", #
                        `8` = "NK cell",`9` = "Macrophage",`10` = "Macrophage",`11` = "Macrophage",`12` = "NK cell",`13` = "Endothelial cell",`14` = "B cell",
                        `15` = "Fibroblast",`16` = "Smooth muscle cell")
sce.all$cellType=Idents(sce.all)
save(sce.all,file="Adrenal_seruat.rda")


#############################################################
######Tissue single-cell data processing--Colon######
#############################################################
###https://figshare.com/articles/dataset/HCL_DGE_Data/7235471#Download data
setwd("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Colon")
count1 <- read.table("AdultAscendingColon1.rmbatchdge.txt.gz", header = TRUE)
anndata1 <- read.csv("Adult-Ascending-Colon1_rmbatchAnno.csv",header = T)
count2 <- read.table("AdultSigmoidColon1.rmbatchdge.txt.gz", header = TRUE)
anndata2 <- read.csv("Adult-Sigmoid-Colon1_rmbatchAnno.csv",header = T)
count3 <- read.table("AdultTransverseColon2.rmbatchdge.txt.gz", header = TRUE)
anndata3 <- read.csv("Adult-Transverse-Colon2_rmbatchAnno.csv",header = T)
inter_gene<- intersect(rownames(count1),intersect(rownames(count2),rownames(count3)))
count1 <- count1[inter_gene,]
count2 <- count2[inter_gene,]
count3 <- count3[inter_gene,]
allcount <- cbind(count1,count2,count3)
allanndata <- rbind(anndata1,anndata2,anndata3)
allanndata <- allanndata[,c(4,14)]
head(colnames(allcount))
table(allanndata$Celltype)
allanndata$Celltype <- apply(as.matrix(allanndata$Celltype),1,function(x){
  x <- strsplit(x,"_")[[1]][1]
  return(x)})
table(allanndata$Celltype)
allanndata$Celltype[which(allanndata$Celltype == "B cell(Plasmocyte)")] <- "Plasmocyte"
allanndata <- allanndata[-which(allanndata$Celltype == "Unknown"),]
allcount <- allcount[,allanndata$Cell_id]
sce.all <- CreateSeuratObject(counts = allcount)
sce.all@meta.data$Celltype  <- allanndata$Celltype
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sce.all=subset(sce.all,subset = nFeature_RNA > 200 & nFeature_RNA <2000  & percent.mt<20)
sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize") 
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst")
all.genes <- rownames(sce.all)
sce.all <- ScaleData(sce.all, features = all.genes)
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
ElbowPlot(sce.all,ndims = 50)
sce.all<- RunUMAP(sce.all, reduction = "pca", dims = 1:20)
sce.all<- FindNeighbors(sce.all, reduction = "pca", dims = 1:20)
sce.all<- FindClusters(sce.all, resolution = 0.7)
pdf("FirstClustering.pdf",width=9)
DimPlot(sce.all,reduction="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
diff.wilcox = FindAllMarkers(sce.all, logfc.threshold = 0.25, min.pct = 0.1, 
                             only.pos = TRUE, test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(sce.all)) 
pdf("topMarkerGeneFeature.pdf",width=12,height=20)
DoHeatmap(sce.all, features = top10, group.by = "seurat_clusters", group.bar = T, size = 4)
dev.off()
#Epithelial cell 1，2，3，6，7，8，11
#B cell 0，4
#Macrophage 9
#Endothelial cell 17
#Smooth muscle cell 15
#Fibroblast 10，14
#T cell 5，18
#Mast cell 12
#Dendritic cell 13
#Enteric nerval cell 16
sce.all <- RenameIdents(sce.all,`0` = "B cell",`1` = "Epithelial cell",`2` = "Epithelial cell",`3` = "Epithelial cell",`4` = "B cell",
                        `5` = "T cell",`6` = "Epithelial cell",`7` = "Epithelial cell",`8` = "Epithelial cell",`9` = "Macrophage",
                        `10` = "Fibroblast",`11` = "Epithelial cell",`12` = "Mast cell",`13` = "Dendritic cell",`14` = "Fibroblast",
                        `15` = "Smooth muscle cell",`16` = "Enteric nerval cell",`17` = "Endothelial cell",`18` = "T cell"
)
sce.all$cellType=Idents(sce.all)
save(sce.all,file="Colon_seruat.rda")

#############################################################
######Tissue single-cell data processing--Esophagus######
#############################################################
###https://figshare.com/articles/dataset/HCL_DGE_Data/7235471#Download data
setwd("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Esophagus")
count1 <- read.table("AdultEsophagus1.rmbatchdge.txt.gz", header = TRUE)
anndata1 <- read.csv("Adult-Esophagus2_rmbatchAnno.csv",header = T)
count2 <- read.table("AdultEsophagus2.rmbatchdge.txt.gz", header = TRUE)
anndata2 <- read.csv("Adult-Esophagus2_rmbatchAnno.csv",header = T)
inter_gene<- intersect(rownames(count1),rownames(count2))
count1 <- count1[inter_gene,]
count2 <- count2[inter_gene,]
allcount <- cbind(count1,count2)
allanndata <- rbind(anndata1,anndata2)
allanndata <- allanndata[,c(4,14)]
head(colnames(allcount))
table(allanndata$Celltype)
allanndata$Celltype <- apply(as.matrix(allanndata$Celltype),1,function(x){
  x <- strsplit(x,"_")[[1]][1]
  return(x)})
table(allanndata$Celltype)
allanndata$Celltype[which(allanndata$Celltype == "Mucosal aquamous Epithelial cell")] <- "Epithelial cell"
allanndata$Celltype[which(allanndata$Celltype == "Neutrophil ")] <- "Neutrophil"
allcount <- allcount[,allanndata$Cell_id]
save(allcount,file = "Esophagus_count.rda")
save(allanndata,file = "Esophagus_allanndata.rda")
sce.all <- CreateSeuratObject(counts = allcount)
sce.all@meta.data$Celltype  <- allanndata$Celltype
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sce.all=subset(sce.all,subset = nFeature_RNA > 200 & nFeature_RNA <2000  & percent.mt<20)
sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize") 
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst")
all.genes <- rownames(sce.all)
sce.all <- ScaleData(sce.all, features = all.genes)
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
ElbowPlot(sce.all,ndims = 50)
sce.all<- RunUMAP(sce.all, reduction = "pca", dims = 1:20)
sce.all<- FindNeighbors(sce.all, reduction = "pca", dims = 1:20)
sce.all<- FindClusters(sce.all, resolution = 0.7)
pdf("FirstClustering.pdf",width=9)
DimPlot(sce.all,reduction="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
diff.wilcox = FindAllMarkers(sce.all, logfc.threshold = 0.25, min.pct = 0.1, 
                             only.pos = TRUE, test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(sce.all)) 
pdf("topMarkerGeneFeature.pdf",width=12,height=20)
DoHeatmap(sce.all, features = top10, group.by = "seurat_clusters", group.bar = T, size = 4)
dev.off()
#Epithelial cell 7,8,13,15，6
#Macrophage 4,9
#Neutrophil 14
#Mast cell 12
#Smooth muscle cell 18
#Endothelial cell 5
#Goblet cell 11
#Fibroblast 0，1，2，10
#B cell 3
#删除16，17
sce.all <- subset(sce.all,subset = seurat_clusters!=16)
sce.all <- subset(sce.all,subset = seurat_clusters!=17)
sce.all <- RenameIdents(sce.all,
                        `0` = "Fibroblast",`1` = "Fibroblast",`2` = "Fibroblast",`3` = "B cell",`4` = "Macrophage",`5` = "Endothelial cell",
                        `6` = "Epithelial cell",`7` = "Epithelial cell",`8` = "Epithelial cell", `9` = "Macrophage",`10` = "Fibroblast",
                        `11` = "Epithelial cell",`12` = "Mast cell",`13` = "Epithelial cell",`14` = "Neutrophil",`15` = "Epithelial cell",`18` = "Smooth muscle cell")
sce.all$cellType=Idents(sce.all)
save(sce.all,file="Esophagus_seruat.rda")



#############################################################
######Tissue single-cell data processing--kidney######
#############################################################
###https://figshare.com/articles/dataset/HCL_DGE_Data/7235471#Download data
setwd("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Kidney")
count1 <- read.table("AdultKidney2.rmbatchdge.txt.gz", header = TRUE)
anndata1 <- read.csv("Adult-Kidney2_rmbatchAnno.csv",header = T)

count2 <- read.table("AdultKidney3.rmbatchdge.txt.gz", header = TRUE)
anndata2 <- read.csv("Adult-Kidney3_rmbatchAnno.csv",header = T)

count3 <- read.table("AdultKidney4.rmbatchdge.txt.gz", header = TRUE)
anndata3 <- read.csv("Adult-Kidney4_rmbatchAnno.csv",header = T)

#合并样本
inter_gene<- intersect(rownames(count1),intersect(rownames(count2),rownames(count3)))
count1 <- count1[inter_gene,]
count2 <- count2[inter_gene,]
count3 <- count3[inter_gene,]
allcount <- cbind(count1,count2,count3)
allanndata <- rbind(anndata1,anndata2,anndata3)
allanndata <- allanndata[,c(4,14)]
head(colnames(allcount))
table(allanndata$Celltype)
allanndata$Celltype <- apply(as.matrix(allanndata$Celltype),1,function(x){
  x <- strsplit(x,"_")[[1]][1]
  return(x)})
table(allanndata$Celltype)
allanndata$Celltype[which(allanndata$Celltype %in% c("Loop of henle ","Loop of Henle (Thick ascending limb)","Loop of Henle(Thick ascending limb)"))] <- "Loop of henle"
allanndata$Celltype[which(allanndata$Celltype %in% c("Kidney Epithelial cell","Ureteric Epithelial cell","Epithelial"))] <- "Epithelial cell"
allanndata$Celltype[which(allanndata$Celltype %in% c("Fenestrated endothelial cell","Glomerular endothelial cell"))] <- "Endothelial cell"
allanndata$Celltype[which(allanndata$Celltype == "Conventional dendritic cell")] <- "Dendritic cell"
allanndata$Celltype[which(allanndata$Celltype == "B cell (Plasmocyte)")] <- "B cell(Plasmocyte)"
allanndata$Celltype[which(allanndata$Celltype == "B cell(Plasmocyte)")] <- "Plasmocyte"
allanndata <- allanndata[-which(allanndata$Celltype == "Unknown"),]
allcount <- allcount[,allanndata$Cell_id]
sce.all <- CreateSeuratObject(counts = allcount)
sce.all@meta.data$Celltype  <- allanndata$Celltype
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sce.all=subset(sce.all,subset = nFeature_RNA > 200 & nFeature_RNA <2000  & percent.mt<20)
sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize") 
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst")
all.genes <- rownames(sce.all)
sce.all <- ScaleData(sce.all, features = all.genes)
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
ElbowPlot(sce.all,ndims = 50)
sce.all<- RunUMAP(sce.all, reduction = "pca", dims = 1:20)
sce.all<- FindNeighbors(sce.all, reduction = "pca", dims = 1:20)
sce.all<- FindClusters(sce.all, resolution = 0.7)
pdf("FirstClustering.pdf",width=9)
DimPlot(sce.all,reduction="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
diff.wilcox = FindAllMarkers(sce.all, logfc.threshold = 0.25, min.pct = 0.1, 
                             only.pos = TRUE, test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(sce.all)) 
pdf("topMarkerGeneFeature.pdf",width=12,height=20)
DoHeatmap(sce.all, features = top10, group.by = "seurat_clusters", group.bar = T, size = 4)
dev.off()
#Intercalated cell 1,3,20
#Macrophage 2
#Neutrophil 13
#T cell 11
#Endothelial cell 7,12,15,18
#Plasmocyte 21
#Smooth muscle cell 19
#Proximal tubule cell 4,5,10,14
#Loop of henle 0,9,16
#Epithelial cell 6
#Principal cell 8
#Distal tubule cell 17
sce.all <- RenameIdents(sce.all,
                        `0` = "Loop of henle",`1` = "Intercalated cell",`2` = "Macrophage",`3` = "Intercalated cell",`4` = "Proximal tubule cell",`5` = "Proximal tubule cell",#
                        `6` = "Epithelial cell",`7` = "Endothelial cell",`8` = "Principal cell",`9` = "Loop of henle",`10` = "Proximal tubule cell",`11` = "T cell",
                        `12` = "Endothelial cell",`13` = "Neutrophil",`14` = "Proximal tubule cell",`15` = "Endothelial cell",`16` = "Loop of henle",`17` = "Distal tubule cell",
                        `18` = "Endothelial cell",`19` = "Smooth muscle cell",`20` = "Intercalated cell",`21` = "Plasmocyte")
pdf("FirstClustering_new.pdf",width=9)
DimPlot(sce.all,reduction="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
sce.all$cellType=Idents(sce.all)
save(sce.all,file="Kidney_seruat.rda")



#############################################################
######Tissue single-cell data processing--Liver######
#############################################################
###https://figshare.com/articles/dataset/HCL_DGE_Data/7235471#Download data
setwd("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Liver")
count1 <- read.table("AdultLiver1.rmbatchdge.txt.gz", header = TRUE)
anndata1 <- read.csv("Adult-Liver1_rmbatchAnno.csv",header = T)
count2 <- read.table("AdultLiver2.rmbatchdge.txt.gz", header = TRUE)
anndata2 <- read.csv("Adult-Liver2_rmbatchAnno.csv",header = T)
count3 <- read.table("AdultLiver4.rmbatchdge.txt.gz", header = TRUE)
anndata3 <- read.csv("Adult-Liver4_rmbatchAnno.csv",header = T)
inter_gene<- intersect(rownames(count1),intersect(rownames(count2),rownames(count3)))
count1 <- count1[inter_gene,]
count2 <- count2[inter_gene,]
count3 <- count3[inter_gene,]
allcount <- cbind(count1,count2,count3)
allanndata <- rbind(anndata1,anndata2,anndata3)
allanndata <- allanndata[,c(4,14)]
head(colnames(allcount))
table(allanndata$Celltype)
allanndata$Celltype <- apply(as.matrix(allanndata$Celltype),1,function(x){
  x <- strsplit(x,"_")[[1]][1]
  return(x)})
table(allanndata$Celltype)
allanndata$Celltype[which(allanndata$Celltype == "Conventional dendritic cell")] <- "Dendritic cell"
allanndata$Celltype[which(allanndata$Celltype == "Motile liver macrophage")] <- "Macrophage"
allanndata$Celltype[which(allanndata$Celltype == "Vascular endothelial cell")] <- "Endothelial cell"
allanndata$Celltype[which(allanndata$Celltype == "Kuppfer Cell")] <- "Kuppfer cell"
allanndata$Celltype[which(allanndata$Celltype == "Activated T cell")] <- "T cell"
allanndata$Celltype[which(allanndata$Celltype == "B cell (Plasmocyte)")] <- "Plasmocyte"
allcount <- allcount[,allanndata$Cell_id]
sce.all <- CreateSeuratObject(counts = allcount)
sce.all@meta.data$Celltype  <- allanndata$Celltype
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sce.all=subset(sce.all,subset = nFeature_RNA > 200 & nFeature_RNA <2000  & percent.mt<20)
sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize") 
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst")
all.genes <- rownames(sce.all)
sce.all <- ScaleData(sce.all, features = all.genes)
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
ElbowPlot(sce.all,ndims = 50)
sce.all<- RunUMAP(sce.all, reduction = "pca", dims = 1:20)
sce.all<- FindNeighbors(sce.all, reduction = "pca", dims = 1:20)
sce.all<- FindClusters(sce.all, resolution = 0.7)
pdf("FirstClustering.pdf",width=9)
DimPlot(sce.all,reduction="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
diff.wilcox = FindAllMarkers(sce.all, logfc.threshold = 0.25, min.pct = 0.1, 
                             only.pos = TRUE, test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(sce.all)) 
pdf("topMarkerGeneFeature.pdf",width=12,height=20)
DoHeatmap(sce.all, features = top10, group.by = "seurat_clusters", group.bar = T, size = 4)
dev.off()
#Endothelial cell 2,11,10
#Macrophage 0,16
#Epithelial cell 8
#Dendritic cell 3，6
#Plasmocyte 7
#T cell 1
#Monocyte 4，5，13
#Neutrophil 12
#Smooth muscle cell 17
#Hepatocyte 14,15
#Mast cell 9
sce.all <- subset(sce.all,subset = seurat_clusters!=17)
sce.all <- RenameIdents(sce.all,
                        `0` = "Macrophage",`1` = "T cell",`2` = "Endothelial cell",`3` = "Dendritic cell",`4` = "Monocyte",`5` = "Monocyte",`6` = "Dendritic cell", #
                        `7` = "Plasmocyte",`8` = "Epithelial cell",`9` = "Mast cell",`10` = "Endothelial cell",`11` = "Endothelial cell",`12` = "Neutrophil",
                        `13` = "Monocyte",`14` = "Hepatocyte",`15` = "Hepatocyte",`16` = "Macrophage")
sce.all$cellType=Idents(sce.all)
save(sce.all,file="Liver_seruat.rda")



#############################################################
######Tissue single-cell data processing--Lung######
#############################################################
###https://figshare.com/articles/dataset/HCL_DGE_Data/7235471#Download data
setwd("/boot3/cjl/Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Lung")
count1 <- read.table("AdultLung1.rmbatchdge.txt.gz", header = TRUE)
anndata1 <- read.csv("Adult-Lung1_rmbatchAnno.csv",header = T)
count2 <- read.table("AdultLung2.rmbatchdge.txt.gz", header = TRUE)
anndata2 <- read.csv("Adult-Lung2_rmbatchAnno.csv",header = T)
count3 <- read.table("AdultLung3.rmbatchdge.txt.gz", header = TRUE)
anndata3 <- read.csv("Adult-Lung3_rmbatchAnno.csv",header = T)
inter_gene<- intersect(rownames(count1),intersect(rownames(count2),rownames(count3)))
count1 <- count1[inter_gene,]
count2 <- count2[inter_gene,]
count3 <- count3[inter_gene,]
allcount <- cbind(count1,count2,count3)
allanndata <- rbind(anndata1,anndata2,anndata3)
allanndata <- allanndata[,c(4,14)]
head(colnames(allcount))
table(allanndata$Celltype)
allanndata$Celltype <- apply(as.matrix(allanndata$Celltype),1,function(x){
  x <- strsplit(x,"_")[[1]][1]
  return(x)})
table(allanndata$Celltype)
allanndata$Celltype[which(allanndata$Celltype %in% c("Arterial endothelial cell","Artry endothelial cell","Lymphatic endothelial cell"))] <- "Endothelial cell"
allanndata$Celltype[which(allanndata$Celltype %in% c("Kidney Epithelial cell","Ureteric Epithelial cell","Epithelial","Basal/Epithelial cell"))] <- "Epithelial cell"
allanndata$Celltype[which(allanndata$Celltype %in% c("Fenestrated endothelial cell","Glomerular endothelial cell"))] <- "Endothelial cell"
allanndata$Celltype[which(allanndata$Celltype == "Conventional dendritic cell")] <- "Dendritic cell"
allanndata$Celltype[which(allanndata$Celltype == "AT1 cell ")] <- "AT1 cell"
allanndata$Celltype[which(allanndata$Celltype == "Proliferating alveolar bipotent progenitor cell")] <- "Proliferating cell"
allanndata$Celltype[which(allanndata$Celltype == "M2 macrophage ")] <- "Macrophage"
allanndata$Celltype[which(allanndata$Celltype == "Alveolar bipotent/intermediate cell")] <- "Alveolar bipotent cell"
allanndata$Celltype[which(allanndata$Celltype == "Actived T cell")] <- "T cell"
allanndata$Celltype[which(allanndata$Celltype == "B cell (Plasmocyte)")] <- "Plasmocyte"
allcount <- allcount[,allanndata$Cell_id]
sce.all <- CreateSeuratObject(counts = allcount)
sce.all@meta.data$Celltype  <- allanndata$Celltype
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sce.all=subset(sce.all,subset = nFeature_RNA > 200 & nFeature_RNA <2000  & percent.mt<15)
sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize") 
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst")
all.genes <- rownames(sce.all)
sce.all <- ScaleData(sce.all, features = all.genes)
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
ElbowPlot(sce.all,ndims = 50)
sce.all<- RunUMAP(sce.all, reduction = "pca", dims = 1:20)
sce.all<- FindNeighbors(sce.all, reduction = "pca", dims = 1:20)
sce.all<- FindClusters(sce.all, resolution = 0.7)
pdf("FirstClustering.pdf",width=9)
DimPlot(sce.all,reduction="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
diff.wilcox = FindAllMarkers(sce.all, logfc.threshold = 0.25, min.pct = 0.1, 
                             only.pos = TRUE, test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(sce.all)) 
pdf("topMarkerGeneFeature.pdf",width=12,height=20)
DoHeatmap(sce.all, features = top10, group.by = "seurat_clusters", group.bar = T, size = 4)
dev.off()
#Macrophage 1，3，17，20
#Dendritic cel 6
#Neutrophil 7
#B cell 14
#T cell 8
#Mast cell 11
#Fibroblast 9,10
#Smooth muscle cell 13
#Endothelial cell 4,5
#Epithelial cell 0,2,12,15,16,18,19
sce.all <- RenameIdents(sce.all,
                        `0` = "Epithelial cell",`1` = "Macrophage",`2` = "Epithelial cell",`3` = "Macrophage",`4` = "Endothelial cell",`5` = "Endothelial cell",`6` = "Dendritic cell",
                        `7` = "Neutrophil",`8` = "T cell",`9` = "Fibroblast",`10` = "Fibroblast",`11` = "Mast cell",`12` = "Epithelial cell",`13` = "Smooth muscle cell",
                        `14` = "B cell",`15` = "Epithelial cell",`16` = "Epithelial cell",`17` = "Macrophage",`18` = "Epithelial cell",`19` = "Epithelial cell",`20` = "Macrophage")
sce.all$cellType=Idents(sce.all)
save(sce.all,file="Lung_seruat.rda")


#############################################################
######Tissue single-cell data processing--Pancreas######
#############################################################
###https://figshare.com/articles/dataset/HCL_DGE_Data/7235471#Download data
setwd("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Pancreas")
allcount <- read.table("AdultPancreas1.rmbatchdge.txt.gz", header = TRUE)
allanndata <- read.csv("Adult-Pancreas1_rmbatchAnno.csv",header = T)
allanndata <- allanndata[,c(4,14)]
head(colnames(allcount))
table(allanndata$Celltype)
allanndata$Celltype <- apply(as.matrix(allanndata$Celltype),1,function(x){
  x <- strsplit(x,"_")[[1]][1]
  return(x)})
table(allanndata$Celltype)
allanndata$Celltype[which(allanndata$Celltype == "Acniar cell")] <- "Acinar cell"
allanndata$Celltype[which(allanndata$Celltype == "M2 Macrophage")] <- "Macrophage"
allcount <- allcount[,allanndata$Cell_id]
sce.all <- CreateSeuratObject(counts = allcount)
sce.all@meta.data$Celltype  <- allanndata$Celltype
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sce.all=subset(sce.all,subset = nFeature_RNA > 200 & nFeature_RNA <2000  & percent.mt<15)
sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize") 
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst")
all.genes <- rownames(sce.all)
sce.all <- ScaleData(sce.all, features = all.genes)
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
ElbowPlot(sce.all,ndims = 50)
sce.all<- RunUMAP(sce.all, reduction = "pca", dims = 1:20)
sce.all<- FindNeighbors(sce.all, reduction = "pca", dims = 1:20)
sce.all<- FindClusters(sce.all, resolution = 0.7)
pdf("FirstClustering.pdf",width=9)
DimPlot(sce.all,reduction="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
diff.wilcox = FindAllMarkers(sce.all, logfc.threshold = 0.25, min.pct = 0.1, 
                             only.pos = TRUE, test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(sce.all)) 
pdf("topMarkerGeneFeature.pdf",width=12,height=20)
DoHeatmap(sce.all, features = top10, group.by = "seurat_clusters", group.bar = T, size = 4)
dev.off()
#Macrophage 8
#Endothelial cell 10
#Fibroblast 11
#	Alpha cell 9
#Eipthelial cell 4
#Smooth muscle cell 6
#Acinar cell 1,2,3
#Exocrine cell 0，5，7
sce.all <- RenameIdents(sce.all,
                        `0` = "Exocrine cell",`1` = "Acinar cell",`2` = "Acinar cell",`3` = "Acinar cell",`4` = "Epithelial cell",`5` = "Exocrine cell",#
                        `6` = "Smooth muscle cell",`7` = "Exocrine cell",`8` = "Macrophage",`9` = "Alpha cell",`10` = "Endothelial cell",`11` = "Fibroblast")

sce.all$cellType=Idents(sce.all)
save(sce.all,file="Pancreas_seruat.rda")




#############################################################
######Tissue single-cell data processing--Prostate######
#############################################################
###https://figshare.com/articles/dataset/HCL_DGE_Data/7235471#Download data
setwd("/boot3/cjl/Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Prostate")
count <- read.table("AdultProstate1.rmbatchdge.txt.gz", header = TRUE)
anndata <- read.csv("Adult-Prostate1_rmbatchAnno.csv",header = T)
anndata <- anndata[,c(4,14)]
head(colnames(count))
table(anndata$Celltype)
anndata$Celltype <- apply(as.matrix(anndata$Celltype),1,function(x){
  x <- strsplit(x,"_")[[1]][1]
  return(x)})
table(anndata$Celltype)
#将Unknow上皮细胞变为上皮细胞
anndata$Celltype[which(anndata$Celltype == "Unknown Epithelial cell")] <- "Epithelial cell"
anndata$Celltype[which(anndata$Celltype == "Intermediate Epithelial cell")] <- "Mesothelial cell"
anndata$Celltype[which(anndata$Celltype == "M1 Macrophage")] <- "Macrophage"
#删除未知细胞
anndata <- anndata[-which(anndata$Celltype == "Unknown"),]
count <- count[,anndata$Cell_id]
identical(colnames(count),anndata$Cell_id)
sce.all <- CreateSeuratObject(counts = count)
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sce.all=subset(sce.all,subset = nFeature_RNA > 200 & nFeature_RNA <2000  & percent.mt<15)
sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize") 
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst")
all.genes <- rownames(sce.all)
sce.all <- ScaleData(sce.all, features = all.genes)
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
ElbowPlot(sce.all,ndims = 50)
sce.all<- RunUMAP(sce.all, reduction = "pca", dims = 1:15)
sce.all<- FindNeighbors(sce.all, reduction = "pca", dims = 1:15)
sce.all<- FindClusters(sce.all, resolution = 0.9)
pdf("FirstClustering.pdf",width=9)
DimPlot(sce.all,reduction="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
diff.wilcox = FindAllMarkers(sce.all, logfc.threshold = 0.25, min.pct = 0.1, 
                             only.pos = TRUE, test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(sce.all)) 
pdf("topMarkerGeneFeature.pdf",width=12,height=20)
DoHeatmap(sce.all, features = top10, group.by = "seurat_clusters", group.bar = T, size = 4)
dev.off()
#Basal cell 6
#Epithelial cell 0，1，2，3，4，5
#Macrophage 7
#Endothelial cell 8
#Fibroblast 9
#T cll 10
sce.all <- RenameIdents(sce.all,
                        `0` = "Epithelial cell",`1` = "Epithelial cell",`2` = "Epithelial cell",`3` = "Epithelial cell",`4` = "Epithelial cell",
                        `5` = "Epithelial cell",`6` = "Basal cell",`7` = "Macrophage",`8` = "Endothelial cell",`9` = "Fibroblast",`10` = "T cell")

sce.all$cellType=Idents(sce.all)
save(sce.all,file="Prostate_seruat.rda")



#############################################################
######Tissue single-cell data processing--Uteurs######
#############################################################
###https://figshare.com/articles/dataset/HCL_DGE_Data/7235471#Download data

setwd("/boot3/cjl/Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Uterus")
count <- read.table("AdultUterus1.rmbatchdge.txt.gz", header = TRUE)
anndata <- read.csv("Adult-Uterus1_rmbatchAnno.csv",header = T)
anndata <- anndata[,c(4,14)]
head(colnames(count))
table(anndata$Celltype)
anndata$Celltype <- apply(as.matrix(anndata$Celltype),1,function(x){
  x <- strsplit(x,"_")[[1]][1]
  return(x)})
table(anndata$Celltype)
anndata$Celltype[which(anndata$Celltype %in% c("Endothelial cell in EMT"))] <- "Endothelial cell"
anndata$Celltype[which(anndata$Celltype == "Luminal epithelium ")] <- "Epithelial cell"
anndata$Celltype[which(anndata$Celltype == "M1 Macrophage")] <- "Macrophage"
anndata <- anndata[match(colnames(count),anndata$Cell_id,nomatch = 0),]
sce.all <- CreateSeuratObject(counts = count)
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sce.all=subset(sce.all,subset = nFeature_RNA > 200 & nFeature_RNA <2000  & percent.mt<15)
sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize") 
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst")
all.genes <- rownames(sce.all)
# 结果存放于 `pbmc[["RNA"]]@scale.data`
sce.all <- ScaleData(sce.all, features = all.genes)
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
ElbowPlot(sce.all,ndims = 50)
sce.all<- RunUMAP(sce.all, reduction = "pca", dims = 1:15)
sce.all<- FindNeighbors(sce.all, reduction = "pca", dims = 1:15)
sce.all<- FindClusters(sce.all, resolution = 0.9)
pdf("FirstClustering.pdf",width=9)
DimPlot(sce.all,reduction="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
diff.wilcox = FindAllMarkers(sce.all, logfc.threshold = 0.25, min.pct = 0.1, 
                             only.pos = TRUE, test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(sce.all)) 
pdf("topMarkerGeneFeature.pdf",width=12,height=20)
DoHeatmap(sce.all, features = top10, group.by = "seurat_clusters", group.bar = T, size = 4)
dev.off()

#Endothelial cell 0，1，2，3，6，9
#Smooth muscle cell 4，5
#Fibroblasts 7,8,10
#Epithelial cell 11
#T cell 13
#Mast cell 14
sce.all <- RenameIdents(sce.all,
                        `0` = "Endothelial cell",`1` = "Endothelial cell",`2` = "Endothelial cell",`3` = "Endothelial cell",`4` = "Smooth muscle cell",
                        `5` = "Smooth muscle cell",`6` = "Endothelial cell",`7` = "Fibroblast",`8` = "Fibroblast",`9` = "Endothelial cell",`10` = "Fibroblast",
                        `11` = "Epithelial cell",`12` = "Endothelial cell",`13` = "T cell",`14` = "Mast cell")
sce.all$cellType=Idents(sce.all)
save(sce.all,file="Uterus_seruat.rda")




#############################################################
######Tissue single-cell data processing--Blood######
#############################################################
###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149938#Download data
setwd("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Blood/GSE149938_RAW")
count=t(fread(file.path("GSE149938_umi_matrix.csv.gz"),data.table = T))
count_colname <- count[1,]
count <- as.data.frame(count[-1,])
colnames(count) <- count_colname
dir="./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Blood/GSE149938_RAW/cell_info"
samples=list.files( dir )
samples
all_cell_info <- data_frame()
for (pro in samples) {
  onecell_info <- fread(file.path(dir,pro),data.table = F,header = F)
  onecell_info$cellAnnto <- gsub(".*?_(.*?)\\..*", "\\1", pro)
  all_cell_info <- rbind(all_cell_info,onecell_info)
}
sce.all <- CreateSeuratObject(counts = count)
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize") 
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst")
all.genes <- rownames(sce.all)
# 结果存放于 `pbmc[["RNA"]]@scale.data`
sce.all <- ScaleData(sce.all, features = all.genes)
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
ElbowPlot(sce.all,ndims = 50)
sce.all<- RunUMAP(sce.all, reduction = "pca", dims = 1:10)
sce.all<- FindNeighbors(sce.all, reduction = "pca", dims = 1:10)
sce.all<- FindClusters(sce.all, resolution = 0.5)
pdf("FirstClustering_BloodCluster.pdf",width=9)
DimPlot(sce.all,reduction="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
diff.wilcox = FindAllMarkers(sce.all, logfc.threshold = 0.25, min.pct = 0.1, 
                             only.pos = TRUE, test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(sce.all)) 
pdf("topMarkerGeneFeature.pdf",width=12,height=20)
DoHeatmap(sce.all, features = top10, group.by = "seurat_clusters", group.bar = T, size = 4)
dev.off()
#erythrocytes 0，14
#B cell 2，7，12，13
#Neutrophils 3，4，6
#Dendritic cells 5
#hematopoietic stem cell 1，8，11
#NK cell 10
#Monocytes 15
#T cell 9
sce.all <- subset(sce.all,subset = seurat_clusters!=16)
sce.all <- RenameIdents(sce.all,
                        `0` = "Erythrocyte",`1` = "Hematopoietic stem cell",`2` = "B cell",`3` = "Neutrophil",`4` = "Neutrophil",
                        `5` = "Dendritic cell",`6` = "Neutrophil",`7` = "B cell",`8` = "Hematopoietic stem cell",`9` = "T cell",`10` = "NK cell",
                        `11` = "Hematopoietic stem cell",`12` = "B cell",`13` = "B cell",`14` = "Erythrocyte",`15` = "Monocyte")
sce.all$cellType=Idents(sce.all)
save(sce.all,file="Blood_seruat.rda")、




#############################################################
######Tissue single-cell data processing--Brain######
#############################################################
###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157827#Download data
setwd("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Brain")
dir='./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Brain/GSE157827_RAW/'     
samples=list.files( dir)
samples 
sceList = lapply(samples,function(pro){ 
  #pro=samples[1]
  folder=file.path(dir,pro) 
  print(pro)
  print(folder)
  print(list.files(folder))
  sce=CreateSeuratObject(counts = Read10X(folder),
                         project =  pro )
  return(sce)
})
names(sceList)  
names(sceList) = samples
sceList
sce.all<- lapply(sceList, FUN = function(x) {
  x[["parcent_mt"]] <- PercentageFeatureSet(x,pattern = "^MT-")#计算线粒体比例
  x<-subset(x,subset = nFeature_RNA>200&nFeature_RNA<20000&parcent_mt<20)#过滤低质量细胞
  x <- NormalizeData(x)##标准化
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = sce.all)
sce.all <- FindIntegrationAnchors(object.list = sce.all,anchor.features = features ,dims=1:20)
sce.all <- IntegrateData(anchorset = sce.all)
sce.all <- ScaleData(sce.all)##缩放
sce.all <- RunPCA(sce.all, npcs = 50, verbose = FALSE)#线性降维
sce.all <- FindNeighbors(sce.all, reduction = "pca",dims = 1:20)#
sce.all <- FindClusters(sce.all, resolution = 1)#找集群
sce.all <- RunUMAP(sce.all, reduction = "pca", dims = 1:20)
pdf("FirstClustering_BrainCluster.pdf",width=9)
DimPlot(sce.all,reduction="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
save(sce.all,file = "frist_seurat")
marker <- read.table("Cell-marker.txt",header = T)
marker <- as.data.frame(marker)
colnames(marker)[2] <- 'Cell Type'
colnames(marker)
table(marker$`Cell Type`)
Astro <- marker[which(marker$`Cell Type`=="Astro"),1]
Endo <- marker[which(marker$`Cell Type`=="Endo"),1]
Excit <- marker[which(marker$`Cell Type`=="Excit"),1]
Inhit <- marker[which(marker$`Cell Type`=="Inhit"),1]
Mic <- marker[which(marker$`Cell Type`=="Mic"),1]
Oligo <- marker[which(marker$`Cell Type`=="Oligo"),1]
set.seed(123)
subobj <- subset(sce.all, downsample = 500)
DoHeatmap(subobj, features = c(Astro[1:100]), disp.min=-2.5, disp.max=2.5)#4,11,23,30,36
DoHeatmap(subobj, features = c(Endo[1:100]), disp.min=-2.5, disp.max=2.5)#26
DoHeatmap(subobj, features = c(Excit ), disp.min=-2.5, disp.max=2.5)#2,6,8,9,13,14,18,19,20,21,22,23,25,28,35      29,34,
DoHeatmap(subobj, features = c(Inhit), disp.min=-2.5, disp.max=2.5)#3,7,10,15,16,27,28,29,32,34,35,36
DoHeatmap(subobj, features = c(Mic[1:100]), disp.min=-2.5, disp.max=2.5)#12,24,31
DoHeatmap(subobj, features = c(Oligo[1:100]), disp.min=-2.5, disp.max=2.5)#0,1,5,22,30,31,32,33
#Astrocytes  4，11，30
#Endothelial cell 26
#Excitatory neurons 2，6，8，9，13，14，18，19，35
#Inhibitory neurons 3，7，10，15，16，17，27，29，34，32
#Microglia 12，24，31
#Oligodendrocytes 0，1，5，33
#删除 20，22，23，25，28，36，
sce.all <- subset(sce.all,subset = seurat_clusters!=20)
sce.all <- subset(sce.all,subset = seurat_clusters!=22)
sce.all <- subset(sce.all,subset = seurat_clusters!=23)
sce.all <- subset(sce.all,subset = seurat_clusters!=25)
sce.all <- subset(sce.all,subset = seurat_clusters!=28)
sce.all <- subset(sce.all,subset = seurat_clusters!=36)
sce.all <- RenameIdents(sce.all,
                        `0` = "Oligodendrocytes",`1` = "Oligodendrocytes",`2` = "Excitatory neurons",`3` = "Inhibitory neurons",`4` = "Astrocytes",`5` = "Oligodendrocytes",
                        `6` = "Excitatory neurons",`7` = "Inhibitory neurons",`8` = "Excitatory neurons",`9` = "Excitatory neurons",`10` = "Inhibitory neurons",
                        `11` = "Astrocytes",`12` = "Microglia",`13` = "Excitatory neurons",`14` = "Excitatory neurons",`15` = "Inhibitory neurons",`16` = "Inhibitory neurons",
                        `17` = "Inhibitory neurons",`18` = "Excitatory neurons",`19` = "Excitatory neurons",`21` = "Excitatory neurons",`24` = "Microglia",
                        `26` = "Endothelial cell", `27` = "Inhibitory neurons",`29` = "Inhibitory neurons",`30` = "Astrocytes",`31` = "Microglia",`32` = "Inhibitory neurons",#Oligo
                        `33` = "Oligodendrocytes",`34` = "Inhibitory neurons",`35` = "Excitatory neurons")
save(sce.all,file = "Brain_seruat.rda")





#############################################################
######Tissue single-cell data processing--Breast######
#############################################################
###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164898#Download data
setwd("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Breast")
dir='./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Breast/data/' 
samples=list.files( dir )
samples
sceList = lapply(samples,function(pro){ 
  print(pro) 
  sce =CreateSeuratObject(counts =  Read10X_h5( file.path(dir,pro)) ,
                          project =   gsub('_filtered_feature_bc_matrix.h5','',gsub('^GSM[0-9]*_','',pro) ) ,
                          min.cells = 5,
                          min.features = 300 )
  return(sce)
})
sceList
samples
#整合数据
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids =   gsub('_filtered_feature_bc_matrix.h5','',gsub('^GSM[0-9]*_','',samples)))
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
sce.all=subset(sce.all,subset = nFeature_RNA > 200)
counts <- GetAssayData(object = sce.all, slot = "counts")
# 根据在每个细胞的计数是否大于0为每个基因输出一个逻辑向量
nonzero <- counts > 0
# 将所有TRUE值相加，如果每个基因的TRUE值超过10个，则返回TRUE。
keep_genes <- Matrix::rowSums(nonzero) >= 5
# 仅保留那些在10个以上细胞中表达的基因
filtered_counts <- counts[keep_genes, ]
# 重新赋值给经过过滤的Seurat对象
sce.all <- CreateSeuratObject(filtered_counts, meta.data = sce.all@meta.data)
sce.all=subset(sce.all,subset = nFeature_RNA > 200 & percent.mt<20)
sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize") 
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst")
all.genes <- rownames(sce.all)
# 结果存放于 `pbmc[["RNA"]]@scale.data`
sce.all <- ScaleData(sce.all, features = all.genes)
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
ElbowPlot(sce.all,ndims = 20)
sce.all<- RunUMAP(sce.all, reduction = "pca", dims = 1:10)
sce.all<- FindNeighbors(sce.all, reduction = "pca", dims = 1:10)
sce.all<- FindClusters(sce.all, resolution = 0.5)
pdf("FirstClustering_BreastCluster_new.pdf",width=9)
DimPlot(sce.all,reduction="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
#细胞注释
diff.wilcox = FindAllMarkers(sce.all, logfc.threshold = 0.25, min.pct = 0.1, 
                             only.pos = TRUE, test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(sce.all)) 
pdf("topMarkerGeneFeature.pdf",width=12,height=20)
DoHeatmap(sce.all, features = top10, group.by = "seurat_clusters", group.bar = T, size = 4)
dev.off()
#NK cell 7
#T cell 0
#Monocyte 4，12，18
#Epithelial cell  1，5，9，10，13，17
#Endothelial cell 2，6，11
#Pericyte 14
#Fibroblast 3，8，16
sce.all <- subset(sce.all,subset = seurat_clusters!=15)
sce.all <- RenameIdents(sce.all,
                        `0` = "T cell",`1` = "Epithelial cell",`2` = "Endothelial cell",`3` = "Fibroblast",`4` = "Monocyte",`5` = "Epithelial cell",`6` = "Endothelial cell",
                        `7` = "NK cell",`8` = "Fibroblast",`9` = "Epithelial cell",`10` = "Epithelial cell",`11` = "Endothelial cell",`12` = "Monocyte",
                        `13` = "Epithelial cell",`14` = "Pericyte",`16` = "Fibroblast",`17` = "Epithelial cell",`18` = "Monocyte")
sce.all$cellType=Idents(sce.all)
save(sce.all,file="Breast_seruat.rda")




#############################################################
######Tissue single-cell data processing--Ovary######
#############################################################
###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118127#Download data
setwd("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Ovary")
dir='./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Ovary/GSE118127_RAW/' 
samples=list.files( dir )
samples
sceList = lapply(samples,function(pro){ 
  print(pro) 
  sce =CreateSeuratObject(counts =  Read10X_h5( file.path(dir,pro)) ,
                          project =   gsub('_filtered_feature_bc_matrix.h5','',gsub('^GSM[0-9]*_','',pro) ) ,
                          min.cells = 5,
                          min.features = 300 )
  return(sce)
})
sceList
samples
#整合数据
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids =   gsub('_filtered_feature_bc_matrix.h5','',gsub('^GSM[0-9]*_','',samples)))
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
sce.all=subset(sce.all,subset = nFeature_RNA > 200 & nFeature_RNA <2500 & nCount_RNA >300 & nCount_RNA < 15000 & percent.mt<10)
sce.all <- NormalizeData(object = sce.all, 
                         normalization.method = "LogNormalize", 
                         scale.factor = 10000)
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst")
all.genes <- rownames(sce.all)
sce.all <- ScaleData(sce.all, features = all.genes)
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
ElbowPlot(sce.all,ndims = 20)
sce.all<- RunUMAP(sce.all, reduction = "pca", dims = 1:15)
sce.all<- FindNeighbors(sce.all, reduction = "pca", dims = 1:15)
sce.all<- FindClusters(sce.all, resolution = 0.9)
pdf("FirstClustering_OvaryCluster.pdf",width=9)
DimPlot(sce.all,reduction="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
mymarker <- c("AMH", "HSD17B1", "SERPINE2", "GSTA1",
              "DCN", "LUM",
              "TAGLN" , "RGS5",
              "VWF" , "CLDN5", 
              "CD53","CXCR4")
pdf("Integrated_cellTypeMarkerDoHeatmap.pdf",width=15,height=20)
DotPlot(sce.all, features = mymarker)
dev.off()
diff.wilcox = FindAllMarkers(sce.all, logfc.threshold = 0.25, min.pct = 0.1, 
                             only.pos = TRUE, test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(sce.all)) 
pdf("topMarkerGeneFeature.pdf",width=12,height=20)
DoHeatmap(sce.all, features = top10, group.by = "seurat_clusters", group.bar = T, size = 4)
dev.off()
#Granulosa cell 13,14,15,21
#Fibroblast 2,3,4,5,9,12,17
#Smooth muscle cell 6,16
#Endothelial cell 7,8,11,18,20,22
#T cell 10
#Dendritic cell 19
#B cell 24
#Theca cell 0，1
sce.all <- subset(sce.all,subset = seurat_clusters!=23)
sce.all <- RenameIdents(sce.all,
                        `0` = "Theca cell",`1` = "Theca cell",`2` = "Fibroblast",`3` = "Fibroblast",`4` = "Fibroblast",`5` = "Fibroblast",`6` = "Smooth muscle cell",
                        `7` = "Endothelial cell",`8` = "Endothelial cell",`9` = "Fibroblast",`10` = "T cell",`11` = "Endothelial cell",`12` = "Fibroblast",
                        `13` = "Granulosa cell",`14` = "Granulosa cell",`15` = "Granulosa cell",`16` = "Smooth muscle cell",`17` = "Fibroblast",`18` = "Endothelial cell",
                        `19` = "Dendritic cell",`20` = "Endothelial cell",`21` = "Granulosa cell",`22` = "Endothelial cell",`24` = "B cell")
sce.all$cellType=Idents(sce.all)
save(sce.all,file="Ovary_seruat.rda")



#############################################################
######Tissue single-cell data processing--Skin######
#############################################################
###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193304 #Download data
setwd("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Skin")
dir='./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Skin/data/' 
samples=list.files( dir )
samples
sceList = lapply(samples,function(pro){ 
  print(pro) 
  sce =CreateSeuratObject(counts =  Read10X_h5( file.path(dir,pro)) ,
                          project =   gsub('_filtered_feature_bc_matrix.h5','',gsub('^GSM[0-9]*_','',pro) ) ,
                          min.cells = 5,
                          min.features = 300 )
  return(sce)
})
sceList
samples
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids =   gsub('_filtered_feature_bc_matrix.h5','',gsub('^GSM[0-9]*_','',samples)))
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
sce.all=subset(sce.all,subset = nFeature_RNA > 200 & percent.mt<10)
sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize") 
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst")
all.genes <- rownames(sce.all)
sce.all <- ScaleData(sce.all, features = all.genes)
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
ElbowPlot(sce.all,ndims = 20)
sce.all<- RunUMAP(sce.all, reduction = "pca", dims = 1:10)
sce.all<- FindNeighbors(sce.all, reduction = "pca", dims = 1:10)
sce.all<- FindClusters(sce.all, resolution = 0.5)
pdf("FirstClustering_SkinCluster.pdf",width=9)
DimPlot(sce.all,reduction="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
MyMarkerGene=c("COL17A1","KRT5","KRT14","KRT1","KRT10","FLG","LOR","MKI67","TOP2A","KRT6B","KRT17","SFRP1",
               "CD207","CD1A","CD3D","PTPRC","PMEL","TYRP1","DCN","COL1A1")
pdf("Integrated_MyMarkerGeneFeature.pdf",width=15,height=10)
DotPlot(sce.all, features = MyMarkerGene)
dev.off()
Basal_cells <- c("COL17A1", "KRT5", "KRT14")#1，2，5，6，11，16
Epithelial_cell <-  c("KRT1","KRT10")#0，3，4，7，10，15
Granule_cell <-  c("FLG","LOR")#12
follicular_cells <- c("KRT6B", "KRT17", "SFRP1")#13
Dendritic_cell <- c("CD207", "CD1A")#8
T_cell <- c("CD3D", "PTPRC")#9
Melanocyte <-  c("PMEL", "TYRP1")#14 
Fibroblast <- c("DCN", "COL1A1")#17
sce.all <- subset(sce.all,subset = seurat_clusters!=18)
sce.all <- RenameIdents(sce.all,
                        `0` = "Epithelial cell",`1` = "Basal cell",`2` = "Basal cell",`3` = "Epithelial cell",`4`  = "Epithelial cell",`5` = "Basal cell",
                        `6` = "Basal cell",`7` = "Epithelial cell",`8` = "Dendritic cell",`9` = "T cell",`10` = "Epithelial cell",`11` = "Basal cell",
                        `12` = "Granule cell",`13` = "Follicular cell",`14` = "Melanocyte",`15` = "Epithelial cell",`16` = "Basal cell",`17` = "Fibroblast")
sce.all$cellType=Idents(sce.all)
save(sce.all,file="Skin_seruat.rda")



#############################################################
######Tissue single-cell data processing--Testis######
#############################################################
###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112013 #Download data
setwd("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Testis")
count1 <- read.table("GSE112013_Combined_UMI_table.txt.gz", header = TRUE,row.names = 1)
sce.all <- CreateSeuratObject(counts=count1)
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
sce.all=subset(sce.all,subset = nFeature_RNA > 500  & percent.mt<20)
sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize", scale.factor = 10000)
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sce.all)
sce.all <- ScaleData(sce.all, features = all.genes)
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
ElbowPlot(sce.all,ndims = 20)
sce.all<- RunUMAP(sce.all, reduction = "pca", dims = 1:10)
sce.all<- FindNeighbors(sce.all, reduction = "pca", dims = 1:10)
sce.all<- FindClusters(sce.all, resolution = 0.5)
pdf("FirstClustering_TesitisCluster.pdf",width=9)
DimPlot(sce.all,reduction="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

diff.wilcox = FindAllMarkers(sce.all, logfc.threshold = 0.25, min.pct = 0.1, 
                             only.pos = TRUE, test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(sce.all)) 
pdf("topMarkerGeneFeature.pdf",width=12,height=20)
DoHeatmap(sce.all, features = top10, group.by = "seurat_clusters", group.bar = T, size = 4)
dev.off()
#Peritubular myoid cells  2
#Germ cell 0,1,3,4,5,8,9,11,13,14
#Endothelial 6
#Dendritic cell 7
#Smooth muscle cells 10
#Pericytes 15
sce.all <- RenameIdents(sce.all,
                        `0` = "Germ cell",`1` = "Germ cell",`2` = "Smooth muscle cell",`3` = "Germ cell",`4` = "Germ cell",`5` = "Germ cell",
                        `6` = "Endothelial cell",`7` = "Dendritic cell",`8` = "Germ cell",`9` = "Germ cell",`10` = "Smooth muscle cell",
                        `11` = "Germ cell",`12` = "Germ cell",`13` = "Germ cell",`14` = "Germ cell",`15` = "Pericytes")
sce.all$cellType=Idents(sce.all)
save(sce.all,file="Tesitis_seruat.rda")


library(ggrepel)
#############################################################
######Cycle through the UMAP charts for each organization####
#############################################################
#---Fig.S3
tissue <- c("Uterus","Prostate","Liver","Breast","Lung","Colon","Kidney","Skin","Ovary","Esophagus","Brain","Blood","Adrenal","Pancreas","Testis")
all_umap <- list()
for(x in 1:15){
  setwd(paste0("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/",tissue[x]))
  load(file = paste0(tissue[x],"_seruat.rda"))
  umap = sce.all@reductions$umap@cell.embeddings %>%
    as.data.frame() %>% 
    cbind(cell_type = sce.all@meta.data$cellType)
  head(umap)
  allcolour <- c("Endothelial cell"="#936CA3","Epithelial cell"="#EE7C6F","Fibroblast"="#BEA2CA","Mast cell"="#72A3C6","Smooth muscle cell"="#89CCC0", 
                 "Macrophage"="#023756","T cell"="#901D43","Dendritic cell"="#136BA5","Hepatocyte"="#A4C9DD","Monocyte"="#085787",              
                 "Neutrophil"="#3989B5","NK cell"="#C11B53","Pericyte"="#F8F5B3","B cell"="#D63081","Enteric nerval cell"="#C6B0D2",    
                 "Plasma cell"="#DF65B0","Principal cell"="#99000D","Basal cell"="#FDB462","Follicular cell"="#FB6A4A","Granule cell"="#EF3B2C",           
                 "Melanocyte"="#FC9272","Granulosa cell"="#FCCDE5","Astrocytes"="#D9F0D3","Excitatory neurons"="#1B7837","Microglia"="#A6DBA0",              
                 "Oligodendrocytes"="#5AAE61","Hematopoietic stem cell"="#A6BDDB","Zona fasciculata cell"="#6A3D9A","Germ cell"="#B15928",
                 "Loop of henle"="#637A9F","Proximal tubule cell"="#C9D7DD","Distal tubule cell"="#AC7D88","Intercalated cell"="#6D2932","Plasmocyte"="#E6A4B4",
                 "Theca cell"="#2D3250","Inhibitory neurons"="#86A789","Platelets"="#B67352","Alpha cell"="#E78895","Acinar cell"="#CDFADB","Exocrine cell"="#B4D4FF",
                 "Pericytes"="#BBAB8C","Erythrocyte"="#FF3EA5")
  p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +
    
    geom_point(size = 0.6, alpha =0.75)+
    
    scale_color_manual(values = allcolour)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = 'white'),
          plot.background=element_rect(fill="white"))+
    theme(legend.title = element_blank(),
          legend.key=element_rect(fill='white'),
          legend.text = element_text(size=10),
          legend.key.size=unit(0.5,'cm') ) + 
    guides(color = guide_legend(override.aes = list(size=5)))
  
  cell_type_med <- umap %>%
    group_by(cell_type) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  
  p1 <- p+geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med,
                           point.padding=unit(0.5, "lines"))+
    theme(legend.position = "none")
  all_umap[[x]] <-p1
  pdf(paste0("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/All_UMAP/",tissue[x],"_UMAP.pdf"),width=8,height = 6)
  print(p1)
  dev.off()
}
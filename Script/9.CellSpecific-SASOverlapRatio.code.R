######Cell type-specific SAS in cancer and normal intersections######

####----Fig.3B--------
load(file="/boot3/cjl/Aging/data/UCSC/TCGA_RNAseq_data/CellPercentage/所有组织的细胞上下调基因.rda")
load(file="/boot3/cjl/Aging/data/UCSC/TCGA_RNAseq_data/CellPercentage/interscertCellList.rda")
load(file="/boot3/cjl/Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/组织统计/Snakey_data.rda")
load(file="/boot3/cjl/Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/组织统计/all_tissue_module_list.rda")
cellmodule_list <- all_tissue_module_list[unique(Sankey_data$pathway)]
#UTUA_info----
UTUA_gene_info <- interscertCellList[c("ESCA_UTUACellInter","LGG_UTUACellInter","GBM_UTUACellInter","TGCT_UTUACellInter")]
names(UTUA_gene_info) <- c("ESCA_Esophagus_UTUACellInter","LGG_Brain_UTUACellInter","GBM_Brain_UTUACellInter","TGCT_Testis_UTUACellInter")
plot_list <- list()
for (i in 1:length(names(UTUA_gene_info))) {
  one_tissue <- str_extract(names(UTUA_gene_info)[i], "(?<=_).*(?=_)")
  oneTissueCellInfo <- allTissue_AllCelllist[[one_tissue]]
  oneTissueCellInfo_Up <- oneTissueCellInfo[grep("Up",names(oneTissueCellInfo))]
  one_tisue_CellmoduleInfo <- Sankey_data[which(Sankey_data$tissue == one_tissue & Sankey_data$moduleType == "Up"),]
  one_tisue_CellmoduleList <- all_tissue_module_list[one_tisue_CellmoduleInfo$pathway]
  cellzhanbi <- as.data.frame(matrix(NA,nrow = length(names(oneTissueCellInfo_Up)),ncol=1))
  one_UTUA_gene_info <- UTUA_gene_info[i]
  cellzhanbiList <- list()
  for (z in 1:length(names(oneTissueCellInfo_Up))) {
    cellzhanbi[z,1]<- length(intersect(oneTissueCellInfo_Up[[z]],unlist(one_UTUA_gene_info)))
    cellzhanbiList[[z]]<- intersect(oneTissueCellInfo_Up[[z]],unlist(one_UTUA_gene_info))
  }
  split_strings <- strsplit(names(oneTissueCellInfo_Up), "_")
  rownames(cellzhanbi) <- lapply(split_strings, function(x) x[[1]])
  colnames(cellzhanbi) <- "UTUA"
  names(cellzhanbiList) <- rownames(cellzhanbi)
  cellzhanbi$module <- 0
  allmoduleInfo <- list()
  for (j in 1:length(rownames(cellzhanbi))) {
    oneTissueCellmoduleInfo <- one_tisue_CellmoduleInfo[which(one_tisue_CellmoduleInfo$cell_name == rownames(cellzhanbi)[j]),]
    oneTissueCellmodule <- one_tisue_CellmoduleList[oneTissueCellmoduleInfo$pathway]
    genes_to_find <- cellzhanbiList[[j]]
    found_lists <- c()
    # 遍历大列表中的每个子列表
    for (lst in oneTissueCellmodule) {
      # 检查子列表是否包含要查找的基因
      if (any(genes_to_find %in% lst)) {
        found_lists <- c(found_lists, names(oneTissueCellmodule)[which(lst %in% genes_to_find)])
      }
    }
    moduleInfo <- oneTissueCellmodule[unique(found_lists)[!is.na(unique(found_lists))]]
    cellzhanbi[j,2] <- length(names(moduleInfo))
    allmoduleInfo <- c(allmoduleInfo,list(moduleInfo))
  }
  names(allmoduleInfo) <- rownames(cellzhanbi)
  save(cellzhanbi,file=paste0("/boot3/cjl/Aging/data/UCSC/TCGA_RNAseq_data/CellPercentage/UTUA/",names(UTUA_gene_info)[i],"_zhanbi.rda"))
  save(allmoduleInfo,file=paste0("/boot3/cjl/Aging/data/UCSC/TCGA_RNAseq_data/CellPercentage/UTUA/",names(UTUA_gene_info)[i],"_moduleInfo.rda"))
  
  data1 <- cellzhanbi
  data1$Cell <- rownames(cellzhanbi)
  data1 <- data1[order(data1$module),]
  data1$Cell <- factor(data1$Cell,levels = data1$Cell)
  p1 <- ggplot(data1, aes(x="", y = module, fill = Cell))+
    geom_bar(width = 1, stat = "identity",color="white")+
    coord_polar('y')+
    theme_void()+
    ggtitle(names(UTUA_gene_info)[i])+
    scale_fill_manual(values=c("Endothelial cell"="#9970AB","Epithelial cell"="#FB8072","Fibroblast"="#C2A5CF","Mast cell"="#74A9CF","Smooth muscle cell"="#8DD3C7", 
                               "Macrophage"="#023858","T cell"="#980043","Dendritic cell"="#0570B0","Hepatocyte"="#A6CEE3","Monocyte"="#045A8D",              
                               "Neutrophil"="#3690C0","NK cell"="#CE1256","Pericyte"="#FFFFB3","B cell"="#E7298A","Enteric nerval cell"="#CAB2D6",    
                               "Plasma cell"="#DF65B0","Principal cell"="#99000D","Basal cell"="#FDB462","Follicular cell"="#FB6A4A","Granule cell"="#EF3B2C",           
                               "Melanocyte"="#FC9272","Granulosa cell"="#FCCDE5","Astrocytes"="#D9F0D3","Excitatory neurons"="#1B7837","Microglia"="#A6DBA0",              
                               "Oligodendrocytes"="#5AAE61","Hematopoietic stem cell"="#A6BDDB","Zona fasciculata cell"="#6A3D9A","Germ cell"="#B15928"))+#自定义颜色
    geom_text(aes(y = sum(module)-cumsum(module)+module/2,
                  label = paste0(Cell,"\n",module)), size=4.5)
  plot_list[[i]] <- p1
}
pdf("/boot3/cjl/Aging/data/UCSC/TCGA_RNAseq_data/CellPercentage/UTUA/UTUA_gene.pdf",
    width = 6,height = 26)
plot(plot_grid(plotlist = plot_list, ncol = 1))
dev.off()


#DTDA_info----
DTDA_gene_info <- interscertCellList[c("UCEC_DTDACellInter","LIHC_DTDACellInter","ESCA_DTDACellInter","TGCT_DTDACellInter")]
names(DTDA_gene_info) <- c("UCEC_Uterus_DTDACellInter","LIHC_Liver_DTDACellInter","ESCA_Esophagus_DTDACellInter","TGCT_Testis_DTDACellInter")

plot_list <- list()
for (i in 1:length(names(DTDA_gene_info))) {
  one_tissue <- str_extract(names(DTDA_gene_info)[i], "(?<=_).*(?=_)")
  oneTissueCellInfo <- allTissue_AllCelllist[[one_tissue]]
  oneTissueCellInfo_Down <- oneTissueCellInfo[grep("Down",names(oneTissueCellInfo))]
  
  one_tisue_CellmoduleInfo <- Sankey_data[which(Sankey_data$tissue == one_tissue & Sankey_data$moduleType == "Down"),]
  one_tisue_CellmoduleList <- all_tissue_module_list[one_tisue_CellmoduleInfo$pathway]
  
  cellzhanbi <- as.data.frame(matrix(NA,nrow = length(names(oneTissueCellInfo_Down)),ncol=1))
  
  one_DTDA_gene_info <- DTDA_gene_info[i]
  
  cellzhanbiList <- list()
  for (z in 1:length(names(oneTissueCellInfo_Down))) {
    cellzhanbi[z,1]<- length(intersect(oneTissueCellInfo_Down[[z]],unlist(one_DTDA_gene_info)))
    cellzhanbiList[[z]]<- intersect(oneTissueCellInfo_Down[[z]],unlist(one_DTDA_gene_info))
  }
  split_strings <- strsplit(names(oneTissueCellInfo_Down), "_")
  rownames(cellzhanbi) <- lapply(split_strings, function(x) x[[1]])
  colnames(cellzhanbi) <- "DTDA"
  names(cellzhanbiList) <- rownames(cellzhanbi)
  
  cellzhanbi$module <- 0
  
  allmoduleInfo <- list()
  for (j in 1:length(rownames(cellzhanbi))) {
    oneTissueCellmoduleInfo <- one_tisue_CellmoduleInfo[which(one_tisue_CellmoduleInfo$cell_name == rownames(cellzhanbi)[j]),]
    oneTissueCellmodule <- one_tisue_CellmoduleList[oneTissueCellmoduleInfo$pathway]
    genes_to_find <- cellzhanbiList[[j]]
    found_lists <- c()
    
    # 遍历大列表中的每个子列表
    for (lst in oneTissueCellmodule) {
      # 检查子列表是否包含要查找的基因
      if (any(genes_to_find %in% lst)) {
        found_lists <- c(found_lists, names(oneTissueCellmodule)[which(lst %in% genes_to_find)])
      }
    }
    
    moduleInfo <- oneTissueCellmodule[unique(found_lists)[!is.na(unique(found_lists))]]
    cellzhanbi[j,2] <- length(names(moduleInfo))
    allmoduleInfo <- c(allmoduleInfo,list(moduleInfo))
  }
  names(allmoduleInfo) <- rownames(cellzhanbi)
  
  save(cellzhanbi,file=paste0("/boot3/cjl/Aging/data/UCSC/TCGA_RNAseq_data/CellPercentage/DTDA/",names(DTDA_gene_info)[i],"_zhanbi.rda"))
  save(allmoduleInfo,file=paste0("/boot3/cjl/Aging/data/UCSC/TCGA_RNAseq_data/CellPercentage/DTDA/",names(DTDA_gene_info)[i],"_moduleInfo.rda"))
  
  data1 <- cellzhanbi
  data1$Cell <- rownames(cellzhanbi)
  data1 <- data1[order(data1$module),]
  data1$Cell <- factor(data1$Cell,levels = data1$Cell)
  p1 <- ggplot(data1, aes(x="", y = module, fill = Cell))+
    geom_bar(width = 1, stat = "identity",color="white")+
    coord_polar('y')+
    theme_void()+
    ggtitle(names(DTDA_gene_info)[i])+
    scale_fill_manual(values=c("Endothelial cell"="#9970AB","Epithelial cell"="#FB8072","Fibroblast"="#C2A5CF","Mast cell"="#74A9CF","Smooth muscle cell"="#8DD3C7", 
                               "Macrophage"="#023858","T cell"="#980043","Dendritic cell"="#0570B0","Hepatocyte"="#A6CEE3","Monocyte"="#045A8D",              
                               "Neutrophil"="#3690C0","NK cell"="#CE1256","Pericyte"="#FFFFB3","B cell"="#E7298A","Enteric nerval cell"="#CAB2D6",    
                               "Plasma cell"="#DF65B0","Principal cell"="#99000D","Basal cell"="#FDB462","Follicular cell"="#FB6A4A","Granule cell"="#EF3B2C",           
                               "Melanocyte"="#FC9272","Granulosa cell"="#FCCDE5","Astrocytes"="#D9F0D3","Excitatory neurons"="#1B7837","Microglia"="#A6DBA0",              
                               "Oligodendrocytes"="#5AAE61","Hematopoietic stem cell"="#A6BDDB","Zona fasciculata cell"="#6A3D9A","Germ cell"="#B15928"))+#自定义颜色
    geom_text(aes(y = sum(module)-cumsum(module)+module/2,
                  label = paste0(Cell,"\n",module)), size=4.5)
  plot_list[[i]] <- p1
}

pdf("/boot3/cjl/Aging/data/UCSC/TCGA_RNAseq_data/CellPercentage/DTDA/DTDA_gene.pdf",
    width = 6,height = 26)
plot(plot_grid(plotlist = plot_list, ncol = 1))
dev.off()



#DTUA_info-----
DTUA_gene_info <- interscertCellList[c("PRAD_DTUACellInter","BRCA_DTUACellInter","LUAD_DTUACellInter","LUSC_DTUACellInter",
                                       "KICH_DTUACellInter","ESCA_DTUACellInter","THYM_DTUACellInter")]
names(DTUA_gene_info) <- c("PRAD_Prostate_DTUACellInter","BRCA_Breast_DTUACellInter","LUAD_Lung_DTUACellInter","LUSC_Lung_DTUACellInter",
                           "KICH_Kidney_DTUACellInter","ESCA_Esophagus_DTUACellInter","THYM_Blood_DTUACellInter")
plot_list <- list()
for (i in 1:length(names(DTUA_gene_info))) {
  one_tissue <- str_extract(names(DTUA_gene_info)[i], "(?<=_).*(?=_)")
  oneTissueCellInfo <- allTissue_AllCelllist[[one_tissue]]
  oneTissueCellInfo_Up <- oneTissueCellInfo[grep("Up",names(oneTissueCellInfo))]
  one_tisue_CellmoduleInfo <- Sankey_data[which(Sankey_data$tissue == one_tissue & Sankey_data$moduleType == "Up"),]
  one_tisue_CellmoduleList <- cellmodule_list[unique(one_tisue_CellmoduleInfo$pathway)]
  cellzhanbi <- as.data.frame(matrix(NA,nrow = length(names(oneTissueCellInfo_Up)),ncol=1))
  one_DTUA_gene_info <- DTUA_gene_info[i]
  cellzhanbiList <- list()
  for (z in 1:length(names(oneTissueCellInfo_Up))) {
    cellzhanbi[z,1]<- length(intersect(oneTissueCellInfo_Up[[z]],unlist(one_DTUA_gene_info)))
    cellzhanbiList[[z]]<- intersect(oneTissueCellInfo_Up[[z]],unlist(one_DTUA_gene_info))
  }
  split_strings <- strsplit(names(oneTissueCellInfo_Up), "_")
  rownames(cellzhanbi) <- lapply(split_strings, function(x) x[[1]])
  colnames(cellzhanbi) <- "DTUA"
  names(cellzhanbiList) <- rownames(cellzhanbi)
  cellzhanbi$module <- 0
  allmoduleInfo <- list()
  for (j in 1:length(rownames(cellzhanbi))) {
    oneTissueCellmoduleInfo <- one_tisue_CellmoduleInfo[which(one_tisue_CellmoduleInfo$cell_name == rownames(cellzhanbi)[j]),]
    oneTissueCellmodule <- one_tisue_CellmoduleList[oneTissueCellmoduleInfo$pathway]
    genes_to_find <- cellzhanbiList[[j]]
    found_lists <- c()
    # 遍历大列表中的每个子列表
    for (lst in oneTissueCellmodule) {
      # 检查子列表是否包含要查找的基因
      if (any(genes_to_find %in% lst)) {
        found_lists <- c(found_lists, names(oneTissueCellmodule)[which(lst %in% genes_to_find)])
      }
    }
    moduleInfo <- oneTissueCellmodule[unique(found_lists)[!is.na(unique(found_lists))]]
    cellzhanbi[j,2] <- length(names(moduleInfo))
    allmoduleInfo <- c(allmoduleInfo,list(moduleInfo))
  }
  names(allmoduleInfo) <- rownames(cellzhanbi)
  save(cellzhanbi,file=paste0(".Aging/data/UCSC/TCGA_RNAseq_data/CellPercentage/DTUA/",names(DTUA_gene_info)[i],"_zhanbi.rda"))
  save(allmoduleInfo,file=paste0("./Aging/data/UCSC/TCGA_RNAseq_data/CellPercentage/DTUA/",names(DTUA_gene_info)[i],"_moduleInfo.rda"))
  data1 <- cellzhanbi
  data1$Cell <- rownames(cellzhanbi)
  data1 <- data1[order(data1$module),]
  data1$Cell <- factor(data1$Cell,levels = data1$Cell)
  p1 <- ggplot(data1, aes(x="", y = module, fill = Cell))+
    geom_bar(width = 1, stat = "identity",color="white")+
    coord_polar('y')+
    theme_void()+
    ggtitle(names(DTUA_gene_info)[i])+
    scale_fill_manual(values=c("Endothelial cell"="#9970AB","Epithelial cell"="#FB8072","Fibroblast"="#C2A5CF","Mast cell"="#74A9CF","Smooth muscle cell"="#8DD3C7", 
                               "Macrophage"="#023858","T cell"="#980043","Dendritic cell"="#0570B0","Hepatocyte"="#A6CEE3","Monocyte"="#045A8D",              
                               "Neutrophil"="#3690C0","NK cell"="#CE1256","Pericyte"="#FFFFB3","B cell"="#E7298A","Enteric nerval cell"="#CAB2D6",    
                               "Plasma cell"="#DF65B0","Principal cell"="#99000D","Basal cell"="#FDB462","Follicular cell"="#FB6A4A","Granule cell"="#EF3B2C",           
                               "Melanocyte"="#FC9272","Granulosa cell"="#FCCDE5","Astrocytes"="#D9F0D3","Excitatory neurons"="#1B7837","Microglia"="#A6DBA0",              
                               "Oligodendrocytes"="#5AAE61","Hematopoietic stem cell"="#A6BDDB","Zona fasciculata cell"="#6A3D9A","Germ cell"="#B15928"))+#自定义颜色
    geom_text(aes(y = sum(module)-cumsum(module)+module/2,
                  label = paste0(Cell,"\n",module)), size=4.5)
  plot_list[[i]] <- p1
}
pdf("./Aging/data/UCSC/TCGA_RNAseq_data/CellPercentage/DTUA/DTUA_gene.pdf",
    width = 6,height = 18)
plot(plot_grid(plotlist = plot_list, ncol = 1))
dev.off()

#UTDA----
UTDA_gene_info <- interscertCellList[c("PRAD_UTDACellInter","LIHC_UTDACellInter","BRCA_UTDACellInter","LUAD_UTDACellInter",
                                       "LUSC_UTDACellInter","COAD_UTDACellInter","READ_UTDACellInter","KIRC_UTDACellInter",
                                       "KIRP_UTDACellInter","OV_UTDACellInter","DLBC_UTDACellInter","PAAD_UTDACellInter",
                                       "THYM_UTDACellInter")]
names(UTDA_gene_info) <- c("PRAD_Prostate_UTDACellInter","LIHC_Liver_UTDACellInter","BRCA_Breast_UTDACellInter","LUAD_Lung_UTDACellInter",
                           "LUSC_Lung_UTDACellInter","COAD_Colon_UTDACellInter","READ_Colon_UTDACellInter","KIRC_Kidney_UTDACellInter",
                           "KIRP_Kidney_UTDACellInter","OV_Ovary_UTDACellInter","DLBC_Blood_UTDACellInter","PAAD_Pancreas_UTDACellInter",
                           "THYM_Blood_UTDACellInter")
plot_list <- list()
for (i in 1:length(names(UTDA_gene_info))) {
  one_tissue <- str_extract(names(UTDA_gene_info)[i], "(?<=_).*(?=_)")
  oneTissueCellInfo <- allTissue_AllCelllist[[one_tissue]]
  oneTissueCellInfo_Down <- oneTissueCellInfo[grep("Down",names(oneTissueCellInfo))]
  one_tisue_CellmoduleInfo <- Sankey_data[which(Sankey_data$tissue == one_tissue & Sankey_data$moduleType == "Down"),]
  one_tisue_CellmoduleList <- all_tissue_module_list[one_tisue_CellmoduleInfo$pathway]
  cellzhanbi <- as.data.frame(matrix(NA,nrow = length(names(oneTissueCellInfo_Down)),ncol=1))
  one_UTDA_gene_info <- UTDA_gene_info[i]
  cellzhanbiList <- list()
  for (z in 1:length(names(oneTissueCellInfo_Down))) {
    cellzhanbi[z,1]<- length(intersect(oneTissueCellInfo_Down[[z]],unlist(one_UTDA_gene_info)))
    cellzhanbiList[[z]]<- intersect(oneTissueCellInfo_Down[[z]],unlist(one_UTDA_gene_info))
  }
  split_strings <- strsplit(names(oneTissueCellInfo_Down), "_")
  rownames(cellzhanbi) <- lapply(split_strings, function(x) x[[1]])
  colnames(cellzhanbi) <- "DTUA"
  names(cellzhanbiList) <- rownames(cellzhanbi)
  cellzhanbi$module <- 0
  allmoduleInfo <- list()
  for (j in 1:length(rownames(cellzhanbi))) {
    oneTissueCellmoduleInfo <- one_tisue_CellmoduleInfo[which(one_tisue_CellmoduleInfo$cell_name == rownames(cellzhanbi)[j]),]
    oneTissueCellmodule <- one_tisue_CellmoduleList[oneTissueCellmoduleInfo$pathway]
    genes_to_find <- cellzhanbiList[[j]]
    found_lists <- c()
    # 遍历大列表中的每个子列表
    for (lst in oneTissueCellmodule) {
      # 检查子列表是否包含要查找的基因
      if (any(genes_to_find %in% lst)) {
        found_lists <- c(found_lists, names(oneTissueCellmodule)[which(lst %in% genes_to_find)])
      }
    }
    moduleInfo <- oneTissueCellmodule[unique(found_lists)[!is.na(unique(found_lists))]]
    cellzhanbi[j,2] <- length(names(moduleInfo))
    allmoduleInfo <- c(allmoduleInfo,list(moduleInfo))
  }
  names(allmoduleInfo) <- rownames(cellzhanbi)
  save(cellzhanbi,file=paste0("./Aging/data/UCSC/TCGA_RNAseq_data/CellPercentage/UTDA/",names(UTDA_gene_info)[i],"_zhanbi.rda"))
  save(allmoduleInfo,file=paste0("./Aging/data/UCSC/TCGA_RNAseq_data/CellPercentage/UTDA/",names(UTDA_gene_info)[i],"_moduleInfo.rda"))
  data1 <- cellzhanbi
  data1$Cell <- rownames(cellzhanbi)
  data1 <- data1[order(data1$module),]
  data1$Cell <- factor(data1$Cell,levels = data1$Cell)
  p1 <- ggplot(data1, aes(x="", y = module, fill = Cell))+
    geom_bar(width = 1, stat = "identity",color="white")+
    coord_polar('y')+
    theme_void()+
    ggtitle(names(UTDA_gene_info)[i])+
    scale_fill_manual(values=c("Endothelial cell"="#9970AB","Epithelial cell"="#FB8072","Fibroblast"="#C2A5CF","Mast cell"="#74A9CF","Smooth muscle cell"="#8DD3C7", 
                               "Macrophage"="#023858","T cell"="#980043","Dendritic cell"="#0570B0","Hepatocyte"="#A6CEE3","Monocyte"="#045A8D",              
                               "Neutrophil"="#3690C0","NK cell"="#CE1256","Pericyte"="#FFFFB3","B cell"="#E7298A","Enteric nerval cell"="#CAB2D6",    
                               "Plasma cell"="#DF65B0","Principal cell"="#99000D","Basal cell"="#FDB462","Follicular cell"="#FB6A4A","Granule cell"="#EF3B2C",           
                               "Melanocyte"="#FC9272","Granulosa cell"="#FCCDE5","Astrocytes"="#D9F0D3","Excitatory neurons"="#1B7837","Microglia"="#A6DBA0",              
                               "Oligodendrocytes"="#5AAE61","Hematopoietic stem cell"="#A6BDDB","Zona fasciculata cell"="#6A3D9A","Germ cell"="#B15928"))+#自定义颜色
    geom_text(aes(y = sum(module)-cumsum(module)+module/2,
                  label = paste0(Cell,"\n",module)), size=4.5)
  plot_list[[i]] <- p1
}

pdf("./Aging/data/UCSC/TCGA_RNAseq_data/CellPercentage/UTDA/UTDA_gene.pdf",
    width = 6,height = 26)
plot(plot_grid(plotlist = plot_list, ncol = 1))
dev.off()
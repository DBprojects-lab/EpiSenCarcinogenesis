#############################################################
######SACMs Cell Type Annotation######
#############################################################
#differential genes were calculated for each cell type using the "wilcoxauc" function in the single-cell data. 
#Then, gene set enrichment analyses were performed using the "fgseaMultilevel" function to identify the cell type of each senescence module.
load(file = "./Aging/data/UCSC/GTEX_RNAseq_data/All_sig_module.Rdata")
tissue <- c("Uterus","Prostate","Liver","Breast","Lung","Colon","Kidney","Skin","Ovary","Esophagus","Brain","Blood","Adrenal","Pancreas","Testis")
for (i in 1:length(tissue)) {
  load(file=paste0(tissue[i],"_seruat.rda"))
  integratedMarker <- wilcoxauc(sce.all, 'cellType')
  table(integratedMarker$group)
  a <- data.frame(module_name = names(All_sig_module[[paste0(tissue[i],"_up")]]),
                  type=rep("Up",length(names(All_sig_module[[paste0(tissue[i],"_up")]]))))
  b <- data.frame(module_name = names(All_sig_module[[paste0(tissue[i],"_down")]]),
                  type=rep("Down",length(names(All_sig_module[[paste0(tissue[i],"_down")]]))))
  one_tissue_module_info <- rbind(a,b)
  one_tissue_module <- c(All_sig_module[[paste0(tissue[i],"_up")]],All_sig_module[[paste0(tissue[i],"_down")]])
  one_tissue_fgseaRes <- data.frame()
  for(cluster in unique(Blood.integratedMarker$group)){
    print (cluster)
    clusterCell<- Blood.integratedMarker %>% dplyr::filter(group == cluster) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
    ranks<- deframe(clusterCell)
    ranks=na.omit(ranks)
    fgseaRes<- fgseaMultilevel(one_tissue_module, stats = ranks,eps=0, nPermSimple = 10000)
    fgseaRes$cell_name <- cluster
    one_tissue_fgseaRes <- rbind(one_tissue_fgseaRes,fgseaRes)
  }
  save(one_tissue_fgseaRes,file = paste0(tissue[i],"_fgseaRes.rda"))
  
  for (module in unique(one_tissue_fgseaRes$pathway)) {
    print (module)
    module_fgseaRes <- one_tissue_fgseaRes %>% dplyr::filter(pathway == module)
    if(length(table(module_fgseaRes$pval<0.05&module_fgseaRes$NES>0))== 1){
      one_tissue_fgseaRes <- one_tissue_fgseaRes[-which(one_tissue_fgseaRes$pathway == module),]
    }
  }
  one_tissue_module_info1 <- one_tissue_module_info[which(one_tissue_module_info$module_name %in% unique(one_tissue_fgseaRes$pathway)),]
  one_tissue_fgseaRes$sig=ifelse(one_tissue_fgseaRes$pval<0.05,"Sig","NonSig")
  one_tissue_fgseaRes$pathway <- factor(one_tissue_fgseaRes$pathway,levels = one_tissue_module_info1$module_name)
  #----Fig.S4
  p1 <- ggplot(one_tissue_fgseaRes,aes(pathway,cell_name,size=-1*log(pval),colour=NES,shape=sig))+geom_point()+
    scale_color_gradient2(low="white",mid="white",high = "red")+
    scale_shape_manual(values=c(3,19))+
    theme_bw()+#theme(legend.position="bottom")+
    theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=rel(1.0),colour = "black",angle=90, vjust = 0.5, hjust=1),axis.text.y = element_text(size=rel(1.0),colour = "black"))
  pdf(paste0("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/fgsea/",tissue[x],"_fgsea.pdf"),width=8,height = 6)
  print(p1)
  dev.off()
}

#############################################################
######Analysis of different cell types######
#############################################################
##Get information about significant overlapping modules
tissue <- c("Uterus","Prostate","Liver","Breast","Lung","Colon","Kidney","Skin","Ovary","Esophagus","Brain","Blood","Adrenal","Pancreas","Testis")
alltissue_fgseaRes <- data_frame()
for (x in 1:length(tissue)) {
  setwd(paste0("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/",tissue[x]))
  load(file = paste0(tissue[x],"_fgseaRes.rda"))
  alltissue_fgseaRes <- rbind(alltissue_fgseaRes,one_tissue_fgseaRes)
}
alltissue_fgseaRes$cell_name[which(alltissue_fgseaRes$cell_name == "Epithelial_cell")] <- "Epithelial cell"
alltissue_fgseaRes$cell_name[which(alltissue_fgseaRes$cell_name == "Pericytes")] <- "Pericyte"
alltissue_fgseaRes$cell_name[which(alltissue_fgseaRes$cell_name == "Plasmocyte")] <- "Plasma cell"
alltissue_fgseaRes$sig=ifelse(alltissue_fgseaRes$pval<0.05,"Sig","NonSig")
cell_table <- as.data.frame(table(alltissue_fgseaRes$cell_name))
colnames(cell_table) <- c("cell_name","Freq")
cell_table <- cell_table[order(cell_table$Freq,decreasing = F),]
Sankey_data <- alltissue_fgseaRes_sig[,c(1,9)]
Sankey_data$tissue <- str_extract(Sankey_data$pathway, ".*(?=_M_)")

tissue1 <- c("Uterus","Prostate","Liver","Breast","Lung","Colon","Kidney","Skin","Ovary","Esophagus","Brain","Blood","Adrenal_Gland","Pancreas","Testis")
load(file = "./Aging/data/UCSC/GTEX_RNAseq_data/All_sig_module.Rdata")
all_tissue_module_info <- data_frame()
for (onetissue in tissue1) {
  a <- data.frame(module_name = names(All_sig_module[[paste0(onetissue,"_up")]]),
                  type=rep("Up",length(names(All_sig_module[[paste0(onetissue,"_up")]]))))
  b <- data.frame(module_name = names(All_sig_module[[paste0(onetissue,"_down")]]),
                  type=rep("Down",length(names(All_sig_module[[paste0(onetissue,"_down")]]))))
  one_tissue_module_info <- rbind(a,b)
  all_tissue_module_info <- rbind(all_tissue_module_info,one_tissue_module_info)
}
all_tissue_module_list <- list()
for (onetissue in tissue1) {
  one_tissue_module <- c(All_sig_module[[paste0(onetissue,"_up")]],All_sig_module[[paste0(onetissue,"_down")]])
  all_tissue_module_list <- c(all_tissue_module_list,one_tissue_module)
}
Sankey_data$moduleType <- all_tissue_module_info$type[match(Sankey_data$pathway,all_tissue_module_info$module_name)]
Sankey_dataNew <- Sankey_data[,c(3,4,2)]
data <- unique(Sankey_dataNew)
data$value <- 1
for (x in 1:length(data$tissue)) {
  one <- Sankey_dataNew[which(Sankey_dataNew$tissue == data$tissue[x]&Sankey_dataNew$moduleType == data$moduleType[x]&Sankey_dataNew$cell_name == data$cell_name[x]),]
  data$value[x] <- length(one$tissue)
}
save(Sankey_data,file="./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Tissuesstatistics/Snakey_data.rda")#Preservation of significantly enriched SACMs results





##--Fig.S5 A-----
cell <- c("T cell","NK cell","B cell","Plasma cell","Macrophage","Monocyte","Dendritic cell","Neutrophil","Mast cell","Hematopoietic stem cell",
          "Endothelial cell","Fibroblast","Epithelial cell","Smooth muscle cell","Basal cell","Pericyte",
          "Astrocytes","Microglia","Oligodendrocytes","Excitatory neurons","Granulosa cell",
          "Melanocyte","Follicular cell","Granule cell",
          "Germ cell","Hepatocyte","Enteric nerval cell",
          "Zona fasciculata cell","Principal cell")
data$tissue <- factor(data$tissue,levels = rev(tissue1))
data$cell_name <- factor(data$cell_name,levels = rev(cell))
data$moduleType <- factor(data$moduleType,levels = rev(c("Up","Down")))
library(ggplot2)
library(ggalluvial)
plot_theme <- theme_minimal() +
  theme(
    text = element_text(family = "sans",
                        face = "bold",
                        size = 14),
    panel.grid = element_blank())
gg1 <- ggplot(data, aes(y = value, axis1 =tissue , axis2 =moduleType ,axis3=cell_name)) +
  geom_alluvium(aes(fill = cell_name), curve_type = "spline", size = 1.2, width = 0.2) +
  geom_stratum(aes(fill = tissue), width = 0.2) +
  geom_stratum(aes(fill = moduleType), width = 0.2) +
  geom_stratum(aes(fill = cell_name), width = 0.2) +
  coord_flip()+
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 3,
            angle = 0) +
  scale_x_continuous(breaks = 1:3,
                     labels = c("ID", "Description", "Category"),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Class", y = "Value")
print(gg1 + plot_theme)

#-----Fig.2C-----top15 cell types-----
cell_tableplot <- cell_table[15:29,]

cell_tableplot$cell_name <- factor(cell_tableplot$cell_name,levels = cell_tableplot$cell_name)
ggplot(cell_tableplot,aes(cell_name,Freq,fill=cell_name))+
  coord_flip()+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c("#e8d738","#c62574","#34a5b1","#d4a2c5","#ac468c","#d1636a","#b31e23",
                               "#f1b3b9","#4558a6","#11398d","#5c308e","#977eb8","#d17613","#6f779f",
                               "#3c97c6"))+
  theme_bw()

#-----Fig.2D--6 cell types-----
paixu <- c("Uterus","Prostate","Liver","Breast","Lung","Colon","Kidney","Skin","Ovary",
           "Esophagus","Brain","Blood","Adrenal_Gland","Pancreas","Testis","Thyroid","Stomach")
load(file="./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Tissuesstatistics/Snakey_data.rda")
load(file="./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Tissuesstatistics/all_tissue_module_info.rda")
load(file="./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Tissuesstatistics/all_tissue_module_list.rda")
cell <- c("Endothelial cell","Epithelial cell","Fibroblast","T cell","Macrophage","Mast cell")
for (i in 1:length(cell)) {
  CellData <- Sankey_data[which(Sankey_data$cell_name == cell[i]),]
  cell_clust <- all_tissue_sigClust[match(CellData$pathway,all_tissue_sigClust$module.id),]
  cell_data <- CellData[,c(1,4,3)]
  tiqu2 <- function(x){
    module <- paste0("M",strsplit(as.character(x),split = "M")[[1]][2])
    return(module)
  }
  colnames(cell_data) <- c("module","type","tissue")
  cell_data$module_name <- apply(as.data.frame(as.character(cell_data$module)),1,tiqu2)
  cell_data$tissue <- factor(cell_data$tissue,levels = paixu)
  cell_data$type <- factor(cell_data$type,levels = c("Up","Down"))
  cell_data <- cell_data[order(cell_data$tissue,cell_data$type),]
  cell_data$module <- factor(cell_data$module,levels = cell_data$module)
  cell_freq <- Endothelial_cell_freq
  cell_freq$Var1 <- factor(cell_freq$Var1,levels = rev(paixu))
  cell_freq <- cell_freq[order(cell_freq$Var1),]
  cell_data$id <- seq(1, nrow(cell_data))
  angle <- 90 - 360 * (cell_data$id-0.5) /nrow(cell_data)
  cell_data$angle <- ifelse(angle < -90, angle+180, angle)
  ##middle ring
  ggplot(cell_data, aes(x = 3, y = module, fill = tissue)) +
    geom_col(position=position_dodge(0.7),width=1,color="black",size=0.25) +
    coord_polar(theta = "y") +
    geom_text(aes(label = module_name),position = position_stack(vjust = 0.5),size=2,angle=cell_data$angle) +
    scale_fill_manual(values = c(Adrenal_Gland="#8DD3C7",Prostate="#FFFFB3",Breast="#BEBADA",Kidney="#FB8072",
                                 Lung="#80B1D3",Ovary="#FDB462",Esophagus="#9fcba6",Liver="#E5C494",Brain="#B3DE69",
                                 Blood="#FCCDE5",Pancreas="#977EB8",Skin="#BC80BD",Uterus="#CCEBC5",Colon="#FFED6F",Testis="#5c308e")
    )+xlim(c(-15,3.5))+
    theme( panel.background = element_blank(),
           axis.title = element_blank(),
           axis.text = element_blank(),
           axis.ticks = element_blank(),
           legend.title = element_blank(),
           legend.position = "right")
  #outer ring
  ggplot(cell_data, aes(x = 3, y = module, fill = type)) +
    geom_col(position=position_dodge(0.7),width=1,color="black",size=0.25) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c(Up="#F20628",Down="#0FB0F7")
    )+xlim(c(-15,3.5))+
    theme( panel.background = element_blank(),
           axis.title = element_blank(),
           axis.text = element_blank(),
           axis.ticks = element_blank(),
           legend.title = element_blank(),
           legend.position = "right")
  #inner ring
  ggplot(cell_freq, aes(x = 3, y = Freq, fill = Var1)) +
    geom_col() +
    geom_text(aes(label = paste0(Var1,"\n","(",Freq,")")),
              position = position_stack(vjust = 0.5),size=3) +
    coord_polar(theta = "y") +
    xlim(c(0.2, 3.5))+
    scale_fill_manual(values = c(Adrenal_Gland="#8DD3C7",Prostate="#FFFFB3",Breast="#BEBADA",Kidney="#FB8072",
                                 Lung="#80B1D3",Ovary="#FDB462",Esophagus="#9fcba6",Liver="#E5C494",Brain="#B3DE69",
                                 Blood="#FCCDE5",Pancreas="#977EB8",Skin="#BC80BD",Uterus="#CCEBC5",Colon="#FFED6F",Testis="#5c308e")
    )+
    theme( panel.background = element_blank(),
           axis.title = element_blank(),
           axis.text = element_blank(),
           axis.ticks = element_blank(),
           legend.title = element_blank(),
           legend.position = "right")
}






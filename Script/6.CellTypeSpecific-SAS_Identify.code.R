library(ggplot2)
#############################################################
######Sharing of SAS from different cells between tissues####
#############################################################
#--Fig.S5B-----
####1.Endothelial cells####
load(file="./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/组织统计/Snakey_data.rda")
Endothelial_cell <- Sankey_data[which(Sankey_data$cell_name == "Endothelial cell"),]
#Up
Endothelial_cell_Up <- Endothelial_cell[which(Endothelial_cell$moduleType== "Up"),]
EndCell_Up <- all_tissue_module_list[Endothelial_cell$pathway[which(Endothelial_cell$moduleType== "Up")]]
Endothelial_cell_all_tissue <- list()
for (onetissue in unique(Endothelial_cell_Up$tissue)) {
  Endothelial_cell_one <- Endothelial_cell_Up[which(Endothelial_cell_Up$tissue == onetissue),]
  Endothelial_cell_one_tissue <- list(unique(unlist(EndCell_Up[Endothelial_cell_one$pathway])))
  names(Endothelial_cell_one_tissue) <- onetissue
  Endothelial_cell_all_tissue <- c(Endothelial_cell_all_tissue,Endothelial_cell_one_tissue)
}
EndCell_Up_total <- as.data.frame(table(as.data.frame(unlist(Endothelial_cell_all_tissue))))
EndCell_Up_More2 <- as.character(EndCell_Up_total[which(EndCell_Up_total$Freq>=2),]$Var1)
EndCell_Up_total_freq <- as.data.frame(table(EndCell_Up_total$Freq))
colnames(EndCell_Up_total_freq) <- c("pinshu","cishu")
#Down
Endothelial_cell_Down <- Endothelial_cell[which(Endothelial_cell$moduleType== "Down"),]
EndCell_Down <- all_tissue_module_list[Endothelial_cell$pathway[which(Endothelial_cell$moduleType== "Down")]]
Endothelial_cell_all_tissue <- list()
for (onetissue in unique(Endothelial_cell_Down$tissue)) {
  Endothelial_cell_one <- Endothelial_cell_Down[which(Endothelial_cell_Down$tissue == onetissue),]
  Endothelial_cell_one_tissue <- list(unique(unlist(EndCell_Down[Endothelial_cell_one$pathway])))
  names(Endothelial_cell_one_tissue) <- onetissue
  Endothelial_cell_all_tissue <- c(Endothelial_cell_all_tissue,Endothelial_cell_one_tissue)
}

EndCell_Down_total <- as.data.frame(table(as.data.frame(unlist(Endothelial_cell_all_tissue))))
EndCell_Down_More2 <- as.character(EndCell_Down_total[which(EndCell_Down_total$Freq>=2),]$Var1)
EndCell_Down_total_freq <- as.data.frame(table(EndCell_Down_total$Freq))
colnames(EndCell_Down_total_freq) <- c("pinshu","cishu")
#两组放在一起画
two_group_data <- cbind(EndCell_Up_total_freq,EndCell_Down_total_freq[,2])
colnames(two_group_data) <- c("tissue_num","Up","Down")
df1 <- reshape2::melt(two_group_data,id.vars = 'tissue_num')
df1$gene <- c(434,152,25,5,274,47,9,0)
ggplot(df1, aes(
  x = factor(tissue_num,levels = unique(tissue_num)),
  y = ifelse(variable == "Up", value, -value),
  fill = variable)) +
  geom_col(width = .85,alpha=0.5)+
  geom_col(
    aes(x = factor(tissue_num,levels = unique(tissue_num)),
        y = ifelse(variable == "Up", gene, -gene),
        fill = variable), width = .5 
  )+
  scale_fill_manual(values = c("#EE3432","#1B64A4"))+   
  geom_text(                                                  
    aes(label=value,                                        
    ),
    size=3)+
  ylab("Number of gene")+xlab("Number of tissues with sharing genes")+
  scale_y_continuous(                                       
    labels = abs,                            
    expand = expansion(mult = c(0.1, 0.1)))+
  theme_bw()+
  theme(panel.grid=element_blank (),legend.position = "none")


#---2.T_cell----
T_cell <- Sankey_data[which(Sankey_data$cell_name == "T cell"),]
#Up
T_cell_Up <- T_cell[which(T_cell$moduleType== "Up"),]
TCell_Up <- all_tissue_module_list[T_cell$pathway[which(T_cell$moduleType== "Up")]]
T_cell_all_tissue <- list()
for (onetissue in unique(T_cell_Up$tissue)) {
  T_cell_one <- T_cell_Up[which(T_cell_Up$tissue == onetissue),]
  T_cell_one_tissue <- list(unique(unlist(TCell_Up[T_cell_one$pathway])))
  names(T_cell_one_tissue) <- onetissue
  T_cell_all_tissue <- c(T_cell_all_tissue,T_cell_one_tissue)
}
TCell_Up_total <- as.data.frame(table(as.data.frame(unlist(T_cell_all_tissue))))
TCell_Up_More2 <- TCell_Up_total[which(TCell_Up_total$Freq>=2),]
TCell_Up_total_freq <- as.data.frame(table(TCell_Up_total$Freq))
colnames(TCell_Up_total_freq) <- c("pinshu","cishu")
#Down
T_cell_Down <- T_cell[which(T_cell$moduleType== "Down"),]
TCell_Down <- all_tissue_module_list[T_cell$pathway[which(T_cell$moduleType== "Down")]]
T_cell_all_tissue <- list()
for (onetissue in unique(T_cell_Down$tissue)) {
  T_cell_one <- T_cell_Down[which(T_cell_Down$tissue == onetissue),]
  T_cell_one_tissue <- list(unique(unlist(TCell_Down[T_cell_one$pathway])))
  names(T_cell_one_tissue) <- onetissue
  T_cell_all_tissue <- c(T_cell_all_tissue,T_cell_one_tissue)
}
TCell_Down_total <- as.data.frame(table(as.data.frame(unlist(T_cell_all_tissue))))
TCell_Down_More2 <- TCell_Down_total[which(TCell_Down_total$Freq>=2),]
TCell_Down_total_freq <- as.data.frame(table(TCell_Down_total$Freq))
colnames(TCell_Down_total_freq) <- c("pinshu","cishu")

two_group_data <- cbind(TCell_Up_total_freq,TCell_Down_total_freq[,2])
colnames(two_group_data) <- c("tissue_num","Up","Down")
df1 <- reshape2::melt(two_group_data,id.vars = 'tissue_num')
ggplot(df1, aes(
  x = factor(tissue_num,levels = unique(tissue_num)), 
  y = ifelse(variable == "Up", value, -value),
  fill = variable)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = c("#EE3432","#1B64A4"))+
  geom_text(
    aes(label=value,
    ),
    size=3)+
  ylab("Number of gene")+xlab("Number of tissues with sharing genes")+
  scale_y_continuous(
    labels = abs,
    expand = expansion(mult = c(0.1, 0.1)))+
  theme_bw()+
  theme(panel.grid=element_blank ())


#---3.Epithelial cell----
Epithelial_cell <- Sankey_data[which(Sankey_data$cell_name == "Epithelial cell"),]
Epithelial_cell_Up <- Epithelial_cell[which(Epithelial_cell$moduleType== "Up"),]
EpiCell_Up <- all_tissue_module_list[Epithelial_cell$pathway[which(Epithelial_cell$moduleType== "Up")]]
Epithelial_cell_all_tissue <- list()
for (onetissue in unique(Epithelial_cell_Up$tissue)) {
  Epithelial_cell_one <- Epithelial_cell_Up[which(Epithelial_cell_Up$tissue == onetissue),]
  Epithelial_cell_one_tissue <- list(unique(unlist(EpiCell_Up[Epithelial_cell_one$pathway])))
  names(Epithelial_cell_one_tissue) <- onetissue
  Epithelial_cell_all_tissue <- c(Epithelial_cell_all_tissue,Epithelial_cell_one_tissue)
}
EpiCell_Up_total <- as.data.frame(table(as.data.frame(unlist(Epithelial_cell_all_tissue))))
EpiCell_Up_total_freq <- as.data.frame(table(EpiCell_Up_total$Freq))
colnames(EpiCell_Up_total_freq) <- c("pinshu","cishu")
#Down
Epithelial_cell_Down <- Epithelial_cell[which(Epithelial_cell$moduleType== "Down"),]
EpiCell_Down <- all_tissue_module_list[Epithelial_cell$pathway[which(Epithelial_cell$moduleType== "Down")]]
Epithelial_cell_all_tissue <- list()
for (onetissue in unique(Epithelial_cell_Down$tissue)) {
  Epithelial_cell_one <- Epithelial_cell_Down[which(Epithelial_cell_Down$tissue == onetissue),]
  Epithelial_cell_one_tissue <- list(unique(unlist(EpiCell_Down[Epithelial_cell_one$pathway])))
  names(Epithelial_cell_one_tissue) <- onetissue
  Epithelial_cell_all_tissue <- c(Epithelial_cell_all_tissue,Epithelial_cell_one_tissue)
}
EpiCell_Down_total <- as.data.frame(table(as.data.frame(unlist(Epithelial_cell_all_tissue))))
EpiCell_Down_More2 <- as.character(EpiCell_Down_total[which(EpiCell_Down_total$Freq>=2),])
EpiCell_Down_total_freq <- as.data.frame(table(EpiCell_Down_total$Freq))
colnames(EpiCell_Down_total_freq) <- c("pinshu","cishu")
#两组放在一起画
EpiCell_Up_total_freq <- rbind(EpiCell_Up_total_freq,c("1",0),c("1",0))
two_group_data <- cbind(EpiCell_Down_total_freq,as.numeric(EpiCell_Up_total_freq[,2]))
colnames(two_group_data) <- c("tissue_num","Down","Up")
df1 <- reshape2::melt(two_group_data,id.vars = 'tissue_num')
ggplot(df1, aes(
  x = factor(tissue_num,levels = unique(tissue_num)),             
  y = ifelse(variable == "Up", value, -value),
  fill = variable)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = c("#1B64A4","#EE3432"))+
  geom_text(
    aes(label=value,
    ),
    size=3)+
  ylab("Number of gene")+xlab("Number of tissues with sharing genes")+
  scale_y_continuous(                              
    labels = abs,                                    
    expand = expansion(mult = c(0.1, 0.1)))+
  theme_bw()+
  theme(panel.grid=element_blank ())

#---4.Fibroblast----
Fibroblast <- Sankey_data[which(Sankey_data$cell_name == "Fibroblast"),]
Fibroblast_Up <- Fibroblast[which(Fibroblast$moduleType== "Up"),]
FibroCell_Up <- all_tissue_module_list[Fibroblast$pathway[which(Fibroblast$moduleType== "Up")]]
Fibroblast_all_tissue <- list()
for (onetissue in unique(Fibroblast_Up$tissue)) {
  Fibroblast_one <- Fibroblast_Up[which(Fibroblast_Up$tissue == onetissue),]
  Fibroblast_one_tissue <- list(unique(unlist(FibroCell_Up[Fibroblast_one$pathway])))
  names(Fibroblast_one_tissue) <- onetissue
  Fibroblast_all_tissue <- c(Fibroblast_all_tissue,Fibroblast_one_tissue)
}
FibroCell_Up_total <- as.data.frame(table(as.data.frame(unlist(Fibroblast_all_tissue))))
FibroCell_Up_More2 <- as.character(FibroCell_Up_total[which(FibroCell_Up_total$Freq>=2),]$Var1)
FibroCell_Up_total_freq <- as.data.frame(table(FibroCell_Up_total$Freq))
colnames(FibroCell_Up_total_freq) <- c("pinshu","cishu")
#Down
Fibroblast_Down <- Fibroblast[which(Fibroblast$moduleType== "Down"),]
FibroCell_Down <- all_tissue_module_list[Fibroblast$pathway[which(Fibroblast$moduleType== "Down")]]
Fibroblast_all_tissue <- list()
for (onetissue in unique(Fibroblast_Down$tissue)) {
  Fibroblast_one <- Fibroblast_Down[which(Fibroblast_Down$tissue == onetissue),]
  Fibroblast_one_tissue <- list(unique(unlist(FibroCell_Down[Fibroblast_one$pathway])))
  names(Fibroblast_one_tissue) <- onetissue
  Fibroblast_all_tissue <- c(Fibroblast_all_tissue,Fibroblast_one_tissue)
}

FibroCell_Down_total <- as.data.frame(table(as.data.frame(unlist(Fibroblast_all_tissue))))
FibroCell_Down_More2 <- as.character(FibroCell_Down_total[which(FibroCell_Down_total$Freq>=2),]$Var1)
FibroCell_Down_total_freq <- as.data.frame(table(FibroCell_Down_total$Freq))
colnames(FibroCell_Down_total_freq) <- c("pinshu","cishu")
#两组放在一起画
FibroCell_Down_total_freq <-rbind(FibroCell_Down_total_freq,c("4",1))
two_group_data <- cbind(FibroCell_Up_total_freq,as.numeric(FibroCell_Down_total_freq[,2]))
colnames(two_group_data) <- c("tissue_num","Up","Down")
df1 <- reshape2::melt(two_group_data,id.vars = 'tissue_num')
ggplot(df1, aes(
  x = factor(tissue_num,levels = unique(tissue_num)),           
  y = ifelse(variable == "Up", value, -value),
  fill = variable)) +
  geom_bar(stat = 'identity')+ 
  scale_fill_manual(values = c("#EE3432","#1B64A4"))+
  geom_text(
    aes(label=value,
    ),
    size=3)+
  ylab("Number of gene")+xlab("Number of tissues with sharing genes")+
  scale_y_continuous(
    labels = abs,                                            
    expand = expansion(mult = c(0.1, 0.1)))+               
  theme_bw()+
  theme(panel.grid=element_blank ())

#---5.Macrophage----
Macrophage <- Sankey_data[which(Sankey_data$cell_name == "Macrophage"),]
#Up
Macrophage_Up <- Macrophage[which(Macrophage$moduleType== "Up"),]
MacroCell_Up <- all_tissue_module_list[Macrophage$pathway[which(Macrophage$moduleType== "Up")]]
Macrophage_all_tissue <- list()
for (onetissue in unique(Macrophage_Up$tissue)) {
  Macrophage_one <- Macrophage_Up[which(Macrophage_Up$tissue == onetissue),]
  Macrophage_one_tissue <- list(unique(unlist(MacroCell_Up[Macrophage_one$pathway])))
  names(Macrophage_one_tissue) <- onetissue
  Macrophage_all_tissue <- c(Macrophage_all_tissue,Macrophage_one_tissue)
}
MacroCell_Up_total <- as.data.frame(table(as.data.frame(unlist(Macrophage_all_tissue))))
MacroCell_Up_More2 <- as.character(MacroCell_Up_total[which(MacroCell_Up_total$Freq>=2),]$Var1)
MacroCell_Up_total_freq <- as.data.frame(table(MacroCell_Up_total$Freq))
colnames(MacroCell_Up_total_freq) <- c("pinshu","cishu")

#Down
Macrophage_Down <- Macrophage[which(Macrophage$moduleType== "Down"),]
MacroCell_Down <- all_tissue_module_list[Macrophage$pathway[which(Macrophage$moduleType== "Down")]]
Macrophage_all_tissue <- list()
for (onetissue in unique(Macrophage_Down$tissue)) {
  Macrophage_one <- Macrophage_Down[which(Macrophage_Down$tissue == onetissue),]
  Macrophage_one_tissue <- list(unique(unlist(MacroCell_Down[Macrophage_one$pathway])))
  names(Macrophage_one_tissue) <- onetissue
  Macrophage_all_tissue <- c(Macrophage_all_tissue,Macrophage_one_tissue)
}
MacroCell_Down_total <- as.data.frame(table(as.data.frame(unlist(Macrophage_all_tissue))))
MacroCell_Down_More2 <- as.character(MacroCell_Down_total[which(MacroCell_Down_total$Freq>=2),]$Var1)
MacroCell_Down_total_freq <- as.data.frame(table(MacroCell_Down_total$Freq))
colnames(MacroCell_Down_total_freq) <- c("pinshu","cishu")
#两组放在一起画
MacroCell_Up_total_freq <- rbind(MacroCell_Up_total_freq,c("1",0))
two_group_data <- cbind(MacroCell_Down_total_freq,as.numeric(MacroCell_Up_total_freq[,2]))
colnames(two_group_data) <- c("tissue_num","Down","Up")
df1 <- reshape2::melt(two_group_data,id.vars = 'tissue_num')
ggplot(df1, aes(
  x = factor(tissue_num,levels = unique(tissue_num)),            
  y = ifelse(variable == "Up", value, -value), 
  fill = variable)) +
  geom_bar(stat = 'identity')+                               
  scale_fill_manual(values = c("#1B64A4","#EE3432"))+  
  geom_text(                                                 
    aes(label=value,                                         
    ),
    size=3)+
  ylab("Number of gene")+xlab("Number of tissues with sharing genes")+
  scale_y_continuous(                                     
    labels = abs,                                           
    expand = expansion(mult = c(0.1, 0.1)))+ 
  theme_bw()+
  theme(panel.grid=element_blank ())

#---6.Mast cell----
Mast_cell <- Sankey_data[which(Sankey_data$cell_name == "Mast cell"),]
#Up
Mast_cell_Up <- Mast_cell[which(Mast_cell$moduleType== "Up"),]
MastCell_Up <- all_tissue_module_list[Mast_cell$pathway[which(Mast_cell$moduleType== "Up")]]
Mast_cell_all_tissue <- list()
for (onetissue in unique(Mast_cell_Up$tissue)) {
  Mast_cell_one <- Mast_cell_Up[which(Mast_cell_Up$tissue == onetissue),]
  Mast_cell_one_tissue <- list(unique(unlist(MastCell_Up[Mast_cell_one$pathway])))
  names(Mast_cell_one_tissue) <- onetissue
  Mast_cell_all_tissue <- c(Mast_cell_all_tissue,Mast_cell_one_tissue)
}
MastCell_Up_total <- as.data.frame(table(as.data.frame(unlist(Mast_cell_all_tissue))))
MastCell_Up_total_freq <- as.data.frame(table(MastCell_Up_total$Freq))
colnames(MastCell_Up_total_freq) <- c("pinshu","cishu")
#Down
Mast_cell_Down <- Mast_cell[which(Mast_cell$moduleType== "Down"),]
MastCell_Down <- all_tissue_module_list[Mast_cell$pathway[which(Mast_cell$moduleType== "Down")]]
Mast_cell_all_tissue <- list()
for (onetissue in unique(Mast_cell_Down$tissue)) {
  Mast_cell_one <- Mast_cell_Down[which(Mast_cell_Down$tissue == onetissue),]
  Mast_cell_one_tissue <- list(unique(unlist(MastCell_Down[Mast_cell_one$pathway])))
  names(Mast_cell_one_tissue) <- onetissue
  Mast_cell_all_tissue <- c(Mast_cell_all_tissue,Mast_cell_one_tissue)
}
MastCell_Down_total <- as.data.frame(table(as.data.frame(unlist(Mast_cell_all_tissue))))
MastCell_Down_More2 <- as.character(MastCell_Down_total[which(MastCell_Down_total$Freq>=2),]$Var1)
MastCell_Down_total_freq <- as.data.frame(table(MastCell_Down_total$Freq))
colnames(MastCell_Down_total_freq) <- c("pinshu","cishu")

#两组放在一起画
MastCell_Up_total_freq <- rbind(MastCell_Up_total_freq,c("1",0),c("1",0))
two_group_data <- cbind(MastCell_Down_total_freq,as.numeric(MastCell_Up_total_freq[,2]))
colnames(two_group_data) <- c("tissue_num","Down","Up")
df1 <- reshape2::melt(two_group_data,id.vars = 'tissue_num')
ggplot(df1, aes(
  x = factor(tissue_num,levels = unique(tissue_num)),           
  y = ifelse(variable == "Up", value, -value),
  fill = variable)) +
  geom_bar(stat = 'identity')+ 
  scale_fill_manual(values = c("#1B64A4","#EE3432"))+ 
  geom_text(                                             
    aes(label=value,                                    
    ),
    size=3)+
  ylab("Number of gene")+xlab("Number of tissues with sharing genes")+
  scale_y_continuous(                                    
    labels = abs,                                            
    expand = expansion(mult = c(0.1, 0.1)))+ 
  theme_bw()+
  theme(panel.grid=element_blank ())

#---Fig.2E-----
setwd("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Tissuesstatistics/CellTypeEnrich")
#Up and Down regulated
for (cell in c("Endothelial","Fibroblast","Tcell")) {
  Cell_up <- read.csv(paste0(cell,"_Up_Enrich.csv"),header = T)
  Cell_down <- read.csv(paste0(cell,"_Down_Enrich.csv"),header = T)
  Cell_data <- rbind(Cell_up[,c(4,5,6)],Cell_down[,c(4,5,6)])
  Cell_data$up_down <- c(rep("Up",length(rownames(Cell_up))),rep("Down",length(rownames(Cell_down))))
  Cell_data <- Cell_data %>% arrange(up_down,desc(LogP))
  Cell_data$LogP <- ifelse(Cell_data$up_down == "Down",Cell_data$LogP,-Cell_data$LogP)
  Cell_data$Description <- factor(Cell_data$Description,Cell_data$Description)
  x_max = max(abs(Cell_data$LogP))
  color_list <- c("#B0CAE6","#E5B3B3")
  
  p <- ggplot(Endothelial,aes(x =LogP, y = Description, fill = up_down)) + 
    geom_col() + 
    scale_x_continuous(limits = c(-x_max,x_max)) +
    scale_fill_manual(values = color_list) +
    geom_text(data = Endothelial[which(Endothelial$up_down == "Up"),],
              aes(x = -0.5, y = Description, label = Description,color = up_down),
              size = 4,
              hjust = 1) + 
    geom_text(data = Endothelial[which(Endothelial$up_down == "Down"),],
              aes(x = 0.5, y = Description, label = Description,color = up_down),
              size = 4,
              hjust = 0)  + 
    scale_colour_manual(values = color_list) +
    labs(x="log(Pvalue)", y=NULL, title="enriched pathway") +
    annotate("text", x = 25, y = 15, label = "Up", size = 10, fontface = "bold", color="#E5B3B3") +
    annotate("text", x = -25, y =5, label = "Down", size =10, fontface = "bold", color="#B0CAE6") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 25,hjust=0.5), 
      axis.text.x = element_text(size = 15),
      axis.title.x = element_text(size = 15),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = 'none'
    )
  pdf(paste0("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Tissuesstatistics/CellTypeEnrich",cell,"_enrichplot.pdf"),width=8,height = 6)
  print(p)
  dev.off()
}
#Down regulated
for (cell in c("Epithelial","MastCell","Macrophage")) {
  Cell_data <- read.csv(paste0(cell,"_Down_Enrich.csv"),header = T)
  Cell_data$up_down <- c(rep("Down",length(rownames(Cell_up))))
  Cell_data <- Cell_data %>% arrange(up_down,desc(LogP))
  Cell_data$LogP <- ifelse(Cell_data$up_down == "Down",Cell_data$LogP,-Cell_data$LogP)
  Cell_data$Description <- factor(Cell_data$Description,Cell_data$Description)
  x_max = max(abs(Cell_data$LogP))
  color_list <- c("#B0CAE6","#E5B3B3")
  
  p <- ggplot(Endothelial,aes(x =LogP, y = Description, fill = up_down)) + 
    geom_col() + 
    scale_x_continuous(limits = c(-x_max,x_max)) +
    scale_fill_manual(values = color_list) +
    geom_text(data = Endothelial[which(Endothelial$up_down == "Up"),],
              aes(x = -0.5, y = Description, label = Description,color = up_down),
              size = 4,
              hjust = 1) + 
    geom_text(data = Endothelial[which(Endothelial$up_down == "Down"),],
              aes(x = 0.5, y = Description, label = Description,color = up_down),
              size = 4,
              hjust = 0)  + 
    scale_colour_manual(values = color_list) +
    labs(x="log(Pvalue)", y=NULL, title="enriched pathway") +
    annotate("text", x = 25, y = 15, label = "Up", size = 10, fontface = "bold", color="#E5B3B3") +
    annotate("text", x = -25, y =5, label = "Down", size =10, fontface = "bold", color="#B0CAE6") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 25,hjust=0.5), 
      axis.text.x = element_text(size = 15),
      axis.title.x = element_text(size = 15),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = 'none'
    )
  pdf(paste0("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Tissuesstatistics/CellTypeEnrich",cell,"_enrichplot.pdf"),width=8,height = 6)
  print(p)
  dev.off()
}
##############################################################################################################
######Identification of SACMs——Overlap and screen for age-related modules and senescence characteristics######
##############################################################################################################
library(GeneOverlap)
Tissues <- c("Adrenal_Gland","Breast","Esophagus","Brain","Kidney","Prostate","Liver","Lung","Ovary","Pancreas","Colon",
             "Testis","Uterus","Skin","Stomach","Thyroid","Blood")
Tissues_proflie <- c("Adrenal_Gland_ACC","Breast_BRCA","Esophagus_ESCA","Brain_LGG_GBM","Kidney_KICH_KIRC_KIRP","Prostate_PRAD",
                     "Liver_LIHC","Lung_LUAD_LUSC","Ovary_OV","Pancreas_PAAD","Colon_COAD_READ","Testis_TGCT","Uterus_UCEC_UCS","Skin_SKCM",
                     "Stomach_STAD","Thyroid_THCA","Blood_DLBC_THYM")
load(file = "./Aging/data/aging_gene/aging_gene_list_new.Rdata")
background_geneNumber<- as.data.frame(read.csv("./Aging\/data/UCSC/GTEX_RNAseq_data/background_geneNumber.csv",header = T,sep = ","))
All_tissue_pvalue <- data.frame()
All_tissue_inter_num <- data.frame()
All_module_info <- data.frame()
All_sig_module <- list()
for (i in 1:17) {
  load(paste0("./Aging/data/UCSC/GTEX_RNAseq_data/",Tissues_proflie[i],"/",Tissues[i],"_AgingRelatedModuleList.rda"))
  All_module <- get(paste0(Tissues[i],"_AgingRelatedModuleList"))
  names(All_module[[paste0(Tissues[i],"_up_modules")]]) <- gsub("_c._", "_M_", names(All_module[[paste0(Tissues[i],"_up_modules")]]))
  names(All_module[[paste0(Tissues[i],"_down_modules")]]) <- gsub("_c._", "_M_", names(All_module[[paste0(Tissues[i],"_down_modules")]]))
  
  tissue <- Tissues[i]
  go.obj <- newGOM(All_module[[paste0(tissue,"_up_modules")]],all_tissue_aging_gene_list[[tissue]],genome.size = background_geneNumber[i,2])
  Pvalue_matrix <- as.data.frame(getMatrix(go.obj, name="pval"))
  Pvalue_matrix_sig <- Pvalue_matrix[apply(Pvalue_matrix, 1, function(x) any(x < 0.05)), ]
  inter.num <- as.data.frame(getMatrix(go.obj, name="intersection"))
  inter.num_sig <- inter.num[rownames(Pvalue_matrix_sig),]
  
  Pvalue_matrix_sig$module <- rownames(Pvalue_matrix_sig)
  Pvalue_matrix_sig_gp <- reshape2::melt(Pvalue_matrix_sig)
  
  inter.num_sig$module <- rownames(inter.num_sig)
  inter.num_sig_gp <- reshape2::melt(inter.num_sig)
  
  sigmodule_Up_info_ggplot <- merge(Pvalue_matrix_sig_gp,inter.num_sig_gp,by = c("module","variable"))
  colnames(sigmodule_Up_info_ggplot) <- c("module","variable","Pvalue","inter_num")
  sigmodule_Up_info_more2 <- sigmodule_Up_info_ggplot[which(sigmodule_Up_info_ggplot$Pvalue <0.05 & sigmodule_Up_info_ggplot$inter_num >=2),]
  sigmodule_Up_name <- unique(sigmodule_Up_info_more2$module)
  if(length(sigmodule_Up_name)!=0){
    Up_module_info <- data.frame(module=sigmodule_Up_name,type = "Up")
  }else{
    Up_module_info <- data.frame()
  }
  tissue_up_sig <- All_module[[paste0(tissue,"_up")]][sigmodule_Up_name]
  Sig_Pvalue_matrix <- Pvalue_matrix[sigmodule_Up_name,]
  Sig_inter_num <- inter.num[sigmodule_Up_name,]
  
  go.obj1 <- newGOM(All_module[[paste0(tissue,"_down_modules")]],all_tissue_aging_gene_list[[tissue]],genome.size = tumor_tissue[i,2])
  Pvalue_matrix1 <- as.data.frame(getMatrix(go.obj1, name="pval"))
  Pvalue_matrix1_sig <- Pvalue_matrix1[apply(Pvalue_matrix1, 1, function(x) any(x < 0.05)), ]
  inter.num1 <- as.data.frame(getMatrix(go.obj1, name="intersection"))
  inter.num1_sig <- inter.num1[rownames(Pvalue_matrix1_sig),]
  
  Pvalue_matrix1_sig$module <- rownames(Pvalue_matrix1_sig)
  Pvalue_matrix1_sig_gp <- reshape2::melt(Pvalue_matrix1_sig)
  
  inter.num1_sig$module <- rownames(inter.num1_sig)
  inter.num1_sig_gp <- reshape2::melt(inter.num1_sig)
  
  sigmodule_Down_info_ggplot <- merge(Pvalue_matrix1_sig_gp,inter.num1_sig_gp,by = c("module","variable"))
  colnames(sigmodule_Down_info_ggplot) <- c("module","variable","Pvalue","inter_num")
  sigmodule_Down_info_more2 <- sigmodule_Down_info_ggplot[which(sigmodule_Down_info_ggplot$Pvalue <0.05 & sigmodule_Down_info_ggplot$inter_num >=2),]
  sigmodule_Down_name <- unique(sigmodule_Down_info_more2$module)
  if(length(sigmodule_Down_name)!=0){
    Down_module_info <- data.frame(module=sigmodule_Down_name,type = "Down")
  }else{
    Down_module_info <- data.frame()
  }
  tissue_down_sig <- All_module[[paste0(tissue,"_down")]][sigmodule_Down_name]
  Sig_Pvalue_matrix1 <- Pvalue_matrix1[sigmodule_Down_name,]
  Sig_inter_num1 <- inter.num1[sigmodule_Down_name,]
  tissue_Pvalue <- rbind(Sig_Pvalue_matrix,Sig_Pvalue_matrix1)
  tissue_inter_num <- rbind(Sig_inter_num,Sig_inter_num1)
  tissue_info <- rbind(Up_module_info,Down_module_info)
  
  
  All_tissue_pvalue <- rbind(All_tissue_pvalue,tissue_Pvalue)
  All_tissue_inter_num <- rbind(All_tissue_inter_num,tissue_inter_num)
  All_module_info <- rbind(All_module_info,tissue_info)
  tissue_up_sig <- list(tissue_up_sig)
  names(tissue_up_sig) <- paste0(tissue,"_up")
  
  tissue_down_sig <- list(tissue_down_sig)
  names(tissue_down_sig) <- paste0(tissue,"_down")
  
  one_tissue_sig_module <- c(tissue_up_sig,tissue_down_sig)
  All_sig_module <- c(All_sig_module,one_tissue_sig_module)
}
save(All_sig_module,file = "./Aging/data/UCSC/GTEX_RNAseq_data/Screening for Senescence Characteristics/All_sig_module.Rdata")
save(All_tissue_pvalue,file = "./Aging/data/UCSC/GTEX_RNAseq_data/Screening for Senescence Characteristics/All_tissue_pvalue.Rdata")
save(All_tissue_inter_num,file = "./Aging/data/UCSC/GTEX_RNAseq_data/Screening for Senescence Characteristics/All_tissue_inter_num.Rdata")
save(All_module_info,file = "./Aging/data/UCSC/GTEX_RNAseq_data/Screening for Senescence Characteristics/All_module_info.Rdata")

#----Figure S2A--------
library(aplot)
library(ggplot2)
library(RColorBrewer)
All_tissue_pvalue$module <- rownames(All_tissue_pvalue)
All_tissue_inter_num$module <- rownames(All_tissue_inter_num)
All_tissue_pvalue_plot <- reshape2::melt(All_tissue_pvalue)
All_tissue_inter_num_plot <- reshape2::melt(All_tissue_inter_num)
All_data <- merge(All_tissue_pvalue_plot,All_tissue_inter_num_plot,by = c("module","variable"))
colnames(All_data) <- c("module","variable","Pvalue","inter_num")
All_module_info$p <- ""
group <- ggplot(All_module_info,aes(p,module,fill=type))+
  geom_tile() +
  scale_fill_manual(values = c("#007cc0","#ce181e"))+
  scale_y_discrete(position="right") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x =element_blank())+
  labs(fill = "module_type")
All_module_info$tissue <- apply(as.data.frame(All_module_info$module),1,function(x) strsplit(x,"_")[[1]][1])
TYPE <- ggplot(All_module_info,aes(p,module,fill=tissue))+
  geom_tile() +
  scale_fill_manual(values = c("Adrenal"="#e8d738","Blood"="#c62574","Brain"="#34a5b1",
                               "Colon"="#d4a2c5","Esophagus"="#ac468c" ,"Kidney"="#d1636a",
                               "Liver"="#b31e23","Lung"="#f1b3b9","Ovary"="#4558a6",
                               "Pancreas"="#11398d","Prostate"="#5c308e","Skin"="#977eb8",
                               "Stomach"="#d17613","Testis"="#3c97c6","Uterus"="#a6cd3f","Breast"="#A6CEE3",Thyroid="#E7298A")) +
  scale_y_discrete(position="right") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x =element_blank())+
  labs(fill = "tissue")
All_data$Pvalue <- -log10(All_data$Pvalue)
All_data$module <- factor(All_data$module,levels = rev(All_module_info$module))
P1 <- ggplot(All_data,aes(x=variable,y=module,fill=Pvalue))+
  geom_raster()+
  theme_minimal()+
  scale_fill_distiller(palette = "OrRd",limits = c(0, 1), breaks = c(0, 0.05,1))+
  ylab(NULL)+xlab(NULL) +
  theme(panel.background = element_blank(),
        legend.key = element_blank(),
        axis.text = element_text(color = "black",size = 10),
        axis.text.x = element_text(angle = 45,hjust = 1),#x轴文字角度
        axis.text.y = element_blank(),
        panel.grid = element_blank()
  )
P1%>%insert_left(group,width = .02)%>%
  insert_left(TYPE,width = .02)




###########################################################
########Bubble plot of SACMs versus age correlation########
###########################################################
load("./Aging/data/UCSC/GTEX_RNAseq_data//Screening for Senescence Characteristics/ALLmodule_SigClust.rda")
load("./Aging/data/UCSC/GTEX_RNAseq_data//Screening for Senescence Characteristics/All_module_info.Rdata")
library(ggrepel)
all_tissue_sigClust <- data.frame()
for (x in 1:17) {
  setwd(paste0("./Aging/data/UCSC/GTEX_RNAseq_data/",Tissues_proflie[x]))
  one_tissue_sigClust <- read.table(file="sigCluster.txt",sep="\t",header = T)
  one_tissue_sigClust$module.id <- paste0(Tissues[i],"_",gsub("c._", "M_", one_tissue_sigClust$module.id))
  all_tissue_sigClust <- rbind(all_tissue_sigClust,one_tissue_sigClust)
}
All_module_info_clust <- all_tissue_sigClust[match(All_module_info$module,all_tissue_sigClust$module.id,nomatch = 0),]
All_module_info_clust$logP <- -log10(All_module_info_clust$Pvalue)
All_module_info_clust$module.size <- sapply(All_module_info_clust$module.id, function(x) length(All_module_list[[x]]) )
All_module_info_clust <- all_tissue_sigClust[match(All_module_info$module,all_tissue_sigClust$module.id,nomatch = 0),]
data <- All_module_info_clust[,c(1,2,6,7)]
data$tissue <- sapply(data$module.id,function(x) strsplit(x, "_")[[1]][1])
data$module.size <- log10(data$module.size)
#----Figure S2B--------
ggplot(data, aes(x =module.size ,y = Cor ,size = logP,color = logP)) +
  geom_point(shape = 16,alpha=0.8) +
  scale_size(range = c(0, 15))+ 
  scale_color_distiller(palette = "YlOrRd", limits = c(0, 12),direction = 1) + 
  geom_text_repel(data = data  %>% 
                    dplyr::arrange(desc(Cor)) %>% head(5),
                  aes(x = module.size, y = Cor, label = module.id),
                  color = "#d6604d",
                  size= 3.5, #字体大小
                  nudge_x = .15,
                  box.padding = 0.5,
                  nudge_y = 0.15,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20) +
  geom_text_repel(data = data  %>% 
                    dplyr::arrange(Cor) %>% head(5),
                  aes(x = module.size, y = Cor, label = module.id),
                  color = "#4393c3",
                  nudge_x = .15,
                  box.padding = 0.5,
                  size= 3.5, #字体大小
                  nudge_y = 0.15,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20) + 
  theme_bw()


###############################
######Statistics on SACMs######
###############################
load(file = "./Aging/data/UCSC/GTEX_RNAseq_data/Screening for Senescence Characteristics/All_sig_module.Rdata")
tumor_tissue<- as.data.frame(read.csv("./Aging/data/UCSC/GTEX_RNAseq_data/tissue.csv",header = T,sep = ","))
paixu <- c("Uterus","Prostate","Liver","Breast","Lung","Colon","Kidney","Skin","Ovary",
           "Esophagus","Brain","Blood","Adrenal_Gland","Pancreas","Testis","Thyroid","Stomach")
tumor_tissue <- tumor_tissue[match(paixu,tumor_tissue$tissue,nomatch = 0),]
aging_module_all_gene <- list()
tissue_module_sig_total <- data.frame(tissue=tumor_tissue$tissue,up_module =0,down_module = 0, all_module = 0,up_gene=0,down_gene=0)
for (x in 1:17) {
  tissue <- tumor_tissue[x,1]
  #up
  Sigupmodule <- All_sig_module[[paste0(tissue,"_up")]]
  Sigupgene <- list(unique(unlist(Sigupmodule)))
  names(Sigupgene) <- paste0(tissue,"_up_gene")
  tissue_module_sig_total[x,5] <- length(unique(unlist(Sigupmodule)))
  tissue_module_sig_total[x,2] <- length(Sigupmodule)
  #down
  Sigdownmodule <- All_sig_module[[paste0(tissue,"_down")]]
  Sigdowngene <- list(unique(unlist(Sigdownmodule)))
  names(Sigdowngene) <- paste0(tissue,"_down_gene")
  tissue_module_sig_total[x,6] <- length(unique(unlist(Sigdownmodule)))
  tissue_module_sig_total[x,3] <- length(Sigdownmodule)
  
  tissue_module_sig_total[x,4] <- length(Sigdownmodule)+length(Sigupmodule)
  a <- c(Sigupgene,Sigdowngene)
  
  aging_module_all_gene <- c(aging_module_all_gene,a)
}
save(aging_module_all_gene,file = "./Aging/data/UCSC/GTEX_RNAseq_data/Screening for Senescence Characteristics/aging_module_all_gene.Rdata")
write.table(tissue_module_sig_total,file ="./Aging/data/UCSC/GTEX_RNAseq_data/Screening for Senescence Characteristics/tissue_module_sig_total.txt",row.names = F,col.names = T,sep="\t")




#----Figure 1D--------
SACMs_moduleNum <- tissue_module_sig_total[,c(1,2,3,4)]
SACMs_moduleNum <- SACMs_moduleNum[order(SACMs_moduleNum$all_module,decreasing = F),]
SACMs_moduleNum$tissue <- factor(SACMs_moduleNum$tissue,levels = SACMs_moduleNum$tissue)
SACMs_moduleNum <- SACMs_moduleNum[,-4]
SACMs_moduleNum <- reshape2::melt(SACMs_moduleNum,id.vars = 'tissue')
colnames(SACMs_moduleNum)<-c('tissue_type','module_type','Pair_number')
ggdotchart(SACMs_moduleNum, x="tissue_type", y="Pair_number", color = "module_type",group = "module_type", 
           palette = c("#B92022", "#3675BB"),
           legend = "right",
           sorting = "none",
           add = "segments",
           add.params = list(color = "lightgray", size = 3),
           dot.size = 7,
           label = round(SACMs_moduleNum$Pair_number), 
           font.label = list(color="white",size=9, vjust=0.35),
           rotate = T,
           ggtheme = theme_pubr(),
)

#----Figure 1E--------
SACMs_geneNum <- tissue_module_sig_total[,c(1,7,8)]
SACMs_geneNum$tissue <- factor(SACMs_geneNum$tissue,levels = SACMs_geneNum$tissue)
SACMs_geneNum <- reshape2::melt(SACMs_geneNum,id.vars = 'tissue')
ggplot(SACMs_geneNum, aes(
  x = factor(tissue,levels = unique(tissue)),             
  y = ifelse(variable == "up_gene", value, -value),
  fill = variable)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = c("#EE3432","#1B64A4"))+
  coord_flip()+       
  geom_text(                                                
    aes(label=value,                                         
        hjust = ifelse(variable == "up_gene", -0.4, 1.1)  
    ),
    size=3)+
  ylab("")+xlab("Tissues")+
  scale_y_continuous(                                         
    labels = abs,                                             
    expand = expansion(mult = c(0.1, 0.1)))+
  theme_bw()+
  theme(panel.grid=element_blank ())



###########################################################
######Sharing of genes in SACMs between organisations######
###########################################################
load(file = "./Aging/data/UCSC/GTEX_RNAseq_data/Screening for Senescence Characteristics/aging_module_all_gene.Rdata")
tissue_up_aging_gene<- aging_module_all_gene[seq(1,34,2)]
tissue_down_aging_gene<- aging_module_all_gene[seq(2,34,2)]
for(y in 1:length(names(tissue_up_aging_gene))){
  names(tissue_up_aging_gene)[y] <- paste0(names(tissue_up_aging_gene[y])," (",length(tissue_up_aging_gene[[y]]),")")
}
for(y in 1:length(names(tissue_down_aging_gene))){
  names(tissue_down_aging_gene)[y] <- paste0(names(tissue_down_aging_gene[y])," (",length(tissue_down_aging_gene[[y]]),")")
}
#down
age_down_total <- as.data.frame(table(as.data.frame(unlist(tissue_down_aging_gene))))
age_down_total <- age_down_total[order(age_down_total$Freq,decreasing = T),]
down_Freq <- as.data.frame(table(age_down_total$Freq))
colnames(down_Freq) <- c("pinshu","cishu")
#up
age_up_total <- as.data.frame(table(as.data.frame(unlist(tissue_up_aging_gene))))
age_up_total <- age_up_total[order(age_up_total$Freq,decreasing = T),]
up_Freq <- as.data.frame(table(age_up_total$Freq))
colnames(up_Freq) <- c("pinshu","cishu")

#----Figure 1F right----
up_Freq <- rbind(up_Freq,c("1",0),c("1",0))
two_group_data <- cbind(down_Freq,as.numeric(up_Freq[,2]))
colnames(two_group_data) <- c("tissue_num","Down","Up")
df1 <- reshape2::melt(two_group_data,id.vars = 'tissue_num')
library(ggplot2)
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


#----Figure 1F left----
age_up_total_top <- age_up_total[which(age_up_total$Freq>=5),]
colnames(age_up_total_top) <- c("gene","freq")
age_up_total_top <- age_up_total_top[order(age_up_total_top$freq,decreasing = F),]

age_up_total_top$freq <- factor(age_up_total_top$freq,levels = c("5","6"))
age_up_total_top$gene <- factor(age_up_total_top$gene,levels = age_up_total_top$gene)
ggplot(age_up_total_top,aes(x=gene,y=freq))+
  geom_bar(aes(fill=freq),stat = "identity",width = 0.8,colour="black",
           size=0.2,alpha=1)+
  scale_fill_manual(values =c(brewer.pal(n=3,name = "Set3")))+
  coord_flip()+
  theme_classic()+theme(
    legend.position = "none"
  )

age_down_total_top <- age_down_total[which(age_down_total$Freq>=6),]
colnames(age_down_total_top) <- c("gene","freq")
age_down_total_top <- age_down_total_top[order(age_down_total_top$freq,decreasing = F),]

age_down_total_top$freq <- factor(age_down_total_top$freq,levels = c("6","7","8"))
age_down_total_top$gene <- factor(age_down_total_top$gene,levels = age_down_total_top$gene)
ggplot(age_down_total_top,aes(x=gene,y=freq))+
  geom_bar(aes(fill=freq),stat = "identity",width = 0.8,colour="black",
           size=0.2,alpha=1)+
  scale_fill_manual(values =c(brewer.pal(n=3,name = "Set3")))+
  coord_flip()+
  theme_classic()+theme(
    legend.position = "none"
  )

#####################################################################################
######cell type over-representation analysis of the senescence-related genes.######
#####################################################################################
#-----Figure 1G
Up_cellEnrich <- as.data.frame(read.table("./Aging/data/UCSC/GTEX_RNAseq_data/enrichr/上调/PanglaoDB_Augmented_2021_table.txt",header = T,sep = "\t"))
Down_cellEnrich <- as.data.frame(read.table("./Aging/data/UCSC/GTEX_RNAseq_data/enrichr/下调/PanglaoDB_Augmented_2021_table (1).txt",header = T,sep = "\t"))
Up_Down <- rbind(Up_cellEnrich[1:15,1:3],Down_cellEnrich[1:15,1:3])
Up_Down$type <- c(rep("up",15),rep("down",15))
Up_Down$LogP <- log10(Up_Down$P.value)
Up_Down <- Up_Down %>% arrange(type,desc(LogP))
Up_Down$LogP <- ifelse(Up_Down$type == "down",Up_Down$LogP,-Up_Down$LogP)
Up_Down$Term <- factor(Up_Down$Term,Up_Down$Term)
x_max = max(abs(Up_Down$LogP))
color_list <- c("#B0CAE6","#E5B3B3")
ggplot(Up_Down,aes(x =LogP, y = Term, fill = type)) + 
  geom_col(width = 0.1) + #添加条形图，收窄柱子为一条线
  geom_point(size=3,aes(color = type)) + #添加散点/气泡
  scale_x_continuous(limits = c(-x_max,x_max)) +
  scale_fill_manual(values = color_list) +
  scale_color_manual(values = color_list)+
  geom_text(data = Up_Down[which(Up_Down$type == "up"),],
            aes(x = -0.5, y = Term, label = Term,color = type),
            size = 4,
            hjust = 1) + 
  geom_text(data = Up_Down[which(Up_Down$type == "down"),],
            aes(x = 0.5, y = Term, label = Term,color = type),
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

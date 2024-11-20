library(DESeq2)
library(data.table)
library(stringr)
library(ggplot2)
library(DEGreport)
library(ggrepel)
library(plyr)
###################################################
######Combining cancer and normal tissue data######
###################################################
tumor_tissue <- as.data.frame(read.csv("./Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/tumor_tissue.csv",header = T,sep = ","))
#去除版本号
ID.ver <- function(id){
  id2 <- unlist(strsplit(id,".",fixed = T))
  id3 <- id2[1]
  return(id3)
}
genetype <- as.data.frame(fread("./Aging/data/UCSC/probeMap_gencode.v23.annotation.gene.probemap",header = T,sep = "\t"))
colnames(genetype)[1:2] <- c("Ensemble","Symbol")
genetype$Ensemble <- apply(as.data.frame(genetype[,1]), 1, ID.ver)
for (x in 1:24) {
  setwd(paste0("./Aging/data/UCSC/TCGA_RNAseq_data/",tumor_tissue[x,1]))
  #Normaltissue_count
  normal_count <- as.data.frame(fread(paste0(tumor_tissue[x,2],"_RNAseq_count.txt"),header = T,sep = "\t"))
  rownames(normal_count) <- normal_count[,1]
  normal_count <- normal_count[,-1]
  #Normaltissue_clinical
  normal_clinical <- as.data.frame(fread(paste0(tumor_tissue[x,2],"_sample_info.txt"),header = T,sep = "\t"))
  normal_clinical <- normal_clinical[,c("SAMPID","SEX","Age")]
  identical(colnames(normal_count),normal_clinical$SAMPID)#TRUE
  #Tumor_count
  cancer_count <- as.data.frame(fread(paste0(tumor_tissue[x,1],"_RNAseq_count.txt"),header = T,sep = "\t"))
  rownames(cancer_count) <- cancer_count[,1]
  cancer_count <- cancer_count[,-1]
  #Tumor_clinical
  cancer_clinical <- as.data.frame(fread(paste0(tumor_tissue[x,1],"_sample_info.txt"),header = T,sep = "\t"))
  cancer_clinical <- cancer_clinical[,c("sample","age_at_initial_pathologic_diagnosis","gender")]
  identical(colnames(cancer_count),cancer_clinical$sample)#TRUE
  #Converting age from cancer sample information to age intervals
  cancer_clinical$age <- ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=19,"10-19",
                                ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=29,"20-29",
                                       ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=39,"30-39",
                                              ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=49,"40-49",
                                                     ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=59,"50-59",
                                                            ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=69,"60-69",
                                                                   ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=79,"70-79",
                                                                          ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=89,"80-89",
                                                                                 ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=99,"90-99")))))))))

  cancer_clinical <- cancer_clinical[complete.cases(cancer_clinical[,"age_at_initial_pathologic_diagnosis"]),]
  cancer_count <- cancer_count[,cancer_clinical$sample]
  #Combining data
  data <- cbind(cancer_count,normal_count)
  data <- 2^data-1
  data <- round(data)
  gene_id <- apply(as.data.frame(rownames(data)), 1, ID.ver)
  rownames(data) <- gene_id
  count_data=data[which(apply(data,1,function(x){return(sum(x>1))})>ncol(data)*0.20),]
  count_data$Ensemble <- rownames(count_data)
  count_data_symbol <- merge(genetype,count_data,by="Ensemble")
  count_data_symbol=count_data_symbol[!duplicated(count_data_symbol$Symbol),]
  rownames(count_data_symbol)=count_data_symbol$Symbol
  count_data_symbol=count_data_symbol[,c(7:ncol(count_data_symbol))]
  count_nc <- count_data_symbol
  length(rownames(count_nc))
  normal_clinical$SEX <- ifelse(normal_clinical$SEX == "1","MALE","FEMALE")
  a <- cancer_clinical[,c("sample","gender","age")]
  b <- normal_clinical[,c("SAMPID","SEX","Age")]
  colnames(b) <- c("sample","gender","age")
  age_clinical <- rbind(a,b)
  identical(colnames(count_nc),age_clinical$sample)
  condition <- factor(ifelse(substr(colnames(count_nc),14,14) == "0" & substr(colnames(count_nc),1,4)=="TCGA","Tumor","Normal"))
  colData <- data.frame(row.names = age_clinical$sample,condition,age_clinical$age,age_clinical$gender)
  colnames(colData)[2] <- "Age"
  colnames(colData)[3] <- "gender"
  colData$Age <- str_replace(colData$Age,"-","_")
  colData$Age <- as.factor(colData$Age)
  colData$gender <- as.factor(colData$gender)
  identical(colnames(count_nc),rownames(colData))
  table(colData$condition)
  write.table(count_nc,paste0(tumor_tissue[x,1],"_normal_count_data.txt"),sep="\t",row.names=TRUE,quote=F)
  write.table(colData,paste0(tumor_tissue[x,1],"_normal_sampleInfomation.txt"),sep="\t",row.names=TRUE,quote=F)
}

######Fig.3B----######
fourtime_gene <- c("HMGA1","KCNK1","MIF","SH3BP1","TALDO1")
#HMGA1   c("Colon","Skin","Esophagus","Pancreas") c(15,16,20,3,14)
#KCNK1   c("Breast","Colon","Skin","Esophagus")   c(2,15,16,20,3)
#MIF     c("Prostate","Breast","Colon","Esophagus") c(9,2,15,16,3)
#SH3BP1  c("Colon","Skin","Esophagus","Pancreas")  c(15,16,20,3,14)
#TALDO1  c("Prostate","Colon","Skin","Esophagus")  c(9,15,16,20,3)
fourtime_gene_tissueList <- list(c(15,16,20,3,14),c(2,15,16,20,3),c(9,2,15,16,3),c(15,16,20,3,14),c(9,15,16,20,3))
names(fourtime_gene_tissueList) <- fourtime_gene
plot_list <- list()
plot_listAll <- list()
for (i in 1:5) {
  one_gene <- fourtime_gene[i]
  onetissue <- tumor_tissue[fourtime_gene_tissueList[[one_gene]],]
  for(x in 1:length(onetissue$tumor)){
    setwd(paste0("./Aging/data/UCSC/TCGA_RNAseq_data/",onetissue[x,1]))
    vsd_count_nc <- read.table(paste0(onetissue[x,1],"_normal_count_data.txt"),header=T,row.names=1,check.names=F)
    colData <- read.table(paste0(onetissue[x,1],"_normal_sampleInfomation.txt"),header=T,row.names=1)
    one_gene_count <- as.data.frame(t(vsd_count_nc[one_gene,]))
    identical(rownames(one_gene_count),rownames(colData))#TRUE
    one_gene_count <- cbind(one_gene_count,colData$Age,colData$condition)
    colnames(one_gene_count) <- c("Expression","Age","condition")
    one_gene_count$Age <- as.numeric(one_gene_count$Age)
    p1 <- ggplot(one_gene_count,aes(x=Age,y=Expression,colour=condition))+
      geom_point(size = 3)+
      labs(title = paste0(one_gene," expression with age in ",onetissue[x,1]),x = "Age",y = "Normalize_Expression",fill = "Condition")+
      scale_color_manual(values = c(Tumor="#ce181e",Normal="#007cc0"))+
      geom_smooth(method = 'lm')+
      scale_x_continuous(breaks = seq(10, 80, 10), limits = c(10, 80))+
      theme_bw()+
      theme(axis.text=element_text(colour='black',size=9))
    plot_list[[x]] <- p1
  }
  plot_listAll <- c(plot_listAll,plot_list)
}
pdf("./Aging/data/UCSC/TCGA_RNAseq_data/five_gene.pdf",
    width = 22,height = 15)
plot(plot_grid(plotlist = plot_listAll, ncol = 5))
dev.off()

########################################################################################################
####Overlap of cancer-suppressing genes, oncogenes, aging-positive genes, and aging-negative genes.#####
########################################################################################################
load(file = "D:\\procedure\\Aging\\data\\UCSC\\GTEX_RNAseq_data\\最新结果文件\\衰老基因集合筛选\\aging_module_all_gene.Rdata")
TSG <- read.table("D:\\procedure\\Aging\\data\\UCSC\\TCGA_RNAseq_data\\Human_TSGs.txt",header = T,sep = "\t")
Cancer_gene <- read.table("D:\\procedure\\Aging\\data\\UCSC\\TCGA_RNAseq_data\\cancerGeneList.csv",header = T,sep = ",")
Onocgene <- Cancer_gene[which(Cancer_gene$Is.Oncogene == "Yes"),]
Aging_Up_gene <- unique(unlist(aging_module_all_gene[grep("up",names(aging_module_all_gene))]))
Aging_Down_gene <- unique(unlist(aging_module_all_gene[grep("down",names(aging_module_all_gene))]))
library(ggvenn)
library(GeneOverlap)
#Aging Up with Onocgene
ggvenn(list(Aging_Up_Allgene=Aging_Up_gene,Onocgene=Onocgene$Hugo.Symbol),fill_color=c("#DE7962", "#EAE5D0"),fill_alpha = .7,stroke_linetype = "longdash",set_name_size = 8,
       text_size=5)
gen.ov <- testGeneOverlap(newGeneOverlap(Aging_Up_gene,Onocgene$Hugo.Symbol,genome.size = 40027))

#Aging Down with Onocgene
ggvenn(list(Aging_Down_Allgene=Aging_Down_gene,Onocgene=Onocgene$Hugo.Symbol),fill_color=c("#6F84AF", "#EAE5D0"),fill_alpha = .7,stroke_linetype = "longdash",set_name_size = 8,
       text_size=5)
gen.ov <- testGeneOverlap(newGeneOverlap(Aging_Down_gene,Onocgene$Hugo.Symbol,genome.size = 40027))

#Aging Up with Tumor Suppressor genes
ggvenn(list(Aging_Up_Allgene=Aging_Up_gene,Onocgene=TSG$GeneSymbol),fill_color=c("#DE7962", "#EDF6FC"),fill_alpha = .7,stroke_linetype = "longdash",set_name_size = 8,
       text_size=5)
gen.ov <- testGeneOverlap(newGeneOverlap(Aging_Up_gene,TSG$GeneSymbol,genome.size = 40027))

#Aging Down with Tumor Suppressor genes
ggvenn(list(Aging_Down_Allgene=Aging_Down_gene,Onocgene=TSG$GeneSymbol),fill_color=c("#6F84AF", "#EDF6FC"),fill_alpha = .7,stroke_linetype = "longdash",set_name_size = 8,
       text_size=5)
gen.ov <- testGeneOverlap(newGeneOverlap(Aging_Down_gene,TSG$GeneSymbol,genome.size = 40027))



###############################################################
######Cycling for DEseq2 differential expression analysis######
###############################################################
tumor_tissue <- as.data.frame(read.csv("./Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/tumor_tissue.csv",header = T,sep = ","))
DEGs_tolal <- data.frame(tumor_up_DEGs=c(1),
                         tumor_down_DEGs=c(1),
                         up_age_gene=c(1),
                         up_age_tumor_up_DEGs=c(1),
                         up_age_tumor_down_DEGs=c(1),
                         down_age_gene=c(1),
                         down_age_tumor_up_DEGs=c(1),
                         down_age_tumor_down_DEGs=c(1),
                         tumor_down_DEGs_potein =c(1),
                         tumor_up_DEGs_potein = c(1),
                         up_age_gene_potein = c(1),
                         down_age_gene_potein = c(1))
ID.ver <- function(id){
  id2 <- unlist(strsplit(id,".",fixed = T))
  id3 <- id2[1]
  return(id3)
}
gene_annotation <- read.table("./Aging/data/UCSC/gene_id_type.txt",sep = "\t",header = T)
procoding_gene <- gene_annotation[which(gene_annotation$type == "protein_coding"),]
load(file = "./Aging/data/UCSC/GTEX_RNAseq_data/age_list_uniq_new.Rdata")
DEGs_list <- list()
for(x in 1:24){
  DEGs_tolal[x,] <- c(0,0,0,0,0,0,0,0,0,0,0,0)
  setwd(paste0("./Aging/data/UCSC/TCGA_RNAseq_data/",tumor_tissue[x,1]))
  #Reading the combined GTEX and TCGA data
  count_nc <- read.table(paste0(tumor_tissue[x,1],"_normal_count_data.txt"),header=T,row.names=1,check.names=F)
  colData <- read.table(paste0(tumor_tissue[x,1],"_normal_sampleInfomation.txt"),header=T,row.names=1)
  
  #DEseq2----Age as a covariate
  if(length(table(colData$gender)) == 1){
    colData$condition <- factor(colData$condition,levels = c("Tumor","Normal"))
    colData$Age <- factor(colData$Age,levels = unique(colData$Age)[order(unique(colData$Age),decreasing = F)])
    dds <- DESeqDataSetFromMatrix(countData = count_nc,colData = colData,design = ~ condition + Age)
    dds1 <- DESeq(dds,test="LRT",reduced = ~Age)
    res <- results(dds1,contrast=c("condition","Tumor","Normal"),
                   independentFiltering=TRUE,alpha=0.05,pAdjustMethod="BH")
  }else{
    #DEseq2----Age and Sex as a covariate
    colData$condition <- factor(colData$condition,levels = c("Tumor","Normal"))
    colData$Age <- factor(colData$Age,levels = unique(colData$Age)[order(unique(colData$Age),decreasing = F)])
    colData$gender <- factor(colData$gender,levels = c("MALE","FEMALE"))
    dds <- DESeqDataSetFromMatrix(countData = count_nc,colData = colData,design = ~ condition + Age + gender )
    dds1 <- DESeq(dds,test="LRT",reduced = ~Age + gender)
    res <- results(dds1,contrast=c("condition","Tumor","Normal"),
                   independentFiltering=TRUE,alpha=0.05,pAdjustMethod="BH")
  }
  res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  write.table(res1,"DEGs_all_info.txt",col.names = T,row.names = T,quote = F,sep = "\t")
  # padj<0.01，|log2FC|>1
  res1_up<- res1[which(res1$log2FoldChange >= 1.5 & res1$padj < 0.01),] 
  res1_down<- res1[which(res1$log2FoldChange <= -1.5 & res1$padj < 0.01),]
  res1_total <- rbind(res1_up,res1_down)
  DEGs_tolal[x,1] <- length(rownames(res1_up))
  DEGs_tolal[x,2] <- length(rownames(res1_down))
  res1_up_potein <- procoding_gene[match(rownames(res1_up),procoding_gene$gene_name,nomatch = 0),]
  res1_down_potein <- procoding_gene[match(rownames(res1_down),procoding_gene$gene_name,nomatch = 0),]
  DEGs_tolal[x,10]<- length(rownames(res1_up_potein))/length(rownames(res1_up))
  DEGs_tolal[x,9] <- length(rownames(res1_down_potein))/length(rownames(res1_down))
  #DEGs_list
  up_list <- list(rownames(res1_up),rownames(res1_down))
  names(up_list) <- c(paste0(tumor_tissue[x,1],"_DEGs_up"),paste0(tumor_tissue[x,1],"_DEGs_down"))
  DEGs_list <- c(DEGs_list,up_list)
  #write DEGs result
  write.table(res1_up,"gene_up.txt",col.names = T,row.names = T,quote = F,sep = "\t")
  write.table(rownames(res1_up),"up_gene.txt",sep="\t",quote = F,row.names = F,col.names = F)
  write.table(res1_down,"gene_down.txt",col.names = T,row.names = T,quote = F,sep = "\t")
  write.table(rownames(res1_down),"down_gene.txt",sep="\t",quote = F,row.names = F,col.names = F)
  age_up_gene <- age_list[[paste0(tumor_tissue[x,2],"_age_up")]]
  age_down_gene <- age_list[[paste0(tumor_tissue[x,2],"_age_down")]]
  DEGs_tolal[x,3] <- length(age_up_gene)
  DEGs_tolal[x,6] <- length(age_down_gene)
  age_up_potein <- procoding_gene[match(age_up_gene,procoding_gene$gene_name,nomatch = 0),]
  age_down_potein <- procoding_gene[match(age_down_gene,procoding_gene$gene_name,nomatch = 0),]
  DEGs_tolal[x,11]<- length(rownames(age_up_potein))/length(age_up_gene)
  DEGs_tolal[x,12] <- length(rownames(age_down_potein))/length(age_down_gene)
  DEGs_tolal[x,4] <- length(intersect(age_up_gene,rownames(res1_up)))
  DEGs_tolal[x,5] <- length(intersect(age_up_gene,rownames(res1_down)))
  DEGs_tolal[x,7] <- length(intersect(age_down_gene,rownames(res1_up)))
  DEGs_tolal[x,8] <- length(intersect(age_down_gene,rownames(res1_down)))
}
rownames(DEGs_tolal) <- tumor_tissue$tumor[1:24]
save(DEGs_list,file = "./Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/DEGs_list.Rdata")
save(DEGs_tolal,file = "./Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/DEGs_tolal_new.Rdata")



##----Fig.3D-------
data <- cbind(rownames(DEGs_tolal),DEGs_tolal[,1:2])
colnames(data)[1] <- "cancer"
data <- data[match(rev(data$cancer),data$cancer),]
data$cancer <- factor(data$cancer,levels = data$cancer)
data <- reshape2::melt(data,id.vars='cancer')
ggplot(data, aes(
  x = cancer,             
  y = ifelse(variable == "tumor_down_DEGs", value, -value),  
  fill = variable)) +
  geom_bar(stat = 'identity')+                         
  scale_fill_manual(values = c("#D6251F","#3778AD"))+
  coord_flip()+                                              
  geom_text(                                                 
    aes(label=value,                                          
        hjust = ifelse(variable == "tumor_down_DEGs", -0.4, 1.1)           
    ),
    size=3)+
  ylab("")+xlab("Cancer")+
  scale_y_continuous(                                   
    labels = abs,                                        
    expand = expansion(mult = c(0.1, 0.1)))+                 
  theme_bw()+
  theme(
    axis.title = element_text(size = 15,face = "plain",color = "black"),
    axis.text = element_text(size = 12,face = "plain",color = "black"),
    legend.title = element_text(size = 14,face = "plain",color = "black"),
    legend.position = "right",
    panel.grid=element_blank ()
  )


####################################################
######Overlap between Age DEGs and cancer DEGs######
####################################################
######Fig.3E-
library(GeneOverlap)
library(ggplot2)
load(file = "./Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/DEGs_list.Rdata")
load(file = "./Aging/data/UCSC/GTEX_RNAseq_data/all_tissue_gene_list.Rdata")
load(file = "./Aging/data/UCSC/GTEX_RNAseq_data/aging_module_all_gene.Rdata")
tumor_tissue <- as.data.frame(read.csv("./Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/tumor_tissue.csv",header = T,sep = ","))
DEGs_tolal <- data.frame(tumor_up_DEGs=c(1),
                         tumor_down_DEGs=c(1),
                         up_age_gene=c(1),
                         down_age_gene=c(1),
                         up_age_tumor_up_DEGs=c(1),
                         up_age_tumor_up_DEGs_P=c(1),
                         up_age_tumor_down_DEGs=c(1),
                         up_age_tumor_down_DEGs_P=c(1),
                         down_age_tumor_up_DEGs=c(1),
                         down_age_tumor_up_DEGs_P=c(1),
                         down_age_tumor_down_DEGs=c(1),
                         down_age_tumor_down_DEGs_P=c(1))

UT_UA_list <- list()
UT_DA_list <- list()
DT_UA_list <- list()
DT_DA_list <- list()

for(x in 1:24){
  DEGs_tolal[x,] <- c(0,0,0,0,0,0,0,0,0,0,0,0)
  rownames(DEGs_tolal)[x] <- tumor_tissue[x,1]
  cancer_DEGs_up <- DEGs_list[[paste0(tumor_tissue[x,1],"_DEGs_up")]]
  DEGs_tolal[x,1] <- length(cancer_DEGs_up)
  
  cancer_DEGs_down <- DEGs_list[[paste0(tumor_tissue[x,1],"_DEGs_down")]]
  DEGs_tolal[x,2] <- length(cancer_DEGs_down)
  
  age_DEGs_up <- aging_module_all_gene[[paste0(tumor_tissue[x,2],"_up_gene")]]
  DEGs_tolal[x,3] <- length(age_DEGs_up)
  
  age_DEGs_down <- aging_module_all_gene[[paste0(tumor_tissue[x,2],"_down_gene")]]
  DEGs_tolal[x,4] <- length(age_DEGs_down)
  
  
  Up_up_go.obj <- testGeneOverlap(newGeneOverlap(cancer_DEGs_up,age_DEGs_up,genome.size = tumor_tissue[x,4]))
  DEGs_tolal[x,5] <- length(getIntersection(Up_up_go.obj))
  DEGs_tolal[x,6] <- getPval(Up_up_go.obj)
  
  UT_UA_list <- c(UT_UA_list,list(getIntersection(Up_up_go.obj)))
  
  Down_up_go.obj <- testGeneOverlap(newGeneOverlap(cancer_DEGs_down,age_DEGs_up,genome.size = tumor_tissue[x,4]))
  DEGs_tolal[x,7] <- length(getIntersection(Down_up_go.obj))
  DEGs_tolal[x,8] <- getPval(Down_up_go.obj)
  
  DT_UA_list <- c(DT_UA_list,list(getIntersection(Down_up_go.obj)))
  
  Up_down_go.obj <- testGeneOverlap(newGeneOverlap(cancer_DEGs_up,age_DEGs_down,genome.size = tumor_tissue[x,4]))
  DEGs_tolal[x,9] <- length(getIntersection(Up_down_go.obj))
  DEGs_tolal[x,10] <- getPval(Up_down_go.obj)
  UT_DA_list <- c(UT_DA_list,list(getIntersection(Up_down_go.obj)))
  
  
  
  Down_down_go.obj <- testGeneOverlap(newGeneOverlap(cancer_DEGs_down,age_DEGs_down,genome.size = tumor_tissue[x,4]))
  DEGs_tolal[x,11] <- length(getIntersection(Down_down_go.obj))
  DEGs_tolal[x,12] <- getPval(Down_down_go.obj)
  DT_DA_list <- c(DT_DA_list,list(getIntersection(Down_down_go.obj)))
}

names(DT_DA_list) <- tumor_tissue$tumor
names(DT_UA_list) <- tumor_tissue$tumor
names(UT_DA_list) <- tumor_tissue$tumor
names(UT_UA_list) <- tumor_tissue$tumor

save(DT_DA_list,file = "./Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/DT_DA_list.Rdata")
save(DT_UA_list,file = "./Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/DT_UA_list.Rdata")
save(UT_DA_list,file = "./Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/UT_DA_list.Rdata")
save(UT_UA_list,file = "./Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/UT_UA_list.Rdata")




paixu <- as.data.frame(read.csv("./Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/paixu.csv",header = T,sep = ","))
DEGs_tolal <- DEGs_tolal[match(paixu$tumor,rownames(DEGs_tolal)),]

#UTUA
library(pheatmap)
library(ggthemes)
Up_up_data <- DEGs_tolal[,c(1,3,5)]
colnames(Up_up_data) <- c("UpInTumor","UpInAge","Overlap")
bk <- c(seq(0,4500,by=100))
pheatmap(as.matrix(Up_up_data),
         scale = "none",
         color = c(colorRampPalette(colors = c("white","#3CCF4E"))(length(bk))),
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 12,
         cellwidth = 28,#小格子的宽度
         cellheight = 20, #小格子的长度
         legend_breaks=seq(0),
         display_numbers = TRUE,
         number_format = "%.0f",
         angle_col = 45,
         legend = FALSE,
         fontsize_col = 10,
         
)
Up_up_Pvalue <- data.frame(cancer=rownames(DEGs_tolal),p_value=DEGs_tolal$up_age_tumor_up_DEGs_P)
Up_up_Pvalue$p_value <- -1*log10(Up_up_Pvalue$p_value)
Up_up_Pvalue$tp <- ifelse(Up_up_Pvalue$p_value < -log10(0.05),"",
                          ifelse(Up_up_Pvalue$p_value >= -log10(0.05) & Up_up_Pvalue$p_value<=-log10(0.01),"*",
                                 ifelse(Up_up_Pvalue$p_value >= -log10(0.01) & Up_up_Pvalue$p_value<=-log10(0.001),"**","***")))
Up_up_Pvalue$p_value <- ifelse(Up_up_Pvalue$p_value > 20,20,Up_up_Pvalue$p_value)
Up_up_Pvalue$cancer <- factor(Up_up_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Up_up_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c("STAD"="#e8d738","THCA"="#c62574","PAAD"="#34a5b1","ACC"="#d4a2c5","THYM"="#ac468c","DLBC"="#d1636a","TGCT"="#b31e23",
                               "GBM"="#f1b3b9","LGG"="#f1b3b9","ESCA"="#11398d","SKCM"="#5c308e","READ"="#977eb8","COAD"="#d17613","OV"="#6f779f",
                               "KIRP"="#3c97c6","KIRC"="#1881b1","KICH"="#b31e23","LUSC"="#e18b1c","LUAD"="#9fcba6","LIHC"="#a6cd3f","BRCA"="#286598",
                               "PRAD"="#ece93d","UCS"="#4558a6","UCEC"="#d1b091"))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())

#DTDA
Down_down_data <- DEGs_tolal[,c(2,4,11)]
colnames(Down_down_data) <- c("DownInTumor","DownInAge","Overlap")
bk <- c(seq(0,4500,by=100))
p2 <- pheatmap(Down_down_data,
               scale = "none",
               color = c(colorRampPalette(colors = c("white","#E55807"))(length(bk))),
               cluster_rows = F,
               cluster_cols = F,
               fontsize = 12,
               cellwidth = 28,#小格子的宽度
               cellheight = 20, #小格子的长度
               legend_breaks=seq(0),
               display_numbers = TRUE,
               number_format = "%.0f",
               angle_col = 45,
               legend = FALSE,
               fontsize_col = 10,
               show_rownames = F,
)

Down_down_Pvalue <- data.frame(cancer=rownames(DEGs_tolal),p_value=DEGs_tolal$down_age_tumor_down_DEGs_P)
Down_down_Pvalue$p_value <- -1*log10(Down_down_Pvalue$p_value)
Down_down_Pvalue$tp <- ifelse(Down_down_Pvalue$p_value < -log10(0.05),"",
                              ifelse(Down_down_Pvalue$p_value >= -log10(0.05) & Down_down_Pvalue$p_value<=-log10(0.01),"*",
                                     ifelse(Down_down_Pvalue$p_value >= -log10(0.01) & Down_down_Pvalue$p_value<=-log10(0.001),"**","***")))
Down_down_Pvalue$p_value <- ifelse(Down_down_Pvalue$p_value > 20,20,Down_down_Pvalue$p_value)
Down_down_Pvalue$cancer <- factor(Down_down_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Down_down_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c("#e8d738","#c62574","#34a5b1","#d4a2c5","#ac468c","#d1636a","#b31e23",
                               "#f1b3b9","#f1b3b9","#11398d","#5c308e","#977eb8","#d17613","#6f779f",
                               "#3c97c6","#1881b1","#b31e23","#e18b1c","#9fcba6","#a6cd3f","#286598",
                               "#ece93d","#4558a6","#d1b091"))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())

#DTUA
Down_up_data <- DEGs_tolal[,c(2,3,7)]
colnames(Down_up_data) <- c("DownInTumor","UpInAge","Overlap")
bk <- c(seq(0,3500,by=100))
p3 <- pheatmap(Down_up_data,
               scale = "none",
               color = c(colorRampPalette(colors = c("white","#7E1717"))(length(bk))),
               cluster_rows = F,
               cluster_cols = F,
               fontsize = 12,
               cellwidth = 28,#小格子的宽度
               cellheight = 20, #小格子的长度
               legend_breaks=seq(0),
               display_numbers = TRUE,
               number_format = "%.0f",
               angle_col = 45,
               legend = FALSE,
               fontsize_col = 10,
               show_rownames = F,
)

Down_up_Pvalue <- data.frame(cancer=rownames(DEGs_tolal),p_value=DEGs_tolal$up_age_tumor_down_DEGs_P)
Down_up_Pvalue$p_value <- -1*log10(Down_up_Pvalue$p_value)
Down_up_Pvalue$tp <- ifelse(Down_up_Pvalue$p_value < -log10(0.05),"",
                            ifelse(Down_up_Pvalue$p_value >= -log10(0.05) & Down_up_Pvalue$p_value<=-log10(0.01),"*",
                                   ifelse(Down_up_Pvalue$p_value >= -log10(0.01) & Down_up_Pvalue$p_value<=-log10(0.001),"**","***")))
Down_up_Pvalue$p_value <- ifelse(Down_up_Pvalue$p_value > 20,20,Down_up_Pvalue$p_value)
Down_up_Pvalue$cancer <- factor(Down_up_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Down_up_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c("#e8d738","#c62574","#34a5b1","#d4a2c5","#ac468c","#d1636a","#b31e23",
                               "#f1b3b9","#f1b3b9","#11398d","#5c308e","#977eb8","#d17613","#6f779f",
                               "#3c97c6","#1881b1","#b31e23","#e18b1c","#9fcba6","#a6cd3f","#286598",
                               "#ece93d","#4558a6","#d1b091"))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())

#UTDA
Up_down_data <- DEGs_tolal[,c(1,4,9)]
colnames(Up_down_data) <- c("UpInTumor","DownInAge","Overlap")
bk <- c(seq(0,3500,by=100))
p4 <- pheatmap(Up_down_data,
               scale = "none",
               color = c(colorRampPalette(colors = c("white","#068DA9"))(length(bk))),
               cluster_rows = F,
               cluster_cols = F,
               fontsize = 12,
               cellwidth = 28,#小格子的宽度
               cellheight = 20, #小格子的长度
               legend_breaks=seq(0),
               display_numbers = TRUE,
               number_format = "%.0f",
               angle_col = 45,
               legend = FALSE,
               fontsize_col = 10,
               show_rownames = F,
)

Up_down_Pvalue <- data.frame(cancer=rownames(DEGs_tolal),p_value=DEGs_tolal$down_age_tumor_up_DEGs_P)
Up_down_Pvalue$p_value <- -1*log10(Up_down_Pvalue$p_value)
Up_down_Pvalue$tp <- ifelse(Up_down_Pvalue$p_value < -log10(0.05),"",
                            ifelse(Up_down_Pvalue$p_value >= -log10(0.05) & Up_down_Pvalue$p_value<=-log10(0.01),"*",
                                   ifelse(Up_down_Pvalue$p_value >= -log10(0.01) & Up_down_Pvalue$p_value<=-log10(0.001),"**","***")))
Up_down_Pvalue$p_value <- ifelse(Up_down_Pvalue$p_value > 20,20,Up_down_Pvalue$p_value)
Up_down_Pvalue$cancer <- factor(Up_down_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Up_down_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c("#e8d738","#c62574","#34a5b1","#d4a2c5","#ac468c","#d1636a","#b31e23",
                               "#f1b3b9","#f1b3b9","#11398d","#5c308e","#977eb8","#d17613","#6f779f",
                               "#3c97c6","#1881b1","#b31e23","#e18b1c","#9fcba6","#a6cd3f","#286598",
                               "#ece93d","#4558a6","#d1b091"))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())






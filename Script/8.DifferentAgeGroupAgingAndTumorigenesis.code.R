library(DESeq2)
library(data.table)
library(stringr)
library(ggplot2)
library(DEGreport)
library(ggrepel)
library(plyr)
######################################################################################################
#####Cancer and normal paired data were categorised into three groups: young, middle-aged and old.####
######################################################################################################
tumor_tissue <- as.data.frame(read.csv("./Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/tumor_tissue.csv",header = T,sep = ","))
group_sample_tolal <- data.frame(young_sample=c(1),
                                 young_cancer=c(1),
                                 young_normal=c(1),
                                 middle_sample=c(1),
                                 middle_cancer=c(1),
                                 middle_normal=c(1),
                                 old_sample=c(1),
                                 old_cancer=c(1),
                                 old_normal=c(1))
for (x in 1:24) {
  group_sample_tolal[x,] <- c(0,0,0,0,0,0,0,0,0)
  setwd(paste0(".Aging/data/UCSC/TCGA_RNAseq_data/",tumor_tissue[x,1]))
  count_nc <- read.table(paste0(tumor_tissue[x,1],"_normal_count_data.txt"),header=T,row.names=1,check.names=F)
  colData <- read.table(paste0(tumor_tissue[x,1],"_normal_sampleInfomation.txt"),header=T,row.names=1)
  #young
  colData$topage <- substr(as.character(colData$Age),4,5)
  young_colData <- colData[which(colData$topage<=40),]
  young_count_nc <- count_nc[,rownames(young_colData)]
  identical(colnames(young_count_nc),rownames(young_colData))
  #TRUE
  young_colData <- young_colData[,-which(colnames(young_colData) == "topage")]
  young <- as.data.frame(table(young_colData$condition))
  # Normal  Tumor 
  # 10    248 
  group_sample_tolal[x,1] <- length(rownames(young_colData))
  group_sample_tolal[x,2] <- young[2,2]
  group_sample_tolal[x,3] <- young[1,2]
  write.table(young_count_nc,paste0(tumor_tissue[x,1],"_young_normal_count_data.txt"),sep="\t",row.names=TRUE,quote=F)
  write.table(young_colData,paste0(tumor_tissue[x,1],"_young_normal_sampleInfomation.txt"),sep="\t",row.names=TRUE,quote=F)
  #midlle
  middle_colData <- colData[which(colData$topage<=60&colData$topage>=40),]
  middle_count_nc <- count_nc[,rownames(middle_colData)]
  identical(colnames(middle_count_nc),rownames(middle_colData))
  #TRUE
  middle_colData <- middle_colData[,-which(colnames(middle_colData) == "topage")]
  middle <- as.data.frame(table(middle_colData$condition))
  # Normal  Tumor 
  # 94    204 
  group_sample_tolal[x,4] <- length(rownames(middle_colData))
  group_sample_tolal[x,5] <- middle[2,2]
  group_sample_tolal[x,6] <- middle[1,2]
  write.table(middle_count_nc,paste0(tumor_tissue[x,1],"_middle_normal_count_data.txt"),sep="\t",row.names=TRUE,quote=F)
  write.table(middle_colData,paste0(tumor_tissue[x,1],"_middle_normal_sampleInfomation.txt"),sep="\t",row.names=TRUE,quote=F)
  #old
  old_colData <- colData[which(colData$topage>=60),]
  old_count_nc <- count_nc[,rownames(old_colData)]
  identical(colnames(old_count_nc),rownames(old_colData))
  #TRUE
  old_colData <- old_colData[,-which(colnames(old_colData) == "topage")]
  old <- as.data.frame(table(old_colData$condition))
  # Normal  Tumor 
  # 101     70 
  group_sample_tolal[x,7] <- length(rownames(old_colData))
  group_sample_tolal[x,8] <- old[2,2]
  group_sample_tolal[x,9] <- old[1,2]
  write.table(old_count_nc,paste0(tumor_tissue[x,1],"_old_normal_count_data.txt"),sep="\t",row.names=TRUE,quote=F)
  write.table(old_colData,paste0(tumor_tissue[x,1],"_old_normal_sampleInfomation.txt"),sep="\t",row.names=TRUE,quote=F)
}
rownames(group_sample_tolal) <- tumor_tissue$tumor
group_sample_tolal[is.na(group_sample_tolal)]=0
write.csv(group_sample_tolal,"./Aging/data/UCSC/TCGA_RNAseq_data/青年分组DEseq2/group_sample_tolal.csv",row.names = T,quote = F)



#####Fig.S7 A top-----
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggprism)
age_group_total <- read.csv("D:/procedure/Aging/data/UCSC/TCGA_RNAseq_data/青年分组DEseq2/group_sample_tolal.csv",header = T,row.names = 1)
all_group_total1 <- age_group_total[,-c(1,4,7)]
all_group_total1$cancer <- rownames(all_group_total1)
all_group_total_need1 <- reshape2::melt(all_group_total1,id.vars = 'cancer')#配合reshape2使用
all_group_total_need1$age <- 1
all_group_total_need1$age[grep("young",all_group_total_need1$variable)] <- "young"
all_group_total_need1$age[grep("middle",all_group_total_need1$variable)] <- "middle"
all_group_total_need1$age[grep("old",all_group_total_need1$variable)] <- "old"
all_group_total_need1$class <- 1
all_group_total_need1$class[grep("cancer$",all_group_total_need1$variable)] <- "cancer"
all_group_total_need1$class[grep("normal$",all_group_total_need1$variable)] <- "normal"
all_group_total_need1$age <- factor(all_group_total_need1$age,levels = c("young","middle","old"))
ggplot(all_group_total_need1,aes(age,value,fill=class))+
  labs(x='time',cex=10,y='sample num')+theme_test(base_size = 20)+
  scale_fill_manual(values = c("#fedc5e","#8fd1e1"))+
  geom_bar(stat="identity",position="stack")+
  facet_grid(~cancer)+
  geom_text(aes(label=value),vjust=2,color="black")+
  theme(axis.text.x=element_text(angle=90, vjust=0.5))

######################################################################################
####Analysis of differences between cancer and normal under different age subgroups###
######################################################################################
library(DESeq2)
library(data.table)
library(stringr)
library(ggplot2)
library(DEGreport)
library(ggrepel)
library(plyr)
tumor_tissue <- as.data.frame(read.csv("D:/procedure/Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/tumor_tissue.csv",header = T,sep = ","))
DEGs_tolal <- data.frame(youth_tumor_up_DEGs=c(1),
                         youth_tumor_down_DEGs=c(1),
                         middle_tumor_up_DEGs=c(1),
                         middle_tumor_down_DEGs=c(1),
                         old_tumor_up_DEGs=c(1),
                         old_tumor_down_DEGs=c(1)
)
#Remove version number
ID.ver <- function(id){
  id2 <- unlist(strsplit(id,".",fixed = T))
  id3 <- id2[1]
  return(id3)
}
##Removal of genes with a percentage of zero expression greater than 25% of all samples
QC <- function(exp){
  num_NA<-length(which(exp==0))
  percentage <- num_NA/length(exp)
  return(percentage)
}
gene_annotation <- read.table("D:/procedure/Aging/data/UCSC/gene_id_type.txt",sep = "\t",header = T)
procoding_gene <- gene_annotation[which(gene_annotation$type == "protein_coding"),]
youth_DEGs_list <- list()
middle_DEGs_list <- list()
old_DEGs_list <- list()
for(x in 1:25){
  x=1
  setwd(paste0("D:/procedure/Aging/data/UCSC/TCGA_RNAseq_data/",tumor_tissue[x,1]))
  cancer_count <- as.data.frame(fread(paste0(tumor_tissue[x,1],"_RNAseq_count.txt"),header = T,sep = "\t"))
  cancer_clinical <- as.data.frame(fread(paste0(tumor_tissue[x,1],"_sample_info.txt"),header = T,sep = "\t"))
  normal_count <- as.data.frame(fread(paste0(tumor_tissue[x,2],"_RNAseq_count.txt"),header = T,sep = "\t"))
  normal_clinical <- as.data.frame(fread(paste0(tumor_tissue[x,2],"_sample_info.txt"),header = T,sep = "\t"))
  cancer_clinical$age <- ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=19,"10-19",
                                ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=29,"20-29",
                                       ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=39,"30-39",
                                              ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=49,"40-49",
                                                     ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=59,"50-59",
                                                            ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=69,"60-69",
                                                                   ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=79,"70-79",
                                                                          ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=89,"80-89",
                                                                                 ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=99,"90-99")))))))))
  
  cancer_clinical$group <- ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=40,"Youth_Age",
                                  ifelse(cancer_clinical$age_at_initial_pathologic_diagnosis<=60,"Middle_Age","Old_Age"
                                  ))
  normal_clinical$Age1 <- as.numeric(substring(normal_clinical$Age, 4))-4
  normal_clinical$group <- ifelse(normal_clinical$Age<=40,"Youth_Age",
                                  ifelse(normal_clinical$Age<=60,"Middle_Age","Old_Age"
                                  ))
  NA_Hd_ID <- cancer_clinical$sample[is.na(cancer_clinical$age_at_initial_pathologic_diagnosis)]
  if(length(NA_Hd_ID) != 0){
    index1 <- match(colnames(cancer_count[,-1]),NA_Hd_ID,nomatch = 0)
    index2 <- which(index1!=0)
    cancer_count <- cancer_count[,-index2]
    index3 <- match(cancer_clinical$sample,NA_Hd_ID,nomatch = 0)
    index4 <- which(index3!=0)
    cancer_clinical <- cancer_clinical[-index4,]
  }
  #Youth----
  youth_nomoral_count <- cbind(normal_count[,1],normal_count[,normal_clinical$SAMPID[which(normal_clinical$group =="Youth_Age")]])
  youth_normal_clinical <- normal_clinical[which(normal_clinical$group =="Youth_Age"),]
  youth_cancer_count <- cbind(cancer_count[,1],cancer_count[,cancer_clinical$sample[which(cancer_clinical$group =="Youth_Age")]])
  youth_cancer_clinical <- cancer_clinical[which(cancer_clinical$group =="Youth_Age"),]
  youth_data <- cbind(youth_cancer_count,youth_nomoral_count[,-1])
  youth_data[,-1] <- 2^youth_data[,-1]-1
  youth_data[,-1] <- round(youth_data[,-1])
  youth_percentage_count <- apply(youth_data[,-1],1,QC)
  youth_Filted_read_count <- youth_data[which(youth_percentage_count < 0.25),]
  gene_id <- apply(as.data.frame(youth_Filted_read_count[,1]), 1, ID.ver)
  youth_count_nc <- cbind(gene_id,youth_Filted_read_count[,-1])
  rownames(youth_count_nc) <- youth_count_nc[,1]
  youth_count_nc <- youth_count_nc[,-1]
  
  if(dim(table(youth_cancer_count$gender)) == 1){
    youth_a <- youth_cancer_clinical[,c(1,35)]
    youth_b <- youth_normal_clinical[,c(1,65)]
    colnames(youth_b) <- c("sample","age")
    youth_age_clinical <- rbind(youth_a,youth_b)
    condition <- factor(ifelse(substr(colnames(youth_count_nc),14,15) == "01","Tumor","Normal"))
    youth_colData <- data.frame(row.names = colnames(youth_count_nc),condition,youth_age_clinical$age)
    colnames(youth_colData)[2] <- "Age"
    youth_colData$Age <- str_replace(youth_colData$Age,"-","_")
    youth_colData$Age <- as.factor(youth_colData$Age)
    youth_dds <- DESeqDataSetFromMatrix(countData = youth_count_nc,colData = youth_colData,design = ~ condition + Age )
    youth_dds1 <- DESeq(youth_dds,test="LRT",reduced = ~Age)
    youth_res <- results(youth_dds1,contrast=c("condition","Tumor","Normal"),
                         independentFiltering=TRUE,alpha=0.05,pAdjustMethod="BH")
    
  }else{
    youth_cancer_clinical$gender <- ifelse(youth_cancer_clinical$gender == "FEMALE","2","1")
    youth_a <- youth_cancer_clinical[,c(1,5,35)]
    youth_b <- youth_normal_clinical[,c(1,64,65)]
    colnames(youth_b) <- c("sample","gender","age")
    youth_age_clinical <- rbind(youth_a,youth_b)
    condition <- factor(ifelse(substr(colnames(youth_count_nc),14,15) == "01","Tumor","Normal"))
    youth_colData <- data.frame(row.names = colnames(youth_count_nc),condition,youth_age_clinical$age,youth_age_clinical$gender)
    colnames(youth_colData)[2] <- "Age"
    colnames(youth_colData)[3] <- "gender"
    youth_colData$Age <- str_replace(youth_colData$Age,"-","_")
    youth_colData$Age <- as.factor(youth_colData$Age)
    youth_colData$gender <- as.factor(youth_colData$gender)
    youth_dds <- DESeqDataSetFromMatrix(countData = youth_count_nc,colData = youth_colData,design = ~ condition + Age + gender )
    youth_dds1 <- DESeq(youth_dds,test="LRT",reduced = ~Age + gender)
    youth_res <- results(youth_dds1,contrast=c("condition","Tumor","Normal"),
                         independentFiltering=TRUE,alpha=0.05,pAdjustMethod="BH")
  }
  youth_res1 <- data.frame(youth_res, stringsAsFactors = FALSE, check.names = FALSE)
  youth_res1_up<- youth_res1[which(youth_res1$log2FoldChange >= 1 & youth_res1$padj < 0.01),]
  youth_res1_down<- youth_res1[which(youth_res1$log2FoldChange <= -1 & youth_res1$padj < 0.01),]
  youth_res1_total <- rbind(youth_res1_up,youth_res1_down)
  write.table(youth_res1_up,"youth_gene_up_infor.txt",col.names = T,row.names = T,quote = F,sep = "\t")
  write.table(rownames(youth_res1_up),"youth_up_gene.txt",sep="\t",quote = F,row.names = F,col.names = F)
  write.table(youth_res1_down,"youth_gene_down_infor.txt",col.names = T,row.names = T,quote = F,sep = "\t")
  write.table(rownames(youth_res1_down),"youth_down_gene.txt",sep="\t",quote = F,row.names = F,col.names = F)
  youth_up_list <- list(rownames(youth_res1_up),rownames(youth_res1_down))
  names(youth_up_list) <- c(paste0(tumor_tissue[x,1],"_youth_DEGs_up"),paste0(tumor_tissue[x,1],"_youth_DEGs_down")) 
  youth_DEGs_list <- c(youth_DEGs_list,youth_up_list)
  
  #middle----
  middle_nomoral_count <- cbind(normal_count[,1],normal_count[,normal_clinical$SAMPID[which(normal_clinical$group =="Middle_Age")]])
  middle_normal_clinical <- normal_clinical[which(normal_clinical$group =="Middle_Age"),]
  middle_cancer_count <- cbind(cancer_count[,1],cancer_count[,cancer_clinical$sample[which(cancer_clinical$group =="Middle_Age")]])
  middle_cancer_clinical <- cancer_clinical[which(cancer_clinical$group =="Middle_Age"),]
  middle_data <- cbind(middle_cancer_count,middle_nomoral_count[,-1])
  middle_data[,-1] <- 2^middle_data[,-1]-1
  middle_data[,-1] <- round(middle_data[,-1])
  middle_percentage_count <- apply(middle_data[,-1],1,QC)
  middle_Filted_read_count <- middle_data[which(middle_percentage_count < 0.25),]
  
  gene_id <- apply(as.data.frame(middle_Filted_read_count[,1]), 1, ID.ver)
  middle_count_nc <- cbind(gene_id,middle_Filted_read_count[,-1])
  rownames(middle_count_nc) <- middle_count_nc[,1]
  middle_count_nc <- middle_count_nc[,-1]
  
  if(dim(table(middle_cancer_count$gender)) == 1){
    middle_a <- middle_cancer_clinical[,c(1,35)]
    middle_b <- middle_normal_clinical[,c(1,65)]
    colnames(middle_b) <- c("sample","age")
    middle_age_clinical <- rbind(middle_a,middle_b)
    condition <- factor(ifelse(substr(colnames(middle_count_nc),14,15) == "01","Tumor","Normal"))
    middle_colData <- data.frame(row.names = colnames(middle_count_nc),condition,middle_age_clinical$age)
    colnames(middle_colData)[2] <- "Age"
    middle_colData$Age <- str_replace(middle_colData$Age,"-","_")
    middle_colData$Age <- as.factor(middle_colData$Age)
    middle_dds <- DESeqDataSetFromMatrix(countData = middle_count_nc,colData = middle_colData,design = ~ condition + Age )
    middle_dds1 <- DESeq(middle_dds,test="LRT",reduced = ~Age)
    middle_res <- results(middle_dds1,contrast=c("condition","Tumor","Normal"),
                          independentFiltering=TRUE,alpha=0.05,pAdjustMethod="BH")
    
  }else{
    middle_cancer_clinical$gender <- ifelse(middle_cancer_clinical$gender == "FEMALE","2","1")
    middle_a <- middle_cancer_clinical[,c(1,5,35)]
    middle_b <- middle_normal_clinical[,c(1,64,65)]
    colnames(middle_b) <- c("sample","gender","age")
    middle_age_clinical <- rbind(middle_a,middle_b)
    condition <- factor(ifelse(substr(colnames(middle_count_nc),14,15) == "01","Tumor","Normal"))
    middle_colData <- data.frame(row.names = colnames(middle_count_nc),condition,middle_age_clinical$age,middle_age_clinical$gender)
    colnames(middle_colData)[2] <- "Age"
    colnames(middle_colData)[3] <- "gender"
    middle_colData$Age <- str_replace(middle_colData$Age,"-","_")
    middle_colData$Age <- as.factor(middle_colData$Age)
    middle_colData$gender <- as.factor(middle_colData$gender)
    middle_dds <- DESeqDataSetFromMatrix(countData = middle_count_nc,colData = middle_colData,design = ~ condition + Age + gender )
    middle_dds1 <- DESeq(middle_dds,test="LRT",reduced = ~Age + gender)
    middle_res <- results(middle_dds1,contrast=c("condition","Tumor","Normal"),
                          independentFiltering=TRUE,alpha=0.05,pAdjustMethod="BH")
  }
  middle_res1 <- data.frame(middle_res, stringsAsFactors = FALSE, check.names = FALSE)
  middle_res1_up<- youth_res1[which(middle_res1$log2FoldChange >= 1 & middle_res1$padj < 0.01),]
  middle_res1_down<- youth_res1[which(middle_res1$log2FoldChange <= -1 & middle_res1$padj < 0.01),]
  middle_res1_total <- rbind(middle_res1_up,middle_res1_down)
  write.table(middle_res1_up,"middle_gene_up_infor.txt",col.names = T,row.names = T,quote = F,sep = "\t")
  write.table(rownames(middle_res1_up),"middle_up_gene.txt",sep="\t",quote = F,row.names = F,col.names = F)
  write.table(middle_res1_down,"middle_gene_down_infor.txt",col.names = T,row.names = T,quote = F,sep = "\t")
  write.table(rownames(middle_res1_down),"middle_down_gene.txt",sep="\t",quote = F,row.names = F,col.names = F)
  middle_up_list <- list(rownames(middle_res1_up),rownames(middle_res1_down))
  names(middle_up_list) <- c(paste0(tumor_tissue[x,1],"_middle_DEGs_up"),paste0(tumor_tissue[x,1],"_middle_DEGs_down")) 
  middle_DEGs_list <- c(middle_DEGs_list,middle_up_list)

  #Old----
  old_nomoral_count <- cbind(normal_count[,1],normal_count[,normal_clinical$SAMPID[which(normal_clinical$group =="Old_Age")]])
  old_normal_clinical <- normal_clinical[which(normal_clinical$group =="Old_Age"),]
  old_cancer_count <- cbind(cancer_count[,1],cancer_count[,cancer_clinical$sample[which(cancer_clinical$group =="Old_Age")]])
  old_cancer_clinical <- cancer_clinical[which(cancer_clinical$group =="Old_Age"),]
  old_data <- cbind(old_cancer_count,old_nomoral_count[,-1])
  old_data[,-1] <- 2^old_data[,-1]-1
  old_data[,-1] <- round(old_data[,-1])
  old_percentage_count <- apply(old_data[,-1],1,QC)
  old_Filted_read_count <- old_data[which(old_percentage_count < 0.25),]
  
  gene_id <- apply(as.data.frame(old_Filted_read_count[,1]), 1, ID.ver)
  old_count_nc <- cbind(gene_id,old_Filted_read_count[,-1])
  rownames(old_count_nc) <- old_count_nc[,1]
  old_count_nc <- old_count_nc[,-1]
  if(dim(table(old_cancer_count$gender)) == 1){
    old_a <- old_cancer_clinical[,c(1,35)]
    old_b <- old_normal_clinical[,c(1,65)]
    colnames(old_b) <- c("sample","age")
    old_age_clinical <- rbind(old_a,old_b)
    condition <- factor(ifelse(substr(colnames(old_count_nc),14,15) == "01","Tumor","Normal"))
    old_colData <- data.frame(row.names = colnames(old_count_nc),condition,old_age_clinical$age)
    colnames(old_colData)[2] <- "Age"
    old_colData$Age <- str_replace(old_colData$Age,"-","_")
    old_colData$Age <- as.factor(old_colData$Age)
    old_dds <- DESeqDataSetFromMatrix(countData = old_count_nc,colData = old_colData,design = ~ condition + Age )
    old_dds1 <- DESeq(old_dds,test="LRT",reduced = ~Age)
    old_res <- results(old_dds1,contrast=c("condition","Tumor","Normal"),
                       independentFiltering=TRUE,alpha=0.05,pAdjustMethod="BH")
  }else{
    old_cancer_clinical$gender <- ifelse(old_cancer_clinical$gender == "FEMALE","2","1")
    old_a <- old_cancer_clinical[,c(1,5,35)]
    old_b <- old_normal_clinical[,c(1,64,65)]
    colnames(old_b) <- c("sample","gender","age")
    old_age_clinical <- rbind(old_a,old_b)
    condition <- factor(ifelse(substr(colnames(old_count_nc),14,15) == "01","Tumor","Normal"))
    old_colData <- data.frame(row.names = colnames(old_count_nc),condition,old_age_clinical$age,old_age_clinical$gender)
    colnames(old_colData)[2] <- "Age"
    colnames(old_colData)[3] <- "gender"
    old_colData$Age <- str_replace(old_colData$Age,"-","_")
    old_colData$Age <- as.factor(old_colData$Age)
    old_colData$gender <- as.factor(old_colData$gender)
    old_dds <- DESeqDataSetFromMatrix(countData = old_count_nc,colData = old_colData,design = ~ condition + Age + gender )
    old_dds1 <- DESeq(old_dds,test="LRT",reduced = ~Age + gender)
    old_res <- results(old_dds1,contrast=c("condition","Tumor","Normal"),
                       independentFiltering=TRUE,alpha=0.05,pAdjustMethod="BH")
  }
  old_res1 <- data.frame(old_res, stringsAsFactors = FALSE, check.names = FALSE)
  old_res1_up<- youth_res1[which(old_res1$log2FoldChange >= 1 & old_res1$padj < 0.01),]
  old_res1_down<- youth_res1[which(old_res1$log2FoldChange <= -1 & old_res1$padj < 0.01),]
  old_res1_total <- rbind(old_res1_up,old_res1_down)
  write.table(old_res1_up,"old_gene_up_infor.txt",col.names = T,row.names = T,quote = F,sep = "\t")
  write.table(rownames(old_res1_up),"old_up_gene.txt",sep="\t",quote = F,row.names = F,col.names = F)
  write.table(old_res1_down,"old_gene_down_infor.txt",col.names = T,row.names = T,quote = F,sep = "\t")
  write.table(rownames(old_res1_down),"old_down_gene.txt",sep="\t",quote = F,row.names = F,col.names = F)
  old_up_list <- list(rownames(old_res1_up),rownames(old_res1_down))
  names(old_up_list) <- c(paste0(tumor_tissue[x,1],"_old_DEGs_up"),paste0(tumor_tissue[x,1],"_old_DEGs_down")) 
  old_DEGs_list <- c(old_DEGs_list,old_up_list)
}
save(youth_DEGs_list,file = "youth_DEGs_list")
save(middle_DEGs_list,file = "middle_DEGs_list")
save(old_DEGs_list,file = "old_DEGs_list")

####Fig S7A bottom
paixu <- as.data.frame(read.csv("D:/procedure/Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/paixu.csv",header = T,sep = ","))
load("./Aging/data/UCSC/TCGA_RNAseq_data/青年分组DEseq2/young/youth_DEGs_total.Rdata")
load("./Aging/data/UCSC/TCGA_RNAseq_data/青年分组DEseq2/middle/middle_DEGs_total.Rdata")
load("./Aging/data/UCSC/TCGA_RNAseq_data/青年分组DEseq2/old/old_DEGs_total.Rdata")
youth_Cacncer_DEGs <- youth_DEGs_total[,c(1,2)]
middle_Cacncer_DEGs <- middle_DEGs_total[,c(1,2)]
old_Cacncer_DEGs <- old_DEGs_total[,c(1,2)]
######---分组柱状图#####
youth_Cacncer_DEGs$cancer <- rownames(youth_Cacncer_DEGs)
middle_Cacncer_DEGs$cancer <- rownames(middle_Cacncer_DEGs)
old_Cacncer_DEGs$cancer <- rownames(old_Cacncer_DEGs)
youth_Cacncer_DEGs <- reshape2::melt(youth_Cacncer_DEGs)
middle_Cacncer_DEGs <- reshape2::melt(middle_Cacncer_DEGs)
old_Cacncer_DEGs <- reshape2::melt(old_Cacncer_DEGs)
youth_Cacncer_DEGs <- rbind(youth_Cacncer_DEGs,c("UCS","tumor_up_DEGs",0),c("UCS","tumor_down_DEGs",0),c("PRAD","tumor_up_DEGs",0),c("PRAD","tumor_down_DEGs",0))
youth_Cacncer_DEGs$group <- "young"
middle_Cacncer_DEGs$group <- "middle"
old_Cacncer_DEGs$group <- "old"
data <- rbind(youth_Cacncer_DEGs,middle_Cacncer_DEGs,old_Cacncer_DEGs)
data$cancer <- factor(data$cancer,levels = paixu$tumor)
data$group <- factor(data$group,levels = c("young","middle","old"))
data$value <- as.numeric(data$value)
ggplot(data,aes(group,value,fill=variable))+
  labs(x='time',cex=10,y='sample num')+theme_test(base_size = 20)+
  scale_fill_manual(values = c("#CD403E","#4D91AE"))+
  geom_bar(stat="identity",position="stack")+### position = stack按照数值填充；fill按照百分比填充  identity
  facet_grid(~cancer)+
  theme(axis.text.x=element_text(angle=90, vjust=0.5))

########################################################
#######Overlap of Youth cancer-DEGs and aging DEGs######
########################################################
#UTUA
library(pheatmap)
library(ggthemes)
Up_up_data <- youth_DEGs_total[,c(1,3,5)]
colnames(Up_up_data) <- c("UpInTumor","UpInAge","Overlap")
bk <- c(seq(0,4500,by=100))
pheatmap(Up_up_data,
         scale = "none",
         color = c(colorRampPalette(colors = c("white","#3CCF4E"))(length(bk))),
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 12,
         cellwidth = 28,
         cellheight = 20,
         legend_breaks=seq(0),
         display_numbers = TRUE,
         number_format = "%.0f",
         angle_col = 45,
         legend = FALSE,
         fontsize_col = 10,
         
)
Up_up_Pvalue <- data.frame(cancer=rownames(youth_DEGs_total),p_value=youth_DEGs_total$up_age_tumor_up_DEGs_P)
Up_up_Pvalue$p_value <- -1*log10(Up_up_Pvalue$p_value)
Up_up_Pvalue$tp <- ifelse(Up_up_Pvalue$p_value < -log10(0.05),"",
                          ifelse(Up_up_Pvalue$p_value >= -log10(0.05) & Up_up_Pvalue$p_value<=-log10(0.01),"*",
                                 ifelse(Up_up_Pvalue$p_value >= -log10(0.01) & Up_up_Pvalue$p_value<=-log10(0.001),"**","***")))
Up_up_Pvalue$p_value <- ifelse(Up_up_Pvalue$p_value > 20,20,Up_up_Pvalue$p_value)
Up_up_Pvalue$cancer <- factor(Up_up_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Up_up_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c(rep("#3CCF4E",22))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())
#DTDA
Down_down_data <- youth_DEGs_total[,c(2,4,11)]
colnames(Down_down_data) <- c("DownInTumor","DownInAge","Overlap")
bk <- c(seq(0,4500,by=100))
p2 <- pheatmap(Down_down_data,
               scale = "none",
               color = c(colorRampPalette(colors = c("white","#E55807"))(length(bk))),
               cluster_rows = F,
               cluster_cols = F,
               fontsize = 12,
               cellwidth = 28,
               cellheight = 20,
               legend_breaks=seq(0),
               display_numbers = TRUE,
               number_format = "%.0f",
               angle_col = 45,
               legend = FALSE,
               fontsize_col = 10,
               show_rownames = F,
)
Down_down_Pvalue <- data.frame(cancer=rownames(youth_DEGs_total),p_value=youth_DEGs_total$down_age_tumor_down_DEGs_P)
Down_down_Pvalue$p_value <- -1*log10(Down_down_Pvalue$p_value)
Down_down_Pvalue$tp <- ifelse(Down_down_Pvalue$p_value < -log10(0.05),"",
                              ifelse(Down_down_Pvalue$p_value >= -log10(0.05) & Down_down_Pvalue$p_value<=-log10(0.01),"*",
                                     ifelse(Down_down_Pvalue$p_value >= -log10(0.01) & Down_down_Pvalue$p_value<=-log10(0.001),"**","***")))
Down_down_Pvalue$p_value <- ifelse(Down_down_Pvalue$p_value > 20,20,Down_down_Pvalue$p_value)
Down_down_Pvalue$cancer <- factor(Down_down_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Down_down_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c(rep("#E55807",22)))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())

#DTUA
Down_up_data <- youth_DEGs_total[,c(2,3,7)]
colnames(Down_up_data) <- c("DownInTumor","UpInAge","Overlap")
bk <- c(seq(0,3500,by=100))
p3 <- pheatmap(Down_up_data,
               scale = "none",
               color = c(colorRampPalette(colors = c("white","#7E1717"))(length(bk))),
               cluster_rows = F,
               cluster_cols = F,
               fontsize = 12,
               cellwidth = 28,
               cellheight = 20,
               legend_breaks=seq(0),
               display_numbers = TRUE,
               number_format = "%.0f",
               angle_col = 45,
               legend = FALSE,
               fontsize_col = 10,
               show_rownames = F,
)
Down_up_Pvalue <- data.frame(cancer=rownames(youth_DEGs_total),p_value=youth_DEGs_total$up_age_tumor_down_DEGs_P)
Down_up_Pvalue$p_value <- -1*log10(Down_up_Pvalue$p_value)
Down_up_Pvalue$tp <- ifelse(Down_up_Pvalue$p_value < -log10(0.05),"",
                            ifelse(Down_up_Pvalue$p_value >= -log10(0.05) & Down_up_Pvalue$p_value<=-log10(0.01),"*",
                                   ifelse(Down_up_Pvalue$p_value >= -log10(0.01) & Down_up_Pvalue$p_value<=-log10(0.001),"**","***")))
Down_up_Pvalue$p_value <- ifelse(Down_up_Pvalue$p_value > 20,20,Down_up_Pvalue$p_value)
Down_up_Pvalue$cancer <- factor(Down_up_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Down_up_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c(rep("#7E1717",22)))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())

#UTDA
Up_down_data <- youth_DEGs_total[,c(1,4,9)]
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
Up_down_Pvalue <- data.frame(cancer=rownames(youth_DEGs_total),p_value=youth_DEGs_total$down_age_tumor_up_DEGs_P)
Up_down_Pvalue$p_value <- -1*log10(Up_down_Pvalue$p_value)
Up_down_Pvalue$tp <- ifelse(Up_down_Pvalue$p_value < -log10(0.05),"",
                            ifelse(Up_down_Pvalue$p_value >= -log10(0.05) & Up_down_Pvalue$p_value<=-log10(0.01),"*",
                                   ifelse(Up_down_Pvalue$p_value >= -log10(0.01) & Up_down_Pvalue$p_value<=-log10(0.001),"**","***")))
Up_down_Pvalue$p_value <- ifelse(Up_down_Pvalue$p_value > 20,20,Up_down_Pvalue$p_value)
Up_down_Pvalue$cancer <- factor(Up_down_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Up_down_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c(rep("#068DA9",22)))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())


########################################################
#######Overlap of Youth cancer-DEGs and aging DEGs######
########################################################
#UTUA
library(pheatmap)
library(ggthemes)
Up_up_data <- youth_DEGs_total[,c(1,3,5)]
colnames(Up_up_data) <- c("UpInTumor","UpInAge","Overlap")
bk <- c(seq(0,4500,by=100))
pheatmap(Up_up_data,
         scale = "none",
         color = c(colorRampPalette(colors = c("white","#3CCF4E"))(length(bk))),
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 12,
         cellwidth = 28,
         cellheight = 20,
         legend_breaks=seq(0),
         display_numbers = TRUE,
         number_format = "%.0f",
         angle_col = 45,
         legend = FALSE,
         fontsize_col = 10,
         
)
Up_up_Pvalue <- data.frame(cancer=rownames(youth_DEGs_total),p_value=youth_DEGs_total$up_age_tumor_up_DEGs_P)
Up_up_Pvalue$p_value <- -1*log10(Up_up_Pvalue$p_value)
Up_up_Pvalue$tp <- ifelse(Up_up_Pvalue$p_value < -log10(0.05),"",
                          ifelse(Up_up_Pvalue$p_value >= -log10(0.05) & Up_up_Pvalue$p_value<=-log10(0.01),"*",
                                 ifelse(Up_up_Pvalue$p_value >= -log10(0.01) & Up_up_Pvalue$p_value<=-log10(0.001),"**","***")))
Up_up_Pvalue$p_value <- ifelse(Up_up_Pvalue$p_value > 20,20,Up_up_Pvalue$p_value)
Up_up_Pvalue$cancer <- factor(Up_up_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Up_up_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c(rep("#3CCF4E",22))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())
#DTDA
Down_down_data <- youth_DEGs_total[,c(2,4,11)]
colnames(Down_down_data) <- c("DownInTumor","DownInAge","Overlap")
bk <- c(seq(0,4500,by=100))
p2 <- pheatmap(Down_down_data,
               scale = "none",
               color = c(colorRampPalette(colors = c("white","#E55807"))(length(bk))),
               cluster_rows = F,
               cluster_cols = F,
               fontsize = 12,
               cellwidth = 28,
               cellheight = 20,
               legend_breaks=seq(0),
               display_numbers = TRUE,
               number_format = "%.0f",
               angle_col = 45,
               legend = FALSE,
               fontsize_col = 10,
               show_rownames = F,
)
Down_down_Pvalue <- data.frame(cancer=rownames(youth_DEGs_total),p_value=youth_DEGs_total$down_age_tumor_down_DEGs_P)
Down_down_Pvalue$p_value <- -1*log10(Down_down_Pvalue$p_value)
Down_down_Pvalue$tp <- ifelse(Down_down_Pvalue$p_value < -log10(0.05),"",
                              ifelse(Down_down_Pvalue$p_value >= -log10(0.05) & Down_down_Pvalue$p_value<=-log10(0.01),"*",
                                     ifelse(Down_down_Pvalue$p_value >= -log10(0.01) & Down_down_Pvalue$p_value<=-log10(0.001),"**","***")))
Down_down_Pvalue$p_value <- ifelse(Down_down_Pvalue$p_value > 20,20,Down_down_Pvalue$p_value)
Down_down_Pvalue$cancer <- factor(Down_down_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Down_down_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c(rep("#E55807",22)))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())

#DTUA
Down_up_data <- youth_DEGs_total[,c(2,3,7)]
colnames(Down_up_data) <- c("DownInTumor","UpInAge","Overlap")
bk <- c(seq(0,3500,by=100))
p3 <- pheatmap(Down_up_data,
               scale = "none",
               color = c(colorRampPalette(colors = c("white","#7E1717"))(length(bk))),
               cluster_rows = F,
               cluster_cols = F,
               fontsize = 12,
               cellwidth = 28,
               cellheight = 20,
               legend_breaks=seq(0),
               display_numbers = TRUE,
               number_format = "%.0f",
               angle_col = 45,
               legend = FALSE,
               fontsize_col = 10,
               show_rownames = F,
)
Down_up_Pvalue <- data.frame(cancer=rownames(youth_DEGs_total),p_value=youth_DEGs_total$up_age_tumor_down_DEGs_P)
Down_up_Pvalue$p_value <- -1*log10(Down_up_Pvalue$p_value)
Down_up_Pvalue$tp <- ifelse(Down_up_Pvalue$p_value < -log10(0.05),"",
                            ifelse(Down_up_Pvalue$p_value >= -log10(0.05) & Down_up_Pvalue$p_value<=-log10(0.01),"*",
                                   ifelse(Down_up_Pvalue$p_value >= -log10(0.01) & Down_up_Pvalue$p_value<=-log10(0.001),"**","***")))
Down_up_Pvalue$p_value <- ifelse(Down_up_Pvalue$p_value > 20,20,Down_up_Pvalue$p_value)
Down_up_Pvalue$cancer <- factor(Down_up_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Down_up_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c(rep("#7E1717",22)))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())

#UTDA
Up_down_data <- youth_DEGs_total[,c(1,4,9)]
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
Up_down_Pvalue <- data.frame(cancer=rownames(youth_DEGs_total),p_value=youth_DEGs_total$down_age_tumor_up_DEGs_P)
Up_down_Pvalue$p_value <- -1*log10(Up_down_Pvalue$p_value)
Up_down_Pvalue$tp <- ifelse(Up_down_Pvalue$p_value < -log10(0.05),"",
                            ifelse(Up_down_Pvalue$p_value >= -log10(0.05) & Up_down_Pvalue$p_value<=-log10(0.01),"*",
                                   ifelse(Up_down_Pvalue$p_value >= -log10(0.01) & Up_down_Pvalue$p_value<=-log10(0.001),"**","***")))
Up_down_Pvalue$p_value <- ifelse(Up_down_Pvalue$p_value > 20,20,Up_down_Pvalue$p_value)
Up_down_Pvalue$cancer <- factor(Up_down_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Up_down_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c(rep("#068DA9",22)))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())

########################################################
#######Overlap of middle cancer-DEGs and aging DEGs######
########################################################
#UTUA
library(pheatmap)
library(ggthemes)
Up_up_data <- middle_DEGs_total[,c(1,3,5)]
colnames(Up_up_data) <- c("UpInTumor","UpInAge","Overlap")
bk <- c(seq(0,4500,by=100))
pheatmap(Up_up_data,
         scale = "none",
         color = c(colorRampPalette(colors = c("white","#3CCF4E"))(length(bk))),
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 12,
         cellwidth = 28,
         cellheight = 20,
         legend_breaks=seq(0),
         display_numbers = TRUE,
         number_format = "%.0f",
         angle_col = 45,
         legend = FALSE,
         fontsize_col = 10,
         
)
Up_up_Pvalue <- data.frame(cancer=rownames(middle_DEGs_total),p_value=middle_DEGs_total$up_age_tumor_up_DEGs_P)
Up_up_Pvalue$p_value <- -1*log10(Up_up_Pvalue$p_value)
Up_up_Pvalue$tp <- ifelse(Up_up_Pvalue$p_value < -log10(0.05),"",
                          ifelse(Up_up_Pvalue$p_value >= -log10(0.05) & Up_up_Pvalue$p_value<=-log10(0.01),"*",
                                 ifelse(Up_up_Pvalue$p_value >= -log10(0.01) & Up_up_Pvalue$p_value<=-log10(0.001),"**","***")))
Up_up_Pvalue$p_value <- ifelse(Up_up_Pvalue$p_value > 20,20,Up_up_Pvalue$p_value)
Up_up_Pvalue$cancer <- factor(Up_up_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Up_up_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c(rep("#3CCF4E",24))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())
#DTDA
Down_down_data <- middle_DEGs_total[,c(2,4,11)]
colnames(Down_down_data) <- c("DownInTumor","DownInAge","Overlap")
bk <- c(seq(0,4500,by=100))
p2 <- pheatmap(Down_down_data,
               scale = "none",
               color = c(colorRampPalette(colors = c("white","#E55807"))(length(bk))),
               cluster_rows = F,
               cluster_cols = F,
               fontsize = 12,
               cellwidth = 28,
               cellheight = 20,
               legend_breaks=seq(0),
               display_numbers = TRUE,
               number_format = "%.0f",
               angle_col = 45,
               legend = FALSE,
               fontsize_col = 10,
               show_rownames = F,
)
Down_down_Pvalue <- data.frame(cancer=rownames(middle_DEGs_total),p_value=middle_DEGs_total$down_age_tumor_down_DEGs_P)
Down_down_Pvalue$p_value <- -1*log10(Down_down_Pvalue$p_value)
Down_down_Pvalue$tp <- ifelse(Down_down_Pvalue$p_value < -log10(0.05),"",
                              ifelse(Down_down_Pvalue$p_value >= -log10(0.05) & Down_down_Pvalue$p_value<=-log10(0.01),"*",
                                     ifelse(Down_down_Pvalue$p_value >= -log10(0.01) & Down_down_Pvalue$p_value<=-log10(0.001),"**","***")))
Down_down_Pvalue$p_value <- ifelse(Down_down_Pvalue$p_value > 20,20,Down_down_Pvalue$p_value)
Down_down_Pvalue$cancer <- factor(Down_down_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Down_down_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c(rep("#E55807",24)))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())

#DTUA
Down_up_data <- middle_DEGs_total[,c(2,3,7)]
colnames(Down_up_data) <- c("DownInTumor","UpInAge","Overlap")
bk <- c(seq(0,3500,by=100))
p3 <- pheatmap(Down_up_data,
               scale = "none",
               color = c(colorRampPalette(colors = c("white","#7E1717"))(length(bk))),
               cluster_rows = F,
               cluster_cols = F,
               fontsize = 12,
               cellwidth = 28,
               cellheight = 20,
               legend_breaks=seq(0),
               display_numbers = TRUE,
               number_format = "%.0f",
               angle_col = 45,
               legend = FALSE,
               fontsize_col = 10,
               show_rownames = F,
)
Down_up_Pvalue <- data.frame(cancer=rownames(middle_DEGs_total),p_value=middle_DEGs_total$up_age_tumor_down_DEGs_P)
Down_up_Pvalue$p_value <- -1*log10(Down_up_Pvalue$p_value)
Down_up_Pvalue$tp <- ifelse(Down_up_Pvalue$p_value < -log10(0.05),"",
                            ifelse(Down_up_Pvalue$p_value >= -log10(0.05) & Down_up_Pvalue$p_value<=-log10(0.01),"*",
                                   ifelse(Down_up_Pvalue$p_value >= -log10(0.01) & Down_up_Pvalue$p_value<=-log10(0.001),"**","***")))
Down_up_Pvalue$p_value <- ifelse(Down_up_Pvalue$p_value > 20,20,Down_up_Pvalue$p_value)
Down_up_Pvalue$cancer <- factor(Down_up_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Down_up_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c(rep("#7E1717",24)))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())

#UTDA
Up_down_data <- middle_DEGs_total[,c(1,4,9)]
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
Up_down_Pvalue <- data.frame(cancer=rownames(middle_DEGs_total),p_value=middle_DEGs_total$down_age_tumor_up_DEGs_P)
Up_down_Pvalue$p_value <- -1*log10(Up_down_Pvalue$p_value)
Up_down_Pvalue$tp <- ifelse(Up_down_Pvalue$p_value < -log10(0.05),"",
                            ifelse(Up_down_Pvalue$p_value >= -log10(0.05) & Up_down_Pvalue$p_value<=-log10(0.01),"*",
                                   ifelse(Up_down_Pvalue$p_value >= -log10(0.01) & Up_down_Pvalue$p_value<=-log10(0.001),"**","***")))
Up_down_Pvalue$p_value <- ifelse(Up_down_Pvalue$p_value > 20,20,Up_down_Pvalue$p_value)
Up_down_Pvalue$cancer <- factor(Up_down_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Up_down_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c(rep("#068DA9",24)))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())





########################################################
#######Overlap of old cancer-DEGs and aging DEGs######
########################################################
#UTUA
library(pheatmap)
library(ggthemes)
Up_up_data <- old_DEGs_total[,c(1,3,5)]
colnames(Up_up_data) <- c("UpInTumor","UpInAge","Overlap")
bk <- c(seq(0,4500,by=100))
pheatmap(Up_up_data,
         scale = "none",
         color = c(colorRampPalette(colors = c("white","#3CCF4E"))(length(bk))),
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 12,
         cellwidth = 28,
         cellheight = 20,
         legend_breaks=seq(0),
         display_numbers = TRUE,
         number_format = "%.0f",
         angle_col = 45,
         legend = FALSE,
         fontsize_col = 10,
         
)
Up_up_Pvalue <- data.frame(cancer=rownames(old_DEGs_total),p_value=old_DEGs_total$up_age_tumor_up_DEGs_P)
Up_up_Pvalue$p_value <- -1*log10(Up_up_Pvalue$p_value)
Up_up_Pvalue$tp <- ifelse(Up_up_Pvalue$p_value < -log10(0.05),"",
                          ifelse(Up_up_Pvalue$p_value >= -log10(0.05) & Up_up_Pvalue$p_value<=-log10(0.01),"*",
                                 ifelse(Up_up_Pvalue$p_value >= -log10(0.01) & Up_up_Pvalue$p_value<=-log10(0.001),"**","***")))
Up_up_Pvalue$p_value <- ifelse(Up_up_Pvalue$p_value > 20,20,Up_up_Pvalue$p_value)
Up_up_Pvalue$cancer <- factor(Up_up_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Up_up_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c(rep("#3CCF4E",24))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())
#DTDA
Down_down_data <- old_DEGs_total[,c(2,4,11)]
colnames(Down_down_data) <- c("DownInTumor","DownInAge","Overlap")
bk <- c(seq(0,4500,by=100))
p2 <- pheatmap(Down_down_data,
               scale = "none",
               color = c(colorRampPalette(colors = c("white","#E55807"))(length(bk))),
               cluster_rows = F,
               cluster_cols = F,
               fontsize = 12,
               cellwidth = 28,
               cellheight = 20,
               legend_breaks=seq(0),
               display_numbers = TRUE,
               number_format = "%.0f",
               angle_col = 45,
               legend = FALSE,
               fontsize_col = 10,
               show_rownames = F,
)
Down_down_Pvalue <- data.frame(cancer=rownames(old_DEGs_total),p_value=old_DEGs_total$down_age_tumor_down_DEGs_P)
Down_down_Pvalue$p_value <- -1*log10(Down_down_Pvalue$p_value)
Down_down_Pvalue$tp <- ifelse(Down_down_Pvalue$p_value < -log10(0.05),"",
                              ifelse(Down_down_Pvalue$p_value >= -log10(0.05) & Down_down_Pvalue$p_value<=-log10(0.01),"*",
                                     ifelse(Down_down_Pvalue$p_value >= -log10(0.01) & Down_down_Pvalue$p_value<=-log10(0.001),"**","***")))
Down_down_Pvalue$p_value <- ifelse(Down_down_Pvalue$p_value > 20,20,Down_down_Pvalue$p_value)
Down_down_Pvalue$cancer <- factor(Down_down_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Down_down_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c(rep("#E55807",24)))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())

#DTUA
Down_up_data <- old_DEGs_total[,c(2,3,7)]
colnames(Down_up_data) <- c("DownInTumor","UpInAge","Overlap")
bk <- c(seq(0,3500,by=100))
p3 <- pheatmap(Down_up_data,
               scale = "none",
               color = c(colorRampPalette(colors = c("white","#7E1717"))(length(bk))),
               cluster_rows = F,
               cluster_cols = F,
               fontsize = 12,
               cellwidth = 28,
               cellheight = 20,
               legend_breaks=seq(0),
               display_numbers = TRUE,
               number_format = "%.0f",
               angle_col = 45,
               legend = FALSE,
               fontsize_col = 10,
               show_rownames = F,
)
Down_up_Pvalue <- data.frame(cancer=rownames(old_DEGs_total),p_value=old_DEGs_total$up_age_tumor_down_DEGs_P)
Down_up_Pvalue$p_value <- -1*log10(Down_up_Pvalue$p_value)
Down_up_Pvalue$tp <- ifelse(Down_up_Pvalue$p_value < -log10(0.05),"",
                            ifelse(Down_up_Pvalue$p_value >= -log10(0.05) & Down_up_Pvalue$p_value<=-log10(0.01),"*",
                                   ifelse(Down_up_Pvalue$p_value >= -log10(0.01) & Down_up_Pvalue$p_value<=-log10(0.001),"**","***")))
Down_up_Pvalue$p_value <- ifelse(Down_up_Pvalue$p_value > 20,20,Down_up_Pvalue$p_value)
Down_up_Pvalue$cancer <- factor(Down_up_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Down_up_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c(rep("#7E1717",24)))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())

#UTDA
Up_down_data <- old_DEGs_total[,c(1,4,9)]
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
Up_down_Pvalue <- data.frame(cancer=rownames(old_DEGs_total),p_value=old_DEGs_total$down_age_tumor_up_DEGs_P)
Up_down_Pvalue$p_value <- -1*log10(Up_down_Pvalue$p_value)
Up_down_Pvalue$tp <- ifelse(Up_down_Pvalue$p_value < -log10(0.05),"",
                            ifelse(Up_down_Pvalue$p_value >= -log10(0.05) & Up_down_Pvalue$p_value<=-log10(0.01),"*",
                                   ifelse(Up_down_Pvalue$p_value >= -log10(0.01) & Up_down_Pvalue$p_value<=-log10(0.001),"**","***")))
Up_down_Pvalue$p_value <- ifelse(Up_down_Pvalue$p_value > 20,20,Up_down_Pvalue$p_value)
Up_down_Pvalue$cancer <- factor(Up_down_Pvalue$cancer,levels = rev(paixu$tumor))
ggplot(Up_down_Pvalue,aes(p_value,cancer,fill=cancer))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values = c(rep("#068DA9",24)))+
  scale_x_continuous(limits = c(0,20))+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme_tufte()+
  theme(panel.background = element_blank())

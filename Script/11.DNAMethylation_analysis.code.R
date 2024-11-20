#################################################################
########Processing DNA methylation data from normal tissues######
#################################################################
setwd("./Aging/data/UCSC/DNA Methy")
library(data.table)
Methy <- as.data.frame(fread("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/GSE213478_methylation_DNAm_noob_final_BMIQ_all_tissues_987.txt",sep = ","))
rownames(Methy) <- Methy$V1
Methy <- Methy[,-1]
dim(Methy)
Methy[1:5,1:5]
Methy_clinic <- as.data.frame(read.table("GTEx_Methylation/GSE213478_series_matrix.txt"),sep="\t")
Methy_clinic <- Methy_clinic[c(1,2,3,4,5,6,9,10,11,12,13,14,15,17,19,20,21,33,46),]
Methy_clinic <- Methy_clinic[,-1]
row.names(Methy_clinic) <- c("Sample_title","Sample_geo_accession","Sample_status",
                             "Sample_submission_date","Sample_last_update_date","Sample_type",
                             "Sample_organism","Sample_plate","Sample_plate_position","sentrix_id",
                             "sentrix_position","participant_id","sample_id","Sex","tissue","age",
                             "smoker_status","Sample_platform_id","ID_REF")
Methy_clinic <- as.data.frame(t(Methy_clinic))
rownames(Methy_clinic) <- Methy_clinic$Sample_title
dim(Methy_clinic)
#改变组织
Methy_clinic$tissue <- unlist(lapply(Methy_clinic$tissue, function(x){
  y <- strsplit(x," -")[[1]][1]
  return(strsplit(y,": ")[[1]][2])
}))
table(Methy_clinic$tissue)
#改变年龄
Methy_clinic$age_group <- unlist(lapply(Methy_clinic$age,function(x){return(strsplit(x,":")[[1]][2])}))
Methy_clinic$age <- unlist(lapply(Methy_clinic$age_group, function(x){
  y <- as.numeric(strsplit(x,"-")[[1]][1])
  return(y + 5)}))
#改变性别
Methy_clinic$Sex <- unlist(lapply(Methy_clinic$Sex, function(x){return(strsplit(x,":")[[1]][2])}))
Methy_clinicData <- Methy_clinic[,c("Sample_title","Sex","tissue","age","age_group")]
write.table(Methy_clinicData,"GTEX_Methy_Clin_data.txt",sep = "\t",quote = F,row.names = F)
##  分组织类型存储
GTEx_cancers <- names(table(Methy_clinicData$tissue))
i = 1
GTEx_Methy_sample_num <- c()
for (i in 1:length(GTEx_cancers)) {
  cancer <- GTEx_cancers[i]
  gtex_samples <- Methy_clinicData[Methy_clinicData$tissue == cancer,"Sample_title"]
  GTEx_Methy_sample_num <- c(GTEx_Methy_sample_num,length(gtex_samples))
  sum(gtex_samples %in% colnames(Methy))
  Methy_cancer <- Methy[,gtex_samples]
  Methy_cancer <- cbind(rownames(Methy_cancer),Methy_cancer)
  colnames(Methy_cancer)[1] <- "cg"
  path <- paste(c("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/",cancer,"_Methy.txt"),collapse = "")
  write.table(Methy_cancer,path,sep = "\t",quote = F,row.names = F)
  print(i)
}

############################################################################################
########Obtaining DNA methylation data for the UTDA gene in epithelial cells in cancer######
############################################################################################
load(file="./Aging/data/UCSC/GTEX_RNAseq_data/细胞类型识别/组织统计/Epithelial_cell/Epithelial_cancer_intersect.rda")
ID.ver <- function(id){
  id2 <- unlist(strsplit(id,".",fixed = T))
  id3 <- id2[1]
  return(id3)
}
gene <- read.table("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/gencode_v42_annotation.txt",sep = "\t",header = T)
gene <- gene[,c(1,2,4,5,6,7)]
gene$id <- apply(as.data.frame(gene[,1]), 1, ID.ver)
gene <- gene[which(gene$chrom != "chrM"),]
x <- gene[3,]
class_site <- function(x){
  if(x[6] == "+"){
    site <- as.numeric(x[4]) 
  }else{
    site <- as.numeric(x[5])
  }
  return(site)
}
g_site <- apply(gene,1, class_site) 
g_site <- cbind(gene[,c("id","gene","chrom")],g_site)
colnames(g_site) <- c("id","gene","chrom","site")
g_site$start <- as.numeric(g_site$site) - 2000
g_site$end <- as.numeric(g_site$site) + 2000
length(table(g_site$chrom))
methy_site <- read.table("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/info_final.txt",sep = "\t",header = T)
cg_site <- methy_site[,c(1,2,3)]
dim(cg_site)
for (cancer in names(intersect_gene)) {
  print(cancer)
  one_gene <- intersect_gene[[cancer]]
  one_g_site <- g_site[match(one_gene,g_site$gene,nomatch = 0),]
  index_list <- c()
  res_pmatch <- data.frame()
  i = 1
  for (i in 1:dim(one_g_site)[1]) {
    site_data <- one_g_site[i,c("chrom","start","end")]
    id <- one_g_site[i,"gene"]
    cg_chr <- cg_site[c(cg_site$chrom_hg38 == site_data$chrom),]
    index <- which((as.numeric(cg_chr$pos_hg38) >= as.numeric(site_data$start)) & (as.numeric(site_data$end) >= as.numeric(cg_chr$pos_hg38)))
    if(length(index) != 0){
      cg_chr_pmatch <- cg_chr[index,]
      re <- cbind(c(rep(id,length(index))),cg_chr_pmatch)
      index_list <- c(index_list,cg_chr_pmatch$id)
      res_pmatch <- rbind(res_pmatch,re)
      print(i)
    }
  }
  colnames(res_pmatch) <- c("gene","cg_id","chrom","cg_pos")
  head(res_pmatch)
  write.table(res_pmatch,paste0("./Aging/data/UCSC/DNA Methy/TCGA_Methylation/UTDA_cg_site/",cancer,"res_pmatch.txt"),sep = "\t",quote = F,row.names = F)
  unique_cg <- unique(index_list)
  Methy <- as.data.frame(fread(paste0("./Aging/data/UCSC/DNA Methy/TCGA_Methylation/UCSC_DNAmethylation_data/TCGA-",cancer,".methylation450.tsv"),sep = "\t"))
  rownames(Methy) <- Methy$`Composite Element REF`
  Methy <- Methy[,-1]
  Methy_need <- Methy[unique_cg,] 
  cg_NA <- which(apply(Methy_need,1,function(x){return(sum(is.na(x) == TRUE))}) > ncol(Methy_need)*0.7)
  Methy_DelNA=Methy_need[-cg_NA,]
  res_pmatch_DelNA <-  res_pmatch[-match(names(cg_NA),res_pmatch$cg_id),]
  Methy_DelNA1 <- Methy_DelNA
  for (i in 1:length(rownames(Methy_DelNA1))) {
    if(sum(is.na(Methy_DelNA1[i,]) == TRUE) != 0){
      pd_NA <- is.na(Methy_DelNA1[i,])
      Methy_DelNA1[i,pd_NA] <- rowSums(Methy_DelNA1[i,!pd_NA])/length(Methy_DelNA1[i,!pd_NA])
    }
  }
  print(data.frame(as.numeric(which(is.na(rowSums(Methy_DelNA1))))))
  cacner_Methy_beta <- c()
  onegene <- unique(res_pmatch_DelNA$gene)[1]
  for (onegene in unique(res_pmatch_DelNA$gene)) {
    cg_str <- res_pmatch_DelNA[which(res_pmatch_DelNA$gene == onegene),"cg_id"]
    methy_data <- Methy_DelNA1[cg_str,]
    methy_res <- colSums(methy_data)/length(cg_str)
    cacner_Methy_beta <- rbind(cacner_Methy_beta,methy_res)
    print(onegene)
  }
  rownames(cacner_Methy_beta) <- unique(res_pmatch_DelNA$gene)
  save_beta <- cbind(rownames(cacner_Methy_beta),cacner_Methy_beta)
  colnames(save_beta)[1] <- "gene_name"
  write.table(save_beta,paste0("./Aging/data/UCSC/DNA Methy/TCGA_Methylation/UTDA_methy/",cancer,"methybeta.txt"),sep = "\t",quote = F,row.names = F,col.names = T)
}




####################################################################################################
########Obtaining DNA methylation data for the UTDA gene in normal tissue epithelial cells######
####################################################################################################
library(broom)
library(biomaRt)
load(file="./Aging/data/UCSC/GTEX_RNAseq_data/细胞类型识别/组织统计/Epithelial_cell/Epithelial_cancer_intersect.rda")
Epithelial_cancer <- intersect_gene[c("BRCA","LUAD","LUSC","COAD","READ","OV","TGCT")]
Normal_tissue_Epithelial <- data.frame(cancaer=c("BRCA","LUAD","LUSC","COAD","READ","OV","TGCT"),
                                       tissue=c("Breast","Lung","Lung","Colon","Colon","Ovary","Testis"))
gene <- read.table("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/gencode_v42_annotation.txt",sep = "\t",header = T)
#get epithelial cell UTDA genes methylation matrices
cancer <- Normal_tissue_Epithelial$cancaer[2]
for (cancer in Normal_tissue_Epithelial$cancaer) {
  tissue <- Normal_tissue_Epithelial[Normal_tissue_Epithelial$cancaer == cancer,2]
  Methy_df <- as.data.frame(fread(paste0("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/",tissue,"_Methy.txt"),sep = "\t"))
  res_pmatch <- read.table(paste0("./Aging/data/UCSC/DNA Methy/TCGA_Methylation/UTDA_cg_site/",cancer,"res_pmatch.txt"),sep = "\t",header = T)
  rownames(Methy_df) <- Methy_df$cg
  Methy_df <- Methy_df[,-1]
  Methy_df_gn <- Methy_df[match(res_pmatch$cg_id,rownames(Methy_df),nomatch = 0),]
  gene_cg <- res_pmatch[match(rownames(Methy_df),res_pmatch$cg_id,nomatch = 0),]
  tissue_Methy_beta <- c()
  onegene <- unique(gene_cg$gene)[1]
  for (onegene in unique(gene_cg$gene)) {
    cg_str <- gene_cg[which(gene_cg$gene == onegene),"cg_id"]
    methy_data <- Methy_df_gn[cg_str,]
    methy_res <- colSums(methy_data)/length(cg_str)
    tissue_Methy_beta <- rbind(tissue_Methy_beta,methy_res)
    print(onegene)
  }
  rownames(tissue_Methy_beta) <- unique(gene_cg$gene)
  Epithelial_down_gene_methy_df <- cbind(rownames(tissue_Methy_beta),tissue_Methy_beta)
  colnames(Epithelial_down_gene_methy_df)[1] <- "gene_name"
  write.table(Epithelial_down_gene_methy_df,paste0("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/UTDA_epithelial_methy/Epithelial_",cancer,"_MethyBeta.txt"),sep = "\t",quote = F,row.names = F)
}

###################################################################
#####Identifying genes that increase with age in normal tissues####
###################################################################
methy_clinic <- read.table("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/GTEX_Methy_Clin_data.txt",sep = "\t",header = T)
for (cancer in Normal_tissue_Epithelial$cancaer) {
  Epithelial_Methy_df <- as.data.frame(fread(paste0("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/UTDA_epithelial_methy/Epithelial_",cancer,"_MethyBeta.txt"),sep = "\t")) 
  rownames(Epithelial_Methy_df) <- Epithelial_Methy_df$V1
  Epithelial_Methy_df <- Epithelial_Methy_df[,-1]
  clinic_methy_epi <- methy_clinic[match(colnames(Epithelial_Methy_df),methy_clinic$Sample_title,nomatch = 0),]
    results <- data.frame(gene=NA,estimate=NA,std.error=NA,statistic=NA,p.value=NA)
  for (gene in rownames(Epithelial_Methy_df)) {
    Epithelial_Methy_tmp <- Epithelial_Methy_df[rownames(Epithelial_Methy_df) == gene,]
    Epithelial_Methy_tmp <- as.data.frame(t(Epithelial_Methy_tmp))
    colnames(Epithelial_Methy_tmp) <- "gene"
    df_tmp <- merge(Epithelial_Methy_tmp, clinic_methy_epi[,c(1,2,4)], by.x = "row.names", by.y = "Sample_title")
    colnames(df_tmp) <- c("patient", colnames(df_tmp)[2:ncol(df_tmp)])
    lm_fit <- lm(formula = as.formula(gene ~ age + Sex), data=df_tmp)  # fit linear model
    summary(lm_fit)
    result <- tidy(lm_fit)
    result_df <- as.data.frame(result[result$term == "age",])
    result_df$term <- gene
    colnames(result_df) <- c("gene", "estimate", "std.error", "statistic", "p.value")
    results <- rbind(results,result_df)
  }
  results <- results[-1,]
  results$q.value <- p.adjust(results$p.value, method = "BH")
  results$Sig <- ifelse(results$q.value < 0.05, TRUE, FALSE)
  write.table(results,paste0("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/UTDA_epithelial_methy/",cancer,"_results.txt"),sep = "\t",quote = F,col.names = T,row.names = F)
}

#statistical results
Age_Methy_DEGs <- data.frame(gene_nums=NA,Up_nums=NA,Down_nums=NA,Nosig=NA)
for (cancer in Normal_tissue_Epithelial$cancaer) {
  df <- read.table(paste0("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/UTDA_epithelial_methy/",cancer,"_results.txt"),sep = "\t",header = T)
  df$type <- NA
  for(i in 1:nrow(df)){
    if(df$q.value[i] < 0.05 & df$estimate[i] < 0){
      df$type[i] <- "down"
    } else if(df$q.value[i] < 0.05 & df$estimate[i] > 0){
      df$type[i] <- "up"
    } else {
      df$type[i] <- "Nosig"
    }
  }
  
  df_total <- data.frame(gene_nums=length(rownames(df)),
                         Up_nums=length(rownames(df)[which(df$type == "up")]),
                         Down_nums=length(rownames(df)[which(df$type == "down")]),
                         Nosig=length(rownames(df)[which(df$type == "Nosig")]))
  rownames(df_total) <- cancer
  Age_Methy_DEGs <- rbind(Age_Methy_DEGs,df_total)
}
###Fig.5A-----
Age_Methy_DEGs <- Age_Methy_DEGs[-1,]
write.table(Age_Methy_DEGs,paste0("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/UTDA_epithelial_methy/Age_Methy_DEGs.csv"),sep = ",",quote = F,row.names = T,col.names = T)
Age_Methy_DEGs <- read.table(paste0("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/UTDA_epithelial_methy/Age_Methy_DEGs.csv"),sep = ",",header = T)
data <- reshape2::melt(Age_Methy_DEGs[,-2])
data$cancer <- factor(data$cancer,levels = rev(Age_Methy_DEGs$cancer))
data$variable <- factor(data$variable,levels = c("Up_nums","Down_nums","Nosig"))
paired=c("#DD5746","#4793AF","#FFEFEF")
ggplot(data,aes(cancer,value,fill=variable))+
  labs(x='tissue',cex=10,y='case')+theme_test(base_size = 20)+
  scale_fill_manual(values = paired)+
  geom_bar(stat="identity",position="fill",width = 0.8,size = 0.25)+
  theme(axis.text.x=element_text(angle=90, vjust=0.5))



###############################################################################################
########Tumor-differential methylation analysis of epithelial cells UTDA genes in cancers######
###############################################################################################
load(file="./Aging/data/UCSC/GTEX_RNAseq_data/细胞类型识别/组织统计/Epithelial_cell/Epithelial_cancer_intersect.rda")
Epithelial_sig <- c("UCEC","BRCA","LUAD","LUSC","COAD","READ","PAAD","THCA")
Epithelial_cancer <- intersect_gene[Epithelial_sig]
for (cancer in Epithelial_sig) {
  TCGA_methy <- as.data.frame(fread(paste0("./Aging/data/UCSC/DNA Methy/TCGA_Methylation/UTDA_methy/",cancer,"methybeta.txt"),sep = "\t"))
  rownames(TCGA_methy) <- TCGA_methy$gene_name
  TCGA_methy <- TCGA_methy[,-1]
  type <- ifelse(substr(colnames(TCGA_methy),14,14) == "0","Tumor","Normal")
  sample_type <- data.frame(sample_id=colnames(TCGA_methy),
                            type= type)
  sample_type$type <- factor(sample_type$type,levels = c("Tumor","Normal"))
  Methy_data=as.data.frame(t(TCGA_methy))
  sample_id <- rownames(Methy_data)
  Methy_data <- as.data.frame(lapply(Methy_data,as.numeric))
  rownames(Methy_data) <- sample_id
  identical(rownames(Methy_data),sample_type$sample_id)
  
  Methy_data <- cbind(Methy_data,sample_type$type)
  colnames(Methy_data)<- gsub("\\.", "-",colnames(Methy_data))
  genes <- rownames(TCGA_methy)
  total<-data.frame(gene_name=genes,p.value=NA,detalR=NA)
  #Rank sum test to obtain p-value
  for (gene in genes) {
    print(gene)
    res<- wilcox.test(unlist(Methy_data[,gene])~type,data = Methy_data)###rank-sum test
    total[which(total$gene_name == gene),2]<- res$p.value
  }
  total$q.value <- p.adjust(total$p.value, method = "BH")
  cancer_methy <- TCGA_methy[,sample_type$sample_id[which(sample_type$type == "Tumor")]]
  normal_methy <- TCGA_methy[,sample_type$sample_id[which(sample_type$type == "Normal")]]
  #deltaR
  for (gene in genes) {
    print(gene)
    Onecancer_Methy <- cancer_methy[gene,]
    Onenormal_Methy <- normal_methy[gene,]
    resDetalR <- rowMeans(Onecancer_Methy)-rowMeans(Onenormal_Methy)
    total[which(total$gene_name == gene),3] <- resDetalR
  }
  total$Sig <- ifelse(total$q.value < 0.05, TRUE, FALSE)
  total$type <- NA
  for(i in 1:nrow(total)){
    if(total$q.value[i] < 0.05 & total$detalR[i] < 0){
      total$type[i] <- "down"
    } else if(total$q.value[i] < 0.05 & total$detalR[i] > 0){
      total$type[i] <- "up"
    } else {
      total$type[i] <- "Nosig"
    }
  }
  total$type <- "Nosige"
  total$p.value <- 0
  total$detalR <- 0
  write.table(total,paste0("./Aging/data/UCSC/DNA Methy/TCGA_Methylation/UTDA_methy/Dif_analysis/",cancer,"_diffRes.txt"),sep = "\t",quote = F,row.names = F,col.names = T)
}

#statistical results
All_diff_Methy <- data.frame(gene_nums=NA,Up_nums=NA,Down_nums=NA,Nosig=NA)
cancer <- tumor_normal[1]
for(cancer in tumor_normal){
  one_cancer_total <-  read.table(paste0("./Aging/data/UCSC/DNA Methy/TCGA_Methylation/UTDA_methy/Dif_analysis/",cancer,"_diffRes.txt"),sep = "\t",header = T)
  diff_methy_total <- data.frame(gene_nums=length(rownames(one_cancer_total)),
                                 Up_nums=length(rownames(one_cancer_total)[which(one_cancer_total$type == "up")]),
                                 Down_nums=length(rownames(one_cancer_total)[which(one_cancer_total$type == "down")]),
                                 Nosig=length(rownames(one_cancer_total)[which(one_cancer_total$type == "Nosig")]))
  rownames(diff_methy_total) <- cancer
  All_diff_Methy <- rbind(All_diff_Methy,diff_methy_total)
}

All_diff_Methy <- All_diff_Methy[-1,]

####Fig.5B------
data <- reshape2::melt(All_diff_Methy[,-2])
data$cancer <- factor(data$cancer,levels = rev(All_diff_Methy$cancer))
data$variable <- factor(data$variable,levels = c("Down_nums","Up_nums","Nosig"))
paired=c("#4793AF","#DD5746","#FFEFEF")
ggplot(data,aes(cancer,value,fill=variable))+
  labs(x='tissue',cex=10,y='case')+theme_test(base_size = 20)+
  scale_fill_manual(values = paired)+
  geom_bar(stat="identity",position="fill",color ="black",width = 0.8,size = 0.25)+
  theme(axis.text.x=element_text(angle=90, vjust=0.5))


##########################################
#####Age DMGs and cancer DMGs overlap####
#########################################
library(GeneOverlap)
load("./Aging/data/UCSC/GTEX_RNAseq_data/aging_gene_list_new.Rdata")
all_inter_list <- list()
for (cancer in c("LUAD","LUSC","COAD","READ")) {
  #Normal
  df <- read.table(paste0("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/UTDA_epithelial_methy/",cancer,"_results.txt"),sep = "\t",header = T)
  df$type <- NA
  for(i in 1:nrow(df)){
    if(df$q.value[i] < 0.05 & df$estimate[i] < 0){
      df$type[i] <- "down"
    } else if(df$q.value[i] < 0.05 & df$estimate[i] > 0){
      df$type[i] <- "up"
    } else {
      df$type[i] <- "Nosig"
    }
  }
  #cancer
  one_cancer_total <-  read.table(paste0("./Aging/data/UCSC/DNA Methy/TCGA_Methylation/UTDA_methy/Dif_analysis/",cancer,"_diffRes.txt"),sep = "\t",header = T)
  UpInAge_methy <- df$gene[which(df$type=="up")]
  DownIncancer_methy <- one_cancer_total$gene_name[which(one_cancer_total$type=="down")]
  genes <- intersect(UpInAge_methy,DownIncancer_methy)
  inter_list <- list(genes)
  names(inter_list) <- cancer
  all_inter_list <- c(all_inter_list,inter_list)
}











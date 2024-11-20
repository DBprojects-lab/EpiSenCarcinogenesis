library(string)
library(DESeq2)
library(limma)
library(data.table)
Tissues <- c("Adrenal_Gland","Breast","Esophagus","Brain","Kidney","Prostate","Liver","Lung","Ovary","Pancreas","Colon",
             "Testis","Uterus","Skin","Stomach","Thyroid","Blood")
Tissues_proflie <- c("Adrenal_Gland_ACC","Breast_BRCA","Esophagus_ESCA","Brain_LGG_GBM","Kidney_KICH_KIRC_KIRP","Prostate_PRAD",
                     "Liver_LIHC","Lung_LUAD_LUSC","Ovary_OV","Pancreas_PAAD","Colon_COAD_READ","Testis_TGCT","Uterus_UCEC_UCS","Skin_SKCM",
                     "Stomach_STAD","Thyroid_THCA","Blood_DLBC_THYM")
#########################################################################################################
######Obtaining normalized expression matrices and sample information for different tissues in GTEX######
#########################################################################################################
for (i in 1:length(Tissues_proflie)) {
  setwd(paste0("./Aging/data/UCSC/GTEX_RNAseq_data/",Tissues_proflie[i]))
  read_count <- as.data.frame(fread(paste0(Tissue[i],"_RNAseq_count.txt"),header = T,sep = "\t"))
  rownames(read_count) <- read_count[,1]
  read_count <- read_count[,-1]
  clinical <- as.data.frame(fread(paste0(Tissue[i],"_sample_info.txt"),header = T,sep = "\t"))
  rownames(clinical) <- clinical$SAMPID
  clinical <- clinical[,c("SEX","Age","DTHHRDY","SMCENTER","SMRIN","SMTSISCH","SMEXNCRT","SMRRNART","SMNTERRT")]
  #Removing samples with NA for covariates
  clinical <- clinical[complete.cases(clinical[,c("DTHHRDY","SMCENTER","SMRIN","SMTSISCH","SMEXNCRT","SMRRNART","SMNTERRT")]),]
  read_count <- read_count[,rownames(clinical)]
  read_count <- 2^read_count-1
  read_count <- round(read_count)
  ID.ver <- function(id){
    id2 <- unlist(strsplit(id,".",fixed = T))
    id3 <- id2[1]
    return(id3)
  }
  gene_id <- apply(as.data.frame(rownames(read_count)), 1, ID.ver)
  rownames(read_count) <- gene_id
  clinical$Gender=ifelse(clinical$SEX==1,"Male","Female")
  clinical$SEX=NULL
  clinical$DTHHRDY = as.character(clinical$DTHHRDY)
  clinical$Age <- str_replace(clinical$Age,"-","_")
  clinical$Age = as.factor(clinical$Age)#Age becomes a factor
  clinical$DTHHRDY <- as.factor(clinical$DTHHRDY)#hardy becomes a factor
  clinical$Gender <- as.factor(clinical$Gender)#Sex becomes a factor
  
  #Genes with raw expression values higher than 1 in more than 20% of the samples were retained
  dds=read_count[which(apply(read_count,1,function(x){return(sum(x>1))})>ncol(read_count)*0.20),]
  
  #DESeq2 standardization - methodology VST 
  vsd = varianceStabilizingTransformation(as.matrix(dds),blind=FALSE)
  
  #Handling covariates with limma's removeBatchEffect
  batch.matrix <- model.matrix(~Gender+DTHHRDY+SMCENTER+SMRIN+SMTSISCH+SMEXNCRT+SMRRNART+SMNTERRT,data = clinical)
  vsd <-  limma::removeBatchEffect(vsd,
                                   design = model.matrix(~Age,data=clinical),
                                   covariates = batch.matrix[,-1]
  )
  vsd.adjust=round(as.data.frame(vsd),digits=3)
  vsd.adjust$Ensemble <- rownames(vsd.adjust)
  genetype <- as.data.frame(fread("Aging/data/UCSC/probeMap_gencode.v23.annotation.gene.probemap",header = T,sep = "\t"))
  colnames(genetype)[1:2] <- c("Ensemble","Symbol")
  genetype$Ensemble <- apply(as.data.frame(genetype[,1]), 1, ID.ver)
  vsd.adjust_symbol <- merge(genetype,vsd.adjust,by="Ensemble")
  vsd.adjust_symbol=vsd.adjust_symbol[!duplicated(vsd.adjust_symbol$Symbol),]
  rownames(vsd.adjust_symbol)=vsd.adjust_symbol$Symbol
  vsd.adjust_symbol=vsd.adjust_symbol[,c(7:ncol(vsd.adjust_symbol))]
  write.table(vsd.adjust_symbol,paste0(Tissue[i],".vsd.adjusted.txt"),sep="\t",row.names=TRUE,quote=F) # vst transformed and adjusted expression matrix
  write.table(clinical,paste0(Tissue[i],"_sampleInfomation.txt"),sep="\t",row.names=TRUE,quote=F) # vst transformed and adjusted expression matrix
}

#########################################################################################################
######Statistics on the age distribution of different organizational samples######
#########################################################################################################


tissue_sample_tolal <- as.data.frame(matrix(1,nrow = 6,ncol=1),row.names = c("20_29",
                                                                             "30_39",
                                                                             "40_49",
                                                                             "50_59",
                                                                             "60_69",
                                                                             "70_79"))
for (x in 1:17) {
  tissue_sample_tolal[,x] <- c(0,0,0,0,0,0)
  colnames(tissue_sample_tolal)[x] <- Tissues[x]
  setwd(paste0("D:\\procedure\\Aging\\data\\UCSC\\GTEX_RNAseq_data\\",Tissues_proflie[x]))
  read_count=read.table(paste0(Tissues[x],".vsd.adjusted.txt"),header=T,row.names=1,check.names=F)
  clinical=read.table(paste0(Tissues[x],"_sampleInfomation.txt"),header=T,row.names=1,sep="\t")
  a <- stringr::str_replace(as.character(as.data.frame(table(clinical$Age))$Var1),"-","_")
  tissue_sample_tolal[match(a,rownames(tissue_sample_tolal),nomatch = 0),x]<- as.data.frame(table(clinical$Age))$Freq
}
data <- t(tissue_sample_tolal)
a <- apply(tissue_sample_tolal,2,sum)
tissue_sample_tolal <- rbind(tissue_sample_tolal,a)
need <- as.data.frame(t(tissue_sample_tolal[-7,]))
need$tissue <- rownames(need)
need1 <- reshape2::melt(need,id.vars = 'tissue')
paired = c("#8E0152", "#C51B7D" ,"#DE77AE", "#F1B6DA" ,"#FFBFBF","#FDE0EF")
need1$variable <- factor(need1$variable,levels = rev(c("20_29","30_39","40_49","50_59","60_69","70_79")))
need1$tissue <- factor(need1$tissue,levels = rev(c("Brain","Thyroid","Esophagus","Lung","Breast","Stomach","Liver","Pancreas","Adrenal_Gland",
                                                   "Kidney","Colon","Ovary","Uterus","Prostate","Testis","Blood","Skin")))
need2 = ddply(need1,'tissue',transform,percent_con=value/sum(value)*100)
#----Figure 1A middle
ggplot(need2,aes(tissue,percent_con,fill=variable))+
  labs(x='tissue',cex=10,y='case')+theme_test(base_size = 20)+
  scale_fill_manual(values = paired)+
  geom_bar(stat="identity",position="stack",color ="black",width = 0.8,size = 0.25)+
  theme(axis.text.x=element_text(angle=90, vjust=0.5))

df <- data.frame(tissue = Tissues,
                 background_gene=0)
for (x in 1:17) {
  setwd(paste0("D:\\procedure\\Aging\\data\\UCSC\\GTEX_RNAseq_data\\",Tissues_proflie[x]))
  read_count=read.table(paste0(Tissues[x],".vsd.adjusted.txt"),header=T,row.names=1,check.names=F)
  clinical=read.table(paste0(Tissues[x],"_sampleInfomation.txt"),header=T,row.names=1,sep="\t")
  df$background_gene[x] <- length(rownames(read_count))
}
df$tissue <- factor(df$tissue,levels = c("Brain","Thyroid","Esophagus","Lung","Breast","Stomach","Liver","Pancreas","Adrenal_Gland",
                                         "Kidney","Colon","Ovary","Uterus","Prostate","Testis","Blood","Skin"))
df$loggene <- log10(df$backgroud_gene)
df$ystart <- 20000
#----Figure 1A right
ggplot(df,aes(x=tissue,y=backgroud_gene,group=2))+geom_bar(stat = "identity",width = 0.8)+geom_line() + geom_point()+coord_cartesian(ylim = c(20000,38000))


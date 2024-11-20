###########################################
########Copy Number Variation analysis#####
###########################################
library(maftools)
library(ComplexHeatmap)
setwd("./Aging/data/UCSC/UCSC_copynumber_data")
load(file="./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/组织统计/Epithelial_cell/Epithelial_cancer_intersect.rda")
Epithelial_cancer <- intersect_gene[c("BRCA","LUAD","LUSC","COAD","READ","OV","TGCT")]
Normal_tissue_Epithelial <- data.frame(cancaer=c("BRCA","LUAD","LUSC","COAD","READ","OV","TGCT"),
                                       tissue=c("Breast","Lung","Lung","Colon","Colon","Ovary","Testis"))
ID.ver <- function(id){
  id2 <- unlist(strsplit(id,".",fixed = T))
  id3 <- id2[1]
  return(id3)
}
genetype <- as.data.frame(fread("./Aging/data/UCSC/probeMap_gencode.v23.annotation.gene.probemap",header = T,sep = "\t"))
colnames(genetype)[1:2] <- c("Ensemble","Symbol")
genetype$Ensemble <- apply(as.data.frame(genetype[,1]), 1, ID.ver)
#Copy Number Variation Percentage Function
copyNum_per <- function(x){
  copy_nums <- as.numeric(length(x[which(x != 0)]))
  persent <- copy_nums / as.numeric(length(x))
  return(persent)
}
#Copy Number Extension Percentage Function
Amp_freq <- function(x){
  copy_nums <- as.numeric(length(x[which(x == 1)]))
  persent <- copy_nums / as.numeric(length(x))
  return(persent)
}

#Copy Number Missing Percentage Function
Del_freq <- function(x){
  copy_nums <- as.numeric(length(x[which(x == -1)]))
  persent <- copy_nums / as.numeric(length(x))
  return(persent)
}

for (cancer in names(intersect_gene)) {
  cancer_gene <- intersect_gene[[cancer]]
  cancer_gene_id <- genetype$Ensemble[match(cancer_gene,genetype$Symbol,nomatch = 0)]
  df <- read.table(paste0(cancer,".txt"),sep = "\t",header = T)
  colnames(df)[1] <- "Ensemble"
  df_name <- merge(genetype,df,by="Ensemble")
  need_copynumb_df <- df_name[match(cancer_gene_id,df_name$Ensemble,nomatch = 0),]
  if(length(rownames(need_copynumb_df)) != 0){
    rownames(need_copynumb_df) <- need_copynumb_df$Symbol
    need_copynumb_df <- need_copynumb_df[,7:ncol(need_copynumb_df)]
    write.table(need_copynumb_df,paste0("./Aging/data/UCSC/UCSC_copynumber_data/上皮细胞UTDT_copy/",cancer,"_copyNumber.txt"),sep="\t",row.names=TRUE,quote=F)
    copyNum_per_total <- data.frame(gene=rownames(need_copynumb_df),
                                    persent_num= apply(need_copynumb_df, 1, copyNum_per))
    copyNum_per_total$Amp_freq <- apply(need_copynumb_df, 1, Amp_freq)
    copyNum_per_total$Del_freq <- apply(need_copynumb_df, 1, Del_freq)
    copyNum_per_total$CNV_Freq <- copyNum_per_total$Amp_freq-copyNum_per_total$Del_freq
    copyNum_per_total$Type <- ifelse(copyNum_per_total$CNV_Freq > 0,"Amp",ifelse(copyNum_per_total$CNV_Freq < 0,"Del","NULL"))
    write.table(copyNum_per_total,paste0("./Aging/data/UCSC/UCSC_copynumber_data/上皮细胞UTDT_copy/",cancer,"_persent.txt"),sep="\t",row.names=F,quote=F)
  }
}

##############################
########Mutation analysis#####
##############################
library(maftools)
setwd("./Aging/data/UCSC/TCGA的突变数据")
load(file="./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/组织统计/Epithelial_cell/Epithelial_cancer_intersect.rda")
Epithelial_cancer <- intersect_gene[c("UCEC","UCS","BRCA","LUAD","LUSC","COAD","READ","OV","PAAD","TGCT","THCA")]
for (cancer in names(Epithelial_cancer)) {
  df <- as.data.frame(fread(paste0(cancer,".txt"),sep = "\t",header = T,quote = ""))
  Epithelial_gene <- Epithelial_cancer[[cancer]]
  all_mut <- read.maf(df)
  a <- all_mut@data %>%
    .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")] %>%
    as.data.frame() %>%
    mutate(Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode,1,16))
  gene <- as.character(unique(a$Hugo_Symbol))
  sample <- as.character(unique(a$Tumor_Sample_Barcode))
  mat <- as.data.frame(matrix("",length(gene),length(sample),
                              dimnames = list(gene,sample)))
  mat_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                  dimnames = list(gene,sample)))
  for (i in 1:nrow(a)){
    mat[as.character(a[i,1]),as.character(a[i,3])] <- as.character(a[i,2])
  }
  for (i in 1:nrow(a)){
    mat_0_1[as.character(a[i,1]),as.character(a[i,3])] <- 1
  }
  gene_count <- data.frame(gene=rownames(mat_0_1),
                           count=as.numeric(apply(mat_0_1,1,sum))) %>%
    arrange(desc(count))
  Epithelial_gene_count <- gene_count[which(gene_count$gene%in%Epithelial_gene),]
  Epithelial_gene_count$Mutation_persent <- Epithelial_gene_count$count/length(sample)
  write.table(Epithelial_gene_count,paste0("./Aging/data/UCSC/TCGA的突变数据/Epithelial_UTDA/",cancer,"_mutation_persent.txt"),sep = "\t",row.names=F,quote=F)
}
##########################
########merge results#####
##########################
library(maftools)
load(file="./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/组织统计/Epithelial_cell/Epithelial_cancer_intersect.rda")
Epithelial_cancer <- intersect_gene[c("UCEC","UCS","BRCA","LUAD","LUSC","COAD","READ","OV","PAAD","TGCT","THCA")]
cancer <- names(Epithelial_cancer)[2]
Epithelial_UTDA_genomic <- list()
for (cancer in names(Epithelial_cancer)) {
  #methylation
  methylation_df <-  read.table(paste0("./Aging/data/UCSC/DNA Methy/TCGA_Methylation/UTDA_methy/Dif_analysis/",cancer,"_diffRes.txt"),sep = "\t",header = T)
  methylation_df$detalR_Sig <- ifelse(methylation_df$type != "Nosig",methylation_df$detalR,0)
  #CNV
  copynumb_df <- read.table(paste0("./Aging/data/UCSC/UCSC_copynumber_data/上皮细胞UTDT_copy/",cancer,"_persent.txt"),header=T,row.names=1,check.names=F)
  copynumb_df$gene <- rownames(copynumb_df)
  #mutation
  mutation_df <- read.table(paste0("./Aging/data/UCSC/TCGA的突变数据/Epithelial_UTDA/",cancer,"_mutation_persent.txt"),sep = "\t",header = T)
  #exp diff
  Expression_df <- read.table(paste0("./Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/表达分析结果total/",cancer,"_DEGs_all_info.txt"),sep = "\t",header = T)
  Expression_df$gene <- rownames(Expression_df)
  Expression_df$type <- ifelse(Expression_df$log2FoldChange >= 1.5 & Expression_df$padj < 0.01,"Up",ifelse(Expression_df$log2FoldChange <= -1.5 & Expression_df$padj < 0.01,"Down","Nosig"))
  Expression_df$log2FCNew <- ifelse(Expression_df$type != "Nosig",Expression_df$log2FoldChange,0)
  Intergene_name <- Reduce(union, list(methylation_df$gene_name,copynumb_df$gene, mutation_df$gene))
  Epithelial_UTDA_genomic[[cancer]] <- Intergene_name
  gene_genenomic <- data.frame(gene = Intergene_name,
                               methylation=NA,
                               meth_type=NA,
                               copynumb=NA,
                               CNV_type=NA,
                               mutation=NA,
                               Muta_type=NA,
                               expression=NA,
                               expre_type=NA)
  #methylation diff detalR
  for (i in 1:nrow(gene_genenomic)) {
    if(gene_genenomic[i,1] %in% methylation_df$gene_name){
      gene_genenomic[i,2] <- methylation_df$detalR_Sig[which(methylation_df$gene_name == gene_genenomic[i,1])]
      gene_genenomic[i,3] <- methylation_df$type[which(methylation_df$gene_name == gene_genenomic[i,1])]
    }else{
      gene_genenomic[i,2] <- 0
      gene_genenomic[i,3] <- "Nosig"
    }
  }
  #Proportion of CNV variants
  for (i in 1:nrow(gene_genenomic)) {
    if(gene_genenomic[i,1] %in% copynumb_df$gene){
      gene_genenomic[i,4] <- copynumb_df$CNV_Freq[which(copynumb_df$gene == gene_genenomic[i,1])]
      gene_genenomic[i,5] <- copynumb_df$Type[which(copynumb_df$gene == gene_genenomic[i,1])]
    }else{
      gene_genenomic[i,4] <- 0
      gene_genenomic[i,5] <- "NULL"
    }
  }
  #mutation rate
  for (i in 1:nrow(gene_genenomic)) {
    if(gene_genenomic[i,1] %in% mutation_df$gene){
      gene_genenomic[i,6] <- mutation_df$Mutation_persent[which(mutation_df$gene == gene_genenomic[i,1])]
    }else{
      gene_genenomic[i,6] <- 0
    }
    gene_genenomic[i,7] <- "NULL"
  }
  #Fold FC in expression
  for (i in 1:nrow(gene_genenomic)) {
    if(gene_genenomic[i,1] %in% Expression_df$gene){
      gene_genenomic[i,8] <- Expression_df$log2FCNew[which(Expression_df$gene == gene_genenomic[i,1])]
      gene_genenomic[i,9] <- Expression_df$type[which(Expression_df$gene == gene_genenomic[i,1])]
    }else{
      gene_genenomic[i,8] <- 0
      gene_genenomic[i,9] <- "Nosig"
    }
  }
  write.table(gene_genenomic,paste0("./Aging/data/UCSC/genomic/",cancer,"_genomic.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
}

save(Epithelial_UTDA_genomic,file="./Aging/data/UCSC/genomic/Epithelial_UTDA_genomic_gene.rda")

a <- as.data.frame(table(unlist(Epithelial_UTDA_genomic)))
a <- a[order(a$Freq,decreasing = T),]
gene <- as.character(a$Var1[which(a$Freq>=9)])

conbind_genomic <- data.frame()

for (cancer in names(Epithelial_cancer)) {
  One_cancerDF<- read.table(paste0("./Aging/data/UCSC/genomic/",cancer,"_genomic.txt"),sep = "\t",header = T)
  #表达差异
  Expression_df <- read.table(paste0("./Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/表达分析结果total/",cancer,"_DEGs_all_info.txt"),sep = "\t",header = T)
  Expression_df$gene <- rownames(Expression_df)
  Expression_df$type <- ifelse(Expression_df$log2FoldChange >= 1.5 & Expression_df$padj < 0.01,"Up",ifelse(Expression_df$log2FoldChange <= -1.5 & Expression_df$padj < 0.01,"Down","Nosig"))
  Expression_df$log2FCNew <- ifelse(Expression_df$type != "Nosig",Expression_df$log2FoldChange,0)
  Topgene_genenomic <- data.frame(cancer=cancer,
                                  gene = gene,
                                  methylation=NA,
                                  meth_type=NA,
                                  copynumb=NA,
                                  CNV_type=NA,
                                  mutation=NA,
                                  Muta_type=NA,
                                  expression=NA,
                                  expre_type=NA)
  for (i in 1:nrow(Topgene_genenomic)) {
    if(Topgene_genenomic[i,2] %in% One_cancerDF$gene){
      Topgene_genenomic[i,3] <- One_cancerDF$methylation[which(One_cancerDF$gene == Topgene_genenomic[i,2])]
      Topgene_genenomic[i,4] <- One_cancerDF$meth_type[which(One_cancerDF$gene == Topgene_genenomic[i,2])]
      Topgene_genenomic[i,5] <- One_cancerDF$copynumb[which(One_cancerDF$gene == Topgene_genenomic[i,2])]
      Topgene_genenomic[i,6] <- One_cancerDF$CNV_type[which(One_cancerDF$gene == Topgene_genenomic[i,2])]
      Topgene_genenomic[i,7] <- One_cancerDF$mutation[which(One_cancerDF$gene == Topgene_genenomic[i,2])]
      Topgene_genenomic[i,8] <- One_cancerDF$Muta_type[which(One_cancerDF$gene == Topgene_genenomic[i,2])]
      Topgene_genenomic[i,9] <- One_cancerDF$expression[which(One_cancerDF$gene == Topgene_genenomic[i,2])]
      Topgene_genenomic[i,10] <- One_cancerDF$expre_type[which(One_cancerDF$gene == Topgene_genenomic[i,2])]
    }else{
      Topgene_genenomic[i,3] <- 0
      Topgene_genenomic[i,4] <- "Nosig"
      Topgene_genenomic[i,5] <- 0
      Topgene_genenomic[i,6] <- "NULL"
      Topgene_genenomic[i,7] <- 0
      Topgene_genenomic[i,8] <- "NULL"
      Topgene_genenomic[i,9] <- ifelse(Topgene_genenomic[i,2] %in% Expression_df$gene ,Expression_df$log2FCNew[which(Expression_df$gene == Topgene_genenomic[i,2])],0)
      Topgene_genenomic[i,10] <- ifelse(Topgene_genenomic[i,2] %in% Expression_df$gene ,Expression_df$type[which(Expression_df$gene == Topgene_genenomic[i,2])],"Nosig")
    }
  }
  write.table(Topgene_genenomic,paste0("./Aging/data/UCSC/genomic/共有top30/",cancer,".txt"),sep = "\t",row.names = F,col.names = T,quote = F)
  conbind_genomic <- rbind(conbind_genomic,Topgene_genenomic)
}

write.table(conbind_genomic,paste0("./Aging/data/UCSC/genomic/conbind_genomic.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

conbind_genomic$meth_type[which(conbind_genomic$meth_type == "Nosige")] <- "Nosig"

#####Fig.5F-------
####cnv color
color9<- colorRampPalette(c("#F9877D","#941B14"))(max(conbind_genomic[which(conbind_genomic[,6]=="Amp"),5])*50000)
color1<- colorRampPalette(c("#78A3CC","#014F9C"))(abs(max(conbind_genomic[which(conbind_genomic[,6]=="Del"),5])*50000))
color10<- colorRampPalette(c("#5f3c23"))(max(conbind_genomic[which(conbind_genomic[,6]=="NULL"),5])*50000)

####exp color
color2<- colorRampPalette(c("#F9877D","#941B14"))(max(conbind_genomic[which(conbind_genomic[,10]=="Up"),9])*50000)
color3<- colorRampPalette(c("#78A3CC","#014F9C"))(abs(max(conbind_genomic[which(conbind_genomic[,10]=="Down"),9])*50000))
color8<- colorRampPalette(c("#778899"))(max(conbind_genomic[which(conbind_genomic[,10]=="Nosig"),9])*50000)


####methylation color
color4<- colorRampPalette(c("#78A3CC","#014F9C"))(max(conbind_genomic[which(conbind_genomic[,4]=="up"),3])*50000)
color5<- colorRampPalette(c("#F9877D","#941B14"))(abs(max(conbind_genomic[which(conbind_genomic[,4]=="down"),3])*50000))
color6<- colorRampPalette(c("#778899"))(max(conbind_genomic[which(conbind_genomic[,4]=="Nosig"),3])*50000)


####mutation color
color7<- colorRampPalette(c("#afb4db","#917ab0"))(max(conbind_genomic[which(conbind_genomic[,8]=="NULL"),7])*50000)

op=par(mar=c(7,7,3,5))
cancer <- unique(conbind_genomic$cancer)
x <- rep(0,length(conbind_genomic$cancer))

for(i in 1:length(cancer)){
  index<-which(conbind_genomic$cancer == cancer[i])
  x[index]<-i
}

FOX <- rev(unique(conbind_genomic$gene))
y <- rep(0,length(conbind_genomic$gene))
for(i in 1:length(FOX)){
  index<-which(conbind_genomic$gene == FOX[i])
  y[index]<-i
}
plot(1,xlim=c(1,length(unique(x))+1),ylim=c(1,length(unique(y))+1),
     type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
###菱形四维图
for(i in 1:nrow(conbind_genomic)) {
  #CNV left-top
  if(conbind_genomic$CNV_type[i]=="Amp") {
    polygon(x[i]+c(0,0.5,0.5),y[i]+c(0.5,1,0.5),col=color9[conbind_genomic$copynumb[i]*40000])
  }
  if(conbind_genomic$CNV_type[i]=="Del") {
    polygon(x[i]+c(0,0.5,0.5),y[i]+c(0.5,1,0.5),col=color1[conbind_genomic$copynumb[i]*40000])
  }
  if(conbind_genomic$CNV_type[i]=="NULL") {
    polygon(x[i]+c(0,0.5,0.5),y[i]+c(0.5,1,0.5),col=color10[conbind_genomic$copynumb[i]*40000])
  }
  
  #methylation left-bottom
  if(conbind_genomic$meth_type[i]=="up"){
    polygon(x[i]+c(0,0.5,0.5),y[i]+c(0.5,0,0.5),col=color4[conbind_genomic$methylation[i]*40000])
  }
  if(conbind_genomic$meth_type[i]=="down"){
    polygon(x[i]+c(0,0.5,0.5),y[i]+c(0.5,0,0.5),col=color5[conbind_genomic$methylation[i]*40000])
  }
  if(conbind_genomic$meth_type[i]=="Nosig"){
    polygon(x[i]+c(0,0.5,0.5),y[i]+c(0.5,0,0.5),col=color6[conbind_genomic$methylation[i]*40000])
  }
  
  
  #exp right-top
  if(conbind_genomic$expre_type[i]=="Up"){
    polygon(x[i]+c(0.5,1,0.5),y[i]+c(1,0.5,0.5),col=color2[conbind_genomic$expression[i]*40000])
  }
  if(conbind_genomic$expre_type[i]=="Down"){
    polygon(x[i]+c(0.5,1,0.5),y[i]+c(1,0.5,0.5),col=color3[conbind_genomic$expression[i]*40000])
  }
  if(conbind_genomic$expre_type[i]=="Nosig"){
    polygon(x[i]+c(0.5,1,0.5),y[i]+c(1,0.5,0.5),col=color8[conbind_genomic$expression[i]*40000])
  }
  
  #mutiation right-bottom
  if(conbind_genomic$Muta_type[i]=="NULL"){
    polygon(x[i]+c(1,0.5,0.5),y[i]+c(0.5,0.5,0),col=color7[conbind_genomic$mutation[i]*40000])
  }
}

axis(1,at=sort(unique(x))+.5,labels=cancer,lty=0,las=2)
axis(2,at=sort(unique(y))+.5,labels=FOX,lty=0,las=2)




##################
######Fig.S10B----
##################
AMP_frq <- data.frame(gene=gene,UCEC=NA,UCS=NA,BRCA=NA,LUAD=NA,LUSC=NA,COAD=NA,READ=NA,OV=NA,PAAD=NA,TGCT=NA,THCA=NA)
DEL_frq <- data.frame(gene=gene,UCEC=NA,UCS=NA,BRCA=NA,LUAD=NA,LUSC=NA,COAD=NA,READ=NA,OV=NA,PAAD=NA,TGCT=NA,THCA=NA)
for (cancer in names(Epithelial_cancer)) {
  copynumb_df <- read.table(paste0("./Aging/data/UCSC/UCSC_copynumber_data/上皮细胞UTDT_copy/",cancer,"_persent.txt"),header=T,row.names=1,check.names=F)
  copynumb_df$gene <- rownames(copynumb_df)
  
  for (onegenes in gene) {
    if(onegenes %in% copynumb_df$gene){
      AMP_frq[match(onegenes,AMP_frq$gene),cancer] <- copynumb_df[onegenes,2]
      DEL_frq[match(onegenes,AMP_frq$gene),cancer] <- copynumb_df[onegenes,3]
    }else{
      AMP_frq[match(onegenes,AMP_frq$gene),cancer] <- 0
      DEL_frq[match(onegenes,AMP_frq$gene),cancer] <- 0
    }
  }
}

rownames(AMP_frq)<- AMP_frq$gene
rownames(DEL_frq)<- DEL_frq$gene
AMP_exp <- AMP_frq[,-1]
library(pheatmap)
levels <- c("S100A11","SDC1","PVRL4","SPTSSB","TLCD1","SLC2A1","F12","FAM83H","MEX3A","EFNA3","HMGA1","KREMEN2","TUBA1C","ESPL1",
            "LAMB3","MMP7","PLEKHN1","KRT19","TNS4","RPLP0P2","LIPM","TACSTD2","LRRC8E","ABCA12","C2orf54","PHLDA2","ANKRD22","PERP",                                          
            "SFN","POF1B","SERPINB5")
DEL_exp <- DEL_frq[levels,]
DEL_exp <- DEL_exp[,-1]
DEL_exp[DEL_exp>0.1]=0.1
bk <- c(seq(0,0.045,by=0.02),seq(0.06,0.1,by=0.02))
pheatmap(DEL_exp,
         scale = "none",fontsize = 8,border_color='white',clustering_method = 'ward.D2',
         legend_breaks = seq(0,0.1,0.02),
         breaks = bk,
         cluster_rows = F
)
AMP_exp[AMP_exp>0.1]=0.1
pheatmap(AMP_exp,
         scale = "none",fontsize = 8,border_color='white',clustering_method = 'ward.D2',
         legend_breaks = seq(0,0.1,0.02),
         breaks = bk
)

##################
######Fig.S10C----
##################
setwd("./Aging/data/UCSC/TCGA的突变数据")
load(file="./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/组织统计/Epithelial_cell/Epithelial_cancer_intersect.rda")
Epithelial_cancer <- intersect_gene[c("UCEC","UCS","BRCA","LUAD","LUSC","COAD","READ","OV","PAAD","TGCT","THCA")]

cancer <- names(Epithelial_cancer)[1]

mutation_All <- data.frame(gene=gene_top,UCEC=NA,UCS=NA,BRCA=NA,LUAD=NA,LUSC=NA,COAD=NA,READ=NA,OV=NA,PAAD=NA,TGCT=NA,THCA=NA)

for(cancer in names(Epithelial_cancer)){
  mutation_df <- read.table(paste0("./Aging/data/UCSC/TCGA的突变数据/Epithelial_UTDA/",cancer,"_mutation_persent.txt"),sep = "\t",header = T)
  
  for (gene in mutation_All$gene) {
    if(gene %in% mutation_df$gene){
      mutation_All[match(gene,mutation_All$gene),cancer] <- mutation_df$Mutation_persent[which(mutation_df$gene == gene)]
    }else{
      mutation_All[match(gene,mutation_All$gene),cancer] <- 0
    }
  }
}
rownames(mutation_All) <- mutation_All$gene
mutation_All <- mutation_All[-24,]
mutation_All <- mutation_All[,-1]

bk <- c(seq(0,0.045,by=0.02),seq(0.06,0.1,by=0.02))
pheatmap(mutation_All,
         scale = "none",fontsize = 8,border_color='white',clustering_method = 'ward.D2',
         legend_breaks = seq(0,0.1,0.02),
         breaks = bk
)

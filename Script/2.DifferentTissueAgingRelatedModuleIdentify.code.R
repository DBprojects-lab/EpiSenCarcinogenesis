library(DESeq2)
library(limma)
library(data.table)
library(MEGENA)
library(WGCNA)
Tissues <- c("Adrenal_Gland","Breast","Esophagus","Brain","Kidney","Prostate","Liver","Lung","Ovary","Pancreas","Colon",
             "Testis","Uterus","Skin","Stomach","Thyroid","Blood")
Tissues_proflie <- c("Adrenal_Gland_ACC","Breast_BRCA","Esophagus_ESCA","Brain_LGG_GBM","Kidney_KICH_KIRC_KIRP","Prostate_PRAD",
                     "Liver_LIHC","Lung_LUAD_LUSC","Ovary_OV","Pancreas_PAAD","Colon_COAD_READ","Testis_TGCT","Uterus_UCEC_UCS","Skin_SKCM",
                     "Stomach_STAD","Thyroid_THCA","Blood_DLBC_THYM")

##################################################################
######MEGENA identifies co-expression modules across tissues######
##################################################################
for (i in 1:length(Tissues_proflie)) {
  setwd(paste0("./aging/data/UCSC/GTEX_RNAseq_data/",Tissues_proflie[i]))
  
  datExpr=read.table(paste0(Tissue[i],".vsd.adjusted.txt"),header=T,row.names=1,check.names=F)
  n.cores <- 5 # Number of cores/threads to invoke PCP
  doPar <-TRUE # Whether to set up parallelism
  method = "pearson" # Correlation calculation method pearson or spearman
  FDR.cutoff = 0.05 # Defining the FDR Threshold
  module.pval = 0.05 # Threshold for assuming that the module is meaningful. Recommended is 0.05.
  hub.pval = 0.05 # Thresholds for meaningful connectivity obtained from random tetrahedral networks
  cor.perm = 5 # FDRs used to calculate all relevant pairs
  hub.perm = 100 # Used to calculate the connectivity significance p-value
  # Comments to be executed downstream
  annot.table=NULL
  id.col = 1
  symbol.col= 2
  
  #Calculating gene-to-gene correlations
  rho.out = calculate.rho.signed(datExpr,n.perm = 10,FDR.cutoff = FDR.cutoff,estimator = method,
                                 use.obs = "na.or.complete",
                                 direction = "absolute",
                                 rho.thresh = NULL,sort.el = TRUE)
  save(rho.out,file="rho.out.rda")
  # calculate PFN
  #In this step, Planar Filtered Network (PFN) is calculated by taking significant correlation pairs, ijw. In the case of utilizing a different similarity measure, one can independently format the results into 3-column data frame with column names c("row","col","weight"), and make sure the weight column ranges within 0 to 1. Using this as an input to calculate.PFN() will work just as fine. 
  #### register multiple cores if needed: note that set.parallel.backend() is deprecated. 
  run.par = doPar & (getDoParWorkers() == 1) 
  if (run.par)
  {
    cl <- parallel::makeCluster(n.cores)
    registerDoParallel(cl)
    # check how many workers are there
    cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
  }
  
  el <- calculate.PFN(rho.out$signif.ijw,doPar = doPar,num.cores = n.cores,keep.track = FALSE)
  save(el,file = "el.rda")
  
  g <- graph.data.frame(el,directed = FALSE)
  
  # Multi-scale cluster analysis by MCA clustering. ‘MEGENA.output’’ is the core output used for downstream analysis summarisation and mapping.
  MEGENA.output <- do.MEGENA(g,
                             mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE,
                             min.size = 10,max.size = vcount(g)/2,
                             doPar = doPar,num.cores = n.cores,n.perm = hub.perm,
                             save.output = FALSE)
  
  save(MEGENA.output,file = "MEGENA.output.rda")
  #Summary of results
  summary.output <- MEGENA.ModuleSummary(MEGENA.output,
                                         mod.pvalue = module.pval,
                                         hub.pvalue = hub.pval,
                                         min.size = 10,
                                         max.size = vcount(g)/2,
                                         annot.table = annot.table,
                                         id.col = id.col,
                                         symbol.col = symbol.col,
                                         output.sig = TRUE)
  save(summary.output,file = "summary.output.rda")
  
  if (!is.null(annot.table))
  {
    V(g)$name <- paste(annot.table[[symbol.col]][match(V(g)$name,annot.table[[id.col]])],V(g)$name,sep = "|")
    summary.output <- output[c("mapped.modules","module.table")]
    names(summary.output)[1] <- "modules"
  }
  #Summary of module information
  module.output <- module_convert_to_table(MEGENA.output,mod.pval = 0.05,
                                           hub.pval = 0.05,min.size = 10,max.size=vcount(g)/2)
  save(module.output,file = "module.output.rda")
}


#########################################################
######Identify age-related modules in organisations######
#########################################################

module_total <- data.frame(tissue = Tissues,
                           AllModuleNum=0,
                           UpModuleNum=0,
                           DownModuleNum=0)

for (i in 1:length(Tissues_proflie)) {
  setwd(paste0("./Aging/data/UCSC/GTEX_RNAseq_data/",Tissues_proflie[i]))
  load(file = "summary.output.Rdata")
  load(file = "module.output.rda")
  module.table <- summary.output$module.table
  ModuleList <- summary.output$modules
  datTraits=read.table(paste0(Tissue[i],"_sampleInfomation.txt"),header=T,row.names=1)
  
  #Taking the median of the age segments
  datTraits$AgeC=as.numeric(substring(datTraits$Age,1,2))+5
  
  moduleLabelList=unique(module.table$module.id)
  nSamples=ncol(datExpr)
  moduleAgeCorResult=matrix(data = NA,nrow = length(moduleLabelList), ncol = 3,dimnames = NULL)
  module_Eigengenes=matrix(data = NA,nrow = nSamples, ncol = length(moduleLabelList),dimnames = NULL)
  for(i in 1:length(moduleLabelList)){
    geneList=module.output[module.output$module==moduleLabelList[i],"id"]
    geneExpr=datExpr[geneList,]
    moduleLabel=module.output[module.output$module==moduleLabelList[i],"module"]
    MEs0 = moduleEigengenes(t(geneExpr),moduleLabel)$eigengenes
    module_Eigengenes[,i] <- MEs0[,1]
    all(rownames(MEs0)==rownames(datTraits)) #TRUE
    moduleAgeCor = cor(MEs0, datTraits$AgeC,method = "spearman");
    moduleAgePvalue = corPvalueStudent(moduleAgeCor, nSamples);
    moduleAgeCorResult[i,1]=moduleLabelList[i]
    moduleAgeCorResult[i,2]=moduleAgeCor
    moduleAgeCorResult[i,3]=moduleAgePvalue
  }
  #Eigengenes value of the module (module eigenvalue)
  module_Eigengenes=data.frame(module_Eigengenes)
  rownames(module_Eigengenes) <- rownames(datTraits)
  colnames(module_Eigengenes) <- moduleLabelList
  #Correlation coefficients between modules and age and p-values
  moduleAgeCorResultTmp=data.frame(moduleAgeCorResult)
  colnames(moduleAgeCorResultTmp)=c("module.id","Cor","Pvalue")
  moduleAgeCorResultTmp=moduleAgeCorResultTmp[order(moduleAgeCorResultTmp$Pvalue),]
  moduleAgeCorResultTmp$Cor <- as.numeric(moduleAgeCorResultTmp$Cor)
  moduleAgeCorResultTmp$Pvalue <- as.numeric(moduleAgeCorResultTmp$Pvalue)
  moduleAgeCorResultTmp$Pattern=ifelse(moduleAgeCorResultTmp$Pvalue<0.01,"Sig","NoSig")
  table(moduleAgeCorResultTmp$Pattern)
  sigModule=moduleAgeCorResultTmp[moduleAgeCorResultTmp$Pattern=="Sig",]##Module and age correlation information
  write.table(sigModule,file="sigCluster.txt",sep="\t",row.names=F,quote=F)
  
  #Significantly Related Module Information
  moduleLabel=module.output[module.output$module%in%sigModule$module,]
  write.table(moduleLabel,file = "moduleLable.txt",sep = "\t",row.names = F,quote = F)
 
  sigupModule <- sigModule[sigModule$type=="Up",]
  sigdownModule <- sigModule[sigModule$type=="Down",]
  #positive module list
  UpModuleList <- ModuleList[sigupModule$module]
  #negative module list
  DownModuleList <- ModuleList[sigdownModule$module]
  
  AgingRelatedModuleList <- c(list(UpModuleList),list(DownModuleList))
  names(AgingRelatedModuleList) <- c(paste0(Tissues[i],"_up_modules"),paste0(Tissues[i],"_down_modules"))
  save(AgingRelatedModuleList,file = paste0(Tissues,"_AgingRelatedModuleList.rda"))
  ###positive_genes
  age_up_gene <- module.output[module.output$module%in%sigupModule$module,]
  ###negative_genes
  age_down_gene <- module.output[module.output$module%in%sigdownModule$module,]
  
  AgingRelatedGeneList <- c(list(unique(age_up_gene$id)),list(unique(age_down_gene$id)))
  names(AgingRelatedGeneList) <- c(paste0(Tissues[i],"_up_gene"),paste0(Tissues[i],"_down_gene"))
  save(AgingRelatedGeneList,file = paste0(Tissues,"_AgingRelatedGeneList.rda"))
  module_total[x,2] <- length(module.table[,1])
  module_total[x,3] <- length(UpModuleList)
  module_total[x,4] <- length(DownModuleList)
}
save(module_total,"./Aging/data/UCSC/GTEX_RNAseq_data/module_total.rda")

#----Figure 1B------
load("./Aging/data/UCSC/GTEX_RNAseq_data/module_total.rda")
need <- module_total
for (i in 1:17) {
  need$tissue[i] <- paste0(need$tissue[i],"(n=",module_total$sample_num[i],")")
}
need <- need[order(need$age_association_module,decreasing = F),]
need$tissue <- factor(need$tissue,levels = need$tissue)
need <- need[,-3]
need1 <- reshape2::melt(need,id.vars = 'tissue')


ggplot(need1,aes(tissue,value,fill=variable))+
  geom_bar(stat="identity",position="stack",color ="black",width = 0.55,size = 0.25)+
  scale_fill_manual(values = c("age_up_module"="#E03833",
                               "age_down_module"="#1B5C93",
                               "no_association_module"="#C7D6DB"))+
  coord_flip()+
  theme_bw()+#经典主题s
  geom_text(aes(label=value,y=value+20),position="stack",color="black")+
  theme(
    axis.title = element_text(size = 15,face = "plain",color = "black"),
    axis.text = element_text(size = 12,face = "plain",color = "black"),
    legend.title = element_text(size = 14,face = "plain",color = "black"),
    legend.position = "none",
    panel.grid=element_blank ()
  )

###########################################################################
#######Cyclic plotting of the Sunburst charts for age-related modules######
###########################################################################
#----Figure S1, Figure 1C
library(ggplot2)
library(ggpubr)
library(DEGreport)
for (i in 1:length(Tissues_proflie)) {
  setwd(paste0("./Aging/data/UCSC/GTEX_RNAseq_data/",Tissues_proflie[i]))
  load(file = "summary.output.Rdata")
  sigModule<- read.table(file="sigCluster.txt",sep="\t",header = T)
  mdf <- summary.output$module.table
  mdf$category <- 0
  mdf[sigModule$module.id,]$category <- sigModule$type
  mdf$category[which(mdf$category == "0")] <- "NO"
  mdf$category <- factor(mdf$category,levels = c("NO","Down","Up"))
  mdf_rnames <- mdf
  mdf_rnames$module.id <-  gsub("c._", "M_", mdf_rnames$module.id)
  mdf_rnames$module.parent <-  gsub("c._", "M_", mdf_rnames$module.parent)
  sbobj = draw_sunburst_wt_fill(module.df = mdf_rnames,feat.col = "category",
                                fill.type = "discrete",
                                fill.scale = scale_fill_manual(values = c("NO"="#E6E6E6","Up" = "#C51605","Down" = "#337CCF")), 
                                id.col = "module.id",parent.col = "module.parent",
                                border.col = "black", # sunburst border color
                                border.width = 0.15, # sunburst border line width
                                theme.adjust = theme_classic2())
  pdf(file = paste0("./Aging/data/UCSC/GTEX_RNAseq_data/Sunburst charts/",Tissues[i],".pdf"),width=10,height=8.5)
  print(sbobj)
  dev.off()
}

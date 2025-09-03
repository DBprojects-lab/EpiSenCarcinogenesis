# EpiSenCarcinogenesis
Analysis codes and documentations for "Integrative Transcriptomic Analysis Decodes the InterplayBetween Aging, Senescence, and Cancer".



![](Figure/Grapg Abstract.svg)

## Codes used to identify Aging-dependent Senescence-associated Coregulated Modules(SACMs) in different tissues
|Script name|Comment| 
|-----------|-------| 
|[1.DifferentTissueBlukRNA-seqDataProcessing.code.R](Script/1.DifferentTissueBlukRNA-seqDataProcessing.code.R)|This script was used to deal with bulkRNA-seq data from different tissues|
|[2.DifferentTissueAgingRelatedModuleIdentify.code.R](Script/2.DifferentTissueAgingRelatedModuleIdentify.code.R)|This script was used to identify age-dependent co-expression modules in different tissues|
|[3.SACMs-IdentifyAndStatisticians.code.R](Script/3.SACMs-IdentifyAndStatisticians.code.R)|This script was used to identity SACMs in different tissues|


## Codes used to identification of the cell type specific identity of SACMs
|Script name|Comment| 
|-----------|-------| 
|[4.Single-CellRNASeq-Processing.code.R](Script/4.Single-CellRNASeq-Processing.code.R)|This script was used to deal with Single-cell RNA-seq data from different tissues|
|[5.SACMsCellType-Annotation.code.R](Script/5.SACMsCellType-Annotation.code.R)|This script was used to map the gene list within each SACM to the cell type specific|
|[6.CellTypeSpecific-SAS_Identify.code.R](Script/6.CellTypeSpecific-SAS_Identify.code.R)|This script was used to analysis tissues-Share cell type specific SAS|

## Codes used to analyse the overlap between aging and oncogenesis
|Script name|Comment| 
|-----------|-------| 
|[7.AgingAndTumorigenesisAnalysis.code.R](Script/7.AgingAndTumorigenesisAnalysis.code.R)|This script was used to deal with Single-cell RNA-seq data from different tissues|
|[8.DifferentAgeGroupAgingAndTumorigenesis.code.R](Script/8.DifferentAgeGroupAgingAndTumorigenesis.code.R)|This script was used to map the gene list within each SACM to the cell type specific|

## Codes used to analyse relationship between cellular senescence and carcinogenesis at the single-cell resolution.
|Script name|Comment| 
|-----------|-------| 
|[9.CellSpecific-SASOverlapRatio.code.R](Script/9.CellSpecific-SASOverlapRatio.code.R)|This script was used to calculate the percentage of cell type-specific SAS between tissue aging and its correspondent cancer|
|[10.CellSpecific-SASAndTumor-DEGsOverlap.code.R](Script/10.CellSpecific-SASAndTumor-DEGsOverlap.code.R)|This script was used to overlap cell type specific-SAS and cancer-DEGs across cancers|

## Codes used to analyse multi-omics alterations in epithelial-associated senescence-associated traits (EPC-SAS).
|Script name|Comment| 
|-----------|-------| 
|[11.DNAMethylation_analysis.code.R](Script/11.DNAMethylation_analysis.code.R)|This script was used to analyse epigenetic alterations in epithelial-associated senescence-associated traits (EPC-SAS)|
|[12.CNVAndMutation_analysis.code.R](Script/12.CNVAndMutation_analysis.code.R)|This script was used to analyse genomic alterations in epithelial-associated senescence-associated traits (EPC-SAS)|

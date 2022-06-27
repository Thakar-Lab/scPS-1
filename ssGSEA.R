#### Set up ####
if(!require(GSEABase))install.packages("GSEABase")
if(!require(genefilter))install.packages("genefilter")
if(!require(Biobase))install.packages("Biobase")
if(!require(stringr))install.packages("stringr")
if(!require(GSVA))BiocManager::install("GSVA")
if(!require(RColorBrewer))install.packages("RColorBrewer")
if(!require(umap))install.packages("umap")
if(!require(forcats))install.packages("forcats")
if(!require(MASS))install.packages("MASS")
if(!require(ggplot2))install.packages("ggplot2")
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))BiocManager::install("mindr")
if(!require(tidyverse))BiocManager::install("tidyverse")
if(!require(AUCell))BiocManager::install("AUCell")
rm(list=ls()) 
gc()


#### load TF binding info ####
# read in TF predicted binding sets - updated v2020
predBind <- GSEABase::getGmt(file.path("TF binding v2020- human motifs and human genes 90pct- unique TF names.GMT"))
numTF <- length(unique(names(predBind)))
# FILTER OUT MOTIFS THAT HAVE EITHER TOO MANY OR TOO FEW PREDICTED SITES
predBindCounts <- sort(setNames(unlist(lapply(predBind, function(x){length(geneIds(x))})), names(predBind)))
predBind_maxMin <- list(max = 10000, min = 10) # SPECIFY CUTOFFS
predBindCounts_postThresh <- predBindCounts[(predBindCounts >= predBind_maxMin$min)&
                                              (predBindCounts <= predBind_maxMin$max)]
predBind <- predBind[names(predBind) %in% names(predBindCounts_postThresh)]

#RunPCA
load(file="/scratch/rwang46/spScore/Participant.integrated.RData")
DefaultAssay(Participant.integrated) <- "RNA"

#### ssGSEA score ####
Participant.integrated <- ScaleData(Participant.integrated, features =  rownames(Participant.integrated)) 
mRNAExpr <- as.matrix(Participant.integrated[['RNA']]@counts)
ssGESA <- gsva(mRNAExpr, predBind, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
write.csv(ssGESA,file = 'Participant.integrated.ssGSEA.csv')

#### AUCell score ####
cell_rankings <- AUCell_buildRankings(mRNAExpr)
cell_AUC <- AUCell_calcAUC(predBind,cell_rankings)
AUCell <- mat.or.vec(nr = length(names(cell_AUC)),nc = ncol(mRNAExpr))
rownames(AUCell) <- names(cell_AUC)
colnames(AUCell) <- colnames(mRNAExpr)
write.csv(AUCell,file = 'Participant.integrated.AUCell.csv')
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
Participant.integrated <- ScaleData(Participant.integrated, features =rownames(Participant.integrated[['RNA']]),  assay="RNA" ) 
cal_cv=function(x){
  y=na.omit(x)
  return(sd(y)/mean(y))
}
CoeffientsVariation <- apply(Participant.integrated@assays$RNA@counts, 1, cal_cv)
# RunPCA for genesize <= 500 
tfTgtPca<- list()
for(currTF in names(predBind)){
  if (predBindCounts[[currTF]] <= 500){
    tfTgtPca[[currTF]] <- RunPCA(Participant.integrated, features = geneIds(predBind[[currTF]]), npcs = 10)
  } 
  else {
    tfTgtPca[[currTF]] <- RunPCA(Participant.integrated, features = names(CoeffientsVariation[order(CoeffientsVariation[names(CoeffientsVariation) %in% geneIds(predBind[[currTF]])],decreasing=TRUE)[1:500]]), npcs = 10)
  }
}
rm(currTF)

# PCA.spScore
# for each case, determine how many PCs it takes to get to x% variance explained
varExplThresh <- 0.5 # target for cumulative max var explained
cumulVarExp <- function(sdevPC){cumsum(sdevPC^2)/sum(sdevPC^2)}
varExp <- function(sdevPC){(sdevPC^2) / sum(sdevPC^2)}
maxCombPC <- unlist(lapply(tfTgtPca, FUN = function(x){
  idx = min(which(cumulVarExp(x[['pca']]@stdev) > varExplThresh));
  if(identical(idx, integer(0))){return(1)}else{return(max(idx))}
}))
table(maxCombPC)

mRNAExpr <- Participant.integrated@assays$RNA@counts
# single person transcription factor score (spTFscore)
aggregPC.unscaled <- mat.or.vec(nr = length(names(maxCombPC)[maxCombPC > 1]),
                                nc = ncol(mRNAExpr))
rownames(aggregPC.unscaled) <- names(maxCombPC)[maxCombPC > 1]
colnames(aggregPC.unscaled) <- colnames(mRNAExpr)

# for each TF after filtering  
for(currTF in rownames(aggregPC.unscaled)){
  # get the var% explained per PC
  currVarExp = varExp(tfTgtPca[[currTF]]@reductions$pca@stdev)
  # weight PC scores
  currWeightedPC = sapply(2:maxCombPC[currTF], FUN =
                            function(i){((tfTgtPca[[currTF]]@reductions$pca@cell.embeddings[,i] ^ 2) * currVarExp[i])})
  if(ncol(currWeightedPC) > 1){
    aggregPC.unscaled[ currTF , rownames(currWeightedPC) ] <-
      sqrt(rowSums(currWeightedPC[ , 1:(maxCombPC[currTF]-1)]))
  }else{
    aggregPC.unscaled[ currTF , rownames(currWeightedPC) ] <-
      sqrt(sum(currWeightedPC[ , 1:(maxCombPC[currTF]-1)]))
  }
  rm(currWeightedPC, currVarExp)
}
write.csv(aggregPC.unscaled,file = 'Participant.integrated.spScore.csv')

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


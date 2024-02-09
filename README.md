# scPS (single-cell Pathway Score) 
## scPS is a single-cell RNA-seq gene set analysis (scGSA) method using the functions in R script (without installation)

##   Input Data: 
### 1. A Seurat object (Seurat V5 are required), 
### 2. A gene matrix transposed (GMT) file format of gene sets/pathways. Each gene set is described by a name, a description, and the genes in the gene set. 

##   Output:
###   1. Single cells with corresponding scores for gene sets/pathways data matrix, where rows correspond to gene sets/pathways and columns correspond to cells.

## How to Apply scPS Function
### Seurat_data = Input Seurat object
### GeneSet = Gene set gmt file
Result  =  scPS(Seurat_data,GeneSet) # calling scPS 

## libraries required
library(Seurat) #Seurat V5
library(GSEABase)
library(genefilter)




##     Authors:  Ruoqiao Wang
##     Email:    RuoQiao_Wang@URMC.Rochester.edu

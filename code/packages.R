if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")
}

if (!requireNamespace("lncRNAtools", quietly = TRUE)){
  
  if (!requireNamespace("edgeR", quietly = TRUE))
    BiocManager::install("edgeR")
  if (!requireNamespace("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")
  if (!requireNamespace("BiocParallel", quietly = TRUE))
    BiocManager::install("BiocParallel")
  if (!requireNamespace("TFEA.ChIP", quietly = TRUE))
    BiocManager::install("TFEA.ChIP")
  if (!requireNamespace("scatterplot3d", quietly = TRUE))
    BiocManager::install("scatterplot3d")
  if (!requireNamespace("DOSE", quietly = TRUE))
    BiocManager::install("DOSE")
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    BiocManager::install("clusterProfiler")
  devtools::install_github('zexllin/lncRNAtools')
}
if (!requireNamespace("corrplot", quietly = TRUE)){
  install.packages("corrplot")
}

if (!requireNamespace("GDCRNATools", quietly = TRUE)){
  BiocManager::install("GDCRNATools")
}

suppressMessages(library(lncRNAtools))
suppressMessages(library(limma))
suppressMessages(library(edgeR))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(clusterProfiler))
suppressMessages(library(corrplot))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(survival))
suppressMessages(library(GDCRNATools))
suppressMessages(library(survminer))
suppressMessages(library(glmnet))
#suppressMessages(library())

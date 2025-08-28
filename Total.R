
# Packages

  library(Seurat)
  library(infercnv)
  library(ROGUE)
  library(harmony)
  library(alakazam)
  library(Startrac)
  library(destiny)
  library(SingleCellExperiment)
  library(CytoTRACE)
  library(SCENT)
  library(SingleR)
  library(slingshot)
  library(scMetabolism)
  library(CellChat)
  
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(rlist)
  library(distdimscr)
  library(Hmisc)
  library(smoother)
  library(openxlsx)
  library(clusterProfiler)
  library(msigdbr)
  
  library(ggplot2)
  library(ggridges)
  library(ggthemes)
  library(ggpubr)
  library(ggrastr)
  library(ggforce)
  library(ggsignif)
  library(ggpmisc)
  library(ggrepel)
  library(plotly)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(viridis)
  library(TreeAndLeaf)
  library(RedeR)
  library(igraph)
  library(corrplot)
  
# Functions

  fucs <- list.files("functions",full.names=T)
  lapply(fucs,function(x) {source(x,encoding="utf-8");return("Yes")})


# colors

  orig.ident <- c(paste0("P",1:4),paste0("M",1:8))
  orig.ident.colors <- c("#1f77b4","#aec7e8","#2ca02c","#98df8a",
                         "#ff7f0e","#ffbb78","#d62728","#ff9896","#9467bd","#c5b0d5","chocolate","pink3")
  names(orig.ident.colors) <- orig.ident
  
  Risk_type.colors <- c(Low_risk="cyan3",High_risk="brown3")
  
  Type.colors <- c(P="deepskyblue2",M="indianred2")
  
  Response.colors <- c(SD="seagreen3",PD="#9467bd")
  
  KRAS_mut.colors <- c("KRAS"="red")
  
  Moffit_type.colors <- c(Classical="seagreen3",`Basal-like`="tomato2")
  
  CellType.colors <- c(Malignant="orangered2",Ductal="#8c564b",Acinar="yellow",Cholangiocyte="#DF7BBC",Hepatocyte="#e1ab8c", 
                       CD4T="seagreen2",CD8T="royalblue1",Prof.T="red4",gdT="darkviolet",NK="salmon",ILC="deeppink",      
                       B="cyan2",Plasma="#1f77b4",  
                       Monocyte="#ff9896",Macrophage="gold",cDC="#9EDAE5",Neutrophil="#b83570",Mast="#e377c2",pDC="#BBA0BD",    
                       TEC="brown3",HSEC="#7DCF72",   
                       myCAF="indianred2",iCAF= "#45C7D6",vCAF= "#A26DB7")

  
  
# 确保图片的可重复性
  
  # 56463 cells
  object <- readRDS("E:/PAAD_Project/Data/Intergrated/object.rds")
  # 28001 cells
  object_l <- readRDS("E:/PAAD_Project/Data/Lymphocyte/object_l.rds")
  # 2435 cells
  object_treg <- readRDS("E:/PAAD_Project/Data/Lymphocyte/object_treg.rds")
  # 4045 cells
  object_tex <- readRDS("E:/PAAD_Project/Data/Lymphocyte/object_tex.rds")
  # 2539 cells
  object_nk <- readRDS("E:/PAAD_Project/Data/Lymphocyte/object_nk.rds")
  # 6945 cells
  object_m <- readRDS("E:/PAAD_Project/Data/Myeloid/object_m.rds")
  # 5156 cells
  object_mm <- readRDS("E:/PAAD_Project/Data/Myeloid/object_mm.rds")
  
  
  # 预后相关基因
  load(file.path("E:/PAAD_Project/Data/MetaSurvGenes/meta_result.rda"))
  # 1966
  adverse_genes <- meta_result[gene_type == "PoorPrognosis",gene_name]   
  # 1424
  favorable_genes <- meta_result[gene_type == "FavorablePrognosis",gene_name]
  

  
  
  
  
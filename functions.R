
# for eg. show_col(ColAssign(letters[1:20]))
ColAssign <- function(Var,palettes="Classic 20"){
  require(ggthemes);require(RColorBrewer)
  pal <- tableau_color_pal(palette = palettes,direction = 1,type="regular")
  if (length(Var) > 20) {
    palOut <- colorRampPalette(pal(20))(length(Var))
    names(palOut) <- Var
  } else if (length(Var) == 20) {
    palOut <- pal(20)
    names(palOut) <- Var
  } else if (length(Var) < 20) {
    palOut <- pal(20)
    palOut <- setdiff(palOut,c("#7f7f7f","#c7c7c7"))# remove grey colors
    #palOut <- sample(palOut)
    #palOut <- c(palOut,c("#7f7f7f","#c7c7c7"))
    palOut <- palOut[1:length(Var)]
    names(palOut) <- Var
  }
  return(palOut)
}


diff_calculate <- function(obj,ident_1="High_risk",group="Risk_type",pct.thres=0,log2FC.thres=0.25,p.thres=0.05){
  
  diff_all <- FindMarkers(obj,ident.1 = ident_1,group.by=group,logfc.threshold = 0)
  diff_all <- diff_all %>%
    mutate(pct.diff=pct.1-pct.2)%>%
    rownames_to_column("gene") %>%
    arrange(desc(avg_log2FC))
  diff_all <- data.table(diff_all)
  
  diff_sig <- diff_all
  diff_sig$type <- "NS"
  diff_sig[pct.diff > pct.thres & avg_log2FC > log2FC.thres & p_val_adj < p.thres, type:="Up"]
  diff_sig[pct.diff < pct.thres & avg_log2FC < -log2FC.thres & p_val_adj < p.thres, type:="Down"]
  
  return(diff_sig) 
}


groupMeans <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}


phenoPropotion <- function(data,x_name,y_name,cols,legend_name,bar_width=0.9,out_path,plot_width=5,plot_height=6.5){
  
  df_stats <- data %>% 
    dplyr::group_by(get(x_name),get(y_name),.drop = FALSE) %>% 
    dplyr::summarise(n=n())
  df_stats <- as.data.frame(df_stats)
  names(df_stats)[1:2] <- c(x_name,y_name)
  df_stats <- na.omit(df_stats)
  df_stats[,y_name] <-  factor(as.character(df_stats[,y_name]),levels=names(cols))
  
  ggplot(df_stats, aes(x = get(x_name), y = n, fill = get(y_name))) + 
    geom_bar(position = "fill",stat = "identity",na.rm=T,width=bar_width) + #position="stack" gives numbers
    scale_fill_manual(values= cols) +
    theme_classic()+ 
    theme(legend.title=element_text(size = 8),
          legend.text=element_text(size = 6),
          legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm"),
          axis.title.y = element_text(size = 9),
          axis.title.x = element_text(size = 8),
          #axis.text.x = element_text(size = 8,angle =30, vjust = 1, hjust=0.75),
          axis.text.x = element_text(size = 7,angle =90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 8))+
    xlab("")+ylab("Cell Fraction (%)")+labs(fill=legend_name)+
    scale_y_continuous(labels = scales::percent)
}



scCluster <- function(obj,nfeature=2500,min.d=0.3,res=1,TCRBCR=F,Ribo_MT=T,harmony=F,seed=1,k_param=30,npcs=30) {
  
  if (TCRBCR) {
    TCR.genes <- grep("^TR[AB][VJ]",rownames(obj),value = T) # keep GD TCR genes for GDT cells
    BCR.genes <- grep("^IG[KHL][VJC]",rownames(obj),value = T)
  } else {TCR.genes <- c();BCR.genes <- c()}
  if(Ribo_MT){
    MT.genes <- grep("^MT-",rownames(obj),value = T)
    RP.genes <- grep("^RP[SL]",rownames(obj),value = T)
  }else {MT.genes <- c();RP.genes <- c()}
  
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj,nfeatures = nfeature)
  var.genes <- VariableFeatures(obj)
  var.genes <- setdiff(var.genes,c(TCR.genes,BCR.genes,MT.genes,RP.genes))
  VariableFeatures(obj) <- var.genes
  obj <- ScaleData(obj,features = var.genes)
  obj <- RunPCA(obj, verbose = FALSE,features = VariableFeatures(obj),npcs=npcs)
  
  if(harmony){
    obj <- RunHarmony(obj, group.by.vars=c("orig.ident"),assay.use ="RNA")
    obj <- FindNeighbors(obj, dims=1:30,reduction = "harmony",k.param = k_param)
    obj <- FindClusters(obj,resolution=res,random.seed=seed,graph.name = 'RNA_snn')
    obj <- RunUMAP(obj,reduction = "harmony",seed.use = seed,dims=1:30,
                   umap.method='uwot',min.dist=min.d,spread=1)
  }else{
    obj <- FindNeighbors(obj, dims=1:30,reduction = "pca",k.param = k_param)
    obj <- FindClusters(obj,resolution=res,random.seed=seed,graph.name = 'RNA_snn')
    obj <- RunUMAP(obj,reduction = "pca",seed.use = seed,dims=1:30,
                   umap.method='uwot',min.dist=min.d,spread=1)
  }
  
  return(obj)
}



# expr: data.frame/matrix, row: gene symbols; col: samples
subtype_calculate <- function(expr,markers_2,markers_1,name_subtype,name_2,name_1){
  
  samples <- colnames(expr)
  
  expr <- expr[rowSums(expr != 0) > 0, ]
  markers_2 <- intersect(markers_2,row.names(expr))
  markers_1 <- intersect(markers_1,row.names(expr))
  expr <- expr[c(markers_2,markers_1),]
  
  expr <- t(apply(expr, 1, scale))
  colnames(expr) <- samples
  
  subtype_ttest <- apply(expr,2,function(x) t.test(x[markers_2],x[markers_1]))
  subtype_t <- unlist(lapply(subtype_ttest,function(x) as.numeric(x$statistic)))
  subtype_p <- unlist(lapply(subtype_ttest,function(x) x$p.value))
  
  subtype_result <- data.table(samples=colnames(expr),subtype_t=subtype_t,subtype_p=subtype_p)
  subtype_result[subtype_t > 0,subtype_result:=name_2]
  subtype_result[subtype_t <= 0,subtype_result:=name_1]
  names(subtype_result) <- c("samples",paste0(name_subtype,"_t"),paste0(name_subtype,"_p"),paste0(name_subtype,"_subtype"))
  
  return(subtype_result)
}
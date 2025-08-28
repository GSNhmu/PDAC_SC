

# Lymphocytes

  object_l <- subset(object,subset=CellType_Main == "Lymphocyte")
  object_l <- subset(object_l,subset=CellType2 %in% c("Prof.T","ILC"),invert=T)
  object_l$CellType2[object_l$CellType == "B"] <- "B"
  CellType2.colors <-  c(B="cyan4",Plasma ="#1f77b4",
                         CD4T_Naive="gold2",CD4T_Mem="seagreen3",CD4T_Tfh="#45C7D6",CD4T_Treg="deeppink2",
                         CD8T_Naive="#8c564b",CD8T_Mem="#e1ab8c",CD8T_MAIT="#b83570",CD8T_Tex="royalblue1",
                         gdT="darkviolet",NK="orangered2") 
  object_l$CellType2 <- factor(object_l$CellType2,levels = names(CellType2.colors))                     
  object_l$Type <- factor(object_l$Type,levels = names(Type.colors))
  object_l$Risk_type <- factor(object_l$Risk_type,names(Risk_type.colors))
  
  object_l <- scCluster(object_l,nfeature=2000,min.d=0.5,res=1.5,TCRBCR=T,Ribo_MT=F,harmony=T,seed=2,k_param=50)
  object_l <- RunTSNE(object_l,reduction = "harmony",dims=1:30,perplexity = 5,seed.use = 4)

  

# Fig3A
  
  DimPlot(object_l,reduction = "tsne",group.by="CellType2",cols=CellType2.colors,label=T,raster=T,pt.size=2,label.size = 1)

  
  
# Fig3B

  colour_bk <- c(colorRampPalette(c("cyan3", "lightblue"))(40),  
                 colorRampPalette(c("lightblue", "linen"))(10),
                 colorRampPalette(c("linen", "lightsalmon2"))(5), 
                 colorRampPalette(c("lightsalmon2", "brown2"))(5),
                 colorRampPalette(c("brown2", "red"))(35))
  FeaturePlot(object_l,"Risk_score",raster=T,pt.size=2,reduction = "tsne")+NoAxes()+
    scale_colour_gradientn(colours = colour_bk)


    
# Fig3C
  
  CellType2_select <- rowSums(table(object_l$CellType2,object_l$Risk_type) > 50)==2
  CellType2_select <- names(CellType2_select[which(CellType2_select)])
  
  diff_CellType2 <- lapply(CellType2_select,function(c){
    print(c)
    diff_c <- FindMarkers(object_l[,object_l$CellType2==c],ident.1="High_risk",group.by="Risk_type",logfc.threshold = 0.25,min.pct = 0.1)
    diff_c$gene <- row.names(diff_c)
    diff_c$CellType2 <- c
    setDT(diff_c)
    setorder(diff_c,"avg_log2FC")
    return(diff_c)
  })
  diff_CellType2 <- rbindlist(diff_CellType2)
  diff_CellType2 <- diff_CellType2[p_val_adj < 0.05]
  
  cellID_DEGs <- lapply(CellType2_select,function(c){
    print(c)
    mat <- LayerData(object_l,"counts")[diff_CellType2[CellType2==c,gene],object_l$cellID[object_l$CellType2 == c]]
    index <- colSums(mat>0) 
    cellID_DEGs <- data.table(cellID=names(index),DEGs=index)
    return(cellID_DEGs)
  })
  cellID_DEGs <- rbindlist(cellID_DEGs)
  
  object_l@meta.data <- left_join(object_l@meta.data,cellID_DEGs)
  row.names(object_l@meta.data) <- object_l$cellID
  
  FeaturePlot(object_l,features = "DEGs",pt.size = 2,raster = T,reduction = "tsne")+
    ggtitle("Number of DEGs (High_risk vs. Low_risk)")+
    scale_color_gradientn(colours = c(colorRampPalette(c("linen","#fdf2b5"))(5),
                                      colorRampPalette(c("#fdf2b5","#F67B51"))(10),
                                      colorRampPalette(c("#F67B51","red"))(15)),
                          na.value = "grey90",)+NoAxes()
  
  
  
# Fig3D
  
  hetero_data <- data.table(orig.ident=object_l$orig.ident,
                            CellType2=object_l$CellType2,
                            Risk_type=object_l$Risk_type)
  hetero_data[,risk_diversity_shannon:=vegan::diversity(table(Risk_type),index = "shannon"),by=c("CellType2","orig.ident")]
  hetero_data <- unique(hetero_data,by=c("orig.ident","CellType2","Risk_type"))
  
  summary_data <- hetero_data %>%
    group_by(CellType2) %>%
    summarise(
      mean_value = mean(risk_diversity_shannon),
      se_value = sd(risk_diversity_shannon) / sqrt(n())) %>%  
    arrange(mean_value)
  
  ggplot(summary_data, aes(x = reorder(CellType2, mean_value), y = mean_value,color=CellType2)) +
    geom_point(size = 2) +  # 绘制点
    geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.2) +  # 添加误差棒
    scale_color_manual(values = CellType2.colors)+
    ylab("Risk Type Diversity (Shannon–Weaver index)")+
    theme_bw() +
    theme(axis.text.x=element_text(size = 6,colour = 'black',angle =90, vjust = 0.5, hjust=1))
  
  pheno_data <- object_l@meta.data
  pheno_data$CellType2 <- factor(pheno_data$CellType2,levels = summary_data$CellType2)
  phenoPropotion(pheno_data,x_name = "CellType2",y_name="Risk_type",
                       cols=Risk_type.colors,
                       legend_name="Risk type",bar_width = 0.98)
  
  
  
# FigS2A
  
  markers_CellType2 <- list(
    Tcell=c("CD3D","CD3E","CD4","CD8A"),
    Naive=c("LEF1","TCF7","CCR7","SELL"),
    Mem=c("IL7R","CXCR4","ANXA1","LTB","CD69","CD44"),
    Tfh=c("CXCL13","TOX2","IL21","CXCR5"),
    Treg=c("FOXP3","IL2RA","CTLA4"),
    MAIT=c("KLRB1","SLC4A10","ZBTB16","ME1"),
    Tex=c("LAG3","TOX","PDCD1"))
  DotPlot(object_l[,object_l$CellType %in% c("CD4T","CD8T")], features = markers_CellType2,group.by = "CellType2",dot.scale=2.5)+
    scale_y_discrete(position="right")+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5,size=6,colour = 'black'))+
    scale_color_gradientn(values = seq(0,1,0.2),colours = c("royalblue3","cyan4","linen","#ffbb78","brown3","brown2"))+
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
  

  
# T cells
  
  object_t <- subset(object_l,subset= CellType %in% c("CD4T","CD8T","gdT"))
  object_t <- scCluster(object_t,nfeature=2000,min.d=0.5,res=1.5,TCRBCR=T,Ribo_MT=F,harmony=T,seed=2,k_param=50)
  object_t <- RunUMAP(object_t,reduction = "harmony",seed.use = 3,dims=1:30,
                      umap.method='uwot',min.dist=0.3,spread=1)
  
  
# FigS2B
  
  DimPlot(object_t,group.by="CellType2",cols=CellType2.colors,label=T,raster=T,pt.size=2,label.size = 3)+NoAxes()
  
  
  
# FigS2C
  
  DimPlot(object_t,label =F,cells.highlight = WhichCells(object_t,expression = clone_stats=="Clonal") ,pt.size = 2,sizes.highlight=2,raster=T,
          cols.highlight="seagreen3",cols="gray90")+NoAxes()
  
  
  clonal_data <- countClones(object_t@meta.data,clone="sample_cloneID")
  names(clonal_data) <- c("sample_cloneID","clone_size","clone_freq")
  clonal_data$clone_stats <- ifelse(clonal_data$clone_size > 1,"Clonal","NonClonal")
  clonal_data$clone_size_log2 <- log2(clonal_data$clone_size)
  clonal_data <- clonal_data %>%
    mutate(clone_rank = dense_rank(-clone_size))
  
  object_t@meta.data <- left_join(object_t@meta.data,clonal_data)
  row.names(object_t@meta.data) <- object_t$cellID
  
  FeaturePlot(object_t,features = 'clone_size_log2',pt.size = 2,raster = T)+
    scale_color_gradientn(colours = c("gray90",colorRampPalette(c("#313695","#abd9e9"))(5),
                                      colorRampPalette(c("#abd9e9","#fee090"))(10),
                                      colorRampPalette(c("#fee090","#a50026"))(20)),na.value="#ebeee8")+
    NoAxes()+ggtitle('Clone size in T cells (log2-transformed)')
  
  

# Fig3E
  
  ggplot(object_t@meta.data[object_t$CellType != "gdT",],aes(CellType2,clone_size_log2))+
    geom_boxplot(aes(fill=Risk_type),width=0.8,color='black',outlier.shape = NA,linewidth  = 0.235)+
    scale_fill_manual(values=Risk_type.colors)+
    stat_compare_means(aes(group = Risk_type,label = paste0(..p.format..)),method = "wilcox.test",hide.ns = F,size=1.5)+
    theme_classic()+ggtitle('')+ylab('Clone size in T cells (log2-transformed)')+xlab('')+
    theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,size = 6,colour = 'black'))
  
  
  object_t$CellType2_risk <- paste0(object_t$CellType2,"-",object_t$Risk_type)
    
  clone_index <- na.omit(data.table(object_t@meta.data[,c("sample_cloneID","clone_rank","CellType2_risk")]))
  clone_index[clone_rank == 1, clone_index:= "Top 1"]
  clone_index[clone_rank %in% 2:10, clone_index:= "Top 2~10"]
  clone_index[clone_rank %in% 11:50, clone_index:= "Top 11~50"]
  clone_index[clone_rank >= 51, clone_index:= "Top >51"]
  clone_index[clone_rank == 60, clone_index:= "NonClonal"]
  clone_index$clone_index <- factor(clone_index$clone_index,levels=c("Top 1","Top 2~10","Top 11~50","Top >51","NonClonal"))
  
  plot_data <- clone_index[,.(count=.N),by=.(CellType2_risk,clone_index)]
  plot_data <- plot_data[,.(ratio=count/sum(count),allcount=sum(count),count,clone_index),by=.(CellType2_risk)]
  plot_data$CellType2_risk <- str_remove_all(plot_data$CellType2_risk,"ow_risk|igh_risk")
  plot_data$CellType2_risk <- factor(plot_data$CellType2_risk,
                                     levels=c("CD4T_Naive-L","CD4T_Naive-H","CD4T_Mem-L","CD4T_Mem-H","CD4T_Tfh-L","CD4T_Tfh-H","CD4T_Treg-L","CD4T_Treg-H",
                                              "CD8T_Naive-L","CD8T_Naive-H","CD8T_Mem-L","CD8T_Mem-H","CD8T_MAIT-L","CD8T_MAIT-H","CD8T_Tex-L","CD8T_Tex-H"))
  
  ggplot(plot_data,aes(x=CellType2_risk,y=ratio,fill=clone_index))+
    geom_bar(stat = 'identity',position = 'fill',width = 0.8)+
    scale_y_continuous(labels = scales::percent)+
    scale_fill_manual(values = c('#D04848','#F3B95F','#FDE767','#6895D2','lightcyan2'))+
    theme_classic()+ggtitle('')+ylab('Cell proportion(%)')+xlab('')+labs(fill = "Clonal Indices")+
    theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size = 6,colour = 'black'))
  

  
# Fig3F
  
  curve <- estimateAbundance(object_t@meta.data[object_t$CellType != "gdT",],group = 'Risk_type',ci=0.95, nboot=100, clone="sample_cloneID")
  plot(curve, colors = Risk_type.colors, legend_title="")+
    scale_x_log10(breaks = c(1, 10, 100), labels = c("1", "10", "100"),limits=c(1,150))+
    ylab('Fraction of T cells')+
    xlab('Rank of top 150 clones')
  
  risk_diversity <- alphaDiversity(curve,group="Risk_type", min_q=0, max_q=2, step_q=1, nboot=200)
  plot(risk_diversity, 2, colors=Risk_type.colors, main_title='',legend_title="")+
    ylab('TCR clonotype diversity \n (1/Simpson\'s index,mean±SD)')+
    scale_x_discrete(limits=names(Risk_type.colors))
  
  
  
# Fig3G
  
  CellType2_order <- c("CD4T_Naive","CD4T_Mem","CD4T_Tfh","CD4T_Treg","CD8T_Naive","CD8T_Mem","CD8T_MAIT","CD8T_Tex")
  object_t_L <- object_t[,(object_t$CellType != "gdT" & object_t$Risk_type == "Low_risk")]
  object_t_L$CellType2 <- factor(object_t_L$CellType2,levels=CellType2_order)
  object_t_H <- object_t[,(object_t$CellType != "gdT" & object_t$Risk_type == "High_risk")]
  object_t_H$CellType2 <- factor(object_t_H$CellType2,levels=CellType2_order)
  
  in.dat_L <- data.frame(Cell_Name=object_t_L$cellID, clone.id=object_t_L$sample_cloneID, patient=object_t_L$orig.ident,
                         majorCluster=object_t_L$CellType2,loc=object_t_L$Type)
  out_L <- Startrac.run(in.dat_L, proj="Low_risk")
  in.dat_H <- data.frame(Cell_Name=object_t_H$cellID, clone.id=object_t_H$sample_cloneID, patient=object_t_H$orig.ident,
                         majorCluster=object_t_H$CellType2,loc=object_t_H$Type)
  out_H <- Startrac.run(in.dat_H, proj="High_risk")
  
  L_lower  <- out_L@pIndex.tran[1:8,3:ncol(out_L@pIndex.tran)]
  row.names(L_lower) <- names(L_lower)
  L_lower <- L_lower[CellType2_order,CellType2_order]
  L_lower[upper.tri(L_lower)] <- 0
  
  H_upper <- out_H@pIndex.tran[1:8,3:ncol(out_H@pIndex.tran)]
  row.names(H_upper) <- names(H_upper)
  H_upper <- H_upper[CellType2_order,CellType2_order]
  H_upper[lower.tri(H_upper)] <- 0
  
  tcr_trans <- L_lower+H_upper
  
  tcr_trans[tcr_trans > 0.1] <- 0.1
  colour_bk <- c(colorRampPalette(c("#3288bd","#abdda4"))(2),
                 colorRampPalette(c("#abdda4","#fee08b"))(20),
                 colorRampPalette(c("#fee08b","#fdae61"))(50),
                 colorRampPalette(c("#fdae61","#9e0142"))(50))
  pheatmap::pheatmap(tcr_trans,cluster_rows=F,cluster_cols=F,angle_col=45,fontsize = 10,
                     color=colour_bk,na_col="white",
                     width=5,height=4)

  

# Tregs
  
  object_treg <- subset(object_l,subset=CellType2 == "CD4T_Treg")
  object_treg <- scCluster(object_treg,harmony=T)
  
  clonal_data <- countClones(object_treg@meta.data,clone="sample_cloneID")
  names(clonal_data) <- c("sample_cloneID","clone_size","clone_freq")
  clonal_data$clone_stats <- ifelse(clonal_data$clone_size > 1,"Clonal","NonClonal")
  clonal_data$clone_size_log2 <- log2(clonal_data$clone_size)
  clonal_data <- clonal_data %>%
    mutate(clone_rank = dense_rank(-clone_size))
  
  object_treg@meta.data <- left_join(object_treg@meta.data,clonal_data)
  row.names(object_treg@meta.data) <- object_treg$cellID
  
  `TNFRSF9+ Treg` <- c("FOXP3", "BATF", "IKZF2", "IKZF4", "ZBTB32","IL32", "NAMPT", "EBI3", "CCL22", "LTB", "LTA", "CSF1", "CD70", "IL1RN", "IL7", 
                       "IL2RA", "CCR8", "IL1R2", "IL2RB", "CD74",
                       "LAYN", "TNFRSF18", "TNFRSF9", "TIGIT", "CTLA4", "TNFRSF4", "CARD16", "LAIR2", "TBC1D4", "PMAIP1")
  object_treg <- AddModuleScore(object_treg,features=list(`TNFRSF9+ Treg`),name="TNFRSF9+ Treg")
  
  
  
# FigS2F
  
  diff_treg <- FindMarkers(object_treg,ident.1 = "High_risk",group.by="Risk_type",logfc.threshold = 0)
  diff_treg <- diff_treg %>%
    mutate(pct.diff=pct.1-pct.2)%>%
    rownames_to_column("gene") %>%
    arrange(desc(avg_log2FC))
  diff_treg <- data.table(diff_treg)
  
  diff_treg$type <- "NS"
  diff_treg[pct.diff > 0 & avg_log2FC > 0.25 & p_val_adj < 0.05, type:="Up"]
  diff_treg[pct.diff < 0 &avg_log2FC < -0.25 & p_val_adj < 0.05, type:="Down"]
  table(diff_treg$type)
  
  treg_markers <- c("CCR8","TNFRSF4","TNFRSF9","TNFRSF18","PDCD1","FOXP3","LYAN","LAG3","HAVCR2","CTLA4","BATF",
                    "IL2RA", "CCR8", "IL1R2", "IL2RB","LTA")
  diff_treg[gene %in% treg_markers,label:=gene]
  
  ggplot(data=diff_treg, aes(x = pct.diff, y = avg_log2FC, color=type,label=label)) +
    geom_point_rast(stroke = 0, alpha = 1,shape = 16,size=0.5) + 
    scale_color_manual(values=c(NS="azure2",Up="#7AB774",Down="#7f91c9"))+
    geom_hline(yintercept=0, linetype="dashed", color = "lightgoldenrod")+
    geom_vline(xintercept=0, linetype="dashed", color = "lightgoldenrod")+
    ggrepel::geom_text_repel(aes(label = label),color="red",
                             size=1.5,segment.size=0.1,segment.color='red',max.overlaps=10000)+
    ylab("Log2-Fold change")+
    xlab("Percentage difference") +
    theme_classic()
  
  
  
# Fig3H & FigS2DE
  
  diff_risk <- FindMarkers(object_treg,ident.1 = "High_risk", group.by = "Risk_type",logfc.threshold = 0.25, min.pct = 0.1)
  # diffusion map
  sce <- SingleCellExperiment(assays = list(logcounts = LayerData(object_treg,"data")[row.names(diff_risk),]),
                              colData = object_treg@meta.data)
  dm <- DiffusionMap(sce,n_pcs = 30)
  dpt <- DPT(dm,tips=c(which(rownames(dm@eigenvectors) == row.names(dm@eigenvectors)[which(dm@eigenvectors[,"DC1"]== max(dm@eigenvectors[,"DC1"]))]),
                       which(rownames(dm@eigenvectors) == row.names(dm@eigenvectors)[which(dm@eigenvectors[,"DC2"]== min(dm@eigenvectors[,"DC2"]))]),
                       which(rownames(dm@eigenvectors) == row.names(dm@eigenvectors)[which(dm@eigenvectors[,"DC2"]== max(dm@eigenvectors[,"DC2"]))])))
  plot(dpt)
  
  diffusion_components <- as.data.frame(dm@eigenvectors)
  colnames(diffusion_components) <- paste0("DC", 1:ncol(diffusion_components))
  object_treg[['diffmap']] <- CreateDimReducObject(embeddings = as.matrix(diffusion_components), key = "DC_")
  object_treg$DC1 <- object_treg$DC2 <- object_treg$DPT <- NULL
  object_treg@meta.data <- cbind(object_treg@meta.data,diffusion_components[,c("DC1", "DC2")])
  object_treg$DPT <- dpt$dpt
  
  
  colour_bk <- c(colorRampPalette(c("purple3", "#4575B4"))(30),
                 colorRampPalette(c("#4575B4", "#91CF60"))(10),
                 colorRampPalette(c("#91CF60", "gold"))(5),
                 colorRampPalette(c("gold", "yellow"))(5))
  ggplot(object_treg@meta.data,aes(x=DC1,y=DC2,color=DPT))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk)+
    theme_classic()+NoAxes()
  
  colour_bk <- c(colorRampPalette(c("cyan3", "lightblue"))(35),  
                 colorRampPalette(c("lightblue", "linen"))(10),
                 colorRampPalette(c("linen", "lightsalmon2"))(8), 
                 colorRampPalette(c("lightsalmon2", "brown2"))(40))
  ggplot(object_treg@meta.data,aes(x=DC1,y=DC2,color=Risk_score))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk)+
    theme_classic()+NoAxes()
  
  
  colour_bk <- c(colorRampPalette(c("seagreen3","lightcyan"))(10),
                 colorRampPalette(c("lightcyan","linen"))(5),
                 colorRampPalette(c("linen","#ffbb78"))(5),
                 colorRampPalette(c("#ffbb78","red3"))(5))
  ggplot(object_treg@meta.data,aes(x=DC1,y=DC2,color=`TNFRSF9+ Treg1`))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk)+
    theme_classic()+NoAxes()
  
  
  object_treg$CCR8 <- LayerData(object_treg,"data")["CCR8",]
  object_treg$TNFRSF9 <- LayerData(object_treg,"data")["TNFRSF9",]
  colour_bk <- c("lightcyan",
                 colorRampPalette(c("#f7f7f7","#fdf2b5"))(10),
                 colorRampPalette(c("#fdf2b5","#F67B51"))(10),
                 colorRampPalette(c("#F67B51","#A30023"))(10))
  ggplot(object_treg@meta.data,aes(x=DC1,y=DC2,color=CCR8))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk)+
    theme_classic()+NoAxes()
  ggplot(object_treg@meta.data,aes(x=DC1,y=DC2,color=TNFRSF9))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk)+
    theme_classic()+NoAxes()
  
  
  
# Fig3J  
  
  colour_bk <- c(colorRampPalette(c("seagreen3","lightcyan"))(10),
                 colorRampPalette(c("lightcyan","linen"))(5),
                 colorRampPalette(c("linen","#ffbb78"))(5),
                 colorRampPalette(c("#ffbb78","red3"))(10))
  ggplot(object_treg@meta.data,aes(x=DPT,y=Risk_score,color=`TNFRSF9+ Treg1`))+
    geom_point_rast(aes(size=clone_size),alpha=0.6)+
    scale_size(range=c(0.1,1.5))+
    ylab("Risk score")+
    scale_colour_gradientn(colours = colour_bk)+
    geom_smooth(fullrange=T,method="loess",span=0.5,color="indianred",size=0.5)+
    facet_zoom(xlim = c(0.2, 0.56))+
    theme_classic()+
    theme(title=element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size=8),
          axis.text.x = element_text(size = 7,angle=0),
          axis.text.y = element_text(size = 6))
  
  
  
# Fig3I
  
  object_hp <- FindVariableFeatures(object_treg,nfeatures = 3000)
  hp_data <- LayerData(object_hp,"data")[VariableFeatures(object_hp),]
  hp_data <- hp_data[rowSums(hp_data != 0) > 10,]
  
  cell_info <- data.frame(cellID=object_hp$cellID,DPT=object_hp$DPT)
  cell_info <- cell_info[order(cell_info$DPT),]
  cell_info$group <- as.numeric(cut2(cell_info$DPT, g=100))
  hp_data <- hp_data[,cell_info$cellID]
  hp_data <- groupMeans(hp_data,groups = cell_info$group)
  
  gene_dpt_cor <- data.frame(gene=rownames(hp_data),
                             cor=cor(1:100,t(hp_data),use="na.or.complete")[1,])
  gene_dpt_cor <- gene_dpt_cor[order(gene_dpt_cor$cor),]
  hp_data <- hp_data[gene_dpt_cor$gene,]
  
  hp_data <- apply(hp_data,1,function(x){
    smth(x,window = 0.1,method = "gaussian")
  })
  hp_data <- t(hp_data)
  hp_data[,1:5] <- hp_data[,6]
  hp_data[,96:100] <- hp_data[,95]
  
  
  row_ha_left <- rowAnnotation(
    Type = anno_simple(case_when(rownames(hp_data) %in% favorable_genes ~ "Favorable",
                                 rownames(hp_data) %in% adverse_genes ~ "Adverse",
                                 .default = "nonAnno"), 
                       col = c("Favorable" = "cyan2", "Adverse" = "brown3","nonAnno"="white")))
  
  
  labels_highlight <- c("CCR8","TNFRSF4","TNFRSF9","TNFRSF18","PDCD1","FOXP3","LYAN","TIGIT","LAG3","HAVCR2","CTLA4","BATF","GITR","FAS",
                        "IL2RA", "CCR8", "IL1R2", "IL2RB","LTA","LAYN","PDCD1","CTLA4","FOXP3")                                 
  `TNFRSF9+ Treg` <- c("FOXP3", "BATF", "IKZF2", "IKZF4", "ZBTB32","IL32", "NAMPT", "EBI3", "CCL22", "LTB", "LTA", "CSF1", "CD70", "IL1RN", "IL7", 
                       "IL2RA", "CCR8", "IL1R2", "IL2RB", #"CD74",
                       "LAYN", "TNFRSF18", "TNFRSF9", "TIGIT", "CTLA4", "TNFRSF4", "CARD16", "LAIR2", "TBC1D4", "PMAIP1")
  highlight_idx <- which(rownames(hp_data) %in% c(labels_highlight,`TNFRSF9+ Treg`))
  row_ha_right <- rowAnnotation(Gene = anno_mark(at = highlight_idx, labels = rownames(hp_data)[highlight_idx],labels_gp = gpar(fontsize = 8)))
  
  bk <- c(seq(-1,-0.1,by=0.02),seq(0,1,by=0.02))
  colour_bk <- c(
    colorRampPalette(c("royalblue3", "darkseagreen2"))(38),  
    colorRampPalette(c("darkseagreen2", "linen"))(10),
    colorRampPalette(c("linen", "lightsalmon"))(10), 
    colorRampPalette(c("lightsalmon", "red3"))(39))
  
  hp_data <- t(scale(t(hp_data)))
  
  Heatmap(hp_data,
          name = "Expression",
          col = colorRamp2(bk, colour_bk), 
          show_row_names = F,
          left_annotation = row_ha_left,  
          right_annotation = row_ha_right,
          cluster_rows = F,
          cluster_columns = F)

  
  
# Exhausted T cells
  
  object_tex <- subset(object_l,subset=CellType2 == "CD8T_Tex")
  object_tex <- scCluster(object_tex,harmony=T)

  clonal_data <- countClones(object_tex@meta.data,clone="sample_cloneID")
  names(clonal_data) <- c("sample_cloneID","clone_size","clone_freq")
  clonal_data$clone_stats <- ifelse(clonal_data$clone_size > 1,"Clonal","NonClonal")
  clonal_data$clone_size_log2 <- log2(clonal_data$clone_size)
  clonal_data <- clonal_data %>%
    mutate(clone_rank = dense_rank(-clone_size))
  
  object_tex@meta.data <- left_join(object_tex@meta.data,clonal_data)
  row.names(object_tex@meta.data) <- object_tex$cellID
  
  t_markers <- data.table(read.xlsx("E:/PAAD_Project/Data/Preprocess/L/T_markers.xlsx"))
  t_markers <- as.list(t_markers)
  object_tex <- AddModuleScore(object_tex,features=t_markers,name=names(t_markers))
  
  
  
# FigS2J
  
  diff_tex <- FindMarkers(object_tex,ident.1 = "High_risk",group.by="Risk_type",logfc.threshold = 0)
  diff_tex <- diff_tex %>%
    mutate(pct.diff=pct.1-pct.2)%>%
    rownames_to_column("gene") %>%
    arrange(desc(avg_log2FC))
  diff_tex <- data.table(diff_tex)
  
  diff_tex$type <- "NS"
  diff_tex[pct.diff > 0 & avg_log2FC > 0.25 & p_val_adj < 0.05, type:="Up"]
  diff_tex[pct.diff < 0 & avg_log2FC < -0.25 & p_val_adj < 0.05, type:="Down"]
  table(diff_tex$type)
  
  down_label <- c("GZMK","CD160","SATB1","FCMR","IL7R","ATM","EIF4B",diff_tex[type=="Down" & gene %like% "^MT"]$gene)
  tex_markers <- c("GZMK","PRF1","TCF7","CD28","EOMES","CCL4",
                   "TNFRSF4","IL6R","IGFBP4",
                   "LAYN","CXCL13",
                   "PDCD1","HAVCR2","LAG3","CTLA4","BATF",
                   "ENTPD1","MYO1E",
                   "MYO7A","TOX","TIGIT","TNFRSF9")
  #diff_tex$label <- NA
  diff_tex[gene %in% c(down_label,tex_markers),label:=gene]
  
  
  ggplot(data=diff_tex, aes(x = pct.diff, y = avg_log2FC, color=type,label=label)) +
    geom_point_rast(stroke = 0, alpha = 1,shape = 16,size=0.5) + 
    scale_color_manual(values=c(NS="azure2",Up="#7AB774",Down="#7f91c9"))+
    geom_hline(yintercept=0, linetype="dashed", color = "lightgoldenrod")+
    geom_vline(xintercept=0, linetype="dashed", color = "lightgoldenrod")+
    ggrepel::geom_text_repel(aes(label = label),color="red",
                             size=1.5,segment.size=0.1,segment.color='red',max.overlaps=10000)+
    xlim(-0.2,0.3)+ylim(-3,3)+
    ylab("Log2-Fold change")+
    xlab("Percentage difference") +
    theme_classic()
  
  
  
# Fig3KL & FigS2GHI
  
  harmony_embeddings <- Embeddings(object_tex, reduction = "harmony")[,1:30]
  set.seed(1)
  dm <- DiffusionMap(harmony_embeddings)
  dpt <- DPT(dm,tips=c(
    which(rownames(dm@eigenvectors) == row.names(dm@eigenvectors)[which(dm@eigenvectors[,"DC2"]== max(dm@eigenvectors[,"DC2"]))]),
    which(rownames(dm@eigenvectors) == row.names(dm@eigenvectors)[which(dm@eigenvectors[,"DC2"]== min(dm@eigenvectors[,"DC2"]))]),
    which(rownames(dm@eigenvectors) == row.names(dm@eigenvectors)[which(dm@eigenvectors[,"DC1"]== min(dm@eigenvectors[,"DC1"]))])))
  plot(dpt)
  
  diffusion_components <- as.data.frame(dm@eigenvectors)
  colnames(diffusion_components) <- paste0("DC", 1:ncol(diffusion_components))
  object_tex[['diffmap']] <- CreateDimReducObject(embeddings = as.matrix(diffusion_components), key = "DC_")
  object_tex$DC1 <- object_tex$DC2 <- object_tex$DPT <- NULL
  object_tex@meta.data <- cbind(object_tex@meta.data,diffusion_components[,c("DC1", "DC2")])
  object_tex$DPT <- dpt$dpt
  
  
  colour_bk <- c(colorRampPalette(c("purple3", "#4575B4"))(20),
                 colorRampPalette(c("#4575B4", "#91CF60"))(10),
                 colorRampPalette(c("#91CF60", "gold"))(10),
                 colorRampPalette(c("gold", "yellow"))(10))
  ggplot(object_tex@meta.data,aes(x=DC1,y=DC2,color=DPT))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk)+
    theme_classic()+NoAxes()
  
  colour_bk <- c(colorRampPalette(c("cyan3", "lightblue"))(35),  
                 colorRampPalette(c("lightblue", "linen"))(20),
                 colorRampPalette(c("linen", "lightsalmon"))(3), 
                 colorRampPalette(c("lightsalmon", "brown2"))(42))
  ggplot(object_tex@meta.data,aes(x=DC1,y=DC2,color=Risk_score))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk)+
    theme_classic()+NoAxes()
  
  object_tex$top_clone <- NA
  object_tex$top_clone[object_tex$sample_cloneID == c("P3_clonotype1")] <- "P3_clonotype1"
  object_tex$top_clone[object_tex$sample_cloneID == c("M1_clonotype1")] <- "M1_clonotype1"
  DimPlot(object_tex,reduction = "diffmap",label =F,group.by = "top_clone",cols=c("indianred1","deepskyblue3"),
          pt.size=3,raster = T,na.value="gray95",order=T)+
    labs(col="Top15 clones")+ggtitle("")+
    theme(legend.title=element_text(size = 8),
          legend.text=element_text(size = 8),
          legend.key.height = unit(0.1,"cm"),
          legend.key.width = unit(0.1,"cm"))+NoAxes()
  
  
  colour_bk <- c(colorRampPalette(c("seagreen3","lightcyan"))(10),
                 colorRampPalette(c("lightcyan","linen"))(5),
                 colorRampPalette(c("linen","#ffbb78"))(10),
                 colorRampPalette(c("#ffbb78","red3"))(10))
  for(i in c("CD8.Tumor-reactivity3","Tumor.Specific6","Virus.Specific8","Influenza.TIL9")){
    p <- ggplot(object_tex@meta.data,aes(x=DC1,y=DC2,color=get(i)))+
      geom_point_rast(size=0.1)+
      scale_colour_gradientn(colours = colour_bk)+
      theme_classic()+NoAxes()+
      ggtitle(i)+theme(title = element_text(size=8))
    print(p)
  }
  
  
  metaData <- object_tex@meta.data
  for (i in unique(metaData$Type)) {
    x <- metaData[metaData$Type == i,]
    dim(x)[1]
    if (dim(x)[1]>1000 ) {
      p <- ggplot(data = x,aes(x = DC1, y = DC2)) +
        stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F)+
        geom_point(data = x[sample(seq_along(x$cellID),1000,replace = F),],
                   aes(x = DC1, y = DC2),color = 'white',size = .005)+
        scale_fill_viridis(option="A")+theme_classic()+NoAxes()+NoLegend()+ggtitle(i)
    }else{
      p <- ggplot(data = x,aes(x = DC1, y = DC2)) +
        stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F)+
        geom_point(data = x,
                   aes(x = DC1, y = DC2),color = 'white',size = .005)+
        scale_fill_viridis(option="A")+theme_classic()+NoAxes()+NoLegend()+ggtitle(i)
    }
    print(p)
  }
  
  
  
# Fig3O
  
  ggplot(object_tex@meta.data, aes(x=DPT,y=as.numeric(Response=="SD")))+
    geom_smooth(aes(x=DPT,y=as.numeric(Response=="SD")),fullrange=T,method="loess",span=1,color="red3",fill="gray90")+
    ylab("%SD (Response)") +xlim(0.15,0.75)+ theme_classic()
  
  
  
# Fig3M
  
  tex_fetures <- c("CD8+.exhaustion2","CD8.Tumor-reactivity3","Tumor.Specific6",
                   "Virus.Specific8","Influenza.TIL9")
  VlnPlot(object_tex,features = tex_fetures,group.by = "Type",stack = T,flip = T,cols=Type.colors,fill.by="ident",adjust=1.5)+
    ylab("")+xlab("")+ggtitle("Activity Scores")+
    theme(title=element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size=8),
          axis.text.x = element_text(size = 7,angle=0),
          axis.text.y = element_text(size = 6))+NoLegend()+
    stat_signif(comparisons = list(c("P","M")),textsize=2.5,y_position=0.5,map_signif_level=T,tip_length=0)
  
  
  
# Fig3N
  
  object_hp <- FindVariableFeatures(object_tex,nfeatures = 3000)
  hp_data <- LayerData(object_hp,"data")[VariableFeatures(object_hp),]
  hp_data <- hp_data[rowSums(hp_data != 0) > 10,]

  cell_info <- data.frame(cellID=object_hp$cellID,DPT=object_hp$DPT)
  cell_info <- cell_info[order(cell_info$DPT),]
  cell_info$group <- as.numeric(cut2(cell_info$DPT, g=100))
  hp_data <- hp_data[,cell_info$cellID]
  hp_data <- groupMeans(hp_data,groups = cell_info$group)
  
  gene_dpt_cor <- data.frame(gene=rownames(hp_data),
                             cor=cor(1:100,t(hp_data),use="na.or.complete")[1,])
  gene_dpt_cor <- gene_dpt_cor[order(gene_dpt_cor$cor),]
  hp_data <- hp_data[gene_dpt_cor$gene,]
  
  hp_data <- apply(hp_data,1,function(x){
    smth(x,window = 0.1,method = "gaussian")
  })
  hp_data <- t(hp_data)
  hp_data[,1:5] <- hp_data[,6]
  hp_data[,96:100] <- hp_data[,95]
  

  row_ha_left <- rowAnnotation(
    Type = anno_simple(case_when(rownames(hp_data) %in% favorable_genes ~ "Favorable",
                                 rownames(hp_data) %in% adverse_genes ~ "Adverse",
                                 .default = "nonAnno"), 
                       col = c("Favorable" = "cyan2", "Adverse" = "brown3","nonAnno"="white")))
  labels_highlight <- c("GZMK","GZMA","TCF7","CD28","EOMES","CCR5","CCL5","CD40",
                        "BTLA","TNFRSF4","SELL","CCR7","IL6R","IGFBP4",
                        "LAYN","CXCL13",
                        "PDCD1","HAVCR2","LAG3","CTLA4","BATF",
                        "ENTPD1","MYO1E",
                        "MYO7A","TOX","TIGIT","TNFRSF9")                                     
  highlight_idx <- which(rownames(hp_data) %in% labels_highlight)
  row_ha_right <- rowAnnotation(
    Gene = anno_mark(at = highlight_idx, labels = rownames(hp_data)[highlight_idx],
                     labels_gp = gpar(fontsize = 8)))
  bk <- c(seq(-1,-0.1,by=0.02),seq(0,1,by=0.02))
  colour_bk <- c(
    colorRampPalette(c("royalblue3", "darkseagreen2"))(38),  
    colorRampPalette(c("darkseagreen2", "linen"))(10),
    colorRampPalette(c("linen", "lightsalmon"))(10), 
    colorRampPalette(c("lightsalmon", "red3"))(39))
  
  hp_data <- t(scale(t(hp_data)))
  
  Heatmap(hp_data,
          name = "Expression",
          col = colorRamp2(bk, colour_bk),  # 颜色映射
          show_row_names = F,
          left_annotation = row_ha_left,  
          right_annotation = row_ha_right,
          cluster_rows = F,
          cluster_columns = F)
  
  
  
# Fig3P
  
  colour_bk <- c(colorRampPalette(c("seagreen3","lightcyan"))(10),
                 colorRampPalette(c("lightcyan","linen"))(5),
                 colorRampPalette(c("linen","#ffbb78"))(5),
                 colorRampPalette(c("#ffbb78","red3"))(10))
  ggplot(object_tex@meta.data,aes(x=DPT,y=Risk_score,color=`CD8+.exhaustion2`))+
    geom_point_rast(aes(size=clone_size),alpha=0.6)+
    scale_size(range=c(0.1,1.5))+
    ylab("Risk score")+
    scale_colour_gradientn(colours = colour_bk)+
    geom_smooth(fullrange=T,method="loess",span=0.5,color="indianred",size=0.5)+
    facet_zoom(xlim = c(0.3, 0.6))+
    theme_classic()
  
  
  
# NK cells
  
  object_nk <- subset(object_l,subset=CellType == "NK")
  object_nk <- scCluster(object_nk,harmony = T,min.d=0.6,res=0.5,TCRBCR=F,Ribo_MT=T,seed=2,k_param=50,npcs=50)
  object_nk <- RunUMAP(object_nk,reduction = "harmony",seed.use = 3,dims=1:30,
                       umap.method='uwot',min.dist=0.5,spread=1)

  GeneCounts <- as.matrix(object_nk@assays$RNA@counts)
  iOrd <- rowSums(GeneCounts>0)
  GeneCounts <- GeneCounts[iOrd>10,]#only keep genes expressing in more than 10 cell
  CytoTRACE.res <- CytoTRACE(GeneCounts, batch = object_nk$orig.ident)
  object_nk$CytoTRACE <- CytoTRACE.res$CytoTRACE[colnames(object_nk)]
  
  scent <- object_nk
  scent <- NormalizeData(object_nk,normalization.method="RC")
  data_scent <- scent@assays$RNA@data
  symbol2id <- bitr(row.names(data_scent),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  data_scent <- data_scent[symbol2id$SYMBOL,]
  row.names(data_scent) <- symbol2id$ENTREZID
  data(net17Jan16)
  data_scent <- log2(data_scent+1)
  ccat.v <- CompCCAT(exp = as.matrix(data_scent), ppiA = net17Jan16.m)
  object_nk$SCENT <- ccat.v
  
  nk_markers <- list(
    NK.activation=c("FCGR3A","KLRD1","KLRC2","KLRK1","NCR3","NCR1","CD226","CD244","LAMP1","IFNG","GZMB","PRF1"),
    NK.exhaustion=c("KLRC1","KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR3DL1","TIGIT","PDCD1","TIM3","HAVCR2","SIGLEC7","SIGLEC9"))
  object_nk <- AddModuleScore(object_nk,features=nk_markers,name=names(nk_markers))
  
  
  
# Fig3M
  
  diff_nk <- FindMarkers(object_nk,ident.1 = "High_risk",group.by="Risk_type",logfc.threshold = 0)
  diff_nk <- diff_nk %>%
    mutate(pct.diff=pct.1-pct.2)%>%
    rownames_to_column("gene") %>%
    arrange(desc(avg_log2FC))
  diff_nk <- data.table(diff_nk)
  
  diff_nk$type <- "NS"
  diff_nk[pct.diff > 0 & avg_log2FC > 0.25 & p_val_adj < 0.05, type:="Up"]
  diff_nk[pct.diff < 0 &avg_log2FC < -0.25 & p_val_adj < 0.05, type:="Down"]
  table(diff_nk$type)
  
  down_genes <- c("XCL2","XCL1","NKG7","GZMH","GZMM","KLF2","KLF3","ICAM2","KLF3","CX3CR1")
  NK.activation=c("FCGR3A","KLRD1","KLRK1","NCR3","NCR1","CD226","IFNG")
  NK.exhaustion=c("KLRC1", "KIR2DL2", "KIR3DL1","PDCD1","TIM3","HAVCR2")
  nk_markers <- c("NCAM1","FCGR3A","KLRD1", 
                  "KLRK1","NCR3","NCR1",
                  "KLRC1","PDCD1","TIM3","HAVCR2","CTLA4",
                  "KIR2DL2", "KIR3DL1",
                  "CAPG","LGALS3","NR4A1","LGALS1")  
  
  diff_nk$label <- ""
  diff_nk[gene %in% c(NK.activation,NK.exhaustion,nk_markers,down_genes),label:=gene]
  
  ggplot(data=diff_nk, aes(x = pct.diff, y = avg_log2FC, color=type,label=label)) +
    geom_point_rast(stroke = 0, alpha = 1,shape = 16,size=0.5) + 
    scale_color_manual(values=c(NS="azure2",Up="#7AB774",Down="#7f91c9"))+
    geom_hline(yintercept=0, linetype="dashed", color = "lightgoldenrod")+
    geom_vline(xintercept=0, linetype="dashed", color = "lightgoldenrod")+
    ggrepel::geom_text_repel(aes(label = label),color="red",
                             size=1.5,segment.size=0.1,segment.color='red',max.overlaps=10000,min.segment.length = 0)+
    ylab("Log2-Fold change")+xlab("Percentage difference")+
    theme_classic()
  
  
  
# Fig3Q & FigS2KLN
  
  harmony_embeddings <- Embeddings(object_nk, reduction = "harmony")[,1:40]
  set.seed(2)
  dm <- DiffusionMap(harmony_embeddings)
  dpt <- DPT(dm,tips=c(which(rownames(dm@eigenvectors) == row.names(dm@eigenvectors)[which(dm@eigenvectors[,"DC1"]== min(dm@eigenvectors[,"DC1"]))]),
                       which(rownames(dm@eigenvectors) == row.names(dm@eigenvectors)[which(dm@eigenvectors[,"DC2"]== max(dm@eigenvectors[,"DC2"]))]),
                       which(rownames(dm@eigenvectors) == row.names(dm@eigenvectors)[which(dm@eigenvectors[,"DC1"]== max(dm@eigenvectors[,"DC1"]))])))
  plot(dpt)
  
  object_nk$DC1 <- object_nk$DC2 <- object_nk$DPT <- NULL
  diffusion_components <- as.data.frame(dm@eigenvectors)
  colnames(diffusion_components) <- paste0("DC", 1:ncol(diffusion_components))
  object_nk[['diffmap']] <- CreateDimReducObject(embeddings = as.matrix(diffusion_components), key = "DC_")
  object_nk@meta.data <- cbind(object_nk@meta.data,diffusion_components[,c("DC1", "DC2")])
  object_nk$DPT <- dpt$dpt
  
  
  colour_bk <- c(colorRampPalette(c("purple3", "#4575B4"))(30),
                 colorRampPalette(c("#4575B4", "#91CF60"))(5),
                 colorRampPalette(c("#91CF60", "gold"))(5),
                 colorRampPalette(c("gold", "yellow"))(12))
  ggplot(object_nk@meta.data,aes(x=DC1,y=DC2,color=DPT))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk)+
    theme_classic()+NoAxes()
  
  
  ggplot(object_nk@meta.data, aes(x=DPT,y=nFeature_RNA))+
    geom_smooth(fullrange=T,method="loess",span=1.5,color="purple3",fill="gray90")+
    ylab("Expressed genes (#)") +
    theme_classic()
  ggplot(object_nk@meta.data, aes(x=DPT,y=SCENT))+
    geom_smooth(fullrange=T,method="loess",span=1.5,color="purple3",fill="gray90")+
    ylab("SCENT") +
    theme_classic()
  ggplot(object_nk@meta.data, aes(x=DPT,y=CytoTRACE))+
    geom_smooth(fullrange=T,method="loess",span=2,color="purple3",fill="gray90")+
    ylab("CytoTRACE") +
    theme_classic()
  
  colour_bk <- c(colorRampPalette(c("cyan3", "#9EDAE5"))(10),
                 colorRampPalette(c("#9EDAE5", "linen"))(10),
                 colorRampPalette(c("linen", "#ff9896"))(10),
                 colorRampPalette(c("#ff9896", "brown3"))(12))
  ggplot(object_nk@meta.data,aes(x=DC1,y=DC2,color=Risk_score))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk)+
    theme_classic()+NoAxes()
  
  object_nk$FCGR3A <- LayerData(object_nk,"data")["FCGR3A",]
  object_nk$LGALS3 <- LayerData(object_nk,"data")["LGALS3",]
  colour_bk <- c("lightcyan",
                 colorRampPalette(c("#f7f7f7","#fdf2b5"))(5),
                 colorRampPalette(c("#fdf2b5","#F67B51"))(5),
                 colorRampPalette(c("#F67B51","#A30023"))(5))
  ggplot(object_nk@meta.data,aes(x=DC1,y=DC2,color=FCGR3A))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk,na.value = "lightcyan")+
    theme_classic()+NoAxes()
  ggplot(object_nk@meta.data,aes(x=DC1,y=DC2,color=LGALS3))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk,na.value = "lightcyan")+
    theme_classic()+NoAxes()
  
  colour_bk <- c(colorRampPalette(c("seagreen3","lightcyan"))(10),
                 colorRampPalette(c("lightcyan","linen"))(5),
                 colorRampPalette(c("linen","#ffbb78"))(5),
                 colorRampPalette(c("#ffbb78","red3"))(5))
  ggplot(object_nk@meta.data,aes(x=DC1,y=DC2,color=NK.activation1))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk)+
    theme_classic()+NoAxes()
  

# Fig3R
  
  object_nk$KLRF1 <- LayerData(object_nk,"data")["KLRF1",]
  object_nk$KLRC1 <- LayerData(object_nk,"data")["KLRC1",]
  
  plot_data <- data.table(object_nk@meta.data)
  
  ggplot(plot_data[Risk_type == "Low_risk"],aes(x=FCGR3A,y=KLRF1))+
    geom_point_rast(color=alpha("cyan3",0.6),size=0.2)+
    geom_hline(yintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    geom_vline(xintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    theme_classic()+xlim(0,4)+ylim(0,3.55)+
    ggtitle('Low risk')+xlab('FCGR3A+')+ylab('KLRF1+')
  ggplot(plot_data[Risk_type == "High_risk"],aes(x=FCGR3A,y=KLRF1))+
    geom_point_rast(color=alpha("brown3",0.6),size=0.2)+
    geom_hline(yintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    geom_vline(xintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    theme_classic()+xlim(0,4)+ylim(0,5)+
    ggtitle('High risk')+xlab('FCGR3A+')+ylab('KLRF1+')
  ggplot(plot_data[Risk_type == "Low_risk"],aes(x=KLRC1,y=LGALS3))+
    geom_point_rast(color=alpha("cyan3",0.6),size=0.2)+
    geom_hline(yintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    geom_vline(xintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    theme_classic()+xlim(0,4)+ylim(0,3.7)+
    ggtitle('Low risk')+xlab('KLRC1+')+ylab('LGALS3+')
  ggplot(plot_data[Risk_type == "High_risk"],aes(x=KLRC1,y=LGALS3))+
    geom_point_rast(color=alpha("brown3",0.6),size=0.2)+
    geom_hline(yintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    geom_vline(xintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    theme_classic()+xlim(0,4)+ylim(0,3.7)+
    ggtitle('High risk')+xlab('KLRC1+')+ylab('LGALS3+')
  
  
# Fig3S
  
  ggplot(object_nk@meta.data,aes(x=Risk_type,y=NK.activation1,fill=Risk_type))+
    geom_boxplot(color='gray40',outlier.shape=NA,size=0.1,width=0.4,linewidth  = 0.235)+
    geom_jitter(fill='white',color='black',shape=21,width =0.1,size=0.4,stroke = 0.1)+
    stat_compare_means(method="wilcox.test",size=2,bracket.size= 0.235)+
    scale_fill_manual(values=Risk_type.colors)+
    ylab("NK activation score")+xlab("")+
    theme_classic()
  
  
  
  
  
  
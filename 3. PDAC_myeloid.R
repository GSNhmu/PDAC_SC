
# Myeloid cells

  if(F){
    object_m <- subset(object,subse=CellType_Main == "Myeloid cell")
    object_m <- scCluster(object_m,nfeature=2500,min.d=0.3,res=1,TCRBCR=F,Ribo_MT=T,harmony=T,seed=1,k_param=30)
    object_m <- subset(object_m,subset=seurat_clusters == 10,invert=T)
    object_m <- RunTSNE(object_m,reduction = "harmony",dims=1:30,perplexity = 10,seed.use = 3)
    
    object_m$CellType2[object_m$seurat_clusters == 12 & object_m$CellType2 == "cDC"] <- "cDC1"
    object_m$CellType2[object_m$CellType2 == "cDC"] <- "cDC2"
    object_m$CellType2[object_m$seurat_clusters == 13 & object_m$CellType2 == "cDC"] <- "cDC3"
    object_m$CellType2[object_m$seurat_clusters == 13 & object_m$CellType2 == "cDC"] <- "cDC3"
  }
  
  
  
  
# Fig4A
  
  CellType2.colors <- c(Monocyte="#ff9896",Macrophage="gold",cDC1="seagreen2",cDC2="royalblue1",cDC3="cyan4",
                        Neutrophil="#b83570",Mast="#e377c2",pDC="#BBA0BD")
  DimPlot(object_m,reduction = "tsne",group.by="CellType2",cols=CellType2.colors,label=T,raster=T,pt.size=3,label.size = 1)
  
  
  
# FigS3A
  
  cDC_markers <- list(
    cDC1=c("CLEC9A","XCR1"),
    cDC2=c("CD1C","CLEC10A"),
    LAMP3_cDC3=c("LAMP3","CCR7","CD80"))
  DotPlot(object_m[,object_m$CellType == "cDC"], features = cDC_markers,group.by = "CellType2",dot.scale=2.5)+
    scale_y_discrete(position="right")+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5,size=6,colour = 'black'))+
    scale_color_gradientn(values = seq(0,1,0.2),colours = c("royalblue3","cyan4","linen","#ffbb78","brown3","brown2"))+
    labs(x=NULL,y=NULL)+
    guides(size=guide_legend(order=3))
  
  
# Fig4B
  
  colour_bk <- c(colorRampPalette(c("cyan3", "lightblue"))(40),  
                 colorRampPalette(c("lightblue", "linen"))(10),
                 colorRampPalette(c("linen", "lightsalmon2"))(10), 
                 colorRampPalette(c("lightsalmon2", "brown2"))(10),
                 colorRampPalette(c("brown2", "red"))(25))
  FeaturePlot(object_m,"Risk_score",raster=T,pt.size=3,reduction = "tsne")+NoAxes()+
    scale_colour_gradientn(colours = colour_bk)
  
# Fig4C
  
  object_m$DEGs <- object_m$cell_diff_num <- NULL
  
  CellType2_select <- rowSums(table(object_m$CellType2,object_m$Risk_type) > 50)==2
  CellType2_select <- names(CellType2_select[which(CellType2_select)])
  
  diff_CellType2 <- lapply(CellType2_select,function(c){
    print(c)
    diff_c <- FindMarkers(object_m[,object_m$CellType2==c],ident.1="High_risk",group.by="Risk_type",logfc.threshold = 0.25,min.pct = 0.1)
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
    mat <- LayerData(object_m,"counts")[diff_CellType2[CellType2==c,gene],object_m$cellID[object_m$CellType2 == c]]
    index <- colSums(mat>0) 
    cellID_DEGs <- data.table(cellID=names(index),DEGs=index)
    return(cellID_DEGs)
  })
  cellID_DEGs <- rbindlist(cellID_DEGs)
  
  object_m@meta.data <- left_join(object_m@meta.data,cellID_DEGs)
  row.names(object_m@meta.data) <- object_m$cellID
  
  FeaturePlot(object_m,features = "DEGs",pt.size = 3,raster = T,reduction = "tsne")+
    ggtitle("Number of DEGs (High_risk vs. Low_risk)")+
    scale_color_gradientn(colours = c(colorRampPalette(c("linen","#fdf2b5"))(5),
                                      colorRampPalette(c("#fdf2b5","#F67B51"))(10),
                                      colorRampPalette(c("#F67B51","red"))(15)),
                          na.value = "grey90",)+NoAxes()
  
  
# Fig4D
  
  hetero_data <- data.table(orig.ident=object_m$orig.ident,
                            CellType2=object_m$CellType2,
                            Risk_type=object_m$Risk_type)
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
    theme_bw() 
  
  pheno_data <- object_m@meta.data
  pheno_data$CellType2 <- factor(pheno_data$CellType2,levels = summary_data$CellType2)
  phenoPropotion(pheno_data,x_name = "CellType2",y_name="Risk_type",
                       cols=Risk_type.colors,
                       legend_name="Risk type",bar_width = 0.9)

  
# Fig4E 

  diff_mono <- diff_calculate(obj=object_m[,object_m$CellType == "Monocyte"])
  diff_mac <- diff_calculate(obj=object_m[,object_m$CellType == "Macrophage"])
  diff_dc2 <- diff_calculate(obj=object_m[,object_m$CellType2 == "cDC2"])
  diff_neu <- diff_calculate(obj=object_m[,object_m$CellType == "Neutrophil"])
  
  #Reduce(intersect,list(diff_mono[type=="Up",gene],diff_mac[type=="Up",gene],diff_dc2[type=="Up",gene],diff_neu[type=="Up",gene]))
  
  Hset <- msigdbr(species = "Homo sapiens",category="H") %>%
    dplyr::select(gs_name, gene_symbol)
  
  enrich_mono <- enricher(diff_mono$gene[diff_mono$type == "Up"],TERM2GENE=Hset)
  enrich_mono <- data.table(enrich_mono@result[enrich_mono@result$ID %in% enrich_mono$ID,],`Cell type`="Monocyte")
  enrich_mac <- enricher(diff_mac$gene[diff_mac$type == "Up"],TERM2GENE=Hset)
  enrich_mac <- data.table(enrich_mac@result[enrich_mac@result$ID %in% enrich_mac$ID,],`Cell type`="Macrophage")
  enrich_dc2 <- enricher(diff_dc2$gene[diff_dc2$type == "Up"],TERM2GENE=Hset)
  enrich_dc2 <- data.table(enrich_dc2@result[enrich_dc2@result$ID %in% enrich_dc2$ID,],`Cell type`="cDC2")
  enrich_neu <- enricher(diff_neu$gene[diff_neu$type == "Up"],TERM2GENE=Hset)
  enrich_neu <- data.table(enrich_neu@result[enrich_neu@result$ID %in% enrich_neu$ID,],`Cell type`="Neutrophil")
  
  common_h <- Reduce(intersect,list(enrich_mono$ID,enrich_mac$ID,enrich_dc2$ID,enrich_neu$ID))
  
  plot_data <- rbind(enrich_mono[ID %in% common_h],enrich_mac[ID %in% common_h],enrich_dc2[ID %in% common_h],enrich_neu[ID %in% common_h])
  plot_data <- dplyr::mutate(plot_data, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
  plot_data$ID <- factor(plot_data$ID,
                         levels=c("HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_HYPOXIA","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_APOPTOSIS"))
  plot_data$`Cell type` <- factor(plot_data$`Cell type`,
                                  levels=rev(c("Monocyte","Macrophage","cDC2","Neutrophil")))
  
  ggplot(plot_data)+
    geom_segment(aes(x=`Cell type`,xend=`Cell type`,y=0,yend=richFactor),color="gray50")+
    geom_point(aes(x=`Cell type`,y=richFactor,color=`Cell type`,size=Count))+
    scale_color_manual(values = c(Monocyte="#ff9896",Macrophage="gold",cDC2="royalblue1",Neutrophil="#b83570"))+
    theme_minimal()+
    theme(panel.border=element_blank(),
          panel.spacing=unit(0.1,"lines"),
          strip.text.x=element_text(size=8))+
    coord_flip()+
    xlab("")+
    ylab("richFactor")+
    facet_wrap(~ID,ncol=1,scale="free_y")
  
  
  
# Fig4F
  
  object_m_subset <- subset(object_m,subset=CellType2 %in% c("Monocyte","Macrophage","cDC2","Neutrophil"))
  object_m_subset$Risk_type <- factor(object_m_subset$Risk_type,levels=c("Low_risk","High_risk"))
  object_m_subset$CellType2 <- factor(object_m_subset$CellType2,levels=c("Monocyte","Macrophage","cDC2","Neutrophil"))
  object_m_subset$VEGFA <- LayerData(object_m_subset,"data")["VEGFA",]
  object_m_subset$HIF1A <- LayerData(object_m_subset,"data")["HIF1A",]
  
  ggplot(object_m_subset@meta.data,aes(CellType2,VEGFA,fill=Risk_type))+
    geom_boxplot(width=0.8,color='black',linewidth  = 0.235,outlier.shape=NA)+
    geom_point(position = position_jitterdodge(),color='black',size=0.3,shape=21,alpha=0.6)+
    scale_fill_manual(values=Risk_type.colors)+
    stat_compare_means(aes(group = Risk_type,label = paste0(..p.format..)),method = "wilcox.test",hide.ns = F,size=1.5)+
    theme_classic()+ggtitle("")+ylab("VEGFA")+xlab("")+
    theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,size = 6,colour = 'black'))
  
  ggplot(object_m_subset@meta.data,aes(CellType2,HIF1A,fill=Risk_type))+
    geom_boxplot(width=0.8,color='black',linewidth  = 0.235,outlier.shape=NA)+
    geom_point(position = position_jitterdodge(),color='black',size=0.3,shape=21,alpha=0.6)+
    scale_fill_manual(values=Risk_type.colors)+
    stat_compare_means(aes(group = Risk_type,label = paste0(..p.format..)),method = "wilcox.test",hide.ns = F,size=1.5)+
    theme_classic()+ggtitle("")+ylab("HIF1A")+xlab("")+
    theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,size = 6,colour = 'black'))
  
  
  
# Monocytes and macrophages
  
  object_mm <- subset(object_m,subset=CellType %in% c("Monocyte","Macrophage"))
  object_mm <- scCluster(object_mm,harmony=T,res = 0.5)
  
  mm_markers <- list(
    Phagocytosis = c("MRC1","CD163","MERTK","C1QB"),
    Antigen_presenting = c("CD74","TAP1","TAP2","TAPBP","CD36","HLA-DRA",
                           "HLA-DRB1","HLA-DPB1","HLA-DOA","HLA-DOB","HLA-DRB5",
                           "HLA-DRB6","HLA-A","HLA-B","HLA-C"),
    Angiogenesis = c("CCND2","CCNE1","CD44","CXCR4","E2F3","EDN1","EZH2",
                     "FGF18","FGFR1","FYN","HEY1","ITGAV","JAG1","JAG2",
                     "MMP9","NOTCH1","PDGFA","PTK2","SPP1","STC1","TNFAIP6",
                     "TYMP","VAV2","VCAN","VEGFA"),
    M1Score = c("IL23","TNF","CXCL9","CXCL10","CXCL11","CD86","IL1A",
                "IL1B","IL6","CCL5","IRF5","IRF1","CD40","IDO1","KYNU","CCR7"),
    M2Score = c("IL4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24",
                "LYVE1","VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC",
                "CTSD","TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A","WNT7B",
                "FASL","TNFSF12","TNFSF8","CD276","VTCN1","MSR1","FN1","IRF4"),
    LAM = c("CD163","SPP1","C1QC","FABP3","FABP4","FABP5","LPL","LIPA","LGALS3","TREM2"))
  object_mm <- AddModuleScore(object_mm,mm_markers,name = names(mm_markers))
  m1m2_data <- subtype_calculate(expr=LayerData(object_mm,"data"),markers_2=mm_markers$M2Score,markers_1=mm_markers$M1Score,
                                 name_subtype="M1M2",name_2="M2",name_1="M1")
  m1m2_data$samples <- NULL
  object_mm@meta.data <- cbind(object_mm@meta.data,m1m2_data)
  
  
  
  object_dc <- FindVariableFeatures(object_mm,nfeatures = 1000)
  sce <- SingleCellExperiment(assays = list(logcounts = LayerData(object_dc,"data")[VariableFeatures(object_dc),]),
                              colData = object_dc@meta.data)
  dm <- DiffusionMap(sce)
  dpt <- DPT(dm,tips=c(which(rownames(dm@eigenvectors) == row.names(dm@eigenvectors)[which(dm@eigenvectors[,"DC1"]== max(dm@eigenvectors[,"DC1"]))]),
                       which(rownames(dm@eigenvectors) == row.names(dm@eigenvectors)[which(dm@eigenvectors[,"DC2"]== min(dm@eigenvectors[,"DC2"]))])),
             which(rownames(dm@eigenvectors) == row.names(dm@eigenvectors)[which(dm@eigenvectors[,"DC2"]== min(dm@eigenvectors[,"DC2"]))]))
  plot(dpt)
  
  diffusion_components <- as.data.frame(dm@eigenvectors)
  colnames(diffusion_components) <- paste0("DC", 1:ncol(diffusion_components))
  object_mm[['diffmap']] <- CreateDimReducObject(embeddings = as.matrix(diffusion_components), key = "DC_")
  object_mm$DC1 <- object_mm$DC2 <- object_mm$DPT <- NULL
  object_mm@meta.data <- cbind(object_mm@meta.data,diffusion_components[,c("DC1", "DC2")])
  object_mm$DPT <- dpt$dpt
  
  
# Fig4G & FigS3B
  
  colour_bk <- c(colorRampPalette(c("purple3", "#4575B4"))(15),
                 colorRampPalette(c("#4575B4", "#91CF60"))(10),
                 colorRampPalette(c("#91CF60", "gold"))(5),
                 colorRampPalette(c("gold", "yellow"))(10))
  ggplot(object_mm@meta.data,aes(x=DC1,y=DC2,color=DPT))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk)+
    theme_classic()+NoAxes()
  
  
  colour_bk <- c(colorRampPalette(c("cyan3", "lightblue"))(40),  
                 colorRampPalette(c("lightblue", "linen"))(10),
                 colorRampPalette(c("linen", "lightsalmon"))(35), 
                 colorRampPalette(c("lightsalmon", "brown2"))(20))
  ggplot(object_mm@meta.data,aes(x=DC1,y=DC2,color=Risk_score))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk)+
    theme_classic()+NoAxes()
  
  
  
# Fig4H
  
  object_mm$CellType_risk <- paste0(object_mm$CellType,'-',object_mm$Risk_type)
  CellType_risk.colors <- c(`Monocyte-Low_risk`="royalblue",`Monocyte-High_risk`="hotpink2",
                            `Macrophage-Low_risk`="seagreen3",`Macrophage-High_risk`="gold")
  object_mm$CellType_risk <- factor(object_mm$CellType_risk,levels = names(CellType_risk.colors))
  ggplot(object_mm@meta.data,aes(x=DC1,y=DC2,color=CellType_risk))+
    geom_point_rast(size=0.1)+
    scale_colour_manual(values = CellType_risk.colors,na.value = "gray90")+
    theme_classic()+NoAxes()
  
  
  m.query <- subset(object_mm,subset=CellType_risk == "Macrophage-High_risk")
  m.reference <- subset(object_mm,subset=CellType_risk == "Macrophage-High_risk",invert=T)   
  
  pred <- SingleR(ref=m.reference@assays$RNA@data,
                  test=m.query@assays$RNA@data,
                  labels=m.reference$CellType_risk)     
  pred_singleR <- data.frame(cellID=row.names(pred),SingleR=pred$pruned.labels)
  
  m.query$predict.SingleR <- m.query$predicted.id <- m.query$SingleR <- NULL
  object_mm$predict.SingleR <- object_mm$predicted.id <- object_mm$SingleR <- NULL
  
  m.query@meta.data <- left_join(m.query@meta.data,pred_singleR)
  row.names(m.query@meta.data) <- m.query$cellID
  
  object_mm@meta.data <- left_join(object_mm@meta.data,m.query@meta.data[,c("cellID","SingleR")])
  rownames(object_mm@meta.data) <- object_mm$cellID
  
  
  DimPlot(object_mm,group.by = "SingleR",reduction = "diffmap",na.value = "grey90",raster=T,pt.size=3,
                cols = c(`Monocyte-Low_risk`="royalblue",`Monocyte-High_risk`="hotpink2",`Macrophage-Low_risk`="seagreen3"))+NoLegend()+NoAxes()
  
  pie_data <- data.frame(table(m.query$SingleR))
  pie_data$Var1 <- factor(pie_data$Var1,levels=c("Monocyte-Low_risk","Monocyte-High_risk","Macrophage-Low_risk"))
  ggplot(pie_data, aes(x = "", y = Freq, fill = Var1)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values=c(`Monocyte-Low_risk`="royalblue",`Monocyte-High_risk`="hotpink2",`Macrophage-Low_risk`="seagreen3"))+
    coord_polar(theta = "y") +  # 使用极坐标
    theme_void() +  # 去除背景
    ggtitle("Potency of transition to Macrophage-High_risk")+
    theme(title = element_text(size=6))

  
# FigS3C
  
  phenoPropotion(object_mm@meta.data,x_name = "CellType_risk",y_name="Response",
                 cols=Response.colors,
                 legend_name="Risk type",bar_width = 0.9)

  
# FigS4I & FigS3DE
  
  colour_bk <- c("lightcyan",
                 colorRampPalette(c("#f7f7f7","#fdf2b5"))(5),
                 colorRampPalette(c("#fdf2b5","#F67B51"))(5),
                 colorRampPalette(c("#F67B51","#A30023"))(10))
  for(k in c("Phagocytosis3","Angiogenesis5")){
    
    p1 <- VlnPlot(object_mm,features=k,
            group.by="Risk_type",pt.size=0,cols=Risk_type.colors,adjust=1)+
      stat_compare_means(size=2)+
      stat_summary(fun.y=mean, geom="point", shape=18,size=1.5, color="white")+
      theme_classic()+ylab(k)+xlab("")+
      theme(axis.text.x=element_text(colour = 'black',angle=45,size=6,vjust=1,hjust=1))
    p2 <- ggplot(object_mm@meta.data,aes(x=DC1,y=DC2,color=get(k)))+
      geom_point_rast(size=0.1)+
      scale_colour_gradientn(colours = colour_bk)+
      theme_classic()+ggtitle(k)+NoAxes()
    p3 <- ggplot(object_mm@meta.data, aes(x=Risk_score,y=get(k)))+
      geom_smooth(fullrange=T,method="loess",span=1.5,color="purple3",fill="gray90")+
      ylab(k)+theme_classic()
    print(p1+p2+p3)
  }
  
  
# FigS4F
  
  object_mm$C1QC <- LayerData(object_mm,"data")["C1QC",]
  object_mm$SPP1 <- LayerData(object_mm,"data")["SPP1",]
  
  C1QC_SPP1_data <- data.table(object_mm@meta.data[,c("DC1","DC2","C1QC","SPP1","Risk_type","CellType")])
  C1QC_SPP1_data$C1QC_SPP1 <- ""
  C1QC_SPP1_data[C1QC > 0 & SPP1 == 0,C1QC_SPP1:="C1QC"]
  C1QC_SPP1_data[C1QC == 0 & SPP1 > 0,C1QC_SPP1:="SPP1"]
  C1QC_SPP1_data[C1QC > 0 & SPP1 > 0,C1QC_SPP1:="SPP1&C1QC"]
  
  ggplot(C1QC_SPP1_data,aes(x=DC1,y=DC2,color=C1QC_SPP1))+
    geom_point_rast(size=0.1)+
    scale_colour_manual(values = c(C1QC="seagreen3",SPP1="#F67B51",`SPP1&C1QC`="yellow"),na.value = "gray90")+
    theme_classic()+NoAxes()
  

    
# FigS4G
  
  ggplot(C1QC_SPP1_data[Risk_type == "Low_risk" & CellType == "Macrophage"],aes(x=C1QC,y=SPP1))+
    geom_point_rast(color=alpha("cyan3",0.6),size=0.2)+
    geom_hline(yintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    geom_vline(xintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    theme_classic()+xlim(0,6)+ylim(0,7)+
    ggtitle('Low risk')+xlab('C1QC+')+ylab('SPP1+')
  ggplot(C1QC_SPP1_data[Risk_type == "High_risk" & CellType == "Macrophage"],aes(x=C1QC,y=SPP1))+
    geom_point_rast(color=alpha("brown3",0.6),size=0.2)+
    geom_hline(yintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    geom_vline(xintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    theme_classic()+xlim(0,6)+ylim(0,7)+
    ggtitle('High risk')+xlab('C1QC+')+ylab('SPP1+')
  
  
# FigS4HI
  
  colour_bk <- c("lightcyan",
                 colorRampPalette(c("#f7f7f7","#fdf2b5"))(10),
                 colorRampPalette(c("#fdf2b5","#F67B51"))(5),
                 colorRampPalette(c("#F67B51","#A30023"))(5))
  for(k in c("LAM8","M1M2_t")){
    
    p1 <- VlnPlot(object_mm[,object_mm$CellType == "Macrophage"],features=k,
                  group.by="Risk_type",pt.size=0,cols=Risk_type.colors,adjust=1)+
      stat_compare_means(size=2)+
      stat_summary(fun.y=mean, geom="point", shape=18,size=1.5, color="white")+
      theme_classic()+ylab(k)+xlab("")+
      theme(axis.text.x=element_text(colour = 'black',angle=45,size=6,vjust=1,hjust=1))
    p2 <- ggplot(object_mm@meta.data,aes(x=DC1,y=DC2,color=get(k)))+
      geom_point_rast(size=0.1)+
      scale_colour_gradientn(colours = colour_bk)+
      theme_classic()+ggtitle(k)+NoAxes()
    print(p1+p2)
  }
  
  
  
# FigS4J
  
  object_mm$TREM1 <- LayerData(object_mm,"data")["TREM1",]
  object_mm$TREM2 <- LayerData(object_mm,"data")["TREM2",]
  
  colour_bk <- c("lightcyan",
                 colorRampPalette(c("#f7f7f7","#fdf2b5"))(10),
                 colorRampPalette(c("#fdf2b5","#F67B51"))(5),
                 colorRampPalette(c("#F67B51","#A30023"))(5))
  ggplot(object_mm@meta.data,aes(x=DC1,y=DC2,color=TREM1))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk,na.value = "lightcyan")+
    theme_classic()+NoAxes()
  ggplot(object_mm@meta.data,aes(x=DC1,y=DC2,color=TREM2))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk,na.value = "lightcyan")+
    theme_classic()+NoAxes()
  
  
  
# Fig4J
  
  plot_data <- pivot_longer(object_mm@meta.data[,c("TREM1","TREM2","Risk_score")], 
                            cols = c("TREM1", "TREM2"), 
                            names_to = "Gene", 
                            values_to = "Expression")
  ggplot(plot_data, aes(x = Risk_score, y = Expression, color = Gene, fill = Gene)) +
    geom_smooth(
      method = "loess", 
      span = 1.5, 
      alpha = 0.2,  # 填充色透明度
      linewidth = 0.8
    ) +
    scale_color_manual(values = c("TREM1" = "purple3", "TREM2" = "orange2")) +
    scale_fill_manual(values = c("TREM1" = "gray90", "TREM2" = "gray90")) +
    labs(x = "Risk Score",  y = "Gene Expression", color = "Gene", fill = "Gene")+
    theme_classic()

  
    
# Fig4K
  
  TREM1_TREM2_data <- data.table(object_mm@meta.data[,c("DC1","DC2","TREM1","TREM2","Risk_type","CellType","M1M2_t")])
  TREM1_TREM2_data$TREM1_TREM2 <- ""
  TREM1_TREM2_data[TREM1 > 0 & TREM2 == 0,TREM1_TREM2:="TREM1"]
  TREM1_TREM2_data[TREM1 == 0 & TREM2 > 0,TREM1_TREM2:="TREM2"]
  TREM1_TREM2_data[TREM1 > 0 & TREM2 > 0,TREM1_TREM2:="TREM1&TREM2"]
  
  ggplot(TREM1_TREM2_data,aes(x=DC1,y=DC2,color=TREM1_TREM2))+
    geom_point_rast(size=0.1)+
    scale_colour_manual(values = c(TREM1="seagreen3",TREM2="#F67B51",`TREM1&TREM2`="yellow"),na.value = "gray90")+
    theme_classic()+NoAxes()
  
  
  
# Fig4L
  
  ggplot(TREM1_TREM2_data[Risk_type == "Low_risk" & CellType == "Macrophage"],aes(x=TREM1,y=TREM2))+
    geom_point_rast(color=alpha("cyan3",0.6),size=0.2)+
    geom_hline(yintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    geom_vline(xintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    theme_classic()+xlim(0,3.6)+ylim(0,4.1)+
    ggtitle('Low risk')+xlab('TREM1+')+ylab('TREM2+')
  ggplot(TREM1_TREM2_data[Risk_type == "High_risk" & CellType == "Macrophage"],aes(x=TREM1,y=TREM2))+
    geom_point_rast(color=alpha("brown3",0.6),size=0.2)+
    geom_hline(yintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    geom_vline(xintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    theme_classic()+xlim(0,3.6)+ylim(0,4.1)+
    ggtitle('High risk')+xlab('TREM1+')+ylab('TREM2+')
  
  
  
# Fig4MN  
  
  object_dc <- subset(object_m,subset=CellType2 == "cDC2") 
  object_dc$Risk_type <- factor(object_dc$Risk_type,levels = c("Low_risk","High_risk"))
  
  dc_markers <- list(Antigen_presenting = c("CD74","TAP1","TAP2","TAPBP","CD36","HLA-DRA",
                                            "HLA-DRB1","HLA-DPB1","HLA-DOA","HLA-DOB","HLA-DRB5",
                                            "HLA-DRB6","HLA-A","HLA-B","HLA-C"),
                     Angiogenesis = c("CCND2","CCNE1","CD44","CXCR4","E2F3","EDN1","EZH2",
                                      "FGF18","FGFR1","FYN","HEY1","ITGAV","JAG1","JAG2",
                                      "MMP9","NOTCH1","PDGFA","PTK2","SPP1","STC1","TNFAIP6",
                                      "TYMP","VAV2","VCAN","VEGFA"),
                     cDC2_mature = c("ADAM19","LSP1","PIM2","SLAMF1","SLAMF7","CXCL11","ETV6","PSME2",
                                     "SELPLG","TAP1","TAPBP","ADTRP","MMP25","SEMA7A","IFIH1","DHX58",
                                     "SDC4","TXNIP","TNFSF4","WIPF1","FSTL3","GSTP1"))
  object_dc <- AddModuleScore(object_dc,dc_markers,name = names(dc_markers))
  
  GeneCounts <- as.matrix(object_dc@assays$RNA@counts)
  iOrd <- rowSums(GeneCounts>0)
  GeneCounts <- GeneCounts[iOrd>10,]#only keep genes expressing in more than 10 cell
  CytoTRACE.res <- CytoTRACE(GeneCounts, batch = object_dc$orig.ident)
  object_dc$CytoTRACE <- CytoTRACE.res$CytoTRACE[colnames(object_dc)]
  
  
  
  ggplot(object_dc@meta.data,aes(x=Risk_type,y=Antigen_presenting1,fill=Risk_type))+
    geom_boxplot(color='gray40',outlier.shape=NA,size=0.1,width=0.4,linewidth  = 0.235)+
    geom_jitter(fill='white',color='black',shape=21,width =0.1,size=0.4,stroke = 0.1)+
    stat_compare_means(method="wilcox.test",size=2,bracket.size= 0.235)+
    scale_fill_manual(values=Risk_type.colors)+
    ylab("Antigen presenting")+xlab("")+
    theme_classic() 
  ggplot(object_dc@meta.data,aes(x=Risk_type,y=Angiogenesis2,fill=Risk_type))+
    geom_boxplot(color='gray40',outlier.shape=NA,size=0.1,width=0.4,linewidth  = 0.235)+
    geom_jitter(fill='white',color='black',shape=21,width =0.1,size=0.4,stroke = 0.1)+
    stat_compare_means(method="wilcox.test",size=2,bracket.size= 0.235)+
    scale_fill_manual(values=Risk_type.colors)+
    ylab("Angiogenesis")+xlab("")+
    theme_classic() 
  
  ggplot(object_dc@meta.data, aes(x=CytoTRACE,y=Risk_type,fill=Risk_type))+
    geom_density_ridges(scale=2) +
    scale_fill_manual(values=Risk_type.colors)+
    theme_ridges()+NoLegend()
  ggplot(object_dc@meta.data, aes(x = Risk_score, y = cDC2_mature3))+
    geom_point(shape =3,size=0.2,color="seagreen")+
    theme_classic()+
    geom_smooth(method ="lm",color ="black",fill ="gray80",size=0.6)+
    stat_cor(method ="spearman",cor.coef.name="rho",label.x =1, label.y=0.4,size =2.5) +
    stat_poly_eq(aes(label = ..eq.label..),
      formula =y ~x,parse =TRUE, geom ="text",label.x =1,label.y =0.5, hjust =0,size =2.5)
  
  

# Fig4O
  
  # 654 cells
  object_neu <- subset(object_m,subset=CellType == "Neutrophil") 
  object_neu$Risk_type <- factor(object_neu$Risk_type,levels = c("Low_risk","High_risk"))
  neu_markers <- list(neu_T3 = c("ATF3","CCL3","CCL4","CD274","CSTB","CXCL3","HCAR2","HILPDA","HK2","HMOX1",
                                 "IER3","JUN","LDHA","MIF","PLIN2","SPP1","TGIF1","TNFRSF10D","VEGFA","ZEB2"))
  object_neu <- AddModuleScore(object_neu,neu_markers,name = names(neu_markers))
  
  GeneCounts <- as.matrix(object_neu@assays$RNA@counts)
  iOrd <- rowSums(GeneCounts>0)
  GeneCounts <- GeneCounts[iOrd>10,]#only keep genes expressing in more than 10 cell
  CytoTRACE.res <- CytoTRACE(GeneCounts, batch = object_neu$orig.ident)
  object_neu$CytoTRACE <- CytoTRACE.res$CytoTRACE[colnames(object_neu)]
  
  ggplot(object_neu@meta.data,aes(x=Risk_type,y=CytoTRACE,fill=Risk_type))+
    geom_boxplot(color='gray40',outlier.shape=NA,size=0.1,width=0.4,linewidth  = 0.235)+
    geom_jitter(fill='white',color='black',shape=21,width =0.1,size=0.4,stroke = 0.1)+
    stat_compare_means(method="wilcox.test",size=2,bracket.size= 0.235)+
    scale_fill_manual(values=Risk_type.colors)+
    ylab("CytoTRACE")+xlab("")+
    theme_classic() 
  ggplot(object_neu@meta.data,aes(x=Risk_type,y=neu_T31,fill=Risk_type))+
    geom_boxplot(color='gray40',outlier.shape=NA,size=0.1,width=0.4,linewidth  = 0.235)+
    geom_jitter(fill='white',color='black',shape=21,width =0.1,size=0.4,stroke = 0.1)+
    stat_compare_means(method="wilcox.test",size=2,bracket.size= 0.235)+
    scale_fill_manual(values=Risk_type.colors)+
    ylab("neu_T3")+xlab("")+
    theme_classic() 

  
  
  
  


# FigS1A
  
  CellType_markers <- list(
    epithelium=unlist(list(Malignant=c("S100P"),
                           Ductal=c("EPCAM","KRT8","KRT18"),
                           Acinar=c("PRSS1","CPA1"),
                           #Cholangiocyte
                           Hepatocytes=c("ARG1","ALB"))),
    Immune=unlist(list(Immune=c("PTPRC"),
                       `T`=c("CD3D","CD3E","CD4","CD8A","CD8B","MKI67","TOP2A","TRDV1","TRGV5"),
                       NK=c("NCAM1","FCGR3A","KLRD1"),
                       ILC=c("RORC","KIT"),
                       B=c("MS4A1","CD19","CD79A"),
                       Plasma=c("JCHAIN","IGHG1","IGHA1"),
                       Monocyte=c("FCN1","VCAN"),
                       Macrophage=c("CD68","CD163"),
                       cDC=c("CLEC9A","CD1C","LAMP3"),
                       Neutrophil=c("CXCL8","CSF3R","S100A8"),
                       Mast=c("MS4A2","CPA3"),
                       pDC=c("LILRA4","SPIB"))),
    Stromal=unlist(list(TEC=c("PECAM1","CD34","PLVAP"),
                        HSEC=c("CLEC4G","OIT3","IL33"),
                        Fibroblast=c("ACTA2","TAGLN","COL1A2"),
                        myCAFs=c("MMP2","IGFBP5","POSTN"),
                        iCAF=c("CXCL12","C3","C7"),
                        vCAF=c("MCAM","RGS5","GJA4"))))
  object$CellType <- factor(object$CellType,levels=rev(names(CellType.colors)))
  DotPlot(object, features = CellType_markers,group.by = "CellType",dot.scale=2.5)+
    scale_y_discrete(position="right")+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5,size=6),
          axis.text.y=element_text(hjust = 1,vjust=0.5,size=8))+
    scale_color_gradientn(values = seq(0,1,0.2),colours = c("royalblue3","cyan3","linen","#ffbb78","brown3"))+
    labs(x=NULL,y=NULL)+
    guides(size=guide_legend(order=3))


  
# FigS1B
  
  DimPlot(object,cells.highlight = list(`TCR`=object$cellID[object$has_TCR == "Yes"]),
          cols.highlight="turquoise3",cols="linen",sizes.highlight=2,raster=T)+
    theme_void()+
    NoLegend()+NoAxes()


  
# FigS1E
  
  DimPlot(object,cells.highlight = list(`KRAS-Mutation`=object$cellID[object$MutGene == "KRAS"]),
          cols.highlight="red2",cols="linen",sizes.highlight=2,raster=T)

  

# Epithelium
  
  object_epi <- subset(object,subset=CellType_Main == "Epithelial cell")
  object_epi <- scCluster(object_epi,seed=3,min.d=0.4)


  
# FigS1C
  
  DimPlot(object_epi,group.by = "orig.ident",cols = orig.ident.colors,raster=T,pt.size=1.8)

  

# FigS1D
  
  colour_bk <- c("lightcyan",
                 colorRampPalette(c("#f7f7f7","#fdf2b5"))(5),
                 colorRampPalette(c("#fdf2b5","#F67B51"))(5),
                 colorRampPalette(c("#F67B51","#A30023"))(10))
  FeaturePlot(object_epi,"S100P")+
    scale_color_gradientn(colours=colour_bk)
  

  
# FigS1F
  
  load("E:/PAAD_Project/Data/scMutation/ScVarscan/cell_mut_n.rda")
  cell_mut_n <- cell_mut_n[MutGene=="KRAS"& Mut_N>1]
  object_epi@meta.data <- left_join(object_epi@meta.data,cell_mut_n[,c("cellID","MutType")])
  row.names(object_epi@meta.data) <- object_epi$cellID
  MutType.colors <- ColAssign(unique(object_epi$MutType))
  MutType.colors["KRAS_p.G12C"] <- "springgreen3"
  MutType.colors["KRAS_p.G12D"] <- "royalblue2"
  MutType.colors["KRAS_p.G12V"] <- "tomato"
  MutType.colors["KRAS_p.Q61H"] <- "darkviolet"
  DimPlot(object_epi,group.by = "MutType",cols = MutType.colors,na.value = "grey95",label=F,pt.size=1.5,raster=T)
  
  
  
# FigS1G
  
  # infercnv_obj
  load("E:/PAAD_Project/Data/Tumor_heterogeneity/inferCNV_all/object_CNV.rda")
  tumor_samples <- paste0("malignant_",orig.ident)
  tumor_cells <- lapply(tumor_samples,function(s){
    return(data.table(celltype=s,cellID=infercnv_obj@tumor_subclusters$hc[[s]]$labels))
  })
  tumor_cells <- rbindlist(tumor_cells)
  tumor_cells[celltype %like% "P", Type:="P"]
  tumor_cells[celltype %like% "M", Type:="M"]
  tumor_cells[, sample:=str_remove(celltype,"malignant_")]
  setDF(tumor_cells,rownames = tumor_cells$cellID)
  tumor_cells$cellID <- NULL
  
  hp_data <- t(infercnv_obj@expr.data[,rownames(tumor_cells)])
  
  show_gene <- c("KRAS","TP53","CDKN2A","SMAD4","MYC")
  
  top_anno <- columnAnnotation(
    show=anno_mark(at= which(colnames(hp_data) %in% show_gene),
                   labels = colnames(hp_data)[which(colnames(hp_data) %in% show_gene)],
                   link_width = unit(1, "mm"),link_height=unit(5, "mm"),labels_gp = gpar(fontsize = 5)))
  
  chromosome.colors <- rep(c("gray95","gray80"),11)
  names(chromosome.colors) <- paste0("chr",1:22)
  bottom_anno <- HeatmapAnnotation(
    chromosome = infercnv_obj@gene_order$chr,
    show_legend = F,
    annotation_name_side = "right",
    annotation_height = unit(0.3, "cm"),
    col = list(chromosome = chromosome.colors))
  
  left_anno <- rowAnnotation(
    df = tumor_cells[,c("Type","sample")],
    col=list(Type = Type.colors,sample = orig.ident.colors))
  
  infercnv_col <- c(
    "#436EEE","#5F84F0","#7C9AF3","#99B0F5","#B6C7F8","#D3DDFB","#F0F3FD",
    "white",
    "#FBEEEE","#F3CCCC","#EBABAB","#E48A8A","#DC6868","#D44747","#CD2626")
  
  Heatmap(hp_data,
          col =colorRamp2(seq(0.8, 1.2, length.out = 15), infercnv_col),
          cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
          column_split = factor(infercnv_obj@gene_order$chr, paste0("chr",1:22)),
          column_gap = unit(0.6, "mm"),
          row_split = factor(tumor_cells$Type,levels = c("P","M")),
          row_gap = unit(1,'mm'),
          heatmap_legend_param = list(title = "Modified expression (inferCNV)",    
                                      direction = "vertical",
                                      title_position = "leftcenter-rot",
                                      at=c(0.8,1,1.2),legend_height = unit(3, "cm")),
          bottom_annotation = bottom_anno,left_annotation = left_anno,top_annotation = top_anno)
 


# FigS1H
  
  load("E:/PAAD_Project/Data/Tumor_heterogeneity/inferCNV_all/object_CNV.rda")
  cnv_data <- infercnv_obj@expr.data
  cnv_score <- apply(cnv_data,2,function(x) mean(sum(abs(x-1))))
  cnv_score <- data.table(cnv_score,value.name=names(cnv_score))
  names(cnv_score)[2] <- "cellID"
  
  object_epi@meta.data <- left_join(object_epi@meta.data,cnv_score)
  row.names(object_epi@meta.data) <- object_epi$cellID
  
  cnv_tumor <- object_epi$cnv_score[object_epi$CellType2 == "Malignant"]
  cnv_nontumor <- object_epi$cnv_score[object_epi$CellType2 != "Malignant"]
  # p-value < 2.2e-16
  wilcox.test(cnv_tumor,cnv_nontumor)
  
  ggplot(object_epi@meta.data,aes(CellType3,cnv_score))+
    geom_boxplot(aes(fill=CellType3),color="gray70",size=0.1,show.legend=F,outlier.color=NA,width=0.6)+
    stat_summary(fun.y=median, geom="point", shape=18,size=1, color="white")+
    geom_segment(aes(x = 0.5, xend = 5.5, y = mean(cnv_nontumor), yend = mean(cnv_nontumor)), color = "cyan3", linetype = "dashed", size = 0.5) +
    geom_segment(aes(x = 5.5, xend = 17.5, y = mean(cnv_tumor), yend = mean(cnv_tumor)),  color = "brown3",linetype = "dashed", size = 0.5) +
    ylab("CNV score")+xlab("") +ylim(0,500)+
    annotate("text",label="p-value < 2.2e-16", x = 5.5, y = 490,size=2.5)+
    theme_base()+ theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))

  

# FigS1I
  
  show_gene <- c("KRAS","TP53","CDKN2A","SMAD4","MYC")
  gene_expr <- as.data.frame(hp_data[,show_gene])
  gene_expr$cellID <- rownames(gene_expr)
  gene_expr <- gene_expr %>%
    pivot_longer(
      cols = -cellID,
      names_to = "gene",  
      values_to = "expression")
  gene_expr$gene <- factor(gene_expr$gene,levels = c("KRAS","MYC","TP53","CDKN2A","SMAD4"))
  
  gene_expr <- gene_expr[gene_expr$expression < 0.995 | gene_expr$expression > 1.005, ]
  
  ggplot(gene_expr, aes(x = expression, y = gene, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
    stat_density_ridges(geom = "density_ridges_gradient", scale = 1.2, calc_ecdf = TRUE) +
    scale_fill_viridis_c(option = "H",direction=1,name = "Tail probability")+
    xlim(c(0.7,1.3))+xlab("Modified expression (inferCNV)")+ylab("")


  
# Fig2B
  
  DimPlot(object,group.by="CellType",label=T,raster=T,cols=CellType.colors)


  
# Fig2C
  
  CellType_select <- rowSums(table(object$CellType,object$Risk_type) > 50)==2
  CellType_select <- names(CellType_select[which(CellType_select)])
  diff_CellType <- lapply(CellType_select,function(c){
    print(c)
    diff_c <- FindMarkers(object[,object$CellType==c],ident.1="High_risk",group.by="Risk_type",logfc.threshold = 0.25,min.pct = 0.1)
    diff_c$gene <- row.names(diff_c)
    diff_c$CellType <- c
    setDT(diff_c)
    setorder(diff_c,"avg_log2FC")
    return(diff_c)
  })
  diff_CellType <- rbindlist(diff_CellType)
  diff_CellType <- diff_CellType[p_val_adj < 0.05]
  
  
  cellID_DEGs <- lapply(CellType_select,function(c){
    print(c)
    mat <- LayerData(object,"counts")[diff_CellType[CellType==c,gene],object$cellID[object$CellType == c]]
    index <- colSums(mat>0) 
    cellID_DEGs <- data.table(cellID=names(index),DEGs=index)
    return(cellID_DEGs)
  })
  cellID_DEGs <- rbindlist(cellID_DEGs)
  
  object@meta.data <- left_join(object@meta.data,cellID_DEGs)
  row.names(object@meta.data) <- object$cellID
  
  FeaturePlot(object,features = 'DEGs',pt.size = 1.5,raster = T)+
    scale_color_gradientn(colours = c(colorRampPalette(c("linen","#fdf2b5"))(5),
                                      colorRampPalette(c("#fdf2b5","#F67B51"))(5),
                                      colorRampPalette(c("#F67B51","#A30023"))(20)),
                          na.value = "grey90")+NoAxes()+
    ggtitle("Number of DEGs (High_risk vs. Low_risk)")
  
  
  
# Fig2D
  
  hp_data_prog <- hp_data[,intersect(colnames(hp_data),c(favorable_genes,adverse_genes))]
  
  gene_expr_prog <- as.data.frame(hp_data_prog)
  gene_expr_prog$cellID <- rownames(gene_expr_prog)
  gene_expr_prog <- gene_expr_prog %>%
    pivot_longer(
      cols = -cellID,
      names_to = "gene",  
      values_to = "expression")
  setDT(gene_expr_prog)
  gene_expr_prog[gene %in% favorable_genes, gene_type:="Favorable_genes"]
  gene_expr_prog[gene %in% adverse_genes, gene_type:="Adverse_genes"]
  
  gene_expr_prog <- gene_expr_prog[gene_expr_prog$expression < 0.995 | gene_expr_prog$expression > 1.005, ]
  ggplot(gene_expr_prog, aes(x = expression, y = gene_type, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
    stat_density_ridges(geom = "density_ridges_gradient", scale = 0.8, calc_ecdf = TRUE) +
    scale_fill_viridis_c(option = "H",direction=1,name = "Tail probability")+
    xlim(c(0.7,1.3))+xlab("Modified expression (inferCNV)")+ylab("")
  
  
  
# Fig2E
  
  colour_bk <- c(colorRampPalette(c("cyan3", alpha("cyan3",0.2)))(20),  
                 colorRampPalette(c(alpha("cyan3",0.2), "linen"))(27),
                 colorRampPalette(c("linen", "lightsalmon"))(30), 
                 colorRampPalette(c("lightsalmon", "brown3"))(45))
  FeaturePlot(object,"Risk_score",raster=T,pt.size=1.5)+NoAxes()+
    scale_colour_gradientn(colours = colour_bk)
  
  DimPlot(object,group.by="Risk_type",label=F,raster=T,pt.size=1.5,cols=Risk_type.colors)+NoAxes()
  
  
  
# Fig2F  
  
  object_risk <- subset(object,subset=CellType %in% c("Acinar","Ductal","Hepatocyte","Cholangiocyte","Prof.T"),invert=T)
  object_risk$CellType_Main <- str_remove_all(object_risk$CellType_Main," cell")
  object_risk$CellType_Main <- str_replace_all(object_risk$CellType_Main,"Epithelial","Malignant")
  CellType_Main.colors <- c(Malignant="orangered3",Myeloid="gold2",Fibroblast="#8c564b",Endothelial="#ff9896",Lymphocyte="seagreen3")
  object_risk$CellType_Main <- factor(object_risk$CellType_Main,levels=names(CellType_Main.colors))
  
  VlnPlot(object_risk,"Risk_score",group.by="CellType_Main",pt.size=0,cols = CellType_Main.colors)+
    stat_summary(fun.y=mean, geom="point", shape=18,size=1, color="white")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "royalblue", size = 0.5)+
    theme_base()+ylab("Risk score")+
    theme(legend.position = 'none',
          axis.text.x=element_text(size = 6,colour = 'black',angle=90,hjust=1,vjust=0.5),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.y = element_text(size = 6,colour = 'black'),
          axis.line = element_line(linewidth  = 0.235), #保证坐标轴描边为0.5
          axis.ticks = element_line(linewidth  = 0.235),
          plot.title = element_text(hjust=0.5,size = 8))
  
  phenoPropotion(data=object_risk@meta.data,x_name="CellType_Main",y_name="Risk_type",cols=Risk_type.colors,legend_name="Risk_type",bar_width=0.9,out_path,plot_width=5,plot_height=6.5)

  
  
# Fig2G
  
  DimPlot(object,group.by="CellType",label=F,raster=T,cols=CellType.colors,split.by = "Risk_type")
  
  metaData <- cbind(object@meta.data,Embeddings(object,reduction = "umap"))
  
  data_high_risk <- metaData[metaData$Risk_type == "High_risk",]
  ggplot(data = data_high_risk,aes(x = umap_1, y = umap_2)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F)+
    geom_point(data = data_high_risk[sample(data_high_risk$cellID,8000,replace = F),],
               aes(x = umap_1, y = umap_2),color = 'white',size = .005)+
    scale_fill_viridis(option="A")+theme_classic()+NoAxes()+NoLegend()
  
  data_low_risk <- metaData[metaData$Risk_type == "Low_risk",]
  ggplot(data = data_low_risk,aes(x = umap_1, y = umap_2)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F)+
    geom_point(data = data_low_risk[sample(data_low_risk$cellID,8000,replace = F),],
               aes(x = umap_1, y = umap_2),color = 'white',size = .005)+
    scale_fill_viridis(option="A")+theme_classic()+NoAxes()+NoLegend()
  
  
  
# Fig2H
  
  cluster_data <- as.matrix(LayerData(object,"data")[intersect(c(adverse_genes,favorable_genes),rownames(object)),])
  cluster_data <- t(groupMeans(cluster_data,groups=object$CellType_risk))
  
  info.v <- as.data.frame(table(object$CellType_risk))
  colnames(info.v) <- c("ID","Freq")
  info.v$Freq_log <- log2(info.v$Freq)
  info.v$Type <- ifelse(info.v$ID %like% "High_risk$","High_risk","Low_risk")
  risk_score <- data.table(ID=object$CellType_risk,Risk_score=object$Risk_score)
  risk_score <- risk_score[,mean(Risk_score),by=ID]
  names(risk_score)[2] <- "Risk_score"
  info.v <- left_join(info.v,risk_score)
  
  info.v$ID <- str_replace_all(info.v$ID,"_High_risk","-H")
  info.v$ID <- str_replace_all(info.v$ID,"_Low_risk","-L")
  rownames(cluster_data) <- str_replace_all(rownames(cluster_data),"_High_risk","-H")
  rownames(cluster_data) <- str_replace_all(rownames(cluster_data),"_Low_risk","-L")
  
  hc <- hclust(dist(cluster_data, method = "manhattan"),method = "ward.D")

  tal <- treeAndLeaf(hc)
  tal <- att.mapv(g = tal, dat = info.v )   
  colour_bk <- c(colorRampPalette(c("cyan3","#abd9e9"))(14),
                 colorRampPalette(c("#abd9e9","#fee090"))(10),
                 colorRampPalette(c("#fee090","brown3"))(18))
  
  tal <- att.setv(g = tal, from = "Risk_score", to = "nodeColor", cols = colour_bk)
  tal <- att.setv(g = tal, from = "Freq_log", to = "nodeSize", xlim = c(2, 40,2))
  tal <- att.addv(tal, "nodeFontSize", value = 12, index = V(tal)$isLeaf)
  tal <- att.adde(tal, "edgeWidth", value = 20)
  
  rdp <- RedPort()
  calld(rdp)
  resetd(rdp)
  addGraph(obj = rdp, g = tal, zoom =100)
  relax(rdp, p1=25, p2=200, p3=10, p4=100, p5=10)
  
  

# Fig2I
  
  object_risk$rogue_type <- paste0(object_risk$CellType_Main,"-",object_risk$Risk_type)    
  expr <- as.matrix(LayerData(object_risk,layer="counts"))
  expr <- matr.filter(expr, min.cells = 10, min.genes = 10)

  cell_type_rogue <- rogue(expr,labels = object_risk$rogue_type, samples = object_risk$orig.ident,
                           platform = "UMI",span = 0.6)
  cell_type_rogue$orig.ident <- row.names(cell_type_rogue)
  cell_type_rogue <- pivot_longer(cell_type_rogue,cols=setdiff(names(cell_type_rogue),"orig.ident"),names_to = "rogue_type",values_to = "ROGUE")
  
  ggplot(cell_type_rogue, aes(x = rogue_type, y = ROGUE)) +
    geom_boxplot(width=0.55,outlier.size = 0) +  
    geom_point(aes(color = orig.ident), size = 0.6) +  
    scale_color_manual(values=orig.ident.colors)+
    ylim(0.35,0.85)+
    theme_base()+
    theme(legend.position = 'none',axis.text.x=element_text(size = 6,colour = 'black',angle=45,hjust=1,vjust=1))+
    stat_compare_means(comparisons=list(c("Malignant-High_risk","Malignant-Low_risk"),c("Myeloid-High_risk","Myeloid-Low_risk")))
                                        
    
# Fig2J
  
  object_bhatt <- object_risk
  object_bhatt <- scCluster(object_bhatt,nfeature=2500,min.d=0.45,res=1,TCRBCR=F,Ribo_MT=F,harmony=T,seed=6,k_param=30)
  
  cell_embed <- Reductions(object_bhatt,"harmony")@cell.embeddings
  
  cell_nums <- as.matrix(table(object_bhatt$CellType,object_bhatt$Risk_type)) 
  CellTypes <- rownames(cell_nums)[rowSums(cell_nums>200) == 2]
  CellTypes <- setdiff(CellTypes,"Ductal")
  
  celltype_distance_list <- lapply(CellTypes,function(c){
    
    lapply(1:100,function(r){
      set.seed(r)
      cell_highrisk <- sample(object_bhatt$cellID[object_bhatt$CellType == c & object_bhatt$Risk_type == "High_risk"],100)
      cell_lowrisk <- sample(object_bhatt$cellID[object_bhatt$CellType == c & object_bhatt$Risk_type == "Low_risk"],100)
      embed_highrisk <- cell_embed[cell_highrisk,]
      embed_lowrisk <- cell_embed[cell_lowrisk,]
      
      set.seed(1)
      d <- dim_dist(embed_mat_x=embed_highrisk,embed_mat_y=embed_lowrisk,dims_use=1:10,num_cells_sample=100,random_sample=FALSE)[1,]
      set.seed(2)
      dr  <- dim_dist(embed_mat_x=embed_highrisk,embed_mat_y=embed_lowrisk,dims_use=1:10,num_cells_sample=100,random_sample=TRUE)[1,]
      
      celltype_distance_c <- data.table(CellType=c,bh_distance=d,bh_distance_random=dr)
      return(celltype_distance_c)
    })
  })
  celltype_distance <- rbindlist(unlist(celltype_distance_list,recursive=F))
  
  ggplot(celltype_distance,aes(x=reorder(CellType, bh_distance, FUN = median),y=bh_distance,fill=CellType))+
    geom_boxplot(width=0.6,outlier.size = 0)+
    stat_summary(fun.y=median, geom="point", shape=18,size=1, color="white")+
    xlab("")+
    scale_fill_manual(values = CellType.colors[unique(celltype_distance$CellType)])+
    theme_base()+ylab("Bhatt distance")+
    theme(legend.position = 'none',
          axis.text.x=element_text(size = 6,colour = 'black',angle=45,hjust=1,vjust=1),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.y = element_text(size = 6,colour = 'black'),
          axis.line = element_line(linewidth  = 0.235), 
          axis.ticks = element_line(linewidth  = 0.235),
          plot.title = element_text(hjust=0.5,size = 8))
  
  
  
  
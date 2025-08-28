

# Fibroblasts

  obj_CAF <- subset(object,subset=CellType%in%c('myCAF','vCAF','iCAF'))
  obj_CAF <- scCluster(obj_CAF,nfeature=2000,min.d=0.4,res=0.7,TCRBCR=F,Ribo_MT=T,harmony=T,seed=1,k_param=30,npcs=30)
  obj_CAF <- subset(obj_CAF,subset=seurat_clusters!=7)
  obj_CAF <- scCluster(obj_CAF,nfeature=2000,min.d=0.4,res=0.7,TCRBCR=F,Ribo_MT=T,harmony=T,seed=1,k_param=30,npcs=30)
  activated_CAF <- c("FAP","POSTN","LRRC15","GREM1")
  obj_CAF <- AddModuleScore(obj_CAF,list(activated_CAF=c("FAP","POSTN","LRRC15","GREM1")),name="activated_CAF")



# Fig5A

  color_CAF <- c(myCAF="indianred2",iCAF= "#45C7D6",vCAF= "#A26DB7")
  DimPlot(obj_CAF, reduction = "umap",cols=color_CAF,group.by = 'CellType2',
          label =T,label.size = 3,pt.size = 3,raster=T,shuffle=T)+ggtitle("")

  
  
# Fig5B 
  
  FeaturePlot(obj_CAF,"Risk_score",raster=T,pt.size=3)+NoAxes()+
    scale_colour_gradientn(colours = c(colorRampPalette(c("cyan3", "lightblue"))(20),  
                                       colorRampPalette(c("lightblue", "linen"))(10),
                                       colorRampPalette(c("linen", "lightsalmon2"))(10), 
                                       colorRampPalette(c("lightsalmon2", "brown2"))(10),
                                       colorRampPalette(c("brown2", "red"))(25)))


# Fig5C
  
  phenoPropotion(obj_CAF@meta.data,x_name = "CellType2",y_name="Risk_type",
                 cols=Risk_type.colors,
                 legend_name="Risk type",bar_width = 0.9)


# Fig5D
  
  DEGs_Risk_myCAF <- diff_calculate(subset(obj_CAF,subset=CellType2=='myCAF'),ident_1="High_risk",group="Risk_type",log2FC.thres=0.25,p.thres=0.05)
  DEGs_Risk_iCAF <- diff_calculate(subset(obj_CAF,subset=CellType2=='iCAF'),ident_1="High_risk",group="Risk_type",log2FC.thres=0.25,p.thres=0.05)
  DEGs_Risk_vCAF <- diff_calculate(subset(obj_CAF,subset=CellType2=='vCAF'),ident_1="High_risk",group="Risk_type",log2FC.thres=0.25,p.thres=0.05)
  colnames(DEGs_Risk_myCAF)[-1] <- paste0(colnames(DEGs_Risk_myCAF)[-1],"-","Risk_myCAF")
  colnames(DEGs_Risk_iCAF)[-1] <- paste0(colnames(DEGs_Risk_iCAF)[-1],"-","Risk_iCAF")
  colnames(DEGs_Risk_vCAF)[-1] <- paste0(colnames(DEGs_Risk_vCAF)[-1],"-","Risk_vCAF")
  
  df <- Reduce("inner_join",list(DEGs_Risk_myCAF,DEGs_Risk_iCAF,DEGs_Risk_vCAF))
  
  up_label <- intersect(df[`type-Risk_myCAF` == "Up" & `type-Risk_iCAF` == "Up" & `type-Risk_vCAF` == "Up", gene],adverse_genes)
  down_label <- intersect(df[`type-Risk_myCAF` == "Down" & `type-Risk_iCAF` == "Down" & `type-Risk_vCAF` == "Down", gene],favorable_genes)
  
  df$label <- df$gene
  df$label[!(df$label %in% c(up_label,down_label))] <- NA
  
  ggplot(data=df, aes(x = `avg_log2FC-Risk_myCAF`, y = `avg_log2FC-Risk_iCAF`, color=`avg_log2FC-Risk_vCAF`,label=label)) +
    geom_point_rast(stroke = 0, alpha = 1,shape = 16,size=0.5) +
    geom_hline(yintercept=0, linetype="dashed", color = "gold")+
    geom_vline(xintercept=0, linetype="dashed", color = "gold")+
    geom_text_repel(color="black",size=1.5,segment.size=0.1,force= 0.1,force_pull= 1,point.padding=NA,max.overlaps=10000) +
    scale_color_gradient2(midpoint=0, low="cyan3", mid="ivory",high="brown3",space ="Lab") +
    ylab("Log2FC of iCAF High risk vs. Low risk")+
    xlab("Log2FC of myCAF High risk vs. Low risk")+
    labs(color = "Log2FC of vCAF High risk vs. Low risk")+
    theme_classic()


# Fig5E
  
  common_up <-  df[`type-Risk_myCAF` == "Up" & `type-Risk_iCAF` == "Up" & `type-Risk_vCAF` == "Up", gene]
  
  #Up - Go
  bp1 <- enrichGO(common_up,keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
  enrich_adverse_simplify <- clusterProfiler::simplify(bp1,cutoff = 0.6)
  
  
  select_path <- c(
    'collagen fibril organization',
    'collagen metabolic process',
    'wound healing',
    'cell-substrate junction assembly',
    'regulation of cell-substrate adhesion',
    'biological process involved in symbiotic interaction',
    'connective tissue development',
    'fibroblast proliferation',
    'positive regulation of epithelial to mesenchymal transition',
    'regulation of angiogenesis',
    'regulation of cell morphogenesis',
    'cellular response to amino acid stimulus',
    'hemostasis'
  )
  
  plot_data <- enrich_adverse_simplify@result %>% data.table()
  plot_data <- plot_data[Description%in%select_path]
  setorder(plot_data,p.adjust)
  
  plot_data$plotp <- -log10(plot_data$p.adjust)
  plot_data$Description <- tools::toTitleCase(plot_data$Description)
  plot_data$Description <- factor(plot_data$Description,levels = rev(plot_data$Description))
  
  ggplot(plot_data, aes(x = plotp, y= Description)) +
    geom_col(fill="lightpink1")+
    geom_text(aes(x = 0.05, label = Description),size = 2,hjust = 0)+
    labs(x = '-log10 (p.adjust)', y = '', title = 'High risk vs. Low risk')+theme_classic()

  
  
# Fig5F
  
  corrLM(data=obj_CAF@meta.data,x_col="Risk_score",y_col="activated_CAF1",color_col = "CellType",
         x_lab="Risk score",y_lab="Activated CAF",title="",x_label=6,y_cor_label=1.8,y_formula_label=1.9)+
    scale_color_manual(values = c(myCAF="indianred2",iCAF= "#45C7D6",vCAF= "#A26DB7"))
  

  
# TECs
  
  obj_TEC <- subset(object,subset=CellType=='TEC')

  tip_marker <- c('ESM1','NID2','RGCC','RAMP3','HSPG2','PLVAP', 'COL4A1','CD93',
                  'HSPG2','APLNR','INSR','KDR','VWA1','COL4A2','ANGPT2','CXCR4')
  obj_TEC <- AddModuleScore(obj_TEC,features = list(score=tip_marker),name =c("Tip"))
  
  GPAG <- c("ANG","C3","CX3CR1","CDC42","TGFBR2","CTNNB1","CCR2","CAV1",
            "PTGS2","SRPK2","IL1B","CYR61","KLF4","CASP8","SCG2","EGR3","JUN",
            "TNFSF12","NR4A1","KRIT1","CYSLTR1","FGF10","HHEX","ROCK1","GPX1",
            "EMCN","S1PR1","ECSCR","PLCD3")
  GPAG <- intersect(GPAG,rownames(obj_TEC))
  PPAG <- c("ADAM15","GPR56","ERBB2","VEGFA","GPI","MYH9","SPHK1","COL4A1","MED1",
            "RBM15","ADAM8","CXCL17","IL18","WNT7B","EPHB3","ROCK2","SHB","ITGA5",
            "COL4A2","ANGPT2","E2F7","ADM2","NAA15","TNFRSF1A","PGF","PLXDC1","TGFBR1")
  PPAG <- intersect(PPAG,rownames(obj_TEC))
  
  TEC_normalization <- apply(LayerData(obj_TEC,"data"), 2, function(x){sum(x[GPAG])-sum(x[PPAG])})
  obj_TEC$vessel_normalization <- TEC_normalization
  
  
  
# Fig5G
  
  diff_tec <- diff_calculate(obj_TEC,ident_1="High_risk",group="Risk_type",pct.thres=0.1,log2FC.thres=0.25,p.thres=0.05)
  diff_tec[gene %in% intersect(diff_tec[type=='Up',gene],tip_marker),label:=gene]

  ggplot(diff_tec,aes(avg_log2FC,-log10(p_val_adj),color=type))+
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    geom_vline(xintercept = log2(1.2), linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    geom_vline(xintercept = -log2(1.2), linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    geom_point(size=0.5)+  
    scale_color_manual(values=c(NS="azure2",Up="#7AB774",Down="#7f91c9"))+
    ggrepel::geom_text_repel(aes(label = label),color="red",
                             size=2,segment.size=0.02,segment.color='red',max.overlaps=10000)+
    theme_classic()+xlab("Log2(Fold Change))")+
    xlim(-2,2.5)+
    ylab("-Log10(adj.P.Value)")+ggtitle('TEC cells High_risk vs. Low_risk')

  
  
  
  
# Fig5H

  ggplot(obj_TEC@meta.data,aes(x=Risk_type,y=vessel_normalization))+
    geom_boxplot(aes(fill=Risk_type),width=0.5,outlier.shape = NA,linewidth  = 0.235)+
    geom_jitter(fill='white',color='black',shape=21,width =0.1,size=0.4,stroke = 0.1)+
    scale_fill_manual(values=Risk_type.colors)+
    ylab('Vessel normalization')+
    xlab('')+ggtitle('TEC cells')+
    theme_classic() +
    stat_compare_means(comparisons=list(c('High_risk','Low_risk')),
                       method ="wilcox.test",size=2,bracket.size= 0.235)
  
  ggplot(obj_TEC@meta.data,aes(x=Risk_type,y=Tip1))+
    geom_boxplot(aes(fill=Risk_type),width=0.5,outlier.shape = NA,linewidth  = 0.235)+
    geom_jitter(fill='white',color='black',shape=21,width =0.1,size=0.4,stroke = 0.1)+
    scale_fill_manual(values=Risk_type.colors)+
    ylab('Tip score')+
    xlab('')+ggtitle('TEC cells')+
    theme_classic() +
    stat_compare_means(comparisons=list(c('High_risk','Low_risk')),
                       method ="wilcox.test",size=2,bracket.size= 0.235)
  
  






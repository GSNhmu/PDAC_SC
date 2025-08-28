

# object_pdac <- readRDS("E:/PAAD_Project/Data/Tumor/object_pdac.rds")


# Fig6AB
  
  colour_bk <- c(colorRampPalette(c("cyan3", alpha("cyan3",0.2)))(15),  
                 colorRampPalette(c(alpha("cyan3",0.2), "linen"))(20),
                 colorRampPalette(c("linen", "lightsalmon"))(40), 
                 colorRampPalette(c("lightsalmon", "brown3"))(40))
  FeaturePlot(object_pdac,"Risk_score",raster=T,pt.size=2)+NoAxes()+
    scale_colour_gradientn(colours = colour_bk)

  colour_bk <- c(colorRampPalette(c("seagreen3", alpha("lightgreen",0.2)))(50),  
                 colorRampPalette(c(alpha("lightgreen",0.2), "linen"))(25),
                 colorRampPalette(c("linen", "lightsalmon"))(20), 
                 colorRampPalette(c("lightsalmon", "tomato2"))(40))
  FeaturePlot(object_pdac,"Moffitt_t",raster=T,pt.size=2)+NoAxes()+
    scale_colour_gradientn(colours = colour_bk)

  phenoPropotion(object_pdac@meta.data,x_name = "Type",y_name="Risk_type",
                 cols=Risk_type.colors,
                 legend_name="Risk type",bar_width = 0.6)+
    coord_polar(theta = "y", start = 0)+NoAxes()
  theme(axis.text.x = element_blank())

  phenoPropotion(object_pdac@meta.data,x_name = "Type",y_name="Moffitt_subtype",
                 cols=Moffit_type.colors,
                 legend_name="Risk type",bar_width = 0.75)+
    coord_polar(theta = "y", start = 0)+NoAxes()
  theme(axis.text.x = element_blank()) 

  
# Fig6C
  
  plot_data <- data.table(object_pdac@meta.data)
  
  cal_plot_data <- function(plot_data){
  
    num_A <- plot_data[Risk_type=='High_risk'& Moffitt_subtype=='Basal-like',.N]
    num_B <- plot_data[Risk_type=='High_risk'& Moffitt_subtype=='Classical',.N]
    num_C <- plot_data[Risk_type=='Low_risk'& Moffitt_subtype=='Classical',.N]
    num_D <- plot_data[Risk_type=='Low_risk'& Moffitt_subtype=='Basal-like',.N]
    
    count_mat <- matrix(c(num_B,num_C,num_A,num_D),nrow = 2)
    rownames(count_mat) <- c('High risk','Low risk')
    colnames(count_mat) <- c('Classical','Basal like')
    
    width_res <- c(num_A,num_B,num_C,num_D)/plot_data[,.N]
    sub_plot_data <- data.frame(
      Category = c("A", "B", "C", "D"),
      X = c(0, 0-width_res[2], 0-width_res[3], 0),  
      Y = c(0, 0, 0-width_res[3], 0-width_res[4]),  
      Width = width_res,  
      Height = width_res  
    )
    return(list(count_mat,sub_plot_data))
  }
  
  res_lit_P <- cal_plot_data(plot_data[Type=='P',])
  res_lit_M <- cal_plot_data(plot_data[Type=='M',])
  
  ggplot(res_lit_P[[2]]) +
    geom_rect(aes(xmin = X, xmax = X + Width, ymin = Y, ymax = Y + Height, fill = Category)) +
    geom_vline(xintercept = 0, color = "black",linewidth  = 1)+  
    geom_hline(yintercept = 0, color = "black",linewidth  = 1)+  
    scale_fill_manual(values = c("A" = "brown2", "B" = "royalblue", "C" = "seagreen3", "D" = "lightpink")) +
    theme_classic()+NoLegend()+NoAxes()+ggtitle('Primary')
  ggplot(res_lit_M[[2]]) +
    geom_rect(aes(xmin = X, xmax = X + Width, ymin = Y, ymax = Y + Height, fill = Category)) +
    geom_vline(xintercept = 0, color = "black",linewidth  = 1)+  
    geom_hline(yintercept = 0, color = "black",linewidth  = 1)+  
    scale_fill_manual(values = c("A" = "brown2", "B" = "royalblue", "C" = "seagreen3", "D" = "lightpink")) +
    theme_classic()+NoLegend()+NoAxes()+ggtitle('Metastatic')
  
  
# Fig6D
  
  DEGs_Type <- diff_calculate(object_pdac,ident_1="M",group="Type",log2FC.thres=0.5,p.thres=0.01)
  DEGs_Risk <- diff_calculate(object_pdac,ident_1="High_risk",group="Risk_type",log2FC.thres=0.5,p.thres=0.01)
  DEGs_Moffitt <- diff_calculate(object_pdac,ident_1="Basal-like",group="Moffitt_subtype",log2FC.thres=0.5,p.thres=0.01)
  
  colnames(DEGs_Type)[-1] <- paste0(colnames(DEGs_Type)[-1],"-","Type")
  colnames(DEGs_Risk)[-1] <- paste0(colnames(DEGs_Risk)[-1],"-","Risk")
  colnames(DEGs_Moffitt)[-1] <- paste0(colnames(DEGs_Moffitt)[-1],"-","Moffitt")
  
  #16489
  df <- Reduce("inner_join",list(DEGs_Type,DEGs_Risk,DEGs_Moffitt))
  
  surv_meta <- meta_result[gene_type %in% c("PoorPrognosis","FavorablePrognosis")]
  setorder(surv_meta,-"Z")
  
  up_label <- c(intersect(df[`type-Type` == "Up" & `type-Risk` == "Up" & `type-Moffitt` == "Up", gene],head(surv_meta$gene_name,50)),"HMGA1")
  down_label <- intersect(df[`type-Type` == "Down" & `type-Risk` == "Down" & `type-Moffitt` == "Down", gene],tail(surv_meta$gene_name,100))
  
  df$label <- df$gene
  df$label[!(df$label %in% c(up_label,down_label))] <- NA
  df$label[df$label %like% "^AC" | df$label %like% "^AL" | df$label %like% "^LINC"] <- NA
  
  
  ggplot(data=df, aes(x = `avg_log2FC-Type`, y = `avg_log2FC-Risk`, color=`avg_log2FC-Moffitt`,label=label)) +
    geom_point_rast(stroke = 0, alpha = 1,shape = 16,size=0.5) + 
    geom_hline(yintercept=0, linetype="dashed", color = "gold")+
    geom_vline(xintercept=0, linetype="dashed", color = "gold")+
    geom_text_repel(color="black",size=1.5,segment.size=0.1,force= 0.1,force_pull= 1,point.padding=NA,max.overlaps=10000) +
    scale_color_gradient2(midpoint=0, low="seagreen3", mid="ivory",high="tomato2",space ="Lab") +
    ylab("Log2FC of Metastasis vs. Primay")+
    xlab("Log2FC of Adverse vs. Favorable")+
    labs(color = "Log2FC of Basal-like vs. Classical")+
    ylim(c(-5,8))+xlim(-10,10)+
    theme_classic()
  
  
  
  common_up <-  df[`type-Type` == "Up" & `type-Risk` == "Up" & `type-Moffitt` == "Up", gene]
  common_down <- df[`type-Type` == "Down" & `type-Risk` == "Down" & `type-Moffitt` == "Down", gene]
  
  # Up - Hallmark
  Hset <- msigdbr(species = "Homo sapiens",category="H") %>%
    dplyr::select(gs_name, gene_symbol)
  Hset$gs_name <- str_remove(Hset$gs_name,"HALLMARK_")
  Hset$gs_name <- str_replace_all(Hset$gs_name,"_"," ")
  Hset$gs_name <- str_to_sentence(Hset$gs_name)
  
  enrich_up <- enricher(common_up,TERM2GENE=Hset)
  enrichplot::dotplot(enrich_up, showCategory = 20) # 观察
  
  plot_data <- data.table(enrich_up@result[1:6,])
  plot_data[,q:=-log10(qvalue)]
  
  ggplot(plot_data,aes(x = q, y = reorder(q,Description)))+
    geom_col(alpha=0.7,fill = "firebrick3",width=0.95)+
    geom_text(aes(x = 0.05, label = Description),size = 2,hjust = 0)+
    labs(x = '-log10(qvalue)', y = ' ')+
    scale_x_continuous(breaks=c(0,10,20,30))+
    theme_classic()
  
  
  # Down - KEGG
  load("E:/PAAD_Project/Data/Data_input/KEGGSets.rds")
  KEGGset <-rbindlist(lapply(names(KEGGSets), function(name) {
    data.table(gs_name = name, gene_symbol = KEGGSets[[name]])}))
  
  enrich_down_k <- enricher(common_down,TERM2GENE=KEGGset)
  plot_data <- data.table(enrich_down_k@result[c(1,2,7),])
  plot_data[,q:=-log10(qvalue)]
  
  ggplot(plot_data,aes(x = q, y = reorder(q,Description)))+
    geom_col(alpha=0.7,fill = "royalblue2",width=0.95)+
    geom_text(aes(x = 0.05, label = Description),size = 2,hjust = 0)+
    labs(x = '-log10(qvalue)', y = ' ')+
    scale_x_continuous(breaks=c(0,1,2,3))+
    theme_classic()
  
  
  
# Fig6E
  
  expr <- as.matrix(LayerData(object_pdac,layer="counts"))
  expr <- matr.filter(expr, 
                      min.cells = 10,
                      min.genes = 10)

  rogue_risk_type <- rogue(expr,
                           labels = object_pdac$Risk_type,
                           samples = object_pdac$orig.ident,
                           platform = "UMI",
                           span = 0.6)
  rogue_risk_type_data <- rogue_risk_type
  rogue_risk_type_data$orig.ident <- row.names(rogue_risk_type_data)
  rogue_risk_type_data <- pivot_longer(rogue_risk_type_data,cols=c("High_risk","Low_risk"),names_to = "Risk_type",values_to = "ROGUE")
  rogue_risk_type_data$Risk_type <- factor(rogue_risk_type_data$Risk_type,levels=c("Low_risk","High_risk"))
  ggplot(rogue_risk_type_data, aes(x = Risk_type, y = ROGUE)) +
    geom_boxplot(width=0.5,color="gray40",fill="ivory") +  
    geom_point(aes(color = orig.ident), size = 2) +  
    geom_line(aes(color = orig.ident,group = orig.ident), size = 0.8, linetype = "dashed") +  
    scale_color_manual(values=orig.ident.colors)+
    theme_classic()+stat_compare_means(size=2.5)
  
  
# Fig6F
  
  ggplot(object_pdac@meta.data, aes(x = SCENT, y = Risk_type, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
    stat_density_ridges(geom = "density_ridges_gradient", scale = 1.2, calc_ecdf = TRUE) +
    scale_fill_viridis_c(option = "H",direction=1,name = "Tail probability")+
    xlab("SCENT")+ylab("")+
    theme_base()
  ggplot(object_pdac@meta.data, aes(x = SCENT, y = Moffitt_subtype, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
    stat_density_ridges(geom = "density_ridges_gradient", scale = 1.2, calc_ecdf = TRUE) +
    scale_fill_viridis_c(option = "H",direction=1,name = "Tail probability")+
    xlab("SCENT")+ylab("")+
    theme_base()
  
  
  
# FigS5A
  
  phenoPropotion(object_pdac@meta.data,x_name = "Risk_type",y_name="Response",
                 cols=Response.colors,
                 legend_name="Risk type",bar_width = 0.9)
  phenoPropotion(object_pdac@meta.data,x_name = "Moffitt_subtype",y_name="Response",
                 cols=Response.colors,
                 legend_name="Risk type",bar_width = 0.9)
       
                  
# Fig6G  
  
  object_pdac$DC3 <-  Embeddings(object_pdac, "diffmap")[,"DC_3"]

  colour_bk <- c(colorRampPalette(c("cyan3", alpha("cyan3",0.2)))(15),  
                 colorRampPalette(c(alpha("cyan3",0.2), "linen"))(20),
                 colorRampPalette(c("linen", "lightsalmon"))(40), 
                 colorRampPalette(c("lightsalmon", "brown3"))(40))
  plot_ly(object_pdac@meta.data,x = ~DC1,y = ~DC2,z = ~DC3,
                      color = ~Risk_score,colors = colour_bk,
                      marker = list(size = 3))
  
  
  rd <- Embeddings(object_pdac, "diffmap")[,1:3]
  cl <- object_pdac$Risk_type
  pto <- slingshot(rd, cl,start.clus="Low_risk")
  sds <- as.SlingshotDataSet(pto)
  Slingshot_curv <- as.data.frame(sds@curves$Lineage1$s)
  Slingshot_curv$order <- 1:150
  plot_ly(x = Slingshot_curv[,1], y = Slingshot_curv[,2], z = Slingshot_curv[,3], 
          mode = "lines",line = list(color = ~Slingshot_curv$order,colorscale = 'Viridis', width=6))
     
                     
  
# FigS5B
  
  cor_data <- object_pdac@meta.data[,c("Pseudotime","nFeature_RNA","CytoTRACE","SCENT")]
  M <- cor(cor_data,method = "spearman")
  testRes <- cor.mtest(cor_data, conf.level = 0.95)
  
  colour_bk <- c(colorRampPalette(c("linen", "lightsalmon"))(15), 
                 colorRampPalette(c("lightsalmon", "brown3"))(10))
  corrplot(M, p.mat = testRes$p, tl.pos = 'd', method ="square",col=colour_bk,col.lim=c(0, 1),
           insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05),
           pch.cex = 0.9, pch.col = 'grey20')
  
  
  
# FigS5C
  
  ggplot(object_pdac@meta.data,aes(x=DC1,y=DC2,color=orig.ident))+
    geom_point_rast(size=0.1)+
    scale_colour_manual(values = orig.ident.colors)
  colour_bk <- c(colorRampPalette(c("seagreen3","lightcyan"))(10),
                 colorRampPalette(c("lightcyan","linen"))(5),
                 colorRampPalette(c("linen","#ffbb78"))(5),
                 colorRampPalette(c("#ffbb78","red3"))(10))
  ggplot(object_pdac@meta.data,aes(x=DC1,y=DC2,color=nFeature_RNA))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk)
  ggplot(object_pdac@meta.data,aes(x=DC1,y=DC2,color=CytoTRACE))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk)
  ggplot(object_pdac@meta.data,aes(x=DC1,y=DC2,color=SCENT))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk)
  
  
  
# Fig6H   
  
  ggplot(object_pdac@meta.data, aes(x=Pseudotime,y=Risk_score))+
    geom_smooth(fullrange=T,method="loess",span=1.5,color="brown3",fill="gray90",linewidth=0.5)+
    ylab("Risk score")
  ggplot(object_pdac@meta.data, aes(x=Pseudotime,y=as.numeric(Moffitt_subtype=="Basal-like")))+
    geom_smooth(fullrange=T,method="loess",span=1.5,color="brown3",fill="gray90",linewidth=0.5)+
    ylab("%Basal-likeness")
  ggplot(object_pdac@meta.data, aes(x=Pseudotime,y=as.numeric(Type=="M")))+
    geom_smooth(fullrange=T,method="loess",span=1.5,color="brown3",fill="gray90",linewidth=0.5)+
    ylab("%Metastasis")
  ggplot(object_pdac@meta.data, aes(x=Pseudotime,y=as.numeric(Response=="PD")))+
    geom_smooth(fullrange=T,method="loess",span=1.5,color="brown3",fill="gray90",linewidth=0.5)+
    ylab("%PD (Response)") +
    theme_classic()
  
  
# Fig6I
  
  Risk_score <- object_pdac$Risk_score
  Regulons <- object_pdac@meta.data[,names(object_pdac@meta.data) %like% "\\(\\+\\)$"]
  
  cor_risk_regulon <- lapply(1:ncol(Regulons),function(i){
    cor_data <- cor.test(Risk_score,Regulons[,i],method="spearman")
    data.table(Regulon=names(Regulons)[i],cor=cor_data$estimate,p=cor_data$p.value)})
  cor_risk_regulon <- list.rbind(cor_risk_regulon)
  
  setorder(cor_risk_regulon,"cor")
  cor_risk_regulon[,rank:=1:nrow(cor_risk_regulon)]
  
  cor_risk_regulon[cor > 0 & p < 0.01,label:="Pos_Cor"] 
  cor_risk_regulon[Regulon %in% tail(cor_risk_regulon$Regulon,10),highlight:=Regulon] 
  
  ggplot(data=cor_risk_regulon,aes(x=rank,y=cor,color=label,label=highlight))+
    scale_color_manual(values=c(Pos_Cor="red3"),na.value = "gray70")+
    geom_point(alpha = 0.8,shape = 16,size=1)+
    theme_bw() +
    geom_label_repel(size=3,max.overlaps = 100000,box.padding = 1, )+
    ylab("Correlation coefficients")+
    NoLegend()+
    theme(title=element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size=8),
          axis.text.x = element_text(size = 7,angle=0,colour = 'black'),
          axis.text.y = element_text(size = 6,colour = 'black'),
          axis.line = element_line(linewidth  = 0.235),
          axis.ticks = element_line(linewidth  = 0.235))
  ggsave("E:/PAAD_Project/Data/Tumor/Risk_SCENIC.pdf",width=8,height=8,units = "cm")
  
  ggrepel::geom_text_repel(aes(label = label),color="red",
                           size=2,segment.size=0.02,segment.color='red',max.overlaps=10000)+
  
  
# Fig6JK
  
  colour_bk <- c(colorRampPalette(c("seagreen3","lightcyan"))(10),
                 colorRampPalette(c("lightcyan","linen"))(5),
                 colorRampPalette(c("linen","#ffbb78"))(5),
                 colorRampPalette(c("#ffbb78","red3"))(10))
  ggplot(object_pdac@meta.data,aes(x=DC1,y=DC2,color=`HMGA1(+)`))+
    geom_point_rast(size=0.1)+
    scale_colour_gradientn(colours = colour_bk)+
    theme_classic()+NoAxes()
  
  ggplot(object_pdac@meta.data, aes(x=Pseudotime,y=`HMGA1(+)`))+
    geom_smooth(fullrange=T,method="loess",span=1.5,color="red3",fill="gray90")+
    ylab("HMGA1(+)") +#xlim(0.15,0.75)+
    theme_classic()
  
  colour_bk <- c(colorRampPalette(c("seagreen3","lightcyan"))(50),
                 colorRampPalette(c("lightcyan","linen"))(10),
                 colorRampPalette(c("linen","#ffbb78"))(20),
                 colorRampPalette(c("#ffbb78","red2"))(20))
  ggplot(object_pdac@meta.data,aes(x=SCENT,y=`HMGA1(+)`,color=Moffitt_t))+
    geom_point_rast(alpha=0.6,size =0.1)+
    scale_colour_gradientn(colours = colour_bk)+
    geom_smooth(fullrange=T,method="loess",span=1,color="purple3",size=1)+
    theme_classic()
  
  
# FigS5D
  
  VlnPlot_gene <- function(obj,feature,group,cols){
    VlnPlot(obj,features=feature,group.by=group,
            cols=cols,pt.siz=0)+xlab("")+ylab(feature)+ggtitle("")+
      theme_classic()+
      stat_summary(fun.y=mean, geom="point", shape=18,size=1.5, color="white")+
      stat_compare_means(size=2)}
  
  VlnPlot_gene(obj=object_pdac,feature="HMGA1",group="Type",cols=Type.colors)
  VlnPlot_gene(obj=object_pdac,feature="HMGA1(+)",group="Type",cols=Type.colors)
  
  VlnPlot_gene(obj=object_pdac,feature="HMGA1",group="Moffitt_subtype",cols=Moffit_type.colors)
  VlnPlot_gene(obj=object_pdac,feature="HMGA1(+)",group="Moffitt_subtype",cols=Moffit_type.colors)
  
  VlnPlot_gene(obj=object_pdac,feature="HMGA1",group="Risk_type",cols=Risk_type.colors)
  VlnPlot_gene(obj=object_pdac,feature="HMGA1(+)",group="Risk_type",cols=Risk_type.colors)
  
  VlnPlot_gene(obj=object_pdac[,!is.na(object_pdac$Response)],feature="HMGA1",group="Response",cols=Response.colors)
  VlnPlot_gene(obj=object_pdac[,!is.na(object_pdac$Response)],feature="HMGA1(+)",group="Response",cols=Response.colors)
  
  
  
# FigS5E
  
  hallmarks <- names(object_pdac@meta.data)[names(object_pdac@meta.data) %like% "HALLMARK_"]
  cor_data <- object_pdac@meta.data[,c("HMGA1(+)",hallmarks)]
  
  M <- cor(cor_data,method = "spearman")
  testRes <- cor.mtest(cor_data, conf.level = 0.95)

  filtered_M <- M[(M[1,] > 0.1 & testRes$p[1,] < 0.05),(M[1,] > 0.1 & testRes$p[1,] < 0.05)]
  filtered_p <- testRes$p[(M[1,] > 0.1 & testRes$p[1,] < 0.05),(M[1,] > 0.1 & testRes$p[1,] < 0.05)]
  

  corrplot(filtered_M, p.mat = filtered_p, tl.pos = 'r', method ="square",col=rev(COL2("RdYlBu", 200)),col.lim=c(-1, 1),type="lower",
           insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05),
           pch.cex = 0.3,tl.cex = 0.7,tl.col = "black", pch.col = 'grey20')

  
# FigS5F
  
  FeaturePlot(object,features = "HMGA1",pt.size = 1.8,raster = T,order=T)+
    scale_color_gradientn(colours = c(
      "lightcyan",
      colorRampPalette(c("#f7f7f7","#fdf2b5"))(10),
      colorRampPalette(c("#fdf2b5","#F67B51"))(10),
      colorRampPalette(c("#F67B51","#A30023"))(20)))+
    theme(title=element_text(size = 8),
          legend.text=element_text(size = 5))
  
  object$CellType_HMGA1 <- paste0(object$CellType_Main)
  object$CellType_HMGA1 <- str_remove_all(object$CellType_HMGA1," cell")
  object$CellType_HMGA1[object$CellType == "Malignant"] <- "Epi_Malignant"
  object$CellType_HMGA1[object$CellType_HMGA1 == "Epithelial" & object$CellType != "Malignant"] <- "Epi_NonMalignant"
  CellType_HMGA1.colors <- c(Epi_NonMalignant="steelblue",Epi_Malignant="orangered3",`Myeloid`="gold2",
                             Fibroblast="#8c564b",`Endothelial`="#ff9896",Lymphocyte="seagreen3")
  object$CellType_HMGA1 <- factor(object$CellType_HMGA1,levels = names(CellType_HMGA1.colors))
  
  VlnPlot(object,"HMGA1",group.by = "CellType_HMGA1",pt.size=0,cols=CellType_HMGA1.colors)
  
  
# Fig6M
  
  countexp.Seurat <- sc.metabolism.Seurat(obj = object_pdac, method = "AUCell", imputation = F, ncores = 2, metabolism.type = "KEGG")
  metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score
  countexp.Seurat$`Retinol metabolism` <- as.numeric(metabolism.matrix["Retinol metabolism",])
  
  cor_data <- countexp.Seurat@meta.data[,c("Retinol metabolism","Risk_score","Moffitt_t", "HMGA1", "HMGA1(+)","CytoTRACE","SCENT","nFeature_RNA")]
  names(cor_data)[2:3] <- c("Risk score","Basal-likeness score")
  M <- cor(cor_data,method = "spearman")
  testRes <- cor.mtest(cor_data, conf.level = 0.95)
  
  colour_bk <- c(colorRampPalette(c("seagreen3","lightcyan"))(10),
                 colorRampPalette(c("#FFBB7899","red3"))(10))
  corrplot(M, p.mat = testRes$p,method ="square",col=colour_bk,col.lim=c(-1, 1),
           tl.pos = 'd',tl.cex = 0.8,tl.col = "black",
           insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05),
           pch.cex = 0.8, pch.col = 'grey20')

  
  
  
  
  
  
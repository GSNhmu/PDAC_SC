
# cellchat

  cell_nums <- as.matrix(table(object$CellType,object$Risk_type)) 
  CellTypes <- rownames(cell_nums)[rowSums(cell_nums>200) == 2]
  CellTypes <- setdiff(CellTypes,"Ductal")
  
  object_cc <- subset(object,subset=CellType %in% CellTypes)
  
  object_High_risk <- object_cc[,object_cc$Risk_type == "High_risk"]
  data.input <- GetAssayData(object_High_risk,"data")
  labels <- object_High_risk$CellType
  meta <- data.frame(labels = labels, row.names = names(labels))
  meta$samples <- factor(object_High_risk$orig.ident,levels = c(paste0("P",1:4),paste0("M",1:8)))
  cellchat_High_risk <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  
  object_Low_risk <- object_cc[,object_cc$Risk_type == "Low_risk"]
  data.input <- GetAssayData(object_Low_risk,"data")
  labels <- object_Low_risk$CellType
  meta <- data.frame(labels = labels, row.names = names(labels))
  meta$samples <- factor(object_Low_risk$orig.ident,levels = c(paste0("P",1:4),paste0("M",1:8)))
  cellchat_Low_risk <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  
  CellChatDB <- CellChatDB.human
  showDatabaseCategory(CellChatDB)
  dplyr::glimpse(CellChatDB$interaction)
  CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key = "annotation")
  
  object.list <- list(Low_risk = cellchat_Low_risk,High_risk = cellchat_High_risk)
  object.list <- lapply(object.list,function(cellchat){
    
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    
    return(cellchat)
  })
  names(object.list) <- c("Low_risk","High_risk")

  
  
# Fig7AB 
  
  load("E:/PAAD_Project/Data/Cellchat/cellchat_object.list.rda")
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))

  compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",color.use = c("cyan3","brown3"))+
    scale_x_discrete(labels=c('Low risk','High risk'))
                           
  cellchat_color <- c(Malignant="orangered2",
                      CD4T="seagreen2",
                      CD8T="royalblue1",
                      NK="salmon",    
                      Macrophage="gold",
                      cDC="#9EDAE5",
                      Neutrophil="#b83570",
                      TEC="brown3",
                      myCAF="indianred2")

  mat_low <- object.list[[1]]@net$weight
  diag(mat_low) <- 0 
  mat_low <- mat_low[names(cellchat_color),names(cellchat_color)]
  
  mat_hi <- object.list[[2]]@net$weight
  diag(mat_hi) <- 0 
  mat_hi <- mat_hi[names(cellchat_color),names(cellchat_color)]
  
  netVisual_circle(mat_low, color.use = cellchat_color,
                   weight.scale = F, edge.width.max = 2,edge.weight.max = max(c(mat_low,mat_hi)),
                   label.edge= F, title.name = "Low risk interaction strength")
  
  netVisual_circle(mat_hi, color.use = cellchat_color,
                   weight.scale = F, edge.width.max = 2,edge.weight.max = max(c(mat_low,mat_hi)),
                   label.edge= F, title.name = "High risk interaction strength")
  




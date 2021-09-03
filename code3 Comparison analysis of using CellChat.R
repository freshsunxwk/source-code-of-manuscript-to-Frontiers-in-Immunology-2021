##############xwk###############
#Comparison analysis of using CellChat
library(CellChat)
library(patchwork)
DimPlot(immune.combined, reduction = "umap",label = F)#带有数字标记的分群,即F1的p2
Idents(immune.combined) <- "stim"
#A_ISO_SP B_ALLO_SP C_ISO_HT D_ALLO_HT
immune.combined1 <- subset(immune.combined, idents =c("C_ISO_HT"))#get subset
Idents(immune.combined1) <- "fenzhu"
cellchat.NL <- createCellChat(object = immune.combined1, group.by = "ident")

immune.combined2 <- subset(immune.combined, idents =c("D_ALLO_HT"))#get subset
Idents(immune.combined2) <- "fenzhu"
cellchat.LS <- createCellChat(object = immune.combined2, group.by = "ident")

cellchat.C_ISO_HT=cellchat.NL
cellchat.D_ALLO_HT=cellchat.LS
#####1#####

cellchat=cellchat.C_ISO_HT
##设置参考数据库
# 选择合适的物种，可选CellChatDB.human, CellChatDB.mouse
CellChatDB <- CellChatDB.mouse  
# 使用"Secreted Signaling"用于细胞通讯分析
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# table(CellChatDB$interaction$annotation)
# Cell-Cell Contact       ECM-Receptor Secreted Signaling 
# 378                432               1211 
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") 
# 将数据库传递给cellchat对象
cellchat@DB <- CellChatDB.use 
##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.mouse)
##推测细胞通讯网络
cellchat <- computeCommunProb(cellchat)################
#画树状图的code
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

#cellchat.NL=cellchat
cellchat.C_ISO_HT=cellchat

########2#####
cellchat=cellchat.D_ALLO_HT
##设置参考数据库
# 选择合适的物种，可选CellChatDB.human, CellChatDB.mouse
CellChatDB <- CellChatDB.mouse  
# 使用"Secreted Signaling"用于细胞通讯分析
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# table(CellChatDB$interaction$annotation)
# Cell-Cell Contact       ECM-Receptor Secreted Signaling 
# 378                432               1211 
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") 
# 将数据库传递给cellchat对象
cellchat@DB <- CellChatDB.use 
##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.mouse)
##推测细胞通讯网络
cellchat <- computeCommunProb(cellchat)
#画树状图的code
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")


cellchat.D_ALLO_HT=cellchat

#################
save(cellchat.D_ALLO_HT,cellchat.C_ISO_HT,file="20210831-2.Rdata")
object.list <- list(C_ISO= cellchat.C_ISO_HT, D_ALLO = cellchat.D_ALLO_HT)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat

#> An object of class CellChat created from a merged object with multiple datasets 
#>  555 signaling genes.
#>  7563 cells.
#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
#Compare the number of interactions and interaction strength among different cell populations
#Differential number of interactions or interaction strength among different cell populations
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

netVisual_heatmap(cellchat)

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#####XXXXX#######
#Differential number of interactions or interaction strength among different cell types
# group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
# group.cellType <- factor(group.cellType, levels = c("FIB", "DC", "TC"))
# table(cellchat@idents$joint)
group.cellType <- c("T cells","B cells","NK cells","Neutrophils","Monocytes","DC")
group.cellType <- factor(group.cellType, levels = c("T cells","B cells","NK cells","Neutrophils","Monocytes","DC"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)
#####XXXXX########
#Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

#Part II: Identify the conserved and context-specific signaling pathways
#Identify signaling networks with larger (or less) difference as well as signaling groups based on their functional/structure similarity

#Identify signaling groups based on their functional similarity

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)###?
#Error in do_one(nmeth) : 外接函数调用时不能有NA/NaN/Inf(arg1)

#> 2D visualization of signaling networks from datasets 1 2
# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
#> 2D visualization of signaling networks from datasets 1 2
#Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(cellchat, type = "functional")
#> Compute the distance of signaling networks between datasets 1 2
#Identify and visualize the conserved and context-specific signaling pathways
#Compare the overall information flow of each signaling pathway
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
#Compare outgoing (or incoming) signaling associated with each cell population
library(ComplexHeatmap)
#> Loading required package: grid
#> ========================================
#> ComplexHeatmap version 2.7.1.1010
#> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
#> Github page: https://github.com/jokergoo/ComplexHeatmap
#> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
#> 
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
#>   genomic data. Bioinformatics 2016.
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================
i = 2
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

i = 1
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+2], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2+ ht3, ht_gap = unit(0.5, "cm"))

i = 1
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+2], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2+ ht3, ht_gap = unit(0.5, "cm"))

i = 1
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+2], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2+ ht3, ht_gap = unit(0.5, "cm"))

table(cellchat@idents$joint)
table(cellchat@idents$NL)
table(cellchat@idents$LS)#
# T cells     B cells Neutrophils          DC    NK cells   Monocytes 
# 1294        2020          39          83         138         111
# groupSize <- as.numeric(table(cellchat@idents)) # 后面有用

#Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:6),  comparison = c(1, 2,3), angle.x = 45)
netVisual_bubble(cellchat, sources.use = c(1:6), targets.use = 1,  comparison = c(1,2,3), angle.x = 45)

#> Comparing communications on a merged object
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:6),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:6),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
#
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
#pos.dataset = "LS"  #最后面names(object.list)[2]))
pos.dataset = "NL"  #最后面names(object.list)[1]))

pos.dataset = "HEART"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "WT",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "HT",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]

#2改成1 了  ，要改回来，最后面names(object.list)[2]))
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1:6), targets.use = c(1:6), comparison = c(1, 2,3),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1:6), targets.use = c(1:6), comparison = c(1, 2,3),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
gg1 + gg2
# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Note: The first link end is drawn out of sector 'MIF'.
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#XXXXX#
#Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
cellchat@netP$HEART$pathways #
pathways.show <- c("ICOS") #
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(3,1), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
# Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
  #> Note: The first link end is drawn out of sector 'Inflam. FIB'.
  # Chord diagram
  group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(object.list[[1]]@idents)
pathways.show <- c("MHC-II")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.
par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 1, targets.use = c(1:6), lab.cex = 0.5, title.name = paste0("Signaling from Inflam.FIB - ", names(object.list)[i]))
}
# show all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1), targets.use = c(1:6),slot.name = "netP", title.name = paste0("Signaling pathways sending from T - ", names(object.list)[i]), legend.pos.x = 10)
}

i=1
netVisual_chord_gene(object.list[[i]], sources.use = c(1), targets.use = c(1:6),slot.name = "netP", title.name = paste0("Signaling pathways sending from T - ", names(object.list)[i]), legend.pos.x = 10)
# [1] "ALCAM"       "APP"         "BST2"        "CD200"      
# [5] "CD22"        "CD226"       "CD39"        "CD45"       
# [9] "CD48"        "CD6"         "CD80"        "CD86"       
# [13] "CD96"        "CDH"         "CDH1"        "CLEC"       
# [17] "ICAM"        "ICOS"        "ITGAL-ITGB2" "JAM"        
# [21] "LAIR1"       "LCK"         "MHC-I"       "MHC-II"     
# [25] "PD-L1"       "PDL2"        "PVR"         "SELPLG"     
# [29] "SEMA4"       "THY1"        "TIGIT"  
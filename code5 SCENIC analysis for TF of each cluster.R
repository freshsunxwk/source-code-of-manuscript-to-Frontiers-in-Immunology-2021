##############xwk###############
#analysis for TF
#SCENIC analysis for TF of each cluster
library(ComplexHeatmap)
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(arrow)
library(SCopeLoomR)
loom <- open_loom(scenicLoomPath)
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 
cellInfo <- cellInfo[,c("sample","nGene","nUMI","cluster","CellType")]
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)
                                     regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
                                     regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                                                          function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
                                     regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
                                     #regulons
                                     my.regulons=as.character(rssPlot$plot$data$Topic)
                                     my.regulons <- gsub('\\(\\+\\)','',my.regulons)
                                     myAUCmatrix <- t(AUCmatrix[rownames(AUCmatrix)%in%my.regulons,])
                                     pheatmap(myAUCmatrix, show_colnames=F, annotation_col=CellType)
                                     
                                     
                                     #rssPlot
                                     rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"])
                                     rssPlot <-plotRSS(rss,labelsToDiscard = NULL,
                                                       zThreshold = 1, cluster_columns = FALSE, order_rows = TRUE,
                                                       trh = 0.01,varName = "cellType",col.low = "grey90",
                                                       col.mid = "darkolivegreen3",
                                                       col.high = "darkgreen",
                                                       revCol = FALSE)
                                     rssPlot$plot 
                                     
                                     
                                     
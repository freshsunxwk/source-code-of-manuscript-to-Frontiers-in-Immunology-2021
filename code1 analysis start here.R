##############xwk###############
#0  library packages
# library(Seurat)
# library(SingleR)
# library(ggplot2)
# library(clusterProfiler)
# library(org.Mm.eg.db)
# library(CellChat)

#1 read data
#1.1
A_ISO_SP <- Read10X("ISO_SP/")
B_ALLO_SP <- Read10X("ALLO_SP/")
C_ISO_HT <- Read10X("GSE142564/")
D_ALLO_HT<- Read10X("ALLO_HT/")
#1.2
A_ISO_SP <- CreateSeuratObject(counts = A_ISO_SP, project = "IMMUNE_stim", min.cells = 5)
A_ISO_SP$stim <- "A_ISO_SP"
A_ISO_SP <- subset(x = A_ISO_SP, subset = nFeature_RNA > 500)
A_ISO_SP <- NormalizeData(object = A_ISO_SP, verbose = FALSE)
A_ISO_SP <- FindVariableFeatures(object = A_ISO_SP, selection.method = "vst", nfeatures = 2000)
#1.3
B_ALLO_SP <- CreateSeuratObject(counts = B_ALLO_SP, project = "IMMUNE_stim", min.cells = 5)
B_ALLO_SP$stim <- "B_ALLO_SP"
B_ALLO_SP <- subset(x = B_ALLO_SP, subset = nFeature_RNA > 500)
B_ALLO_SP <- NormalizeData(object = B_ALLO_SP, verbose = FALSE)
B_ALLO_SP <- FindVariableFeatures(object = B_ALLO_SP, selection.method = "vst", nfeatures = 2000)
#1.4
C_ISO_HT <- CreateSeuratObject(counts = C_ISO_HT, project = "IMMUNE_stim", min.cells = 5)
C_ISO_HT$stim <- "C_ISO_HT"
C_ISO_HT <- subset(x = C_ISO_HT, subset = nFeature_RNA > 500)
C_ISO_HT <- NormalizeData(object = C_ISO_HT, verbose = FALSE)
C_ISO_HT <- FindVariableFeatures(object = C_ISO_HT, selection.method = "vst", nfeatures = 2000)
#1.5
D_ALLO_HT <- CreateSeuratObject(counts = D_ALLO_HT, project = "IMMUNE_stim", min.cells = 5)
D_ALLO_HT$stim <- "D_ALLO_HT"
D_ALLO_HT <- subset(x = D_ALLO_HT, subset = nFeature_RNA > 500)
D_ALLO_HT <- NormalizeData(object = D_ALLO_HT, verbose = FALSE)
D_ALLO_HT <- FindVariableFeatures(object = D_ALLO_HT, selection.method = "vst", nfeatures = 2000)

#2. Perform integration
#2.1find integration anchors
immune.anchors <- FindIntegrationAnchors(object.list = list(A_ISO_SP,B_ALLO_SP,C_ISO_HT,D_ALLO_HT), dims = 1:20)
#2.2 integration data
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

#3.1 not integrated 
immune.combined.not <- merge(A_ISO_SP, y = c(B_ALLO_SP,C_ISO_HT,D_ALLO_HT))
immune.combined.not <- FindVariableFeatures(object = immune.combined.not, selection.method = "vst", nfeatures = 2000)
immune.combined.not <- ScaleData(immune.combined.not, verbose = FALSE)
immune.combined.not <- RunPCA(immune.combined.not, npcs = 20, verbose = FALSE)
immune.combined.not <- RunUMAP(immune.combined.not, reduction = "pca",dims = 1:20,seed.use = 1 )
table(immune.combined.not@active.ident)
Idents(immune.combined.not) <- immune.combined$stim
immune.combined.not[["percent.mt"]] <- PercentageFeatureSet(immune.combined.not, pattern = "^mt-")
DimPlot(immune.combined.not, reduction = "umap",group.by ="stim",cols = col)
VlnPlot(immune.combined.not, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3, split.by = "stim",cols = col)
VlnPlot(immune.combined.not, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),cols = col,pt.size = 0)
#3.2. Perform an integrated analysis  
immune.combined <- FindVariableFeatures(object = immune.combined, selection.method = "vst", nfeatures = 2000)
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 20, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca",dims = 1:20,seed.use =3)
DefaultAssay(immune.combined) <- "integrated"

#4.Visualize QC metrics as a violin plot
immune.combined[["percent.mt"]] <- PercentageFeatureSet(immune.combined, pattern = "^mt-")
VlnPlot(immune.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "stim",cols = col,pt.size = 0)
immune.combined <- subset(immune.combined,subset = nFeature_RNA > 500 &nFeature_RNA < 5000 &percent.mt < 5)

#5.DimPlot
# col <-c( '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3',
#          '#E95C59', '#E59CC4', '#AB3282','#BD956A', '#9FA3A8', 
#          '#E0D4CA',  '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
#          '#AA9A59',  '#E39A35',  '#6778AE','#91D0BE')
DimPlot(immune.combined, reduction = "umap",label = F,cols = col)
DimPlot(immune.combined, reduction = "umap", split.by = "stim",label = F,cols = col)

#6.DotPlot
immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top <- immune.combined.markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC)
DotPlot(immune.combined, features = markers.to.plot, cols = c("AliceBlue","#DC143C")
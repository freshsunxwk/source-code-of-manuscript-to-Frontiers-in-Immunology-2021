##############xwk###############
# analysis for myeloid cells for example
# analysis for myeloid T cells  B cells is similar
#1
immune.combined <- FindVariableFeatures(object = immune.combined, selection.method = "vst", nfeatures = 2000)
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 20, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca",dims = 1:20,seed.use =3)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
Idents(immune.combined) <- immune.combined$seurat_clusters#

#2
col <-c( '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3',
         '#E95C59', '#E59CC4', '#AB3282','#BD956A', '#9FA3A8' )
DimPlot(immune.combined, reduction = "umap",label = T,cols = col)
DimPlot(immune.combined,reduction = "umap",label=T, cols = col,split.by = "stim")

#3
immune.combined_for_SingleR <- GetAssayData(immune.combined, slot="data")
clusters=immune.combined@meta.data$seurat_clusters
mouseImmu <- ref_IGD
pred.mouseImmu <- SingleR(test = immune.combined_for_SingleR, ref = mouseImmu, labels = mouseImmu$label.fine,method = "cluster", clusters = clusters,assay.type.test = "logcounts", assay.type.ref = "logcounts")
immune.combined@meta.data$singleR=cellType[match(clusters,cellType$ClusterID),'mouseImmu']

#4
immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'roc')
top10 <- immune.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object = immune.combined, features = top10$gene,group.colors=col)+NoLegend()
table(immune.combined@meta.data$seurat_clusters)

#5
# markers.to.plot=c("Cd14","Fcgr3","Cd83","H2-DMb2","Cd83","H2-Ab1",
#                   "H2-Eb1","H2-Aa","Cd74","H2-DMa","H2-Aa","Cd74",
#                   "Cd40","Ccr2","Ifitm2","Ifitm3","Ifitm6","Ifit1",
#                   "Ifit2","Ifit3","Irf7","Jun","Fos","Egr1","Klf2",
#                   "Atf3","Cd14","Mrc1","Hspa1a","Hspa1b","Ccl5",
#                   "Rsad2","Flt3","Cst3","Cd209a","Clec9a","Ly6i","Fcgr4")
DotPlot(object = immune.combined, features =markers.to.plot, dot.scale =6)+ scale_color_viridis()+coord_flip()

#6
VlnPlot(immune.combined, features ="Cst3",log =F, pt.size = 0,cols =col)+ NoLegend()

#7
df1=immune.combined@meta.data
colnames(df1)
ggplot(df1,aes(x=stim,fill=singleR_2)#,col=col
)+geom_bar(stat="count",position="dodge")+scale_fill_manual(values=col)+theme(legend.key = element_blank())
ggplot(df1,aes(x=stim,fill=singleR_2),col=col)+geom_bar(stat="count",position="fill")+scale_fill_manual(values=col)+theme(legend.key = element_blank())
huatu=table(immune.combined@meta.data$stim,immune.combined@active.ident)
tab.1=table(immune.combined@meta.data$stim,immune.combined@active.ident) 
balloonplot(tab.1)
tab.1=table(immune.combined@meta.data$seurat_clusters,immune.combined@active.ident) 
balloonplot(tab.1)

#8
# immune.combined <- RenameIdents(
#   object = immune.combined,
#   "0" = "Resident",
#   "1" = "Inflammtory",
#   "2" = "Resident",
#   "3" = "Resident",
#   "4" = "Resident",
#   "5" = "Resident",
#   "6" = "Resident",
#   "7" = "DC",
#   "8" = "DC",
#   "9" = "Resident",
#   "10" = "Monocyte")

#9
b.interferon.response <- FindMarkers(immune.combined,
                                     ident.1 = "Inflammtory",
                                     ident.2 = "Resident",
                                     logfc.threshold = 0.00,
                                     verbose = FALSE, min.pct = 0.25)
b.interferon.response <- FindMarkers(immune.combined,
                                     ident.1 = "Inflammtory",
                                     ident.2 = "Monocyte",
                                     logfc.threshold = 0.00,
                                     verbose = FALSE, min.pct = 0.25)
nrDEG=b.interferon.response
nrDEG$adj.P.Val <- ifelse(nrDEG$adj.P.Val < 1e-300,1e-300,nrDEG$adj.P.Val)
nrDEG$logP=-log10(nrDEG$adj.P.Val)
nrDEG$gene_name=rownames(nrDEG)
logFC_cutoff=0.25#
PValue_cutoff=0.05#
yintercept_LINE=-log10(PValue_cutoff)
nrDEG$group='not-significant'
nrDEG$group[which((nrDEG$adj.P.Val<PValue_cutoff)&(nrDEG$logFC>logFC_cutoff))]='up-regulated'
nrDEG$group[which((nrDEG$adj.P.Val<PValue_cutoff)&(nrDEG$logFC< -logFC_cutoff))]='down-regulated'
table(nrDEG$group)
nrDEG <- subset(nrDEG, nrDEG$group!='not-significant')
which(nrDEG$group!='not-significant')
ggscatter(nrDEG,x='logFC',y='logP')


#clusterProfiler

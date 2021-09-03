##############xwk###############
#GSEA for DEG between the scRNA cluster
geneList = sort(geneList, decreasing = TRUE)
# Molecular Signatures Database
library(msigdbr)
library(ggplot2)
library(ggnewscale)
H (hallmark gene sets, 50 gene sets)
C1 (positional gene sets, 278 gene sets)
by chromosome: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT
C2 (curated gene sets, 6290 gene sets)
CGP (chemical and genetic perturbations, 3368 gene sets)
CP (canonical pathways, 2922 gene sets)
CP:BIOCARTA (BioCarta gene sets, 292 gene sets)
CP:KEGG (KEGG gene sets, 186 gene sets)
CP:PID (PID gene sets, 196 gene sets)
CP:REACTOME (Reactome gene sets, 1604 gene sets)
CP:WIKIPATHWAYS (WikiPathways gene sets, 615 gene sets)
C3 (regulatory target gene sets, 3731 gene sets)
MIR (microRNA targets, 2598 gene sets)
MIR:MIR_Legacy (legacy microRNA targets, 221 gene sets)
MIR:MIRDB (MIRDB microRNA targets, 2377 gene sets)
TFT (all transcription factor targets, 1133 gene sets)
TFT:GTRD (GTRD transcription factor targets, 523 gene sets)
TFT:TFT_Legacy (legacy transcription factor targets, 610 gene sets)
C4 (computational gene sets, 858 gene sets)
CGN (cancer gene neighborhoods, 427 gene sets)
CM (cancer modules, 431 gene sets)
C5 (ontology gene sets, 14998 gene sets)
GO (Gene Ontology, 10185 gene sets)
GO:BP (GO biological process, 7481 gene sets)
GO:CC (GO cellular component, 996 gene sets)
GO:MF (GO molecular function, 1708 gene sets)
HPO (Human Phenotye Ontology, 4813 gene sets)
C6 (oncogenic signature gene sets, 189 gene sets)
C7 (immunologic signature gene sets, 5219 gene sets)
IMMUNESIGDB (ImmuneSigDB gene sets, 4872 gene sets)
VAX (vaccine response gene sets, 347 gene sets)
C8 (cell type signature gene sets, 671 gene sets)
# msigdbr_show_species()
# #"Homo sapiens"             "Mus musculus"            
# Dm_msigdbr <- msigdbr(species="Mus musculus")
load("F:/XU/20210329LAST/for-GSEA.Rdata")
DmGOlist = split(x = DmGO$gene_symbol, f = DmGO$gs_name)
cc=names(DmGOlist)
write.csv(cc,file = 'DmGOlist--names.csv')
DmGO <- msigdbr(species="Mus musculus") %>% dplyr::select(gs_name, entrez_gene, gene_symbol)
DmGO <- msigdbr(species="Mus musculus",category="C5", subcategory = "GO:BP") %>% dplyr::select(gs_name, entrez_gene, gene_symbol)
# GSEA(
#   geneList,
#   exponent = 1,
#   minGSSize = 10,
#   maxGSSize = 500,
#   eps = 1e-10,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   TERM2GENE,
#   TERM2NAME = NA,
#   verbose = TRUE,
#   seed = FALSE,
#   by = "fgsea",
#   ...
# )
KEGG <- GSEA(geneList,TERM2GENE=DmGO[,c(1,3)])
dotplot(KEGG,color="pvalue",showCategory = 12)
dotplot(KEGG,split=".sign")+facet_grid(~.sign) 
gseaplot2(KEGG,1,color="red",pvalue_table = T)
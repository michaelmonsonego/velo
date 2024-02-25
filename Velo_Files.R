library(Seurat)
library(tidyverse)

S.obj <- readRDS(r"(D:\Michael\velo\ICOS.Tregsubset_merged.RDS)")
S.obj$Sample

############### by https://github.com/basilkhuder/Seurat-to-RNA-Velocity ##################
write.csv(Embeddings(S.obj, reduction = "umap"), file = r"(D:\Michael\velo\cell_embeddings_umap.csv)")
write.csv(Idents(S.obj), file = r"(D:\Michael\velo\clusters.csv)")
write.csv(S.obj$Sample, file = r"(D:\Michael\velo\Treatment.csv)")

#M# just checking consistency
DimPlot(S.obj, reduction = "tsne", label= TRUE, pt.size=2, label.size = 8)
DimPlot(S.obj, reduction = "umap", label= TRUE, pt.size=2, label.size = 8)

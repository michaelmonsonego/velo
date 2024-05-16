library(Seurat)
library(tidyverse)

S.obj <- readRDS(r"(D:\Michael\velo\<whatever>)")
S.obj$Sample

############### by https://github.com/basilkhuder/Seurat-to-RNA-Velocity ##################
write.csv(Embeddings(S.obj, reduction = "umap"), file = r"(D:\Michael\icos_velo\cell_embeddings_umap.csv)")
write.csv(Idents(S.obj), file = r"(D:\Michael\velo\clusters.csv)")
write.csv(S.obj$Sample, file = r"(D:\Michael\velo\Treatment.csv)")

#M# just checking consistency
DimPlot(S.obj, reduction = "tsne", label= TRUE, pt.size=2, label.size = 8)
DimPlot(S.obj, reduction = "umap", label= TRUE, pt.size=2, label.size = 8)

setwd(r"(D:\Michael\velo)")
#M# Control induced dynamic
FeaturePlot(S.obj, features = c('Tgtp2', 'Pycard', 'Mrpl20', 'Kcna4', 'Reep5'), order=TRUE,pt.size=1, reduction="umap", ncol=3)
ggsave(file = "Rplots/Control induced dynamic.png", dpi=300, width=12, height=10)

#M# Treatment induced dynamic
FeaturePlot(S.obj, features = c('Ldha', 'Gimap4', 'Polr2i', 'Ly6e', 'Ccn4'), order=TRUE,pt.size=1, reduction="umap", ncol=3)
ggsave(file = "Rplots/Treatment induced dynamic.png", dpi=300, width=12, height=10)


S.obj_markers = FindAllMarkers(S.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.05)
Top50Markers =
  S.obj_markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) %>%
  as.data.frame %>% 
  arrange(cluster,-avg_log2FC)
write_csv(Top50Markers, "Rplots/50_DE genes.csv")



expression_matrix <- GetAssayData(object = S.obj, slot = "counts")
gene_expression <- expression_matrix['Pdcd1', ]
FeaturePlot(S.obj, features = c("Ptprc"), reduction="tsne")





library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load matrices
ctrl.data <- read.table("data/immune_control_expression_matrix.txt", header=TRUE, row.names=1, sep="\t")
stim.data <- read.table("data/immune_stimulated_expression_matrix.txt", header=TRUE, row.names=1, sep="\t")

# Initialize Seurat objects
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "IMMUNE_CTRL", min.cells = 3, min.features = 200)
stim <- CreateSeuratObject(counts = stim.data, project = "IMMUNE_STIM", min.cells = 3, min.features = 200)

ctrl$stim <- "CTRL"
stim$stim <- "STIM"

# Normalize and find variable features
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
stim <- NormalizeData(stim, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Integration via CCA
immune.list <- list(ctrl = ctrl, stim = stim)
immune.features <- SelectIntegrationFeatures(object.list = immune.list, nfeatures = 2000)
immune.anchors <- FindIntegrationAnchors(object.list = immune.list, anchor.features = immune.features)
immune.combined <- IntegrateData(anchorset = immune.anchors)

# Downstream analysis
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Save split UMAP
p_split <- DimPlot(immune.combined, reduction = "umap", split.by = "stim")
ggsave("output/umap_split.png", plot = p_split, width = 14, height = 6)

# Switch to RNA assay and Join layers for Seurat v5 compatibility
DefaultAssay(immune.combined) <- "RNA"
immune.combined <- JoinLayers(immune.combined)

# Conserved markers
nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
write.csv(head(nk.markers, 15), file = "output/cluster7_nk_markers.csv")

# Feature plot for canonical markers
canonical_genes <- intersect(c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", "CCL2", "PPBP"), rownames(immune.combined))
p_feat <- FeaturePlot(immune.combined, features = canonical_genes, min.cutoff = "q9", combine = TRUE)
ggsave("output/canonical_markers_featureplot.png", plot = p_feat, width = 12, height = 12)

# Cell type annotation
immune.combined <- RenameIdents(immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", 
    `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "T activated", `7` = "NK", `8` = "DC", `9` = "B Activated", 
    `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets")

p_labeled <- DimPlot(immune.combined, label = TRUE)
ggsave("output/umap_celltypes.png", plot = p_labeled, width = 8, height = 7)

# DotPlot
Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("Mono/Mk Doublets", "pDC", 
    "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated", 
    "CD4 Naive T", "CD4 Memory T"))
markers.to.plot <- intersect(c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5", 
    "CD8A", "MS4A1", "CD79A", "MIR1HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", 
    "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ"), rownames(immune.combined))

p_dot <- DotPlot(immune.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
    split.by = "stim") + RotatedAxis()
ggsave("output/celltype_markers_dotplot.png", plot = p_dot, width = 14, height = 8)

# Scatter plots
t.cells <- subset(immune.combined, idents = "CD4 Naive T")
Idents(t.cells) <- "stim"
avg.t.cells <- as.data.frame(log1p(AverageExpression(t.cells, verbose = FALSE)$RNA))
avg.t.cells$gene <- rownames(avg.t.cells)

cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono) <- "stim"
avg.cd14.mono <- as.data.frame(log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA))
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

genes.to.label <- c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
ggsave("output/avg_expression_scatter.png", plot = (p1 + p2), width = 12, height = 6)

# Differential response
immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"
b.interferon.response <- FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
write.csv(head(b.interferon.response, n = 15), file = "output/b_interferon_response.csv")

# Interferon Feature and Violin plots
p_resp <- FeaturePlot(immune.combined, features = intersect(c("CD3D", "GNLY", "IFI6"), rownames(immune.combined)), split.by = "stim", max.cutoff = 3, cols = c("grey", "red"), combine = TRUE)
ggsave("output/interferon_response_featureplot.png", plot = p_resp, width = 14, height = 5)

plots <- VlnPlot(immune.combined, features = intersect(c("LYZ", "ISG15", "CXCL10"), rownames(immune.combined)), split.by = "stim", group.by = "celltype", pt.size = 0, combine = FALSE)
ggsave("output/interferon_response_violin.png", plot = wrap_plots(plots = plots, ncol = 1), width = 10, height = 12)

cat("R Analysis Workflow Completed Successfully.\n")

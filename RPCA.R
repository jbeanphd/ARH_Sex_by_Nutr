## Integration with Campbell et al. 2017 dataset, neuronal cells only
library(Seurat)
library(Matrix)
library(data.table)
library(ggplot2)
library(patchwork)

ref <- readRDS('ARH_NN_neurons_integrated_RPCA.rds')

arh <- readRDS('GSE282955_ARH_Sex_by_Nutr_neurons_integrated_RPCA.rds')
arh <- subset(arh,subset = cell_type == 'Neurons')

DefaultAssay(arh) <- "RNA"
arh <- NormalizeData(arh)
arh  <- FindVariableFeatures(arh , selection.method = "vst", nfeatures = 2000)
arh  <- ScaleData(arh)
arh <- RunPCA(arh, npcs = 30, verbose = FALSE)
arh <- RunUMAP(arh, reduction = "pca", dims = 1:30)
arh <- RunTSNE(arh, reduction = "pca", dims = 1:30)
arh <- FindNeighbors(arh, reduction = "pca", dims = 1:30)
arh <- FindClusters(arh, resolution = 0.5)

##Integration
ref@meta.data$ident <- "Campbell2017"
arh@meta.data$ident <- "Self_generated"
Idents(ref) <- 'ident'
Idents(arh) <- 'ident'
DefaultAssay(ref)
DefaultAssay(arh)

obj_list <- list("ref" = ref, "self" = arh)

features <- rownames(arh)
features <- intersect(features, rownames(ref))

coembed.anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features, reduction = "rpca")
coembed.combined <- IntegrateData(anchorset = coembed.anchors, new.assay.name = "coembed")

DefaultAssay(coembed.combined) <- "coembed"

# Run the standard workflow for visualization and clustering
coembed.combined <- ScaleData(coembed.combined, verbose = FALSE)
coembed.combined <- RunPCA(coembed.combined, npcs = 30, verbose = FALSE)
coembed.combined <- RunUMAP(coembed.combined, reduction = "pca", dims = 1:30)
coembed.combined <- RunTSNE(coembed.combined, reduction = "pca", dims = 1:30)
coembed.combined <- FindNeighbors(coembed.combined, reduction = "pca", dims = 1:30)
coembed.combined <- FindClusters(coembed.combined, resolution = 0.5, cluster.name = "Rp5_clusters")

saveRDS(coembed.combined, "integration_campbell_neurons_v1.rds")

# Plot DimPlot
DimPlot(coembed.combined, group.by = "ident", order = c("Campbell2017", "Self_generated"))+ labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))

DimPlot(coembed.combined, group.by = "ident", order = c("Campbell2017", "Self_generated"))

DimPlot(coembed.combined, reduction = 'tsne',group.by = "ident", order = c("Campbell2017", "Self_generated"))

p1 <- DimPlot(
  coembed.combined,
  cells = WhichCells(coembed.combined, expression = ident == "Campbell2017"),
  group.by = "ident",
  pt.size = 0.5,
  cols = c("#14B4B8")
) + labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))

p2 <- DimPlot(
  coembed.combined,
  cells = WhichCells(coembed.combined, expression = ident == "Self_generated"),
  group.by = "ident",
  pt.size = 0.5,
  cols = c("#EB7369")
) + labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))



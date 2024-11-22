libs <- c( 'Seurat','dplyr','tidyr','ggplot2','scDblFinder','SoupX')

lapply(libs, require, character.only = TRUE)

#the purpose of this script is to filter out ambient RNA using SoupX program

#create basic clustering for sample
F_Fed_B1_flt <- Read10X('../STARsolo/F_Fed_B1Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fed_B1_flt <- CreateSeuratObject(F_Fed_B1_flt)

F_Fed_B1_flt <- SCTransform(F_Fed_B1_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fed_B1_flt <- FindVariableFeatures(F_Fed_B1_flt)

yx_chrom_genes <- c('Eif2s3y','Sry','Zfy','Rps4y1','Amely','Tbl1y','Pcdh11y','Tgif2ly',
                    'Tspy1','Tspy2','Azfa','Usp9y','Ddx3y','Uty','Tb4y','Azfb',
                    'Cyorf15','Rps4y2','Eif1ay','Kdm5d','Xkry','Hsfy1','Hsfy2',
                    'Pry','Pry2','Rbmy1a1','Azfc','Daz1','Daz2','Daz3','Daz4',
                    'Cdy1','Cdy2','Vcy1','Vcy2','Xist','Tsix'
)

mito.genes <- grep(pattern = "^mt-", x = rownames(x = F_Fed_B1_flt@assays$RNA@data), value = TRUE)
hemoglobin_genes <- c('Hbq1a','Hbb-y','Hbb-bt','Hba-a2','Hba-a1')
ribo.genes <- grep(pattern = '^Rp[sl]', x = rownames(x = F_Fed_B1_flt@assays$RNA@data), value = TRUE)


F_Fed_B1_flt@assays$SCT@var.features <- F_Fed_B1_flt@assays$SCT@var.features[(!F_Fed_B1_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

F_Fed_B1_flt <- RunPCA(F_Fed_B1_flt, features = VariableFeatures(object = F_Fed_B1_flt))
F_Fed_B1_flt <- RunUMAP(F_Fed_B1_flt, reduction = "pca", dims = 1:30)

F_Fed_B1_flt <- FindNeighbors(F_Fed_B1_flt, dims = 1:30)
F_Fed_B1_flt <- FindClusters(F_Fed_B1_flt, resolution = 0.8)


DimPlot(F_Fed_B1_flt, label = TRUE, label.size = 2) + NoLegend()

#FeaturePlot(F_Fed_B1_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B1_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

#FeaturePlot(F_Fed_B1_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B1_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fed_B1_markers_b4_soupx <- FindAllMarkers(F_Fed_B1_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

F_Fed_B1_UMAP <- F_Fed_B1_flt@reductions$umap@cell.embeddings |> as.data.frame()
F_Fed_B1_meta <- F_Fed_B1_flt@meta.data[,c(2,3,7)] |> cbind(F_Fed_B1_UMAP)



#start SoupX in earnest 
F_Fed_B1_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fed_B1Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fed_B1_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fed_B1Solo.out/GeneFull_Ex50pAS/raw/')
F_Fed_B1_soup <-SoupChannel(F_Fed_B1_raw, F_Fed_B1_flt)


F_Fed_B1_soup = setDR(F_Fed_B1_soup, F_Fed_B1_meta[colnames(F_Fed_B1_flt), c("UMAP_1", "UMAP_2")])
F_Fed_B1_soup = setClusters(F_Fed_B1_soup, setNames(F_Fed_B1_meta$seurat_clusters, colnames(F_Fed_B1_flt)))


#ggplot(F_Fed_B1_soup$metaData, aes(UMAP_1, UMAP_2)) +
#  geom_point(aes(color = clusters))

F_Fed_B1_soup = autoEstCont(F_Fed_B1_soup)
# rho estimated at 0.04; however, 0.12 producing cleaner looking clusters

F_Fed_B1_soup <- setContaminationFraction(F_Fed_B1_soup, 0.12)
F_Fed_B1out <- adjustCounts(F_Fed_B1_soup, method = 'multinomial')

# SoupX finished, check results
F_Fed_B1out <- CreateSeuratObject(F_Fed_B1out)

F_Fed_B1out <- SCTransform(F_Fed_B1out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fed_B1out <- FindVariableFeatures(F_Fed_B1out)
F_Fed_B1out@assays$SCT@var.features <- F_Fed_B1out@assays$SCT@var.features[(!F_Fed_B1out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fed_B1out <- RunPCA(F_Fed_B1out, features = VariableFeatures(object = F_Fed_B1out))
F_Fed_B1out <- RunUMAP(F_Fed_B1out, reduction = "pca", dims = 1:30)

F_Fed_B1out <- FindNeighbors(F_Fed_B1out, dims = 1:30)
F_Fed_B1out <- FindClusters(F_Fed_B1out, resolution = 0.8)

DimPlot(F_Fed_B1out, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(F_Fed_B1out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B1_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fed_B1out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B1_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fed_B1_markers_post_soupx <- FindAllMarkers(F_Fed_B1out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#mito.genes <- grep(pattern = "^mt-", x = rownames(x = F_Fed_B1out@assays$RNA@data), value = TRUE)
mitoPercent <- Matrix::colSums(F_Fed_B1out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(F_Fed_B1out@assays$RNA@counts)*100
F_Fed_B1out$mitoPercent <- mitoPercent


hemoglobin <- Matrix::colSums(F_Fed_B1out@assays$RNA@counts[hemoglobin_genes, ])
F_Fed_B1out <- AddMetaData(object = F_Fed_B1out, metadata = hemoglobin, col.name = "hemoglobin")


#ribo.genes <- grep(pattern = '^Rp[sl]', x = rownames(x = F_Fed_B1out@assays$RNA@data), value = TRUE)
riboPercent <- Matrix::colSums(F_Fed_B1out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(F_Fed_B1out@assays$RNA@counts)*100
F_Fed_B1out$riboPercent <- riboPercent

VlnPlot(F_Fed_B1out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fed_B1_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#filter out extreme values, that may be low quality or potential doublet
F_Fed_B1out_filtered <- subset(F_Fed_B1out, 
                               subset = nCount_RNA < (mean(F_Fed_B1out@meta.data$nCount_RNA) + 4*sd(F_Fed_B1out@meta.data$nCount_RNA))
                               & nFeature_RNA < (mean(F_Fed_B1out@meta.data$nFeature_RNA) + 4*sd(F_Fed_B1out@meta.data$nFeature_RNA))
                               & mitoPercent < (mean(F_Fed_B1out@meta.data$mitoPercent) + 4*sd(F_Fed_B1out@meta.data$mitoPercent))
                               & riboPercent < (mean(F_Fed_B1out@meta.data$riboPercent) + 4*sd(F_Fed_B1out@meta.data$riboPercent))
                               & hemoglobin < (mean(F_Fed_B1out@meta.data$hemoglobin) + 4*sd(F_Fed_B1out@meta.data$hemoglobin)))


VlnPlot(F_Fed_B1out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fed_B1_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
DefaultAssay(F_Fed_B1out_filtered) <- 'RNA'
F_Fed_B1out_filtered <- DietSeurat(F_Fed_B1out_filtered, assays = 'RNA')
F_Fed_B1out_filtered <- as.SingleCellExperiment(F_Fed_B1out_filtered)

F_Fed_B1_DBR <- F_Fed_B1out_filtered@colData@nrows/50000
F_Fed_B1_scDbl <- scDblFinder(F_Fed_B1out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = F_Fed_B1_DBR, propRandom = 0.2, iter = 5)
F_Fed_B1_scDbl <- as.Seurat(F_Fed_B1_scDbl)

F_Fed_B1_scDbl <- SCTransform(F_Fed_B1_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fed_B1_scDbl <- FindVariableFeatures(F_Fed_B1_scDbl)
F_Fed_B1_scDbl@assays$SCT@var.features <- F_Fed_B1_scDbl@assays$SCT@var.features[(!F_Fed_B1_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fed_B1_scDbl <- RunPCA(F_Fed_B1_scDbl, features = VariableFeatures(object = F_Fed_B1_scDbl))
F_Fed_B1_scDbl <- RunUMAP(F_Fed_B1_scDbl, reduction = "pca", dims = 1:30)



F_Fed_B1_scDbl <- FindNeighbors(F_Fed_B1_scDbl, dims = 1:30)
F_Fed_B1_scDbl <- FindClusters(F_Fed_B1_scDbl, resolution = 0.8)



DimPlot(F_Fed_B1_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/F_Fed_B1_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)
#F_Fed_B1_scDbl <- DietSeurat(F_Fed_B1_scDbl, dimreducs = NULL)

F_Fed_B1_scDbl <- subset(F_Fed_B1_scDbl, scDblFinder.class == 'singlet')


F_Fed_B1_scDbl <- SCTransform(F_Fed_B1_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fed_B1_scDbl <- FindVariableFeatures(F_Fed_B1_scDbl)
F_Fed_B1_scDbl@assays$SCT@var.features <- F_Fed_B1_scDbl@assays$SCT@var.features[(!F_Fed_B1_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fed_B1_scDbl <- RunPCA(F_Fed_B1_scDbl, features = VariableFeatures(object = F_Fed_B1_scDbl))
F_Fed_B1_scDbl <- RunUMAP(F_Fed_B1_scDbl, reduction = "pca", dims = 1:30)



F_Fed_B1_scDbl <- FindNeighbors(F_Fed_B1_scDbl, dims = 1:30)
F_Fed_B1_scDbl <- FindClusters(F_Fed_B1_scDbl, resolution = 0.8)

#check purity of clusters
FeaturePlot(F_Fed_B1_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B1_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fed_B1_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B1_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fed_B1_markers_post_scDblFnd <- FindAllMarkers(F_Fed_B1_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#FindMarkers(F_Fed_B1_scDbl, ident.1 = c('4','8'), features = c('Agrp','Npy'), only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save F_Fed_B1_scDbl, markers
rm(F_Fed_B1_flt, F_Fed_B1_meta, F_Fed_B1_raw, F_Fed_B1_soup, F_Fed_B1_UMAP, 
   F_Fed_B1out, F_Fed_B1out_filtered, F_Fed_B1_singlets)




#F_Fed_B2
#create basic clustering for sample
F_Fed_B2_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fed_B2Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fed_B2_flt <- CreateSeuratObject(F_Fed_B2_flt)

F_Fed_B2_flt <- SCTransform(F_Fed_B2_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fed_B2_flt <- FindVariableFeatures(F_Fed_B2_flt)
F_Fed_B2_flt@assays$SCT@var.features <- F_Fed_B2_flt@assays$SCT@var.features[(!F_Fed_B2_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

F_Fed_B2_flt <- RunPCA(F_Fed_B2_flt, features = VariableFeatures(object = F_Fed_B2_flt))
F_Fed_B2_flt <- RunUMAP(F_Fed_B2_flt, reduction = "pca", dims = 1:30)

F_Fed_B2_flt <- FindNeighbors(F_Fed_B2_flt, dims = 1:30)
F_Fed_B2_flt <- FindClusters(F_Fed_B2_flt, resolution = 0.8)


DimPlot(F_Fed_B2_flt, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(F_Fed_B2_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B2_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fed_B2_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B2_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fed_B2_markers_b4_soupx <- FindAllMarkers(F_Fed_B2_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

F_Fed_B2_UMAP <- F_Fed_B2_flt@reductions$umap@cell.embeddings |> as.data.frame()
F_Fed_B2_meta <- F_Fed_B2_flt@meta.data[,c(2,3,7)] |> cbind(F_Fed_B2_UMAP)



#start SoupX in earnest 
F_Fed_B2_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fed_B2Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fed_B2_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fed_B2Solo.out/GeneFull_Ex50pAS/raw/')
F_Fed_B2_soup <-SoupChannel(F_Fed_B2_raw, F_Fed_B2_flt)


F_Fed_B2_soup = setDR(F_Fed_B2_soup, F_Fed_B2_meta[colnames(F_Fed_B2_flt), c("UMAP_1", "UMAP_2")])
F_Fed_B2_soup = setClusters(F_Fed_B2_soup, setNames(F_Fed_B2_meta$seurat_clusters, colnames(F_Fed_B2_flt)))


ggplot(F_Fed_B2_soup$metaData, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = clusters))

F_Fed_B2_soup = autoEstCont(F_Fed_B2_soup)
# rho estimated at 0.08; however, 0.20 producing cleaner looking clusters

F_Fed_B2_soup <- setContaminationFraction(F_Fed_B2_soup, 0.2)
F_Fed_B2out <- adjustCounts(F_Fed_B2_soup, method = 'multinomial')

# SoupX finished, check results
F_Fed_B2out <- CreateSeuratObject(F_Fed_B2out)

F_Fed_B2out <- SCTransform(F_Fed_B2out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fed_B2out <- FindVariableFeatures(F_Fed_B2out)
F_Fed_B2out@assays$SCT@var.features <- F_Fed_B2out@assays$SCT@var.features[(!F_Fed_B2out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fed_B2out <- RunPCA(F_Fed_B2out, features = VariableFeatures(object = F_Fed_B2out))
F_Fed_B2out <- RunUMAP(F_Fed_B2out, reduction = "pca", dims = 1:30)

F_Fed_B2out <- FindNeighbors(F_Fed_B2out, dims = 1:30)
F_Fed_B2out <- FindClusters(F_Fed_B2out, resolution = 0.8)

DimPlot(F_Fed_B2out, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(F_Fed_B2out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B2_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fed_B2out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B2_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fed_B2_markers_post_soupx <- FindAllMarkers(F_Fed_B2out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#mito.genes <- grep(pattern = "^mt-", x = rownames(x = F_Fed_B2out@assays$RNA@data), value = TRUE)
mitoPercent <- Matrix::colSums(F_Fed_B2out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(F_Fed_B2out@assays$RNA@counts)*100
F_Fed_B2out$mitoPercent <- mitoPercent


hemoglobin <- Matrix::colSums(F_Fed_B2out@assays$RNA@counts[hemoglobin_genes, ])
F_Fed_B2out <- AddMetaData(object = F_Fed_B2out, metadata = hemoglobin, col.name = "hemoglobin")


#ribo.genes <- grep(pattern = '^Rp[sl]', x = rownames(x = F_Fed_B2out@assays$RNA@data), value = TRUE)
riboPercent <- Matrix::colSums(F_Fed_B2out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(F_Fed_B2out@assays$RNA@counts)*100
F_Fed_B2out$riboPercent <- riboPercent

VlnPlot(F_Fed_B2out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fed_B2_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fed_B2out_filtered <- subset(F_Fed_B2out, 
                               subset = nCount_RNA < (mean(F_Fed_B2out@meta.data$nCount_RNA) + 4*sd(F_Fed_B2out@meta.data$nCount_RNA))
                               & nFeature_RNA < (mean(F_Fed_B2out@meta.data$nFeature_RNA) + 4*sd(F_Fed_B2out@meta.data$nFeature_RNA))
                               & mitoPercent < (mean(F_Fed_B2out@meta.data$mitoPercent) + 4*sd(F_Fed_B2out@meta.data$mitoPercent))
                               & riboPercent < (mean(F_Fed_B2out@meta.data$riboPercent) + 4*sd(F_Fed_B2out@meta.data$riboPercent))
                               & hemoglobin < (mean(F_Fed_B2out@meta.data$hemoglobin) + 4*sd(F_Fed_B2out@meta.data$hemoglobin)))


VlnPlot(F_Fed_B2out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fed_B2_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
DefaultAssay(F_Fed_B2out_filtered) <- 'RNA'
F_Fed_B2out_filtered <- DietSeurat(F_Fed_B2out_filtered, assays = 'RNA')
F_Fed_B2out_filtered <- as.SingleCellExperiment(F_Fed_B2out_filtered)

F_Fed_B2_DBR <- F_Fed_B2out_filtered@colData@nrows/50000
F_Fed_B2_scDbl <- scDblFinder(F_Fed_B2out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = F_Fed_B2_DBR, propRandom = 0.2, iter = 5)
F_Fed_B2_scDbl <- as.Seurat(F_Fed_B2_scDbl)

F_Fed_B2_scDbl <- SCTransform(F_Fed_B2_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fed_B2_scDbl <- FindVariableFeatures(F_Fed_B2_scDbl)
F_Fed_B2_scDbl@assays$SCT@var.features <- F_Fed_B2_scDbl@assays$SCT@var.features[(!F_Fed_B2_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fed_B2_scDbl <- RunPCA(F_Fed_B2_scDbl, features = VariableFeatures(object = F_Fed_B2_scDbl))
F_Fed_B2_scDbl <- RunUMAP(F_Fed_B2_scDbl, reduction = "pca", dims = 1:30)



F_Fed_B2_scDbl <- FindNeighbors(F_Fed_B2_scDbl, dims = 1:30)
F_Fed_B2_scDbl <- FindClusters(F_Fed_B2_scDbl, resolution = 0.8)



DimPlot(F_Fed_B2_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/F_Fed_B2_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)
#F_Fed_B2_scDbl <- DietSeurat(F_Fed_B2_scDbl, dimreducs = NULL)

F_Fed_B2_scDbl <- subset(F_Fed_B2_scDbl, scDblFinder.class == 'singlet')


F_Fed_B2_scDbl <- SCTransform(F_Fed_B2_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fed_B2_scDbl <- FindVariableFeatures(F_Fed_B2_scDbl)
F_Fed_B2_scDbl@assays$SCT@var.features <- F_Fed_B2_scDbl@assays$SCT@var.features[(!F_Fed_B2_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fed_B2_scDbl <- RunPCA(F_Fed_B2_scDbl, features = VariableFeatures(object = F_Fed_B2_scDbl))
F_Fed_B2_scDbl <- RunUMAP(F_Fed_B2_scDbl, reduction = "pca", dims = 1:30)



F_Fed_B2_scDbl <- FindNeighbors(F_Fed_B2_scDbl, dims = 1:30)
F_Fed_B2_scDbl <- FindClusters(F_Fed_B2_scDbl, resolution = 0.8)

#check purity of clusters
FeaturePlot(F_Fed_B2_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B2_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fed_B2_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B2_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

DimPlot(F_Fed_B2_scDbl, label = TRUE) + NoLegend()
F_Fed_B2_markers_post_scDblFnd <- FindAllMarkers(F_Fed_B2_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#FindMarkers(F_Fed_B2_scDbl, ident.1 = c('4','8'), features = c('Agrp','Npy'), only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save F_Fed_B2_scDbl, markers
rm(F_Fed_B2_flt, F_Fed_B2_meta, F_Fed_B2_raw, F_Fed_B2_soup, F_Fed_B2_UMAP, 
   F_Fed_B2out, F_Fed_B2out_filtered, F_Fed_B2_singlets)


#F_Fed_B3
#create basic clustering for sample
F_Fed_B3_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fed_B3Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fed_B3_flt <- CreateSeuratObject(F_Fed_B3_flt)

F_Fed_B3_flt <- SCTransform(F_Fed_B3_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fed_B3_flt <- FindVariableFeatures(F_Fed_B3_flt)
F_Fed_B3_flt@assays$SCT@var.features <- F_Fed_B3_flt@assays$SCT@var.features[(!F_Fed_B3_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

F_Fed_B3_flt <- RunPCA(F_Fed_B3_flt, features = VariableFeatures(object = F_Fed_B3_flt))
F_Fed_B3_flt <- RunUMAP(F_Fed_B3_flt, reduction = "pca", dims = 1:30)

F_Fed_B3_flt <- FindNeighbors(F_Fed_B3_flt, dims = 1:30)
F_Fed_B3_flt <- FindClusters(F_Fed_B3_flt, resolution = 0.8)


DimPlot(F_Fed_B3_flt, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(F_Fed_B3_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B3_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fed_B3_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B3_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fed_B3_markers_b4_soupx <- FindAllMarkers(F_Fed_B3_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

F_Fed_B3_UMAP <- F_Fed_B3_flt@reductions$umap@cell.embeddings |> as.data.frame()
F_Fed_B3_meta <- F_Fed_B3_flt@meta.data[,c(2,3,7)] |> cbind(F_Fed_B3_UMAP)



#start SoupX in earnest 
F_Fed_B3_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fed_B3Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fed_B3_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fed_B3Solo.out/GeneFull_Ex50pAS/raw/')
F_Fed_B3_soup <-SoupChannel(F_Fed_B3_raw, F_Fed_B3_flt)


F_Fed_B3_soup = setDR(F_Fed_B3_soup, F_Fed_B3_meta[colnames(F_Fed_B3_flt), c("UMAP_1", "UMAP_2")])
F_Fed_B3_soup = setClusters(F_Fed_B3_soup, setNames(F_Fed_B3_meta$seurat_clusters, colnames(F_Fed_B3_flt)))


F_Fed_B3_soup = autoEstCont(F_Fed_B3_soup)
# rho estimated at 0.06; however, 0.15 producing cleaner looking clusters

F_Fed_B3_soup <- setContaminationFraction(F_Fed_B3_soup, 0.15)
F_Fed_B3out <- adjustCounts(F_Fed_B3_soup, method = 'multinomial')

# SoupX finished, check results
F_Fed_B3out <- CreateSeuratObject(F_Fed_B3out)

F_Fed_B3out <- SCTransform(F_Fed_B3out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fed_B3out <- FindVariableFeatures(F_Fed_B3out)
F_Fed_B3out@assays$SCT@var.features <- F_Fed_B3out@assays$SCT@var.features[(!F_Fed_B3out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fed_B3out <- RunPCA(F_Fed_B3out, features = VariableFeatures(object = F_Fed_B3out))
F_Fed_B3out <- RunUMAP(F_Fed_B3out, reduction = "pca", dims = 1:30)

F_Fed_B3out <- FindNeighbors(F_Fed_B3out, dims = 1:30)
F_Fed_B3out <- FindClusters(F_Fed_B3out, resolution = 0.8)

DimPlot(F_Fed_B3out, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(F_Fed_B3out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B3_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fed_B3out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B3_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fed_B3_markers_post_soupx <- FindAllMarkers(F_Fed_B3out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#mito.genes <- grep(pattern = "^mt-", x = rownames(x = F_Fed_B3out@assays$RNA@data), value = TRUE)
mitoPercent <- Matrix::colSums(F_Fed_B3out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(F_Fed_B3out@assays$RNA@counts)*100
F_Fed_B3out$mitoPercent <- mitoPercent


hemoglobin <- Matrix::colSums(F_Fed_B3out@assays$RNA@counts[hemoglobin_genes, ])
F_Fed_B3out <- AddMetaData(object = F_Fed_B3out, metadata = hemoglobin, col.name = "hemoglobin")


#ribo.genes <- grep(pattern = '^Rp[sl]', x = rownames(x = F_Fed_B3out@assays$RNA@data), value = TRUE)
riboPercent <- Matrix::colSums(F_Fed_B3out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(F_Fed_B3out@assays$RNA@counts)*100
F_Fed_B3out$riboPercent <- riboPercent

VlnPlot(F_Fed_B3out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fed_B3_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fed_B3out_filtered <- subset(F_Fed_B3out, 
                               subset = nCount_RNA < (mean(F_Fed_B3out@meta.data$nCount_RNA) + 4*sd(F_Fed_B3out@meta.data$nCount_RNA))
                               & nFeature_RNA < (mean(F_Fed_B3out@meta.data$nFeature_RNA) + 4*sd(F_Fed_B3out@meta.data$nFeature_RNA))
                               & mitoPercent < (mean(F_Fed_B3out@meta.data$mitoPercent) + 4*sd(F_Fed_B3out@meta.data$mitoPercent))
                               & riboPercent < (mean(F_Fed_B3out@meta.data$riboPercent) + 4*sd(F_Fed_B3out@meta.data$riboPercent))
                               & hemoglobin < (mean(F_Fed_B3out@meta.data$hemoglobin) + 4*sd(F_Fed_B3out@meta.data$hemoglobin)))


VlnPlot(F_Fed_B3out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fed_B3_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
DefaultAssay(F_Fed_B3out_filtered) <- 'RNA'
F_Fed_B3out_filtered <- DietSeurat(F_Fed_B3out_filtered, assays = 'RNA')
F_Fed_B3out_filtered <- as.SingleCellExperiment(F_Fed_B3out_filtered)

F_Fed_B3_DBR <- F_Fed_B3out_filtered@colData@nrows/50000
F_Fed_B3_scDbl <- scDblFinder(F_Fed_B3out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = F_Fed_B3_DBR, propRandom = 0.2, iter = 5)
F_Fed_B3_scDbl <- as.Seurat(F_Fed_B3_scDbl)

F_Fed_B3_scDbl <- SCTransform(F_Fed_B3_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fed_B3_scDbl <- FindVariableFeatures(F_Fed_B3_scDbl)
F_Fed_B3_scDbl@assays$SCT@var.features <- F_Fed_B3_scDbl@assays$SCT@var.features[(!F_Fed_B3_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fed_B3_scDbl <- RunPCA(F_Fed_B3_scDbl, features = VariableFeatures(object = F_Fed_B3_scDbl))
F_Fed_B3_scDbl <- RunUMAP(F_Fed_B3_scDbl, reduction = "pca", dims = 1:30)



F_Fed_B3_scDbl <- FindNeighbors(F_Fed_B3_scDbl, dims = 1:30)
F_Fed_B3_scDbl <- FindClusters(F_Fed_B3_scDbl, resolution = 0.8)



DimPlot(F_Fed_B3_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/F_Fed_B3_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)


F_Fed_B3_scDbl <- subset(F_Fed_B3_scDbl, scDblFinder.class == 'singlet')


F_Fed_B3_scDbl <- SCTransform(F_Fed_B3_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fed_B3_scDbl <- FindVariableFeatures(F_Fed_B3_scDbl)
F_Fed_B3_scDbl@assays$SCT@var.features <- F_Fed_B3_scDbl@assays$SCT@var.features[(!F_Fed_B3_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fed_B3_scDbl <- RunPCA(F_Fed_B3_scDbl, features = VariableFeatures(object = F_Fed_B3_scDbl))
F_Fed_B3_scDbl <- RunUMAP(F_Fed_B3_scDbl, reduction = "pca", dims = 1:30)



F_Fed_B3_scDbl <- FindNeighbors(F_Fed_B3_scDbl, dims = 1:30)
F_Fed_B3_scDbl <- FindClusters(F_Fed_B3_scDbl, resolution = 0.8)

#check purity of clusters
FeaturePlot(F_Fed_B3_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B3_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fed_B3_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B3_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fed_B3_markers_post_scDblFnd <- FindAllMarkers(F_Fed_B3_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#FindMarkers(F_Fed_B3_scDbl, ident.1 = c('4','8'), features = c('Agrp','Npy'), only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save F_Fed_B3_scDbl, markers
rm(F_Fed_B3_flt, F_Fed_B3_meta, F_Fed_B3_raw, F_Fed_B3_soup, F_Fed_B3_UMAP, 
   F_Fed_B3out, F_Fed_B3out_filtered, F_Fed_B3_singlets)


#F_Fast_B1
#create basic clustering for sample
F_Fast_B1_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B1Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fast_B1_flt <- CreateSeuratObject(F_Fast_B1_flt)

F_Fast_B1_flt <- SCTransform(F_Fast_B1_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fast_B1_flt <- FindVariableFeatures(F_Fast_B1_flt)
F_Fast_B1_flt@assays$SCT@var.features <- F_Fast_B1_flt@assays$SCT@var.features[(!F_Fast_B1_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

F_Fast_B1_flt <- RunPCA(F_Fast_B1_flt, features = VariableFeatures(object = F_Fast_B1_flt))
F_Fast_B1_flt <- RunUMAP(F_Fast_B1_flt, reduction = "pca", dims = 1:30)

F_Fast_B1_flt <- FindNeighbors(F_Fast_B1_flt, dims = 1:30)
F_Fast_B1_flt <- FindClusters(F_Fast_B1_flt, resolution = 0.8)


DimPlot(F_Fast_B1_flt, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(F_Fast_B1_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B1_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fast_B1_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B1_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fast_B1_markers_b4_soupx <- FindAllMarkers(F_Fast_B1_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

F_Fast_B1_UMAP <- F_Fast_B1_flt@reductions$umap@cell.embeddings |> as.data.frame()
F_Fast_B1_meta <- F_Fast_B1_flt@meta.data[,c(2,3,7)] |> cbind(F_Fast_B1_UMAP)



#start SoupX in earnest 
F_Fast_B1_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B1Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fast_B1_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B1Solo.out/GeneFull_Ex50pAS/raw/')
F_Fast_B1_soup <-SoupChannel(F_Fast_B1_raw, F_Fast_B1_flt)


F_Fast_B1_soup = setDR(F_Fast_B1_soup, F_Fast_B1_meta[colnames(F_Fast_B1_flt), c("UMAP_1", "UMAP_2")])
F_Fast_B1_soup = setClusters(F_Fast_B1_soup, setNames(F_Fast_B1_meta$seurat_clusters, colnames(F_Fast_B1_flt)))


F_Fast_B1_soup = autoEstCont(F_Fast_B1_soup)
# rho estimated at 0.03; however, 0.06 producing cleaner looking clusters

F_Fast_B1_soup <- setContaminationFraction(F_Fast_B1_soup, 0.06)
F_Fast_B1out <- adjustCounts(F_Fast_B1_soup, method = 'multinomial')

# SoupX finished, check results
F_Fast_B1out <- CreateSeuratObject(F_Fast_B1out)

F_Fast_B1out <- SCTransform(F_Fast_B1out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fast_B1out <- FindVariableFeatures(F_Fast_B1out)
F_Fast_B1out@assays$SCT@var.features <- F_Fast_B1out@assays$SCT@var.features[(!F_Fast_B1out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fast_B1out <- RunPCA(F_Fast_B1out, features = VariableFeatures(object = F_Fast_B1out))
F_Fast_B1out <- RunUMAP(F_Fast_B1out, reduction = "pca", dims = 1:30)

F_Fast_B1out <- FindNeighbors(F_Fast_B1out, dims = 1:30)
F_Fast_B1out <- FindClusters(F_Fast_B1out, resolution = 0.8)

DimPlot(F_Fast_B1out, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(F_Fast_B1out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B1_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fast_B1out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B1_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fast_B1_markers_post_soupx <- FindAllMarkers(F_Fast_B1out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#mito.genes <- grep(pattern = "^mt-", x = rownames(x = F_Fast_B1out@assays$RNA@data), value = TRUE)
mitoPercent <- Matrix::colSums(F_Fast_B1out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(F_Fast_B1out@assays$RNA@counts)*100
F_Fast_B1out$mitoPercent <- mitoPercent


hemoglobin <- Matrix::colSums(F_Fast_B1out@assays$RNA@counts[hemoglobin_genes, ])
F_Fast_B1out <- AddMetaData(object = F_Fast_B1out, metadata = hemoglobin, col.name = "hemoglobin")


ribo.genes <- grep(pattern = '^Rp[sl]', x = rownames(x = F_Fast_B1out@assays$RNA@data), value = TRUE)
riboPercent <- Matrix::colSums(F_Fast_B1out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(F_Fast_B1out@assays$RNA@counts)*100
F_Fast_B1out$riboPercent <- riboPercent

VlnPlot(F_Fast_B1out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fast_B1_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fast_B1out_filtered <- subset(F_Fast_B1out, 
                                subset = nCount_RNA < (mean(F_Fast_B1out@meta.data$nCount_RNA) + 4*sd(F_Fast_B1out@meta.data$nCount_RNA))
                                & nFeature_RNA < (mean(F_Fast_B1out@meta.data$nFeature_RNA) + 4*sd(F_Fast_B1out@meta.data$nFeature_RNA))
                                & mitoPercent < (mean(F_Fast_B1out@meta.data$mitoPercent) + 4*sd(F_Fast_B1out@meta.data$mitoPercent))
                                & riboPercent < (mean(F_Fast_B1out@meta.data$riboPercent) + 4*sd(F_Fast_B1out@meta.data$riboPercent))
                                & hemoglobin < (mean(F_Fast_B1out@meta.data$hemoglobin) + 4*sd(F_Fast_B1out@meta.data$hemoglobin)))


VlnPlot(F_Fast_B1out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fast_B1_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
DefaultAssay(F_Fast_B1out_filtered) <- 'RNA'
F_Fast_B1out_filtered <- DietSeurat(F_Fast_B1out_filtered, assays = 'RNA')
F_Fast_B1out_filtered <- as.SingleCellExperiment(F_Fast_B1out_filtered)

F_Fast_B1_DBR <- F_Fast_B1out_filtered@colData@nrows/50000
F_Fast_B1_scDbl <- scDblFinder(F_Fast_B1out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = F_Fast_B1_DBR, propRandom = 0.2, iter = 5)
F_Fast_B1_scDbl <- as.Seurat(F_Fast_B1_scDbl)

F_Fast_B1_scDbl <- SCTransform(F_Fast_B1_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fast_B1_scDbl <- FindVariableFeatures(F_Fast_B1_scDbl)
F_Fast_B1_scDbl@assays$SCT@var.features <- F_Fast_B1_scDbl@assays$SCT@var.features[(!F_Fast_B1_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fast_B1_scDbl <- RunPCA(F_Fast_B1_scDbl, features = VariableFeatures(object = F_Fast_B1_scDbl))
F_Fast_B1_scDbl <- RunUMAP(F_Fast_B1_scDbl, reduction = "pca", dims = 1:30)



F_Fast_B1_scDbl <- FindNeighbors(F_Fast_B1_scDbl, dims = 1:30)
F_Fast_B1_scDbl <- FindClusters(F_Fast_B1_scDbl, resolution = 0.8)



DimPlot(F_Fast_B1_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/F_Fast_B1_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)


F_Fast_B1_scDbl <- subset(F_Fast_B1_scDbl, scDblFinder.class == 'singlet')


F_Fast_B1_scDbl <- SCTransform(F_Fast_B1_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fast_B1_scDbl <- FindVariableFeatures(F_Fast_B1_scDbl)
F_Fast_B1_scDbl@assays$SCT@var.features <- F_Fast_B1_scDbl@assays$SCT@var.features[(!F_Fast_B1_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fast_B1_scDbl <- RunPCA(F_Fast_B1_scDbl, features = VariableFeatures(object = F_Fast_B1_scDbl))
F_Fast_B1_scDbl <- RunUMAP(F_Fast_B1_scDbl, reduction = "pca", dims = 1:30)



F_Fast_B1_scDbl <- FindNeighbors(F_Fast_B1_scDbl, dims = 1:30)
F_Fast_B1_scDbl <- FindClusters(F_Fast_B1_scDbl, resolution = 0.8)

#check purity of clusters
FeaturePlot(F_Fast_B1_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B1_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fast_B1_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B1_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fast_B1_markers_post_scDblFnd <- FindAllMarkers(F_Fast_B1_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#FindMarkers(F_Fast_B1_scDbl, ident.1 = c('4','8'), features = c('Agrp','Npy'), only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save F_Fast_B1_scDbl, markers
rm(F_Fast_B1_flt, F_Fast_B1_meta, F_Fast_B1_raw, F_Fast_B1_soup, F_Fast_B1_UMAP, 
   F_Fast_B1out, F_Fast_B1out_filtered, F_Fast_B1_singlets)


#F_Fast_B2
#create basic clustering for sample
F_Fast_B2_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B2Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fast_B2_flt <- CreateSeuratObject(F_Fast_B2_flt)

F_Fast_B2_flt <- SCTransform(F_Fast_B2_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fast_B2_flt <- FindVariableFeatures(F_Fast_B2_flt)
F_Fast_B2_flt@assays$SCT@var.features <- F_Fast_B2_flt@assays$SCT@var.features[(!F_Fast_B2_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

F_Fast_B2_flt <- RunPCA(F_Fast_B2_flt, features = VariableFeatures(object = F_Fast_B2_flt))
F_Fast_B2_flt <- RunUMAP(F_Fast_B2_flt, reduction = "pca", dims = 1:30)

F_Fast_B2_flt <- FindNeighbors(F_Fast_B2_flt, dims = 1:30)
F_Fast_B2_flt <- FindClusters(F_Fast_B2_flt, resolution = 0.8)


DimPlot(F_Fast_B2_flt, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(F_Fast_B2_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B2_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fast_B2_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B2_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fast_B2_markers_b4_soupx <- FindAllMarkers(F_Fast_B2_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

F_Fast_B2_UMAP <- F_Fast_B2_flt@reductions$umap@cell.embeddings |> as.data.frame()
F_Fast_B2_meta <- F_Fast_B2_flt@meta.data[,c(2,3,7)] |> cbind(F_Fast_B2_UMAP)



#start SoupX in earnest 
F_Fast_B2_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B2Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fast_B2_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B2Solo.out/GeneFull_Ex50pAS/raw/')
F_Fast_B2_soup <-SoupChannel(F_Fast_B2_raw, F_Fast_B2_flt)


F_Fast_B2_soup = setDR(F_Fast_B2_soup, F_Fast_B2_meta[colnames(F_Fast_B2_flt), c("UMAP_1", "UMAP_2")])
F_Fast_B2_soup = setClusters(F_Fast_B2_soup, setNames(F_Fast_B2_meta$seurat_clusters, colnames(F_Fast_B2_flt)))


F_Fast_B2_soup = autoEstCont(F_Fast_B2_soup)
# rho estimated at 0.06; however, 0.12 producing cleaner looking clusters

F_Fast_B2_soup <- setContaminationFraction(F_Fast_B2_soup, 0.12)
F_Fast_B2out <- adjustCounts(F_Fast_B2_soup, method = 'multinomial')

# SoupX finished, check results
F_Fast_B2out <- CreateSeuratObject(F_Fast_B2out)

F_Fast_B2out <- SCTransform(F_Fast_B2out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fast_B2out <- FindVariableFeatures(F_Fast_B2out)
F_Fast_B2out@assays$SCT@var.features <- F_Fast_B2out@assays$SCT@var.features[(!F_Fast_B2out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fast_B2out <- RunPCA(F_Fast_B2out, features = VariableFeatures(object = F_Fast_B2out))
F_Fast_B2out <- RunUMAP(F_Fast_B2out, reduction = "pca", dims = 1:30)

F_Fast_B2out <- FindNeighbors(F_Fast_B2out, dims = 1:30)
F_Fast_B2out <- FindClusters(F_Fast_B2out, resolution = 0.8)

DimPlot(F_Fast_B2out, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(F_Fast_B2out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B2_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fast_B2out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B2_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fast_B2_markers_post_soupx <- FindAllMarkers(F_Fast_B2out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#mito.genes <- grep(pattern = "^mt-", x = rownames(x = F_Fast_B2out@assays$RNA@data), value = TRUE)
mitoPercent <- Matrix::colSums(F_Fast_B2out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(F_Fast_B2out@assays$RNA@counts)*100
F_Fast_B2out$mitoPercent <- mitoPercent


hemoglobin <- Matrix::colSums(F_Fast_B2out@assays$RNA@counts[hemoglobin_genes, ])
F_Fast_B2out <- AddMetaData(object = F_Fast_B2out, metadata = hemoglobin, col.name = "hemoglobin")


ribo.genes <- grep(pattern = '^Rp[sl]', x = rownames(x = F_Fast_B2out@assays$RNA@data), value = TRUE)
riboPercent <- Matrix::colSums(F_Fast_B2out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(F_Fast_B2out@assays$RNA@counts)*100
F_Fast_B2out$riboPercent <- riboPercent

VlnPlot(F_Fast_B2out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fast_B2_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fast_B2out_filtered <- subset(F_Fast_B2out, 
                                subset = nCount_RNA < (mean(F_Fast_B2out@meta.data$nCount_RNA) + 4*sd(F_Fast_B2out@meta.data$nCount_RNA))
                                & nFeature_RNA < (mean(F_Fast_B2out@meta.data$nFeature_RNA) + 4*sd(F_Fast_B2out@meta.data$nFeature_RNA))
                                & mitoPercent < (mean(F_Fast_B2out@meta.data$mitoPercent) + 4*sd(F_Fast_B2out@meta.data$mitoPercent))
                                & riboPercent < (mean(F_Fast_B2out@meta.data$riboPercent) + 4*sd(F_Fast_B2out@meta.data$riboPercent))
                                & hemoglobin < (mean(F_Fast_B2out@meta.data$hemoglobin) + 4*sd(F_Fast_B2out@meta.data$hemoglobin)))


VlnPlot(F_Fast_B2out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fast_B2_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
DefaultAssay(F_Fast_B2out_filtered) <- 'RNA'
F_Fast_B2out_filtered <- DietSeurat(F_Fast_B2out_filtered, assays = 'RNA')
F_Fast_B2out_filtered <- as.SingleCellExperiment(F_Fast_B2out_filtered)

F_Fast_B2_DBR <- F_Fast_B2out_filtered@colData@nrows/50000
F_Fast_B2_scDbl <- scDblFinder(F_Fast_B2out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = F_Fast_B2_DBR, propRandom = 0.2, iter = 5)
F_Fast_B2_scDbl <- as.Seurat(F_Fast_B2_scDbl)

F_Fast_B2_scDbl <- SCTransform(F_Fast_B2_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fast_B2_scDbl <- FindVariableFeatures(F_Fast_B2_scDbl)
F_Fast_B2_scDbl@assays$SCT@var.features <- F_Fast_B2_scDbl@assays$SCT@var.features[(!F_Fast_B2_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fast_B2_scDbl <- RunPCA(F_Fast_B2_scDbl, features = VariableFeatures(object = F_Fast_B2_scDbl))
F_Fast_B2_scDbl <- RunUMAP(F_Fast_B2_scDbl, reduction = "pca", dims = 1:30)



F_Fast_B2_scDbl <- FindNeighbors(F_Fast_B2_scDbl, dims = 1:30)
F_Fast_B2_scDbl <- FindClusters(F_Fast_B2_scDbl, resolution = 0.8)



DimPlot(F_Fast_B2_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/F_Fast_B2_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)


F_Fast_B2_scDbl <- subset(F_Fast_B2_scDbl, scDblFinder.class == 'singlet')


F_Fast_B2_scDbl <- SCTransform(F_Fast_B2_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fast_B2_scDbl <- FindVariableFeatures(F_Fast_B2_scDbl)
F_Fast_B2_scDbl@assays$SCT@var.features <- F_Fast_B2_scDbl@assays$SCT@var.features[(!F_Fast_B2_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fast_B2_scDbl <- RunPCA(F_Fast_B2_scDbl, features = VariableFeatures(object = F_Fast_B2_scDbl))
F_Fast_B2_scDbl <- RunUMAP(F_Fast_B2_scDbl, reduction = "pca", dims = 1:30)



F_Fast_B2_scDbl <- FindNeighbors(F_Fast_B2_scDbl, dims = 1:30)
F_Fast_B2_scDbl <- FindClusters(F_Fast_B2_scDbl, resolution = 0.8)

#check purity of clusters
FeaturePlot(F_Fast_B2_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B2_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fast_B2_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B2_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fast_B2_markers_post_scDblFnd <- FindAllMarkers(F_Fast_B2_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#FindMarkers(F_Fast_B2_scDbl, ident.1 = c('4','8'), features = c('Agrp','Npy'), only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save F_Fast_B2_scDbl, markers
rm(F_Fast_B2_flt, F_Fast_B2_meta, F_Fast_B2_raw, F_Fast_B2_soup, F_Fast_B2_UMAP, 
   F_Fast_B2out, F_Fast_B2out_filtered, F_Fast_B2_singlets)



#F_Fast_B3
#create basic clustering for sample
F_Fast_B3_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B3Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fast_B3_flt <- CreateSeuratObject(F_Fast_B3_flt)

F_Fast_B3_flt <- SCTransform(F_Fast_B3_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fast_B3_flt <- FindVariableFeatures(F_Fast_B3_flt)
F_Fast_B3_flt@assays$SCT@var.features <- F_Fast_B3_flt@assays$SCT@var.features[(!F_Fast_B3_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

F_Fast_B3_flt <- RunPCA(F_Fast_B3_flt, features = VariableFeatures(object = F_Fast_B3_flt))
F_Fast_B3_flt <- RunUMAP(F_Fast_B3_flt, reduction = "pca", dims = 1:30)

F_Fast_B3_flt <- FindNeighbors(F_Fast_B3_flt, dims = 1:30)
F_Fast_B3_flt <- FindClusters(F_Fast_B3_flt, resolution = 0.8)


DimPlot(F_Fast_B3_flt, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(F_Fast_B3_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B3_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fast_B3_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B3_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fast_B3_markers_b4_soupx <- FindAllMarkers(F_Fast_B3_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

F_Fast_B3_UMAP <- F_Fast_B3_flt@reductions$umap@cell.embeddings |> as.data.frame()
F_Fast_B3_meta <- F_Fast_B3_flt@meta.data[,c(2,3,7)] |> cbind(F_Fast_B3_UMAP)



#start SoupX in earnest 
F_Fast_B3_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B3Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fast_B3_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B3Solo.out/GeneFull_Ex50pAS/raw/')
F_Fast_B3_soup <-SoupChannel(F_Fast_B3_raw, F_Fast_B3_flt)


F_Fast_B3_soup = setDR(F_Fast_B3_soup, F_Fast_B3_meta[colnames(F_Fast_B3_flt), c("UMAP_1", "UMAP_2")])
F_Fast_B3_soup = setClusters(F_Fast_B3_soup, setNames(F_Fast_B3_meta$seurat_clusters, colnames(F_Fast_B3_flt)))


F_Fast_B3_soup = autoEstCont(F_Fast_B3_soup)
# rho estimated at 0.04; however, 0.10 producing cleaner looking clusters

F_Fast_B3_soup <- setContaminationFraction(F_Fast_B3_soup, 0.10)
F_Fast_B3out <- adjustCounts(F_Fast_B3_soup, method = 'multinomial')

# SoupX finished, check results
F_Fast_B3out <- CreateSeuratObject(F_Fast_B3out)

F_Fast_B3out <- SCTransform(F_Fast_B3out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fast_B3out <- FindVariableFeatures(F_Fast_B3out)
F_Fast_B3out@assays$SCT@var.features <- F_Fast_B3out@assays$SCT@var.features[(!F_Fast_B3out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fast_B3out <- RunPCA(F_Fast_B3out, features = VariableFeatures(object = F_Fast_B3out))
F_Fast_B3out <- RunUMAP(F_Fast_B3out, reduction = "pca", dims = 1:30)

F_Fast_B3out <- FindNeighbors(F_Fast_B3out, dims = 1:30)
F_Fast_B3out <- FindClusters(F_Fast_B3out, resolution = 0.8)

DimPlot(F_Fast_B3out, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(F_Fast_B3out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B3_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fast_B3out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B3_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fast_B3_markers_post_soupx <- FindAllMarkers(F_Fast_B3out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#mito.genes <- grep(pattern = "^mt-", x = rownames(x = F_Fast_B3out@assays$RNA@data), value = TRUE)
mitoPercent <- Matrix::colSums(F_Fast_B3out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(F_Fast_B3out@assays$RNA@counts)*100
F_Fast_B3out$mitoPercent <- mitoPercent


hemoglobin <- Matrix::colSums(F_Fast_B3out@assays$RNA@counts[hemoglobin_genes, ])
F_Fast_B3out <- AddMetaData(object = F_Fast_B3out, metadata = hemoglobin, col.name = "hemoglobin")


ribo.genes <- grep(pattern = '^Rp[sl]', x = rownames(x = F_Fast_B3out@assays$RNA@data), value = TRUE)
riboPercent <- Matrix::colSums(F_Fast_B3out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(F_Fast_B3out@assays$RNA@counts)*100
F_Fast_B3out$riboPercent <- riboPercent

VlnPlot(F_Fast_B3out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fast_B3_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fast_B3out_filtered <- subset(F_Fast_B3out, 
                                subset = nCount_RNA < (mean(F_Fast_B3out@meta.data$nCount_RNA) + 4*sd(F_Fast_B3out@meta.data$nCount_RNA))
                                & nFeature_RNA < (mean(F_Fast_B3out@meta.data$nFeature_RNA) + 4*sd(F_Fast_B3out@meta.data$nFeature_RNA))
                                & mitoPercent < (mean(F_Fast_B3out@meta.data$mitoPercent) + 4*sd(F_Fast_B3out@meta.data$mitoPercent))
                                & riboPercent < (mean(F_Fast_B3out@meta.data$riboPercent) + 4*sd(F_Fast_B3out@meta.data$riboPercent))
                                & hemoglobin < (mean(F_Fast_B3out@meta.data$hemoglobin) + 4*sd(F_Fast_B3out@meta.data$hemoglobin)))


VlnPlot(F_Fast_B3out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fast_B3_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
DefaultAssay(F_Fast_B3out_filtered) <- 'RNA'
F_Fast_B3out_filtered <- DietSeurat(F_Fast_B3out_filtered, assays = 'RNA')
F_Fast_B3out_filtered <- as.SingleCellExperiment(F_Fast_B3out_filtered)

F_Fast_B3_DBR <- F_Fast_B3out_filtered@colData@nrows/50000
F_Fast_B3_scDbl <- scDblFinder(F_Fast_B3out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = F_Fast_B3_DBR, propRandom = 0.2, iter = 5)
F_Fast_B3_scDbl <- as.Seurat(F_Fast_B3_scDbl)

F_Fast_B3_scDbl <- SCTransform(F_Fast_B3_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fast_B3_scDbl <- FindVariableFeatures(F_Fast_B3_scDbl)
F_Fast_B3_scDbl@assays$SCT@var.features <- F_Fast_B3_scDbl@assays$SCT@var.features[(!F_Fast_B3_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fast_B3_scDbl <- RunPCA(F_Fast_B3_scDbl, features = VariableFeatures(object = F_Fast_B3_scDbl))
F_Fast_B3_scDbl <- RunUMAP(F_Fast_B3_scDbl, reduction = "pca", dims = 1:30)



F_Fast_B3_scDbl <- FindNeighbors(F_Fast_B3_scDbl, dims = 1:30)
F_Fast_B3_scDbl <- FindClusters(F_Fast_B3_scDbl, resolution = 0.8)



DimPlot(F_Fast_B3_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/F_Fast_B3_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)


F_Fast_B3_scDbl <- subset(F_Fast_B3_scDbl, scDblFinder.class == 'singlet')


F_Fast_B3_scDbl <- SCTransform(F_Fast_B3_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#F_Fast_B3_scDbl <- FindVariableFeatures(F_Fast_B3_scDbl)
F_Fast_B3_scDbl@assays$SCT@var.features <- F_Fast_B3_scDbl@assays$SCT@var.features[(!F_Fast_B3_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



F_Fast_B3_scDbl <- RunPCA(F_Fast_B3_scDbl, features = VariableFeatures(object = F_Fast_B3_scDbl))
F_Fast_B3_scDbl <- RunUMAP(F_Fast_B3_scDbl, reduction = "pca", dims = 1:30)



F_Fast_B3_scDbl <- FindNeighbors(F_Fast_B3_scDbl, dims = 1:30)
F_Fast_B3_scDbl <- FindClusters(F_Fast_B3_scDbl, resolution = 1.5)

DimPlot(F_Fast_B3_scDbl, label = TRUE) + NoLegend()
#check purity of clusters
FeaturePlot(F_Fast_B3_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B3_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(F_Fast_B3_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B3_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

F_Fast_B3_markers_post_scDblFnd <- FindAllMarkers(F_Fast_B3_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#FindMarkers(F_Fast_B3_scDbl, ident.1 = c('4','8'), features = c('Agrp','Npy'), only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save F_Fast_B3_scDbl, markers
rm(F_Fast_B3_flt, F_Fast_B3_meta, F_Fast_B3_raw, F_Fast_B3_soup, F_Fast_B3_UMAP, 
   F_Fast_B3out, F_Fast_B3out_filtered, F_Fast_B3_singlets)


#M_Fed_B1
#create basic clustering for sample
M_Fed_B1_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B1Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fed_B1_flt <- CreateSeuratObject(M_Fed_B1_flt)

M_Fed_B1_flt <- SCTransform(M_Fed_B1_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fed_B1_flt <- FindVariableFeatures(M_Fed_B1_flt)
M_Fed_B1_flt@assays$SCT@var.features <- M_Fed_B1_flt@assays$SCT@var.features[(!M_Fed_B1_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

M_Fed_B1_flt <- RunPCA(M_Fed_B1_flt, features = VariableFeatures(object = M_Fed_B1_flt))
M_Fed_B1_flt <- RunUMAP(M_Fed_B1_flt, reduction = "pca", dims = 1:30)

M_Fed_B1_flt <- FindNeighbors(M_Fed_B1_flt, dims = 1:30)
M_Fed_B1_flt <- FindClusters(M_Fed_B1_flt, resolution = 0.8)


DimPlot(M_Fed_B1_flt, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(M_Fed_B1_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B1_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fed_B1_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B1_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fed_B1_markers_b4_soupx <- FindAllMarkers(M_Fed_B1_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

M_Fed_B1_UMAP <- M_Fed_B1_flt@reductions$umap@cell.embeddings |> as.data.frame()
M_Fed_B1_meta <- M_Fed_B1_flt@meta.data[,c(2,3,7)] |> cbind(M_Fed_B1_UMAP)



#start SoupX in earnest 
M_Fed_B1_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B1Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fed_B1_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B1Solo.out/GeneFull_Ex50pAS/raw/')
M_Fed_B1_soup <-SoupChannel(M_Fed_B1_raw, M_Fed_B1_flt)


M_Fed_B1_soup = setDR(M_Fed_B1_soup, M_Fed_B1_meta[colnames(M_Fed_B1_flt), c("UMAP_1", "UMAP_2")])
M_Fed_B1_soup = setClusters(M_Fed_B1_soup, setNames(M_Fed_B1_meta$seurat_clusters, colnames(M_Fed_B1_flt)))


M_Fed_B1_soup = autoEstCont(M_Fed_B1_soup)
# rho estimated at 0.06; however, 0.15 producing cleaner looking clusters

M_Fed_B1_soup <- setContaminationFraction(M_Fed_B1_soup, 0.15)
M_Fed_B1out <- adjustCounts(M_Fed_B1_soup, method = 'multinomial')

# SoupX finished, check results
M_Fed_B1out <- CreateSeuratObject(M_Fed_B1out)

M_Fed_B1out <- SCTransform(M_Fed_B1out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fed_B1out <- FindVariableFeatures(M_Fed_B1out)
M_Fed_B1out@assays$SCT@var.features <- M_Fed_B1out@assays$SCT@var.features[(!M_Fed_B1out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fed_B1out <- RunPCA(M_Fed_B1out, features = VariableFeatures(object = M_Fed_B1out))
M_Fed_B1out <- RunUMAP(M_Fed_B1out, reduction = "pca", dims = 1:30)

M_Fed_B1out <- FindNeighbors(M_Fed_B1out, dims = 1:30)
M_Fed_B1out <- FindClusters(M_Fed_B1out, resolution = 0.8)

DimPlot(M_Fed_B1out, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(M_Fed_B1out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B1_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fed_B1out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B1_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fed_B1_markers_post_soupx <- FindAllMarkers(M_Fed_B1out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#mito.genes <- grep(pattern = "^mt-", x = rownames(x = M_Fed_B1out@assays$RNA@data), value = TRUE)
mitoPercent <- Matrix::colSums(M_Fed_B1out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(M_Fed_B1out@assays$RNA@counts)*100
M_Fed_B1out$mitoPercent <- mitoPercent


hemoglobin <- Matrix::colSums(M_Fed_B1out@assays$RNA@counts[hemoglobin_genes, ])
M_Fed_B1out <- AddMetaData(object = M_Fed_B1out, metadata = hemoglobin, col.name = "hemoglobin")


ribo.genes <- grep(pattern = '^Rp[sl]', x = rownames(x = M_Fed_B1out@assays$RNA@data), value = TRUE)
riboPercent <- Matrix::colSums(M_Fed_B1out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(M_Fed_B1out@assays$RNA@counts)*100
M_Fed_B1out$riboPercent <- riboPercent

VlnPlot(M_Fed_B1out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fed_B1_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fed_B1out_filtered <- subset(M_Fed_B1out, 
                               subset = nCount_RNA < (mean(M_Fed_B1out@meta.data$nCount_RNA) + 4*sd(M_Fed_B1out@meta.data$nCount_RNA))
                               & nFeature_RNA < (mean(M_Fed_B1out@meta.data$nFeature_RNA) + 4*sd(M_Fed_B1out@meta.data$nFeature_RNA))
                               & mitoPercent < (mean(M_Fed_B1out@meta.data$mitoPercent) + 4*sd(M_Fed_B1out@meta.data$mitoPercent))
                               & riboPercent < (mean(M_Fed_B1out@meta.data$riboPercent) + 4*sd(M_Fed_B1out@meta.data$riboPercent))
                               & hemoglobin < (mean(M_Fed_B1out@meta.data$hemoglobin) + 4*sd(M_Fed_B1out@meta.data$hemoglobin)))


VlnPlot(M_Fed_B1out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fed_B1_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
DefaultAssay(M_Fed_B1out_filtered) <- 'RNA'
M_Fed_B1out_filtered <- DietSeurat(M_Fed_B1out_filtered, assays = 'RNA')
M_Fed_B1out_filtered <- as.SingleCellExperiment(M_Fed_B1out_filtered)

M_Fed_B1_DBR <- M_Fed_B1out_filtered@colData@nrows/50000
M_Fed_B1_scDbl <- scDblFinder(M_Fed_B1out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = M_Fed_B1_DBR, propRandom = 0.2, iter = 5)
M_Fed_B1_scDbl <- as.Seurat(M_Fed_B1_scDbl)

M_Fed_B1_scDbl <- SCTransform(M_Fed_B1_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fed_B1_scDbl <- FindVariableFeatures(M_Fed_B1_scDbl)
M_Fed_B1_scDbl@assays$SCT@var.features <- M_Fed_B1_scDbl@assays$SCT@var.features[(!M_Fed_B1_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fed_B1_scDbl <- RunPCA(M_Fed_B1_scDbl, features = VariableFeatures(object = M_Fed_B1_scDbl))
M_Fed_B1_scDbl <- RunUMAP(M_Fed_B1_scDbl, reduction = "pca", dims = 1:30)



M_Fed_B1_scDbl <- FindNeighbors(M_Fed_B1_scDbl, dims = 1:30)
M_Fed_B1_scDbl <- FindClusters(M_Fed_B1_scDbl, resolution = 0.8)



DimPlot(M_Fed_B1_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/M_Fed_B1_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)


M_Fed_B1_scDbl <- subset(M_Fed_B1_scDbl, scDblFinder.class == 'singlet')


M_Fed_B1_scDbl <- SCTransform(M_Fed_B1_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fed_B1_scDbl <- FindVariableFeatures(M_Fed_B1_scDbl)
M_Fed_B1_scDbl@assays$SCT@var.features <- M_Fed_B1_scDbl@assays$SCT@var.features[(!M_Fed_B1_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fed_B1_scDbl <- RunPCA(M_Fed_B1_scDbl, features = VariableFeatures(object = M_Fed_B1_scDbl))
M_Fed_B1_scDbl <- RunUMAP(M_Fed_B1_scDbl, reduction = "pca", dims = 1:30)



M_Fed_B1_scDbl <- FindNeighbors(M_Fed_B1_scDbl, dims = 1:30)
M_Fed_B1_scDbl <- FindClusters(M_Fed_B1_scDbl, resolution = 0.8)

#check purity of clusters
FeaturePlot(M_Fed_B1_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B1_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fed_B1_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B1_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fed_B1_markers_post_scDblFnd <- FindAllMarkers(M_Fed_B1_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#FindMarkers(M_Fed_B1_scDbl, ident.1 = c('2','3','5', '10'), features = c('Slc1a3','Slc1a2','Gfap'), only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)
#FindMarkers(M_Fed_B1out, ident.1 = c('0','2','5'), features = c('Slc1a3','Slc1a2','Gfap'), only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# remove unneeded objects, save M_Fed_B1_scDbl, markers
rm(M_Fed_B1_flt, M_Fed_B1_meta, M_Fed_B1_raw, M_Fed_B1_soup, M_Fed_B1_UMAP, 
   M_Fed_B1out, M_Fed_B1out_filtered, M_Fed_B1_singlets)


#M_Fed_B2
#create basic clustering for sample
M_Fed_B2_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B2Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fed_B2_flt <- CreateSeuratObject(M_Fed_B2_flt)

M_Fed_B2_flt <- SCTransform(M_Fed_B2_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fed_B2_flt <- FindVariableFeatures(M_Fed_B2_flt)
M_Fed_B2_flt@assays$SCT@var.features <- M_Fed_B2_flt@assays$SCT@var.features[(!M_Fed_B2_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

M_Fed_B2_flt <- RunPCA(M_Fed_B2_flt, features = VariableFeatures(object = M_Fed_B2_flt))
M_Fed_B2_flt <- RunUMAP(M_Fed_B2_flt, reduction = "pca", dims = 1:30)

M_Fed_B2_flt <- FindNeighbors(M_Fed_B2_flt, dims = 1:30)
M_Fed_B2_flt <- FindClusters(M_Fed_B2_flt, resolution = 0.8)


DimPlot(M_Fed_B2_flt, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(M_Fed_B2_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B2_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fed_B2_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B2_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fed_B2_markers_b4_soupx <- FindAllMarkers(M_Fed_B2_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

M_Fed_B2_UMAP <- M_Fed_B2_flt@reductions$umap@cell.embeddings |> as.data.frame()
M_Fed_B2_meta <- M_Fed_B2_flt@meta.data[,c(2,3,7)] |> cbind(M_Fed_B2_UMAP)



#start SoupX in earnest 
M_Fed_B2_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B2Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fed_B2_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B2Solo.out/GeneFull_Ex50pAS/raw/')
M_Fed_B2_soup <-SoupChannel(M_Fed_B2_raw, M_Fed_B2_flt)


M_Fed_B2_soup = setDR(M_Fed_B2_soup, M_Fed_B2_meta[colnames(M_Fed_B2_flt), c("UMAP_1", "UMAP_2")])
M_Fed_B2_soup = setClusters(M_Fed_B2_soup, setNames(M_Fed_B2_meta$seurat_clusters, colnames(M_Fed_B2_flt)))


M_Fed_B2_soup = autoEstCont(M_Fed_B2_soup)
# rho estimated at 0.03; however, 0.105 producing cleaner looking clusters

M_Fed_B2_soup <- setContaminationFraction(M_Fed_B2_soup, 0.105)
M_Fed_B2out <- adjustCounts(M_Fed_B2_soup, method = 'multinomial')

# SoupX finished, check results
M_Fed_B2out <- CreateSeuratObject(M_Fed_B2out)

M_Fed_B2out <- SCTransform(M_Fed_B2out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fed_B2out <- FindVariableFeatures(M_Fed_B2out)
M_Fed_B2out@assays$SCT@var.features <- M_Fed_B2out@assays$SCT@var.features[(!M_Fed_B2out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fed_B2out <- RunPCA(M_Fed_B2out, features = VariableFeatures(object = M_Fed_B2out))
M_Fed_B2out <- RunUMAP(M_Fed_B2out, reduction = "pca", dims = 1:30)

M_Fed_B2out <- FindNeighbors(M_Fed_B2out, dims = 1:30)
M_Fed_B2out <- FindClusters(M_Fed_B2out, resolution = 1.5)

DimPlot(M_Fed_B2out, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(M_Fed_B2out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B2_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fed_B2out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B2_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fed_B2_markers_post_soupx <- FindAllMarkers(M_Fed_B2out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#mito.genes <- grep(pattern = "^mt-", x = rownames(x = M_Fed_B2out@assays$RNA@data), value = TRUE)
mitoPercent <- Matrix::colSums(M_Fed_B2out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(M_Fed_B2out@assays$RNA@counts)*100
M_Fed_B2out$mitoPercent <- mitoPercent


hemoglobin <- Matrix::colSums(M_Fed_B2out@assays$RNA@counts[hemoglobin_genes, ])
M_Fed_B2out <- AddMetaData(object = M_Fed_B2out, metadata = hemoglobin, col.name = "hemoglobin")


ribo.genes <- grep(pattern = '^Rp[sl]', x = rownames(x = M_Fed_B2out@assays$RNA@data), value = TRUE)
riboPercent <- Matrix::colSums(M_Fed_B2out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(M_Fed_B2out@assays$RNA@counts)*100
M_Fed_B2out$riboPercent <- riboPercent

VlnPlot(M_Fed_B2out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fed_B2_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fed_B2out_filtered <- subset(M_Fed_B2out, 
                               subset = nCount_RNA < (mean(M_Fed_B2out@meta.data$nCount_RNA) + 4*sd(M_Fed_B2out@meta.data$nCount_RNA))
                               & nFeature_RNA < (mean(M_Fed_B2out@meta.data$nFeature_RNA) + 4*sd(M_Fed_B2out@meta.data$nFeature_RNA))
                               & mitoPercent < (mean(M_Fed_B2out@meta.data$mitoPercent) + 4*sd(M_Fed_B2out@meta.data$mitoPercent))
                               & riboPercent < (mean(M_Fed_B2out@meta.data$riboPercent) + 4*sd(M_Fed_B2out@meta.data$riboPercent))
                               & hemoglobin < (mean(M_Fed_B2out@meta.data$hemoglobin) + 4*sd(M_Fed_B2out@meta.data$hemoglobin)))


VlnPlot(M_Fed_B2out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fed_B2_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
DefaultAssay(M_Fed_B2out_filtered) <- 'RNA'
M_Fed_B2out_filtered <- DietSeurat(M_Fed_B2out_filtered, assays = 'RNA')
M_Fed_B2out_filtered <- as.SingleCellExperiment(M_Fed_B2out_filtered)

M_Fed_B2_DBR <- M_Fed_B2out_filtered@colData@nrows/50000
M_Fed_B2_scDbl <- scDblFinder(M_Fed_B2out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = M_Fed_B2_DBR, propRandom = 0.2, iter = 5)
M_Fed_B2_scDbl <- as.Seurat(M_Fed_B2_scDbl)

M_Fed_B2_scDbl <- SCTransform(M_Fed_B2_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fed_B2_scDbl <- FindVariableFeatures(M_Fed_B2_scDbl)
M_Fed_B2_scDbl@assays$SCT@var.features <- M_Fed_B2_scDbl@assays$SCT@var.features[(!M_Fed_B2_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fed_B2_scDbl <- RunPCA(M_Fed_B2_scDbl, features = VariableFeatures(object = M_Fed_B2_scDbl))
M_Fed_B2_scDbl <- RunUMAP(M_Fed_B2_scDbl, reduction = "pca", dims = 1:30)



M_Fed_B2_scDbl <- FindNeighbors(M_Fed_B2_scDbl, dims = 1:30)
M_Fed_B2_scDbl <- FindClusters(M_Fed_B2_scDbl, resolution = 0.8)



DimPlot(M_Fed_B2_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/M_Fed_B2_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)


M_Fed_B2_scDbl <- subset(M_Fed_B2_scDbl, scDblFinder.class == 'singlet')


M_Fed_B2_scDbl <- SCTransform(M_Fed_B2_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fed_B2_scDbl <- FindVariableFeatures(M_Fed_B2_scDbl)
M_Fed_B2_scDbl@assays$SCT@var.features <- M_Fed_B2_scDbl@assays$SCT@var.features[(!M_Fed_B2_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fed_B2_scDbl <- RunPCA(M_Fed_B2_scDbl, features = VariableFeatures(object = M_Fed_B2_scDbl))
M_Fed_B2_scDbl <- RunUMAP(M_Fed_B2_scDbl, reduction = "pca", dims = 1:30)



M_Fed_B2_scDbl <- FindNeighbors(M_Fed_B2_scDbl, dims = 1:30)
M_Fed_B2_scDbl <- FindClusters(M_Fed_B2_scDbl, resolution = 1.2)
DimPlot(M_Fed_B2_scDbl, label = TRUE) + NoLegend()
#check purity of clusters
FeaturePlot(M_Fed_B2_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B2_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fed_B2_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B2_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fed_B2_markers_post_scDblFnd <- FindAllMarkers(M_Fed_B2_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#FindMarkers(M_Fed_B2_scDbl, ident.1 = c('0','2'), features = c('Slc1a2','Gfap'), only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)
#FindMarkers(M_Fed_B2out, ident.1 = c('0','2','26'), features = c('Slc1a2','Gfap'), only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#FindMarkers(M_Fed_B2_scDbl, ident.1 = c('4'), features = c('Mbp'), only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)
#FindMarkers(M_Fed_B2out, ident.1 = c('4','34','26'), features = c('Mbp'), only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save M_Fed_B2_scDbl, markers
rm(M_Fed_B2_flt, M_Fed_B2_meta, M_Fed_B2_raw, M_Fed_B2_soup, M_Fed_B2_UMAP, 
   M_Fed_B2out, M_Fed_B2out_filtered, M_Fed_B2_singlets)


#M_Fed_B3
#create basic clustering for sample
M_Fed_B3_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B3Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fed_B3_flt <- CreateSeuratObject(M_Fed_B3_flt)

M_Fed_B3_flt <- SCTransform(M_Fed_B3_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fed_B3_flt <- FindVariableFeatures(M_Fed_B3_flt)
M_Fed_B3_flt@assays$SCT@var.features <- M_Fed_B3_flt@assays$SCT@var.features[(!M_Fed_B3_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

M_Fed_B3_flt <- RunPCA(M_Fed_B3_flt, features = VariableFeatures(object = M_Fed_B3_flt))
M_Fed_B3_flt <- RunUMAP(M_Fed_B3_flt, reduction = "pca", dims = 1:30)

M_Fed_B3_flt <- FindNeighbors(M_Fed_B3_flt, dims = 1:30)
M_Fed_B3_flt <- FindClusters(M_Fed_B3_flt, resolution = 0.8)


DimPlot(M_Fed_B3_flt, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(M_Fed_B3_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B3_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fed_B3_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B3_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fed_B3_markers_b4_soupx <- FindAllMarkers(M_Fed_B3_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

M_Fed_B3_UMAP <- M_Fed_B3_flt@reductions$umap@cell.embeddings |> as.data.frame()
M_Fed_B3_meta <- M_Fed_B3_flt@meta.data[,c(2,3,7)] |> cbind(M_Fed_B3_UMAP)



#start SoupX in earnest 
M_Fed_B3_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B3Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fed_B3_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B3Solo.out/GeneFull_Ex50pAS/raw/')
M_Fed_B3_soup <-SoupChannel(M_Fed_B3_raw, M_Fed_B3_flt)


M_Fed_B3_soup = setDR(M_Fed_B3_soup, M_Fed_B3_meta[colnames(M_Fed_B3_flt), c("UMAP_1", "UMAP_2")])
M_Fed_B3_soup = setClusters(M_Fed_B3_soup, setNames(M_Fed_B3_meta$seurat_clusters, colnames(M_Fed_B3_flt)))


M_Fed_B3_soup = autoEstCont(M_Fed_B3_soup)
# rho estimated at 0.07; however, 0.20 producing cleaner looking clusters

M_Fed_B3_soup <- setContaminationFraction(M_Fed_B3_soup, 0.20)
M_Fed_B3out <- adjustCounts(M_Fed_B3_soup, method = 'multinomial')

# SoupX finished, check results
M_Fed_B3out <- CreateSeuratObject(M_Fed_B3out)

M_Fed_B3out <- SCTransform(M_Fed_B3out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fed_B3out <- FindVariableFeatures(M_Fed_B3out)
M_Fed_B3out@assays$SCT@var.features <- M_Fed_B3out@assays$SCT@var.features[(!M_Fed_B3out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fed_B3out <- RunPCA(M_Fed_B3out, features = VariableFeatures(object = M_Fed_B3out))
M_Fed_B3out <- RunUMAP(M_Fed_B3out, reduction = "pca", dims = 1:30)

M_Fed_B3out <- FindNeighbors(M_Fed_B3out, dims = 1:30)
M_Fed_B3out <- FindClusters(M_Fed_B3out, resolution = 1.5)

DimPlot(M_Fed_B3out, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(M_Fed_B3out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B3_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fed_B3out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B3_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fed_B3_markers_post_soupx <- FindAllMarkers(M_Fed_B3out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#mito.genes <- grep(pattern = "^mt-", x = rownames(x = M_Fed_B3out@assays$RNA@data), value = TRUE)
mitoPercent <- Matrix::colSums(M_Fed_B3out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(M_Fed_B3out@assays$RNA@counts)*100
M_Fed_B3out$mitoPercent <- mitoPercent


hemoglobin <- Matrix::colSums(M_Fed_B3out@assays$RNA@counts[hemoglobin_genes, ])
M_Fed_B3out <- AddMetaData(object = M_Fed_B3out, metadata = hemoglobin, col.name = "hemoglobin")


ribo.genes <- grep(pattern = '^Rp[sl]', x = rownames(x = M_Fed_B3out@assays$RNA@data), value = TRUE)
riboPercent <- Matrix::colSums(M_Fed_B3out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(M_Fed_B3out@assays$RNA@counts)*100
M_Fed_B3out$riboPercent <- riboPercent

VlnPlot(M_Fed_B3out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fed_B3_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fed_B3out_filtered <- subset(M_Fed_B3out, 
                               subset = nCount_RNA < (mean(M_Fed_B3out@meta.data$nCount_RNA) + 4*sd(M_Fed_B3out@meta.data$nCount_RNA))
                               & nFeature_RNA < (mean(M_Fed_B3out@meta.data$nFeature_RNA) + 4*sd(M_Fed_B3out@meta.data$nFeature_RNA))
                               & mitoPercent < (mean(M_Fed_B3out@meta.data$mitoPercent) + 4*sd(M_Fed_B3out@meta.data$mitoPercent))
                               & riboPercent < (mean(M_Fed_B3out@meta.data$riboPercent) + 4*sd(M_Fed_B3out@meta.data$riboPercent))
                               & hemoglobin < (mean(M_Fed_B3out@meta.data$hemoglobin) + 4*sd(M_Fed_B3out@meta.data$hemoglobin)))


VlnPlot(M_Fed_B3out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fed_B3_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
DefaultAssay(M_Fed_B3out_filtered) <- 'RNA'
M_Fed_B3out_filtered <- DietSeurat(M_Fed_B3out_filtered, assays = 'RNA')
M_Fed_B3out_filtered <- as.SingleCellExperiment(M_Fed_B3out_filtered)

M_Fed_B3_DBR <- M_Fed_B3out_filtered@colData@nrows/50000
M_Fed_B3_scDbl <- scDblFinder(M_Fed_B3out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = M_Fed_B3_DBR, propRandom = 0.2, iter = 5)
M_Fed_B3_scDbl <- as.Seurat(M_Fed_B3_scDbl)

M_Fed_B3_scDbl <- SCTransform(M_Fed_B3_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fed_B3_scDbl <- FindVariableFeatures(M_Fed_B3_scDbl)
M_Fed_B3_scDbl@assays$SCT@var.features <- M_Fed_B3_scDbl@assays$SCT@var.features[(!M_Fed_B3_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fed_B3_scDbl <- RunPCA(M_Fed_B3_scDbl, features = VariableFeatures(object = M_Fed_B3_scDbl))
M_Fed_B3_scDbl <- RunUMAP(M_Fed_B3_scDbl, reduction = "pca", dims = 1:30)



#M_Fed_B3_scDbl <- FindNeighbors(M_Fed_B3_scDbl, dims = 1:30)
#M_Fed_B3_scDbl <- FindClusters(M_Fed_B3_scDbl, resolution = 0.8)



DimPlot(M_Fed_B3_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/M_Fed_B3_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)


M_Fed_B3_scDbl <- subset(M_Fed_B3_scDbl, scDblFinder.class == 'singlet')


M_Fed_B3_scDbl <- SCTransform(M_Fed_B3_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fed_B3_scDbl <- FindVariableFeatures(M_Fed_B3_scDbl)
M_Fed_B3_scDbl@assays$SCT@var.features <- M_Fed_B3_scDbl@assays$SCT@var.features[(!M_Fed_B3_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fed_B3_scDbl <- RunPCA(M_Fed_B3_scDbl, features = VariableFeatures(object = M_Fed_B3_scDbl))
M_Fed_B3_scDbl <- RunUMAP(M_Fed_B3_scDbl, reduction = "pca", dims = 1:30)



M_Fed_B3_scDbl <- FindNeighbors(M_Fed_B3_scDbl, dims = 1:30)
M_Fed_B3_scDbl <- FindClusters(M_Fed_B3_scDbl, resolution = 1.2)
DimPlot(M_Fed_B3_scDbl, label = TRUE) + NoLegend()
#check purity of clusters
FeaturePlot(M_Fed_B3_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B3_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fed_B3_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B3_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fed_B3_markers_post_scDblFnd <- FindAllMarkers(M_Fed_B3_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#FindMarkers(M_Fed_B3_scDbl, ident.1 = c('4','8'), features = c('Agrp','Npy'), only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save M_Fed_B3_scDbl, markers
rm(M_Fed_B3_flt, M_Fed_B3_meta, M_Fed_B3_raw, M_Fed_B3_soup, M_Fed_B3_UMAP, 
   M_Fed_B3out, M_Fed_B3out_filtered, M_Fed_B3_singlets)




#M_Fast_B1
#create basic clustering for sample
M_Fast_B1_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B1Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fast_B1_flt <- CreateSeuratObject(M_Fast_B1_flt)

M_Fast_B1_flt <- SCTransform(M_Fast_B1_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fast_B1_flt <- FindVariableFeatures(M_Fast_B1_flt)
M_Fast_B1_flt@assays$SCT@var.features <- M_Fast_B1_flt@assays$SCT@var.features[(!M_Fast_B1_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

M_Fast_B1_flt <- RunPCA(M_Fast_B1_flt, features = VariableFeatures(object = M_Fast_B1_flt))
M_Fast_B1_flt <- RunUMAP(M_Fast_B1_flt, reduction = "pca", dims = 1:30)

M_Fast_B1_flt <- FindNeighbors(M_Fast_B1_flt, dims = 1:30)
M_Fast_B1_flt <- FindClusters(M_Fast_B1_flt, resolution = 0.8)


DimPlot(M_Fast_B1_flt, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(M_Fast_B1_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B1_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fast_B1_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B1_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fast_B1_markers_b4_soupx <- FindAllMarkers(M_Fast_B1_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

M_Fast_B1_UMAP <- M_Fast_B1_flt@reductions$umap@cell.embeddings |> as.data.frame()
M_Fast_B1_meta <- M_Fast_B1_flt@meta.data[,c(2,3,7)] |> cbind(M_Fast_B1_UMAP)



#start SoupX in earnest 
M_Fast_B1_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B1Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fast_B1_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B1Solo.out/GeneFull_Ex50pAS/raw/')
M_Fast_B1_soup <-SoupChannel(M_Fast_B1_raw, M_Fast_B1_flt)


M_Fast_B1_soup = setDR(M_Fast_B1_soup, M_Fast_B1_meta[colnames(M_Fast_B1_flt), c("UMAP_1", "UMAP_2")])
M_Fast_B1_soup = setClusters(M_Fast_B1_soup, setNames(M_Fast_B1_meta$seurat_clusters, colnames(M_Fast_B1_flt)))


M_Fast_B1_soup = autoEstCont(M_Fast_B1_soup)
# rho estimated at 0.07; however, 0.175 producing cleaner looking clusters

M_Fast_B1_soup <- setContaminationFraction(M_Fast_B1_soup, 0.175)
M_Fast_B1out <- adjustCounts(M_Fast_B1_soup, method = 'multinomial')

# SoupX finished, check results
M_Fast_B1out <- CreateSeuratObject(M_Fast_B1out)

M_Fast_B1out <- SCTransform(M_Fast_B1out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fast_B1out <- FindVariableFeatures(M_Fast_B1out)
M_Fast_B1out@assays$SCT@var.features <- M_Fast_B1out@assays$SCT@var.features[(!M_Fast_B1out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fast_B1out <- RunPCA(M_Fast_B1out, features = VariableFeatures(object = M_Fast_B1out))
M_Fast_B1out <- RunUMAP(M_Fast_B1out, reduction = "pca", dims = 1:30)

M_Fast_B1out <- FindNeighbors(M_Fast_B1out, dims = 1:30)
M_Fast_B1out <- FindClusters(M_Fast_B1out, resolution = 1.5)

DimPlot(M_Fast_B1out, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(M_Fast_B1out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B1_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fast_B1out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B1_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fast_B1_markers_post_soupx <- FindAllMarkers(M_Fast_B1out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#mito.genes <- grep(pattern = "^mt-", x = rownames(x = M_Fast_B1out@assays$RNA@data), value = TRUE)
mitoPercent <- Matrix::colSums(M_Fast_B1out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(M_Fast_B1out@assays$RNA@counts)*100
M_Fast_B1out$mitoPercent <- mitoPercent


hemoglobin <- Matrix::colSums(M_Fast_B1out@assays$RNA@counts[hemoglobin_genes, ])
M_Fast_B1out <- AddMetaData(object = M_Fast_B1out, metadata = hemoglobin, col.name = "hemoglobin")


ribo.genes <- grep(pattern = '^Rp[sl]', x = rownames(x = M_Fast_B1out@assays$RNA@data), value = TRUE)
riboPercent <- Matrix::colSums(M_Fast_B1out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(M_Fast_B1out@assays$RNA@counts)*100
M_Fast_B1out$riboPercent <- riboPercent

VlnPlot(M_Fast_B1out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fast_B1_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fast_B1out_filtered <- subset(M_Fast_B1out, 
                                subset = nCount_RNA < (mean(M_Fast_B1out@meta.data$nCount_RNA) + 4*sd(M_Fast_B1out@meta.data$nCount_RNA))
                                & nFeature_RNA < (mean(M_Fast_B1out@meta.data$nFeature_RNA) + 4*sd(M_Fast_B1out@meta.data$nFeature_RNA))
                                & mitoPercent < (mean(M_Fast_B1out@meta.data$mitoPercent) + 4*sd(M_Fast_B1out@meta.data$mitoPercent))
                                & riboPercent < (mean(M_Fast_B1out@meta.data$riboPercent) + 4*sd(M_Fast_B1out@meta.data$riboPercent))
                                & hemoglobin < (mean(M_Fast_B1out@meta.data$hemoglobin) + 4*sd(M_Fast_B1out@meta.data$hemoglobin)))


VlnPlot(M_Fast_B1out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fast_B1_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
DefaultAssay(M_Fast_B1out_filtered) <- 'RNA'
M_Fast_B1out_filtered <- DietSeurat(M_Fast_B1out_filtered, assays = 'RNA')
M_Fast_B1out_filtered <- as.SingleCellExperiment(M_Fast_B1out_filtered)

M_Fast_B1_DBR <- M_Fast_B1out_filtered@colData@nrows/50000
#M_Fast_B1_scDbl <- scDblFinder(M_Fast_B1out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = 0.1, propRandom = 0.2, iter = 5)
M_Fast_B1_scDbl <- scDblFinder(M_Fast_B1out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = M_Fast_B1_DBR, propRandom = 0.2, iter = 5)
M_Fast_B1_scDbl <- as.Seurat(M_Fast_B1_scDbl)

M_Fast_B1_scDbl <- SCTransform(M_Fast_B1_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fast_B1_scDbl <- FindVariableFeatures(M_Fast_B1_scDbl)
M_Fast_B1_scDbl@assays$SCT@var.features <- M_Fast_B1_scDbl@assays$SCT@var.features[(!M_Fast_B1_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fast_B1_scDbl <- RunPCA(M_Fast_B1_scDbl, features = VariableFeatures(object = M_Fast_B1_scDbl))
M_Fast_B1_scDbl <- RunUMAP(M_Fast_B1_scDbl, reduction = "pca", dims = 1:30)



#M_Fast_B1_scDbl <- FindNeighbors(M_Fast_B1_scDbl, dims = 1:30)
#M_Fast_B1_scDbl <- FindClusters(M_Fast_B1_scDbl, resolution = 0.8)



DimPlot(M_Fast_B1_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/M_Fast_B1_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)


M_Fast_B1_scDbl <- subset(M_Fast_B1_scDbl, scDblFinder.class == 'singlet')


M_Fast_B1_scDbl <- SCTransform(M_Fast_B1_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fast_B1_scDbl <- FindVariableFeatures(M_Fast_B1_scDbl)
M_Fast_B1_scDbl@assays$SCT@var.features <- M_Fast_B1_scDbl@assays$SCT@var.features[(!M_Fast_B1_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fast_B1_scDbl <- RunPCA(M_Fast_B1_scDbl, features = VariableFeatures(object = M_Fast_B1_scDbl))
M_Fast_B1_scDbl <- RunUMAP(M_Fast_B1_scDbl, reduction = "pca", dims = 1:30)



M_Fast_B1_scDbl <- FindNeighbors(M_Fast_B1_scDbl, dims = 1:30)
M_Fast_B1_scDbl <- FindClusters(M_Fast_B1_scDbl, resolution = 1.5)
DimPlot(M_Fast_B1_scDbl, label = TRUE) + NoLegend()
#check purity of clusters
FeaturePlot(M_Fast_B1_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B1_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fast_B1_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B1_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fast_B1_markers_post_scDblFnd <- FindAllMarkers(M_Fast_B1_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#FindMarkers(M_Fast_B1_scDbl, ident.1 = c('4','8'), features = c('Agrp','Npy'), only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save M_Fast_B1_scDbl, markers
rm(M_Fast_B1_flt, M_Fast_B1_meta, M_Fast_B1_raw, M_Fast_B1_soup, M_Fast_B1_UMAP, 
   M_Fast_B1out, M_Fast_B1out_filtered, M_Fast_B1_singlets)


#M_Fast_B2
#create basic clustering for sample
M_Fast_B2_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B2Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fast_B2_flt <- CreateSeuratObject(M_Fast_B2_flt)

M_Fast_B2_flt <- SCTransform(M_Fast_B2_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fast_B2_flt <- FindVariableFeatures(M_Fast_B2_flt)
M_Fast_B2_flt@assays$SCT@var.features <- M_Fast_B2_flt@assays$SCT@var.features[(!M_Fast_B2_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

M_Fast_B2_flt <- RunPCA(M_Fast_B2_flt, features = VariableFeatures(object = M_Fast_B2_flt))
M_Fast_B2_flt <- RunUMAP(M_Fast_B2_flt, reduction = "pca", dims = 1:30)

M_Fast_B2_flt <- FindNeighbors(M_Fast_B2_flt, dims = 1:30)
M_Fast_B2_flt <- FindClusters(M_Fast_B2_flt, resolution = 0.8)


DimPlot(M_Fast_B2_flt, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(M_Fast_B2_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B2_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fast_B2_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B2_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fast_B2_markers_b4_soupx <- FindAllMarkers(M_Fast_B2_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

M_Fast_B2_UMAP <- M_Fast_B2_flt@reductions$umap@cell.embeddings |> as.data.frame()
M_Fast_B2_meta <- M_Fast_B2_flt@meta.data[,c(2,3,7)] |> cbind(M_Fast_B2_UMAP)



#start SoupX in earnest 
M_Fast_B2_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B2Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fast_B2_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B2Solo.out/GeneFull_Ex50pAS/raw/')
M_Fast_B2_soup <-SoupChannel(M_Fast_B2_raw, M_Fast_B2_flt)


M_Fast_B2_soup = setDR(M_Fast_B2_soup, M_Fast_B2_meta[colnames(M_Fast_B2_flt), c("UMAP_1", "UMAP_2")])
M_Fast_B2_soup = setClusters(M_Fast_B2_soup, setNames(M_Fast_B2_meta$seurat_clusters, colnames(M_Fast_B2_flt)))


M_Fast_B2_soup = autoEstCont(M_Fast_B2_soup)
# rho estimated at 0.05; however, 0.175 producing cleaner looking clusters

M_Fast_B2_soup <- setContaminationFraction(M_Fast_B2_soup, 0.175)
M_Fast_B2out <- adjustCounts(M_Fast_B2_soup, method = 'multinomial')

# SoupX finished, check results
M_Fast_B2out <- CreateSeuratObject(M_Fast_B2out)

M_Fast_B2out <- SCTransform(M_Fast_B2out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fast_B2out <- FindVariableFeatures(M_Fast_B2out)
M_Fast_B2out@assays$SCT@var.features <- M_Fast_B2out@assays$SCT@var.features[(!M_Fast_B2out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fast_B2out <- RunPCA(M_Fast_B2out, features = VariableFeatures(object = M_Fast_B2out))
M_Fast_B2out <- RunUMAP(M_Fast_B2out, reduction = "pca", dims = 1:30)

M_Fast_B2out <- FindNeighbors(M_Fast_B2out, dims = 1:30)
M_Fast_B2out <- FindClusters(M_Fast_B2out, resolution = 1.5)

DimPlot(M_Fast_B2out, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(M_Fast_B2out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B2_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fast_B2out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B2_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fast_B2_markers_post_soupx <- FindAllMarkers(M_Fast_B2out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)
#FindMarkers(M_Fast_B2out, ident.1 = c('8'), features = c('Pomc','Lepr','Esr1','Ar'), only.pos = TRUE)

#FindMarkers(M_Fast_B2out, ident.1 = c('2','5','9','21','32'), features = c('Slc1a3'), only.pos = TRUE)
#mito.genes <- grep(pattern = "^mt-", x = rownames(x = M_Fast_B2out@assays$RNA@data), value = TRUE)
mitoPercent <- Matrix::colSums(M_Fast_B2out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(M_Fast_B2out@assays$RNA@counts)*100
M_Fast_B2out$mitoPercent <- mitoPercent


hemoglobin <- Matrix::colSums(M_Fast_B2out@assays$RNA@counts[hemoglobin_genes, ])
M_Fast_B2out <- AddMetaData(object = M_Fast_B2out, metadata = hemoglobin, col.name = "hemoglobin")


ribo.genes <- grep(pattern = '^Rp[sl]', x = rownames(x = M_Fast_B2out@assays$RNA@data), value = TRUE)
riboPercent <- Matrix::colSums(M_Fast_B2out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(M_Fast_B2out@assays$RNA@counts)*100
M_Fast_B2out$riboPercent <- riboPercent

VlnPlot(M_Fast_B2out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fast_B2_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fast_B2out_filtered <- subset(M_Fast_B2out, 
                                subset = nCount_RNA < (mean(M_Fast_B2out@meta.data$nCount_RNA) + 4*sd(M_Fast_B2out@meta.data$nCount_RNA))
                                & nFeature_RNA < (mean(M_Fast_B2out@meta.data$nFeature_RNA) + 4*sd(M_Fast_B2out@meta.data$nFeature_RNA))
                                & mitoPercent < (mean(M_Fast_B2out@meta.data$mitoPercent) + 4*sd(M_Fast_B2out@meta.data$mitoPercent))
                                & riboPercent < (mean(M_Fast_B2out@meta.data$riboPercent) + 4*sd(M_Fast_B2out@meta.data$riboPercent))
                                & hemoglobin < (mean(M_Fast_B2out@meta.data$hemoglobin) + 4*sd(M_Fast_B2out@meta.data$hemoglobin)))


VlnPlot(M_Fast_B2out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fast_B2_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
DefaultAssay(M_Fast_B2out_filtered) <- 'RNA'
M_Fast_B2out_filtered <- DietSeurat(M_Fast_B2out_filtered, assays = 'RNA')
M_Fast_B2out_filtered <- as.SingleCellExperiment(M_Fast_B2out_filtered)

M_Fast_B2_DBR <- M_Fast_B2out_filtered@colData@nrows/50000
M_Fast_B2_scDbl <- scDblFinder(M_Fast_B2out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = M_Fast_B2_DBR, propRandom = 0.2, iter = 5)
M_Fast_B2_scDbl <- as.Seurat(M_Fast_B2_scDbl)

M_Fast_B2_scDbl <- SCTransform(M_Fast_B2_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fast_B2_scDbl <- FindVariableFeatures(M_Fast_B2_scDbl)
M_Fast_B2_scDbl@assays$SCT@var.features <- M_Fast_B2_scDbl@assays$SCT@var.features[(!M_Fast_B2_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fast_B2_scDbl <- RunPCA(M_Fast_B2_scDbl, features = VariableFeatures(object = M_Fast_B2_scDbl))
M_Fast_B2_scDbl <- RunUMAP(M_Fast_B2_scDbl, reduction = "pca", dims = 1:30)



#M_Fast_B2_scDbl <- FindNeighbors(M_Fast_B2_scDbl, dims = 1:30)
#M_Fast_B2_scDbl <- FindClusters(M_Fast_B2_scDbl, resolution = 0.8)



DimPlot(M_Fast_B2_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/M_Fast_B2_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)


M_Fast_B2_scDbl <- subset(M_Fast_B2_scDbl, scDblFinder.class == 'singlet')


M_Fast_B2_scDbl <- SCTransform(M_Fast_B2_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fast_B2_scDbl <- FindVariableFeatures(M_Fast_B2_scDbl)
M_Fast_B2_scDbl@assays$SCT@var.features <- M_Fast_B2_scDbl@assays$SCT@var.features[(!M_Fast_B2_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fast_B2_scDbl <- RunPCA(M_Fast_B2_scDbl, features = VariableFeatures(object = M_Fast_B2_scDbl))
M_Fast_B2_scDbl <- RunUMAP(M_Fast_B2_scDbl, reduction = "pca", dims = 1:30)



M_Fast_B2_scDbl <- FindNeighbors(M_Fast_B2_scDbl, dims = 1:30)
M_Fast_B2_scDbl <- FindClusters(M_Fast_B2_scDbl, resolution = 1.5)
DimPlot(M_Fast_B2_scDbl, label = TRUE) + NoLegend()
#check purity of clusters
FeaturePlot(M_Fast_B2_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B2_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fast_B2_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B2_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fast_B2_markers_post_scDblFnd <- FindAllMarkers(M_Fast_B2_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#FindMarkers(M_Fast_B2_scDbl, ident.1 = c('7'), features = c('Pomc','Lepr'), only.pos = TRUE)


# remove unneeded objects, save M_Fast_B2_scDbl, markers
rm(M_Fast_B2_flt, M_Fast_B2_meta, M_Fast_B2_raw, M_Fast_B2_soup, M_Fast_B2_UMAP, 
   M_Fast_B2out, M_Fast_B2out_filtered, M_Fast_B2_singlets)


#M_Fast_B3
#create basic clustering for sample
M_Fast_B3_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B3Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fast_B3_flt <- CreateSeuratObject(M_Fast_B3_flt)

M_Fast_B3_flt <- SCTransform(M_Fast_B3_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fast_B3_flt <- FindVariableFeatures(M_Fast_B3_flt)
M_Fast_B3_flt@assays$SCT@var.features <- M_Fast_B3_flt@assays$SCT@var.features[(!M_Fast_B3_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

M_Fast_B3_flt <- RunPCA(M_Fast_B3_flt, features = VariableFeatures(object = M_Fast_B3_flt))
M_Fast_B3_flt <- RunUMAP(M_Fast_B3_flt, reduction = "pca", dims = 1:30)

M_Fast_B3_flt <- FindNeighbors(M_Fast_B3_flt, dims = 1:30)
M_Fast_B3_flt <- FindClusters(M_Fast_B3_flt, resolution = 0.8)


DimPlot(M_Fast_B3_flt, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(M_Fast_B3_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B3_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fast_B3_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B3_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fast_B3_markers_b4_soupx <- FindAllMarkers(M_Fast_B3_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

M_Fast_B3_UMAP <- M_Fast_B3_flt@reductions$umap@cell.embeddings |> as.data.frame()
M_Fast_B3_meta <- M_Fast_B3_flt@meta.data[,c(2,3,7)] |> cbind(M_Fast_B3_UMAP)



#start SoupX in earnest 
M_Fast_B3_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B3Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fast_B3_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B3Solo.out/GeneFull_Ex50pAS/raw/')
M_Fast_B3_soup <-SoupChannel(M_Fast_B3_raw, M_Fast_B3_flt)


M_Fast_B3_soup = setDR(M_Fast_B3_soup, M_Fast_B3_meta[colnames(M_Fast_B3_flt), c("UMAP_1", "UMAP_2")])
M_Fast_B3_soup = setClusters(M_Fast_B3_soup, setNames(M_Fast_B3_meta$seurat_clusters, colnames(M_Fast_B3_flt)))


M_Fast_B3_soup = autoEstCont(M_Fast_B3_soup)
# rho estimated at 0.08; however, 0.16 producing cleaner looking clusters

M_Fast_B3_soup <- setContaminationFraction(M_Fast_B3_soup, 0.16)
M_Fast_B3out <- adjustCounts(M_Fast_B3_soup, method = 'multinomial')

# SoupX finished, check results
M_Fast_B3out <- CreateSeuratObject(M_Fast_B3out)

M_Fast_B3out <- SCTransform(M_Fast_B3out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fast_B3out <- FindVariableFeatures(M_Fast_B3out)
M_Fast_B3out@assays$SCT@var.features <- M_Fast_B3out@assays$SCT@var.features[(!M_Fast_B3out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fast_B3out <- RunPCA(M_Fast_B3out, features = VariableFeatures(object = M_Fast_B3out))
M_Fast_B3out <- RunUMAP(M_Fast_B3out, reduction = "pca", dims = 1:30)

M_Fast_B3out <- FindNeighbors(M_Fast_B3out, dims = 1:30)
M_Fast_B3out <- FindClusters(M_Fast_B3out, resolution = 0.8)

DimPlot(M_Fast_B3out, label = TRUE, label.size = 2) + NoLegend()

FeaturePlot(M_Fast_B3out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B3_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fast_B3out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B3_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fast_B3_markers_post_soupx <- FindAllMarkers(M_Fast_B3out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)
#FindMarkers(M_Fast_B3out, ident.1 = c('0','2','20','26'), features = c('Slc1a3'), only.pos = TRUE)
#mito.genes <- grep(pattern = "^mt-", x = rownames(x = M_Fast_B3out@assays$RNA@data), value = TRUE)
mitoPercent <- Matrix::colSums(M_Fast_B3out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(M_Fast_B3out@assays$RNA@counts)*100
M_Fast_B3out$mitoPercent <- mitoPercent


hemoglobin <- Matrix::colSums(M_Fast_B3out@assays$RNA@counts[hemoglobin_genes, ])
M_Fast_B3out <- AddMetaData(object = M_Fast_B3out, metadata = hemoglobin, col.name = "hemoglobin")


ribo.genes <- grep(pattern = '^Rp[sl]', x = rownames(x = M_Fast_B3out@assays$RNA@data), value = TRUE)
riboPercent <- Matrix::colSums(M_Fast_B3out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(M_Fast_B3out@assays$RNA@counts)*100
M_Fast_B3out$riboPercent <- riboPercent

VlnPlot(M_Fast_B3out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fast_B3_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fast_B3out_filtered <- subset(M_Fast_B3out, 
                                subset = nCount_RNA < (mean(M_Fast_B3out@meta.data$nCount_RNA) + 4*sd(M_Fast_B3out@meta.data$nCount_RNA))
                                & nFeature_RNA < (mean(M_Fast_B3out@meta.data$nFeature_RNA) + 4*sd(M_Fast_B3out@meta.data$nFeature_RNA))
                                & mitoPercent < (mean(M_Fast_B3out@meta.data$mitoPercent) + 4*sd(M_Fast_B3out@meta.data$mitoPercent))
                                & riboPercent < (mean(M_Fast_B3out@meta.data$riboPercent) + 4*sd(M_Fast_B3out@meta.data$riboPercent))
                                & hemoglobin < (mean(M_Fast_B3out@meta.data$hemoglobin) + 4*sd(M_Fast_B3out@meta.data$hemoglobin)))


VlnPlot(M_Fast_B3out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fast_B3_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
DefaultAssay(M_Fast_B3out_filtered) <- 'RNA'
M_Fast_B3out_filtered <- DietSeurat(M_Fast_B3out_filtered, assays = 'RNA')
M_Fast_B3out_filtered <- as.SingleCellExperiment(M_Fast_B3out_filtered)

M_Fast_B3_DBR <- M_Fast_B3out_filtered@colData@nrows/50000
M_Fast_B3_scDbl <- scDblFinder(M_Fast_B3out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = M_Fast_B3_DBR, propRandom = 0.2, iter = 5)
M_Fast_B3_scDbl <- as.Seurat(M_Fast_B3_scDbl)

M_Fast_B3_scDbl <- SCTransform(M_Fast_B3_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fast_B3_scDbl <- FindVariableFeatures(M_Fast_B3_scDbl)
M_Fast_B3_scDbl@assays$SCT@var.features <- M_Fast_B3_scDbl@assays$SCT@var.features[(!M_Fast_B3_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fast_B3_scDbl <- RunPCA(M_Fast_B3_scDbl, features = VariableFeatures(object = M_Fast_B3_scDbl))
M_Fast_B3_scDbl <- RunUMAP(M_Fast_B3_scDbl, reduction = "pca", dims = 1:30)



#M_Fast_B3_scDbl <- FindNeighbors(M_Fast_B3_scDbl, dims = 1:30)
#M_Fast_B3_scDbl <- FindClusters(M_Fast_B3_scDbl, resolution = 0.8)



DimPlot(M_Fast_B3_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/M_Fast_B3_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)


M_Fast_B3_scDbl <- subset(M_Fast_B3_scDbl, scDblFinder.class == 'singlet')


M_Fast_B3_scDbl <- SCTransform(M_Fast_B3_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
#M_Fast_B3_scDbl <- FindVariableFeatures(M_Fast_B3_scDbl)
M_Fast_B3_scDbl@assays$SCT@var.features <- M_Fast_B3_scDbl@assays$SCT@var.features[(!M_Fast_B3_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]



M_Fast_B3_scDbl <- RunPCA(M_Fast_B3_scDbl, features = VariableFeatures(object = M_Fast_B3_scDbl))
M_Fast_B3_scDbl <- RunUMAP(M_Fast_B3_scDbl, reduction = "pca", dims = 1:30)



M_Fast_B3_scDbl <- FindNeighbors(M_Fast_B3_scDbl, dims = 1:30)
M_Fast_B3_scDbl <- FindClusters(M_Fast_B3_scDbl, resolution = 0.8)
DimPlot(M_Fast_B3_scDbl, label = TRUE) + NoLegend()
#check purity of clusters
FeaturePlot(M_Fast_B3_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B3_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(M_Fast_B3_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B3_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fast_B3_markers_post_scDblFnd <- FindAllMarkers(M_Fast_B3_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#FindMarkers(M_Fast_B3_scDbl, ident.1 = c('0','2','20','27'), features = c('Slc1a3'), only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save M_Fast_B3_scDbl, markers
rm(M_Fast_B3_flt, M_Fast_B3_meta, M_Fast_B3_raw, M_Fast_B3_soup, M_Fast_B3_UMAP, 
   M_Fast_B3out, M_Fast_B3out_filtered, M_Fast_B3_singlets)


#FindMarkers(F_Fast_B1_scDbl, ident.1 = c('0','2','24'), features = c('Slc1a3'))
#FindMarkers(F_Fast_B2_scDbl, ident.1 = c('0','5','19'), features = c('Slc1a3'))

#soupx_scdblfnd_diff <- readODS::read_ods('arh_sex_by_nutr_soupx_scdblfnd.ods')

#lm(soupx_scdblfnd_diff$Agrp_out ~ soupx_scdblfnd_diff$condition1 + soupx_scdblfnd_diff$condition2) |> summary()
#lm(data = soupx_scdblfnd_diff, formula = Agrp_out ~ condition1 + condition2 + condition1*condition2) |> summary()
#lm(data = soupx_scdblfnd_diff, formula = Agrp_out ~ SoupX + scDblFnd + SoupX*scDblFnd) |> summary()

#soupx_scdblfnd_diff |> filter(condition1 == 'SoupX') |> lm(formula = Agrp_out ~ condition2) |> summary()
#soupx_scdblfnd_diff |> filter(condition1 == 'SoupX') |> lm(formula = Agrp_out ~ scDblFnd) |> summary()

#soupx_scdblfnd_diff |> filter(condition1 == 'SoupX') |> lm(formula = Mbp_out ~ scDblFnd) |> summary()
#soupx_scdblfnd_diff |> filter(condition1 == 'SoupX') |> lm(formula = Slc1a3_proxy_out ~ scDblFnd) |> summary()

#soupx_scdblfnd_diff |> filter(condition2 == 'All') |> lm(formula = Slc1a3_proxy_out ~ condition1) |> summary()
#soupx_scdblfnd_diff |> filter(condition1 == 'SoupX') |> lm(formula = Agrp_out ~ condition2) |> summary()

#t.test(Mbp_out ~ condition2, paired = TRUE, data = soupx_scdblfnd_diff[c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30,32,33,35,36),])
#t.test(Pomc_out ~ condition1, paired = TRUE, data = soupx_scdblfnd_diff[c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32,34,35),])
#t.test(Pomc_out ~ condition1, paired = TRUE, data = soupx_scdblfnd_diff[c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,25,26,28,29,34,35),])
#t.test(Pomc_out ~ condition2, paired = TRUE, data = soupx_scdblfnd_diff[c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30,32,33,35,36),])

#t.test(Agrp_out ~ condition1, paired = TRUE, data = soupx_scdblfnd_diff[c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32,34,35),])
#t.test(Agrp_out ~ condition2, paired = TRUE, data = soupx_scdblfnd_diff[c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30,32,33,35,36),])


#t.test(Tac2_out ~ condition1, paired = TRUE, data = soupx_scdblfnd_diff[c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32,34,35),])
#t.test(Tac2_out ~ condition2, paired = TRUE, data = soupx_scdblfnd_diff[c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30,32,33,35,36),])


#t.test(Ptprc_out ~ condition1, paired = TRUE, data = soupx_scdblfnd_diff[c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32,34,35),])
#t.test(Ptprc_out ~ condition2, paired = TRUE, data = soupx_scdblfnd_diff[c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30,32,33,35,36),])

#t.test(Mbp_out ~ condition1, paired = TRUE, data = soupx_scdblfnd_diff[c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32,34,35),])
#t.test(Mbp_out ~ condition2, paired = TRUE, data = soupx_scdblfnd_diff[c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30,32,33,35,36),])

#t.test(Slc1a3_proxy_out ~ condition1, paired = TRUE, data = soupx_scdblfnd_diff[c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32,34,35),])
#t.test(Slc1a3_proxy_out ~ condition2, paired = TRUE, data = soupx_scdblfnd_diff[c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30,32,33,35,36),])


#soupx_scdblfnd_diff |> 
#  filter(condition2 == 'All') |> 
#  ggplot(aes(condition1, Pomc_out)) + 
#  stat_summary(geom = 'bar', fun.data = 'mean_se', width = 0.5, aes(fill = condition1)) +
#  stat_summary(geom = 'errorbar', fun.data = 'mean_se', width = 0.25) +
#  geom_point() +
#  theme_classic() +
#  scale_fill_manual(values = c('grey','red'), name = '') +
#  coord_cartesian(ylim = c(0,30)) +
#  theme_classic() +
#  labs(y = "POMC outside of POMC cluster %",
#       x = '',
#       title = bquote('Ambient RNA Removal')) +
#  theme(plot.title = element_text(hjust = 0.5,family = "Arial", size = 12, color = 'black'),
#        text = element_text(family = "Arial", size = 8, color = 'black'),
#        axis.text.x= element_text(family = 'Arial', color = 'black', size = 8),
#        axis.text.y= element_text(family = 'Arial', color = 'black', size = 8)) +
#  theme(legend.position='none')
#ggsave('soupx_figures/soupX_pomc.tiff', device = 'tiff', units = 'in', width = 3, height = 3, dpi = 300)  



#soupx_scdblfnd_diff |> 
#  filter(condition1 == 'SoupX') |> 
#  ggplot(aes(condition2, Pomc_out)) + 
#  stat_summary(geom = 'bar', fun.data = 'mean_se', width = 0.5, aes(fill = condition2)) +
#  stat_summary(geom = 'errorbar', fun.data = 'mean_se', width = 0.25) +
#  geom_point() +
#  theme_classic() +
#  scale_fill_manual(values = c('grey','red'), name = '') +
#  coord_cartesian(ylim = c(0,1)) +
#  theme_classic() +
#  labs(y = "POMC outside of POMC cluster %",
#       x = '',
#       title = bquote('Doublet Removal')) +
#  theme(plot.title = element_text(hjust = 0.5,family = "Arial", size = 12, color = 'black'),
#        text = element_text(family = "Arial", size = 8, color = 'black'),
#        axis.text.x= element_text(family = 'Arial', color = 'black', size = 8),
#        axis.text.y= element_text(family = 'Arial', color = 'black', size = 8)) +
#  theme(legend.position='none')
#ggsave('soupx_figures/scDblFnd_pomc.tiff', device = 'tiff', units = 'in', width = 3, height = 3, dpi = 300)  


#soupx_scdblfnd_diff |> 
#  filter(condition2 == 'All') |> 
#  ggplot(aes(condition1, Agrp_out)) + 
#  stat_summary(geom = 'bar', fun.data = 'mean_se', width = 0.5, aes(fill = condition1)) +
#  stat_summary(geom = 'errorbar', fun.data = 'mean_se', width = 0.25) +
#  geom_point() +
#  theme_classic() +
#  scale_fill_manual(values = c('grey','red'), name = '') +
#  coord_cartesian(ylim = c(0,45)) +
#  theme_classic() +
#  labs(y = "AgRP outside of AgRP cluster %",
#       x = '',
#       title = bquote('Ambient RNA Removal')) +
#  theme(plot.title = element_text(hjust = 0.5,family = "Arial", size = 12, color = 'black'),
#        text = element_text(family = "Arial", size = 8, color = 'black'),
#        axis.text.x= element_text(family = 'Arial', color = 'black', size = 8),
#        axis.text.y= element_text(family = 'Arial', color = 'black', size = 8)) +
#  theme(legend.position='none')
#ggsave('soupx_figures/soupX_agrp.tiff', device = 'tiff', units = 'in', width = 3, height = 3, dpi = 300)  



#soupx_scdblfnd_diff |> 
#  filter(condition1 == 'SoupX') |> 
#  ggplot(aes(condition2, Agrp_out)) + 
#  stat_summary(geom = 'bar', fun.data = 'mean_se', width = 0.5, aes(fill = condition2)) +
#  stat_summary(geom = 'errorbar', fun.data = 'mean_se', width = 0.25) +
#  geom_line(aes(condition2, Agrp_out, group = SAMPLE)) +
#  geom_point() +
#  theme_classic() +
#  scale_fill_manual(values = c('grey','red'), name = '') +
#  coord_cartesian(ylim = c(0,10)) +
#  theme_classic() +
#  labs(y = "AgRP outside of AgRP cluster %",
#       x = '',
#       title = bquote('Doublet Removal')) +
#  theme(plot.title = element_text(hjust = 0.5,family = "Arial", size = 12, color = 'black'),
#        text = element_text(family = "Arial", size = 8, color = 'black'),
#        axis.text.x= element_text(family = 'Arial', color = 'black', size = 8),
#        axis.text.y= element_text(family = 'Arial', color = 'black', size = 8)) +
#  theme(legend.position='none')
#ggsave('soupx_figures/scDblFnd_agrp.tiff', device = 'tiff', units = 'in', width = 3, height = 3, dpi = 300)  



#soupx_scdblfnd_diff |> 
#  filter(condition2 == 'All') |> 
#  ggplot(aes(condition1, Tac2_out)) + 
#  stat_summary(geom = 'bar', fun.data = 'mean_se', width = 0.5, aes(fill = condition1)) +
#  stat_summary(geom = 'errorbar', fun.data = 'mean_se', width = 0.25) +
#  geom_point() +
#  theme_classic() +
#  scale_fill_manual(values = c('grey','red'), name = '') +
#  coord_cartesian(ylim = c(0,7)) +
#  theme_classic() +
#  labs(y = "Tac2 outside of Tac2 cluster %",
#       x = '',
#       title = bquote('Ambient RNA Removal')) +
#  theme(plot.title = element_text(hjust = 0.5,family = "Arial", size = 12, color = 'black'),
#        text = element_text(family = "Arial", size = 8, color = 'black'),
#        axis.text.x= element_text(family = 'Arial', color = 'black', size = 8),
#        axis.text.y= element_text(family = 'Arial', color = 'black', size = 8)) +
#  theme(legend.position='none')
#ggsave('soupx_figures/soupX_tac2.tiff', device = 'tiff', units = 'in', width = 3, height = 3, dpi = 300)  



#soupx_scdblfnd_diff |> 
#  filter(condition1 == 'SoupX') |> 
#  ggplot(aes(condition2, Tac2_out)) + 
#  stat_summary(geom = 'bar', fun.data = 'mean_se', width = 0.5, aes(fill = condition2)) +
#  stat_summary(geom = 'errorbar', fun.data = 'mean_se', width = 0.25) +
#  geom_point() +
#  theme_classic() +
#  scale_fill_manual(values = c('grey','red'), name = '') +
#  coord_cartesian(ylim = c(0,1)) +
#  theme_classic() +
#  labs(y = "Tac2 outside of Tac2 cluster %",
#       x = '',
#       title = bquote('Doublet Removal')) +
#  theme(plot.title = element_text(hjust = 0.5,family = "Arial", size = 12, color = 'black'),
#        text = element_text(family = "Arial", size = 8, color = 'black'),
#        axis.text.x= element_text(family = 'Arial', color = 'black', size = 8),
#        axis.text.y= element_text(family = 'Arial', color = 'black', size = 8)) +
#  theme(legend.position='none')
#ggsave('soupx_figures/scDblFnd_tac2.tiff', device = 'tiff', units = 'in', width = 3, height = 3, dpi = 300)  






#soupx_scdblfnd_diff |> 
#  filter(condition2 == 'All') |> 
#  ggplot(aes(condition1, Ptprc_out)) + 
#  stat_summary(geom = 'bar', fun.data = 'mean_se', width = 0.5, aes(fill = condition1)) +
#  stat_summary(geom = 'errorbar', fun.data = 'mean_se', width = 0.25) +
#  geom_jitter(width = 0.05, height = 0) +
#  theme_classic() +
#  scale_fill_manual(values = c('grey','red'), name = '') +
#  coord_cartesian(ylim = c(0,0.5)) +
#  theme_classic() +
#  labs(y = "Ptprc outside of Ptprc cluster %",
#       x = '',
#       title = bquote('Ambient RNA Removal')) +
#  theme(plot.title = element_text(hjust = 0.5,family = "Arial", size = 12, color = 'black'),
#        text = element_text(family = "Arial", size = 8, color = 'black'),
#        axis.text.x= element_text(family = 'Arial', color = 'black', size = 8),
#        axis.text.y= element_text(family = 'Arial', color = 'black', size = 8)) +
#  theme(legend.position='none')
#ggsave('soupx_figures/soupX_Ptprc.tiff', device = 'tiff', units = 'in', width = 3, height = 3, dpi = 300)  



#soupx_scdblfnd_diff |> 
#  filter(condition1 == 'SoupX') |> 
#  ggplot(aes(condition2, Ptprc_out)) + 
#  stat_summary(geom = 'bar', fun.data = 'mean_se', width = 0.5, aes(fill = condition2)) +
#  stat_summary(geom = 'errorbar', fun.data = 'mean_se', width = 0.25) +
#  geom_jitter(width = 0.05, height = 0) +
#  geom_line(aes(condition2, Ptprc_out, group = SAMPLE)) +
#  theme_classic() +
#  scale_fill_manual(values = c('grey','red'), name = '') +
#  coord_cartesian(ylim = c(0,0.5)) +
#  theme_classic() +
#  labs(y = "Ptprc outside of Ptprc cluster %",
#       x = '',
#       title = bquote('Doublet Removal')) +
#  theme(plot.title = element_text(hjust = 0.5,family = "Arial", size = 12, color = 'black'),
#        text = element_text(family = "Arial", size = 8, color = 'black'),
#        axis.text.x= element_text(family = 'Arial', color = 'black', size = 8),
#        axis.text.y= element_text(family = 'Arial', color = 'black', size = 8)) +
#  theme(legend.position='none')
#ggsave('soupx_figures/scDblFnd_ptprc.tiff', device = 'tiff', units = 'in', width = 3, height = 3, dpi = 300)  




#soupx_scdblfnd_diff |> 
#  filter(condition2 == 'All') |> 
#  ggplot(aes(condition1, Mbp_out)) + 
#  stat_summary(geom = 'bar', fun.data = 'mean_se', width = 0.5, aes(fill = condition1)) +
#  stat_summary(geom = 'errorbar', fun.data = 'mean_se', width = 0.25) +
#  geom_point() +
#  theme_classic() +
#  scale_fill_manual(values = c('grey','red'), name = '') +
#  coord_cartesian(ylim = c(0,40)) +
#  theme_classic() +
#  labs(y = "Mbp outside of Mbp cluster %",
#       x = '',
#       title = bquote('Ambient RNA Removal')) +
#  theme(plot.title = element_text(hjust = 0.5,family = "Arial", size = 12, color = 'black'),
#        text = element_text(family = "Arial", size = 8, color = 'black'),
#        axis.text.x= element_text(family = 'Arial', color = 'black', size = 8),
#        axis.text.y= element_text(family = 'Arial', color = 'black', size = 8)) +
#  theme(legend.position='none')
#ggsave('soupx_figures/soupX_mbp.tiff', device = 'tiff', units = 'in', width = 3, height = 3, dpi = 300)  



#soupx_scdblfnd_diff |> 
#  filter(condition1 == 'SoupX') |> 
#  ggplot(aes(condition2, Mbp_out)) + 
#  stat_summary(geom = 'bar', fun.data = 'mean_se', width = 0.5, aes(fill = condition2)) +
#  stat_summary(geom = 'errorbar', fun.data = 'mean_se', width = 0.25) +
#  geom_point() +
#  geom_line(aes(condition2, Mbp_out, group = SAMPLE)) +
#  theme_classic() +
#  scale_fill_manual(values = c('grey','red'), name = '') +
#  coord_cartesian(ylim = c(0,25)) +
#  theme_classic() +
#  labs(y = "Mbp outside of Mbp cluster %",
#       x = '',
#       title = bquote('Doublet Removal')) +
#  theme(plot.title = element_text(hjust = 0.5,family = "Arial", size = 12, color = 'black'),
#        text = element_text(family = "Arial", size = 8, color = 'black'),
#        axis.text.x= element_text(family = 'Arial', color = 'black', size = 8),
#        axis.text.y= element_text(family = 'Arial', color = 'black', size = 8)) +
#  theme(legend.position='none')
#ggsave('soupx_figures/scDblFnd_mbp.tiff', device = 'tiff', units = 'in', width = 3, height = 3, dpi = 300)  



#soupx_scdblfnd_diff |> 
#  filter(condition2 == 'All') |> 
#  ggplot(aes(condition1, Slc1a3_proxy_out)) + 
#  stat_summary(geom = 'bar', fun.data = 'mean_se', width = 0.5, aes(fill = condition1)) +
#  stat_summary(geom = 'errorbar', fun.data = 'mean_se', width = 0.25) +
#  geom_point() +
#  theme_classic() +
#  scale_fill_manual(values = c('grey','red'), name = '') +
#  coord_cartesian(ylim = c(0,30)) +
#  theme_classic() +
#  labs(y = "Slc1a3 outside of Slc1a3 cluster %",
#       x = '',
#       title = bquote('Ambient RNA Removal')) +
#  theme(plot.title = element_text(hjust = 0.5,family = "Arial", size = 12, color = 'black'),
#        text = element_text(family = "Arial", size = 8, color = 'black'),
#        axis.text.x= element_text(family = 'Arial', color = 'black', size = 8),
#        axis.text.y= element_text(family = 'Arial', color = 'black', size = 8)) +
#  theme(legend.position='none')
#ggsave('soupx_figures/soupX_slc1a3_proxy.tiff', device = 'tiff', units = 'in', width = 3, height = 3, dpi = 300)  



#soupx_scdblfnd_diff |> 
#  filter(condition1 == 'SoupX') |> 
#  ggplot(aes(condition2, Slc1a3_proxy_out)) + 
#  stat_summary(geom = 'bar', fun.data = 'mean_se', width = 0.5, aes(fill = condition2)) +
#  stat_summary(geom = 'errorbar', fun.data = 'mean_se', width = 0.25) +
#  geom_point() +
#  geom_line(aes(condition2, Slc1a3_proxy_out, group = SAMPLE)) +
#  theme_classic() +
#  scale_fill_manual(values = c('grey','red'), name = '') +
#  coord_cartesian(ylim = c(0,30)) +
#  theme_classic() +
#  labs(y = "Slc1a3 outside of Slc1a3 cluster %",
#       x = '',
#       title = bquote('Doublet Removal')) +
#  theme(plot.title = element_text(hjust = 0.5,family = "Arial", size = 12, color = 'black'),
#        text = element_text(family = "Arial", size = 8, color = 'black'),
#        axis.text.x= element_text(family = 'Arial', color = 'black', size = 8),
#        axis.text.y= element_text(family = 'Arial', color = 'black', size = 8)) +
#  theme(legend.position='none')
#ggsave('soupx_figures/scDblFnd_slc1a3_proxy.tiff', device = 'tiff', units = 'in', width = 3, height = 3, dpi = 300)  


#soupx_scdblfnd_diff |> 
#  filter(condition1 == 'SoupX', condition2 == 'scDblFnd') |> 
#  ggplot(aes(condition2, Slc1a3_out)) + 
#  stat_summary(geom = 'bar', fun.data = 'mean_se', width = 0.5, aes(fill = condition2)) +
#  stat_summary(geom = 'errorbar', fun.data = 'mean_se', width = 0.25) +
#  geom_point() +
#  theme_classic() +
#  scale_fill_manual(values = c('red'), name = '') +
#  coord_cartesian(ylim = c(0,7)) +
#  theme_classic() +
#  labs(y = "Slc1a3 outside of Slc1a3 cluster %",
#       x = '',
#       title = bquote('Doublet Removal')) +
#  theme(plot.title = element_text(hjust = 0.5,family = "Arial", size = 12, color = 'black'),
#        text = element_text(family = "Arial", size = 8, color = 'black'),
#        axis.text.x= element_text(family = 'Arial', color = 'black', size = 8),
#        axis.text.y= element_text(family = 'Arial', color = 'black', size = 8)) +
#  theme(legend.position='none')
#ggsave('soupx_figures/scDblFnd_slc1a3.tiff', device = 'tiff', units = 'in', width = 1.5, height = 3, dpi = 300) 




#add metadata
F_Fed_B1_scDbl@meta.data$Sex = 'F'
F_Fed_B1_scDbl@meta.data$Sample_ID = 'F_Fed_B1'
F_Fed_B1_scDbl@meta.data$Nutr_State = 'Fed'
F_Fed_B1_scDbl@meta.data$sexXnutr = 'F_Fed'
F_Fed_B1_scDbl@meta.data$Batch = 'B1'

#double check columns
F_Fed_B1_scDbl@meta.data <- F_Fed_B1_scDbl@meta.data[,c(21,20,22:24,2:5,8:10)]


F_Fed_B2_scDbl@meta.data$Sex = 'F'
F_Fed_B2_scDbl@meta.data$Sample_ID = 'F_Fed_B2'
F_Fed_B2_scDbl@meta.data$Nutr_State = 'Fed'
F_Fed_B2_scDbl@meta.data$sexXnutr = 'F_Fed'
F_Fed_B2_scDbl@meta.data$Batch = 'B2'

#double check columns
F_Fed_B2_scDbl@meta.data <- F_Fed_B2_scDbl@meta.data[,c(21,20,22:24,2:5,8:10)]


F_Fed_B3_scDbl@meta.data$Sample_ID = 'F_Fed_B3'
F_Fed_B3_scDbl@meta.data$Sex = 'F'
F_Fed_B3_scDbl@meta.data$Nutr_State = 'Fed'
F_Fed_B3_scDbl@meta.data$sexXnutr = 'F_Fed'
F_Fed_B3_scDbl@meta.data$Batch = 'B3'

#double check columns
F_Fed_B3_scDbl@meta.data <- F_Fed_B3_scDbl@meta.data[,c(20:24,2:5,8:10)]



F_Fast_B1_scDbl@meta.data$Sample_ID = 'F_Fast_B1'
F_Fast_B1_scDbl@meta.data$Sex = 'F'
F_Fast_B1_scDbl@meta.data$Nutr_State = 'Fast'
F_Fast_B1_scDbl@meta.data$sexXnutr = 'F_Fast'
F_Fast_B1_scDbl@meta.data$Batch = 'B1'

#double check columns
F_Fast_B1_scDbl@meta.data <- F_Fast_B1_scDbl@meta.data[,c(20:24,2:5,8:10)]


F_Fast_B2_scDbl@meta.data$Sample_ID = 'F_Fast_B2'
F_Fast_B2_scDbl@meta.data$Sex = 'F'
F_Fast_B2_scDbl@meta.data$Nutr_State = 'Fast'
F_Fast_B2_scDbl@meta.data$sexXnutr = 'F_Fast'
F_Fast_B2_scDbl@meta.data$Batch = 'B2'

#double check columns
F_Fast_B2_scDbl@meta.data <- F_Fast_B2_scDbl@meta.data[,c(20:24,2:5,8:10)]


F_Fast_B3_scDbl@meta.data$Sample_ID = 'F_Fast_B3'
F_Fast_B3_scDbl@meta.data$Sex = 'F'
F_Fast_B3_scDbl@meta.data$Nutr_State = 'Fast'
F_Fast_B3_scDbl@meta.data$sexXnutr = 'F_Fast'
F_Fast_B3_scDbl@meta.data$Batch = 'B3'

#double check columns
F_Fast_B3_scDbl@meta.data <- F_Fast_B3_scDbl@meta.data[,c(21:25,2:5,8:10)]



M_Fed_B1_scDbl@meta.data$Sample_ID = 'M_Fed_B1'
M_Fed_B1_scDbl@meta.data$Sex = 'M'
M_Fed_B1_scDbl@meta.data$Nutr_State = 'Fed'
M_Fed_B1_scDbl@meta.data$sexXnutr = 'M_Fed'
M_Fed_B1_scDbl@meta.data$Batch = 'B1'

#double check columns
M_Fed_B1_scDbl@meta.data <- M_Fed_B1_scDbl@meta.data[,c(20:24,2:5,8:10)]


M_Fed_B2_scDbl@meta.data$Sample_ID = 'M_Fed_B2'
M_Fed_B2_scDbl@meta.data$Sex = 'M'
M_Fed_B2_scDbl@meta.data$Nutr_State = 'Fed'
M_Fed_B2_scDbl@meta.data$sexXnutr = 'M_Fed'
M_Fed_B2_scDbl@meta.data$Batch = 'B2'

#double check columns
M_Fed_B2_scDbl@meta.data <- M_Fed_B2_scDbl@meta.data[,c(22:26,2:5,8:10)]


M_Fed_B3_scDbl@meta.data$Sample_ID = 'M_Fed_B3'
M_Fed_B3_scDbl@meta.data$Sex = 'M'
M_Fed_B3_scDbl@meta.data$Nutr_State = 'Fed'
M_Fed_B3_scDbl@meta.data$sexXnutr = 'M_Fed'
M_Fed_B3_scDbl@meta.data$Batch = 'B3'

#double check columns
M_Fed_B3_scDbl@meta.data <- M_Fed_B3_scDbl@meta.data[,c(21:25,2:5,8:10)]



M_Fast_B1_scDbl@meta.data$Sample_ID = 'M_Fast_B1'
M_Fast_B1_scDbl@meta.data$Sex = 'M'
M_Fast_B1_scDbl@meta.data$Nutr_State = 'Fast'
M_Fast_B1_scDbl@meta.data$sexXnutr = 'M_Fast'
M_Fast_B1_scDbl@meta.data$Batch = 'B1'

#double check columns
M_Fast_B1_scDbl@meta.data <- M_Fast_B1_scDbl@meta.data[,c(21:25,2:5,8:10)]


M_Fast_B2_scDbl@meta.data$Sample_ID = 'M_Fast_B2'
M_Fast_B2_scDbl@meta.data$Sex = 'M'
M_Fast_B2_scDbl@meta.data$Nutr_State = 'Fast'
M_Fast_B2_scDbl@meta.data$sexXnutr = 'M_Fast'
M_Fast_B2_scDbl@meta.data$Batch = 'B2'

#double check columns
M_Fast_B2_scDbl@meta.data <- M_Fast_B2_scDbl@meta.data[,c(20:24,2:5,8:10)]


M_Fast_B3_scDbl@meta.data$Sample_ID = 'M_Fast_B3'
M_Fast_B3_scDbl@meta.data$Sex = 'M'
M_Fast_B3_scDbl@meta.data$Nutr_State = 'Fast'
M_Fast_B3_scDbl@meta.data$sexXnutr = 'M_Fast'
M_Fast_B3_scDbl@meta.data$Batch = 'B3'

#double check columns
M_Fast_B3_scDbl@meta.data <- M_Fast_B3_scDbl@meta.data[,c(20:24,2:5,8:10)]


#M_Fed_B2_scDbl@meta.data <- M_Fed_B2_scDbl@meta.data[,c(1:11)]
#M_Fed_B2_scDbl@meta.data <- M_Fed_B2_scDbl@meta.data |> mutate(mitoPercent = 0, .before = 10)


DefaultAssay(F_Fed_B1_scDbl) <- 'RNA'
Idents(F_Fed_B1_scDbl) <- 'Sample_ID'
F_Fed_B1_scDbl <- DietSeurat(F_Fed_B1_scDbl, assays = 'RNA')

DefaultAssay(F_Fed_B2_scDbl) <- 'RNA'
Idents(F_Fed_B2_scDbl) <- 'Sample_ID'
F_Fed_B2_scDbl <- DietSeurat(F_Fed_B2_scDbl, assays = 'RNA')

DefaultAssay(F_Fed_B3_scDbl) <- 'RNA'
Idents(F_Fed_B3_scDbl) <- 'Sample_ID'
F_Fed_B3_scDbl <- DietSeurat(F_Fed_B3_scDbl, assays = 'RNA')


DefaultAssay(F_Fast_B1_scDbl) <- 'RNA'
Idents(F_Fast_B1_scDbl) <- 'Sample_ID'
F_Fast_B1_scDbl <- DietSeurat(F_Fast_B1_scDbl, assays = 'RNA')

DefaultAssay(F_Fast_B2_scDbl) <- 'RNA'
Idents(F_Fast_B2_scDbl) <- 'Sample_ID'
F_Fast_B2_scDbl <- DietSeurat(F_Fast_B2_scDbl, assays = 'RNA')

DefaultAssay(F_Fast_B3_scDbl) <- 'RNA'
Idents(F_Fast_B3_scDbl) <- 'Sample_ID'
F_Fast_B3_scDbl <- DietSeurat(F_Fast_B3_scDbl, assays = 'RNA')


DefaultAssay(M_Fed_B1_scDbl) <- 'RNA'
Idents(M_Fed_B1_scDbl) <- 'Sample_ID'
M_Fed_B1_scDbl <- DietSeurat(M_Fed_B1_scDbl, assays = 'RNA')

DefaultAssay(M_Fed_B2_scDbl) <- 'RNA'
Idents(M_Fed_B2_scDbl) <- 'Sample_ID'
M_Fed_B2_scDbl <- DietSeurat(M_Fed_B2_scDbl, assays = 'RNA')

DefaultAssay(M_Fed_B3_scDbl) <- 'RNA'
Idents(M_Fed_B3_scDbl) <- 'Sample_ID'
M_Fed_B3_scDbl <- DietSeurat(M_Fed_B3_scDbl, assays = 'RNA')


DefaultAssay(M_Fast_B1_scDbl) <- 'RNA'
Idents(M_Fast_B1_scDbl) <- 'Sample_ID'
M_Fast_B1_scDbl <- DietSeurat(M_Fast_B1_scDbl, assays = 'RNA')

DefaultAssay(M_Fast_B2_scDbl) <- 'RNA'
Idents(M_Fast_B2_scDbl) <- 'Sample_ID'
M_Fast_B2_scDbl <- DietSeurat(M_Fast_B2_scDbl, assays = 'RNA')

DefaultAssay(M_Fast_B3_scDbl) <- 'RNA'
Idents(M_Fast_B3_scDbl) <- 'Sample_ID'
M_Fast_B3_scDbl <- DietSeurat(M_Fast_B3_scDbl, assays = 'RNA')


ARH_Sex_by_Nutr <- merge(F_Fed_B1_scDbl, y= c(F_Fed_B2_scDbl, F_Fed_B3_scDbl, 
                                              F_Fast_B1_scDbl, F_Fast_B2_scDbl, F_Fast_B3_scDbl, 
                                              M_Fed_B1_scDbl, M_Fed_B2_scDbl, M_Fed_B3_scDbl, 
                                              M_Fast_B1_scDbl, M_Fast_B2_scDbl, M_Fast_B3_scDbl), 
                         add.cell.ids = c('ffdb1','ffdb2','ffdb3',
                                          'fftb1','fftb2','fftb3',
                                          'mfdb1','mfdb2','mfdb3',
                                          'mftb1','mftb2','mftb3'), 
                         project = 'ARH_Sex_by_Nutr')


saveRDS(ARH_Sex_by_Nutr, file = 'data/ARH_Sex_by_Nutr.rds')

#remove all objects, libs
rm(cells2HL, cluster_celltype_scores_exon_only, dr10_eo_arh_cell_markers, dr10_eo_arh_cell_markers250pc, dr10_eo_arh_neurons_markers, F_Fast_B1_DBR, F_Fast_B1_markers,
   F_Fast_B1_scDbl, F_Fast_B2_DBR, F_Fast_B1_markers, F_Fast_B2_scDbl, F_Fast_B3_DBR, F_Fast_B3_markers, F_Fast_B3_scDbl, F_Fed_B1_DBR, F_Fed_B1_markers, 
   F_Fed_B1_scDbl, F_Fed_B1_scDbl2, F_Fed_B2_DBR, F_Fed_B2_markers, F_Fed_B2_scDbl, F_Fed_B3_DBR, F_Fed_B3_markers, F_Fed_B3_scDbl, M_Fast_B1_DBR, M_Fast_B1_markers,
   M_Fast_B1_scDbl, hemoglobin, arh_features, libs_min, M_Fast_B2_DBR, M_Fast_B2_markers, M_Fast_B2_scDbl, M_Fast_B3_DBR, M_Fast_B3_scDbl, M_Fed_B1_DBR,
   M_Fed_B1_markers, M_Fed_B1_scDbl, M_Fed_B1_markers, M_Fed_B1_scDbl, M_Fed_B2_DBR, M_Fed_B2_markers, M_Fed_B2_scDbl, M_Fed_B3_DBR, M_Fed_B3_markers, M_Fed_B3_scDbl,
   markers, mito.genes2, mitoPercent, riboPercent, select_arh_features, sn_cell_types, sn_conditions)
# remove any remaining objects

#ARH_Sex_by_Nutr <- readRDS('data/ARH_Sex_by_Nutr.rds')

#ribo.genes <- grep(pattern = '^Rp[sl]', x = rownames(x = ARH_Sex_by_Nutr@assays$RNA@data), value = TRUE)
#riboPercent <- Matrix::colSums(ARH_Sex_by_Nutr@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(ARH_Sex_by_Nutr@assays$RNA@counts)*100
#ARH_Sex_by_Nutr$riboPercent <- riboPercent

#saveRDS(ARH_Sex_by_Nutr, file = 'data/ARH_Sex_by_Nutr.rds')

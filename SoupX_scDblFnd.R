#the purpose of this script is to filter out ambient RNA using SoupX program and exclude doublets using scDblFinder

#load required libraries
libs <- c( 'Seurat','dplyr','tidyr','ggplot2','scDblFinder','SoupX')

lapply(libs, require, character.only = TRUE)



#create basic clustering for sample
#load filtered sample
F_Fed_B1_flt <- Read10X('../STARsolo/F_Fed_B1Solo.out/GeneFull_Ex50pAS/filtered/')
#create seurat object
F_Fed_B1_flt <- CreateSeuratObject(F_Fed_B1_flt)

#normalize expression
F_Fed_B1_flt <- SCTransform(F_Fed_B1_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')

# genes to ignore on y chromosome and Xist Tsix
yx_chrom_genes <- c('Eif2s3y','Sry','Zfy','Rps4y1','Amely','Tbl1y','Pcdh11y','Tgif2ly',
                    'Tspy1','Tspy2','Azfa','Usp9y','Ddx3y','Uty','Tb4y','Azfb',
                    'Cyorf15','Rps4y2','Eif1ay','Kdm5d','Xkry','Hsfy1','Hsfy2',
                    'Pry','Pry2','Rbmy1a1','Azfc','Daz1','Daz2','Daz3','Daz4',
                    'Cdy1','Cdy2','Vcy1','Vcy2','Xist','Tsix'
)

# genes to ignore that are mitochondrial genes
mito.genes <- grep(pattern = "^mt-", x = rownames(x = F_Fed_B1_flt@assays$RNA@data), value = TRUE)
# genes to ignore that are hemoglobin genes
hemoglobin_genes <- c('Hbq1a','Hbb-y','Hbb-bt','Hba-a2','Hba-a1')
# genes to ignore that are ribosomal binding protein genes
ribo.genes <- grep(pattern = '^Rp[sl]', x = rownames(x = F_Fed_B1_flt@assays$RNA@data), value = TRUE)

# find variable features ignoring genes from list above
F_Fed_B1_flt@assays$SCT@var.features <- F_Fed_B1_flt@assays$SCT@var.features[(!F_Fed_B1_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

# find principal components using variable features
F_Fed_B1_flt <- RunPCA(F_Fed_B1_flt, features = VariableFeatures(object = F_Fed_B1_flt))
# run umap reduction using first 30 PCs
F_Fed_B1_flt <- RunUMAP(F_Fed_B1_flt, reduction = "pca", dims = 1:30)

# find neighbors using first 30 PCs, find clusters at a resolution of 0.8
F_Fed_B1_flt <- FindNeighbors(F_Fed_B1_flt, dims = 1:30)
F_Fed_B1_flt <- FindClusters(F_Fed_B1_flt, resolution = 0.8)

# plot UMAP projection
DimPlot(F_Fed_B1_flt, label = TRUE, label.size = 2) + NoLegend()

# check expression of canonical ARH neurons genes
#FeaturePlot(F_Fed_B1_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B1_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# check expression of neuron and glia marker genes
#FeaturePlot(F_Fed_B1_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B1_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes before running soupX
F_Fed_B1_markers_b4_soupx <- FindAllMarkers(F_Fed_B1_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# save embeddings for later use
F_Fed_B1_UMAP <- F_Fed_B1_flt@reductions$umap@cell.embeddings |> as.data.frame()
F_Fed_B1_meta <- F_Fed_B1_flt@meta.data[,c(2,3,7)] |> cbind(F_Fed_B1_UMAP)



#start SoupX in earnest
# read filtered and raw samples in
F_Fed_B1_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fed_B1Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fed_B1_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fed_B1Solo.out/GeneFull_Ex50pAS/raw/')
# create soupX object
F_Fed_B1_soup <-SoupChannel(F_Fed_B1_raw, F_Fed_B1_flt)

# set demension reduction based on saved UMAP projections 
F_Fed_B1_soup = setDR(F_Fed_B1_soup, F_Fed_B1_meta[colnames(F_Fed_B1_flt), c("UMAP_1", "UMAP_2")])
# set clusters based on saved clusters from before
F_Fed_B1_soup = setClusters(F_Fed_B1_soup, setNames(F_Fed_B1_meta$seurat_clusters, colnames(F_Fed_B1_flt)))

# plot UMAP projection
#ggplot(F_Fed_B1_soup$metaData, aes(UMAP_1, UMAP_2)) +
#  geom_point(aes(color = clusters))

# estimate rho
F_Fed_B1_soup = autoEstCont(F_Fed_B1_soup)
# rho estimated at 0.04; however, 0.12 producing cleaner looking clusters based on Agrp expression
# set rho and adjust counts accordingly
F_Fed_B1_soup <- setContaminationFraction(F_Fed_B1_soup, 0.12)
F_Fed_B1out <- adjustCounts(F_Fed_B1_soup, method = 'multinomial')

# SoupX finished, check results
F_Fed_B1out <- CreateSeuratObject(F_Fed_B1out)

#normalize expression, and find varibale genes ignoring select genes from lists above
F_Fed_B1out <- SCTransform(F_Fed_B1out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')

F_Fed_B1out@assays$SCT@var.features <- F_Fed_B1out@assays$SCT@var.features[(!F_Fed_B1out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# find PCs in new object using varibale features
F_Fed_B1out <- RunPCA(F_Fed_B1out, features = VariableFeatures(object = F_Fed_B1out))
# find UMAP projections based on first 30 PCs
F_Fed_B1out <- RunUMAP(F_Fed_B1out, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs and clusters at resolution of 0.8 in new object 
F_Fed_B1out <- FindNeighbors(F_Fed_B1out, dims = 1:30)
F_Fed_B1out <- FindClusters(F_Fed_B1out, resolution = 0.8)

# plot UMAP projections
DimPlot(F_Fed_B1out, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes 
FeaturePlot(F_Fed_B1out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B1_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot neuron and glia marker genes
FeaturePlot(F_Fed_B1out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B1_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find marker genes after running soupX
F_Fed_B1_markers_post_soupx <- FindAllMarkers(F_Fed_B1out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# find and save percent of mitochondria genes expression
mitoPercent <- Matrix::colSums(F_Fed_B1out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(F_Fed_B1out@assays$RNA@counts)*100
F_Fed_B1out$mitoPercent <- mitoPercent

# find and save count of hemoglobin gene expression
hemoglobin <- Matrix::colSums(F_Fed_B1out@assays$RNA@counts[hemoglobin_genes, ])
F_Fed_B1out <- AddMetaData(object = F_Fed_B1out, metadata = hemoglobin, col.name = "hemoglobin")


# find and save percent of ribosomal binding protein genes expression
riboPercent <- Matrix::colSums(F_Fed_B1out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(F_Fed_B1out@assays$RNA@counts)*100
F_Fed_B1out$riboPercent <- riboPercent

# plot quality markers saved above
VlnPlot(F_Fed_B1out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fed_B1_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#filter out extreme values, that may be low quality or potential doublet (mean + 4sd)
F_Fed_B1out_filtered <- subset(F_Fed_B1out, 
                               subset = nCount_RNA < (mean(F_Fed_B1out@meta.data$nCount_RNA) + 4*sd(F_Fed_B1out@meta.data$nCount_RNA))
                               & nFeature_RNA < (mean(F_Fed_B1out@meta.data$nFeature_RNA) + 4*sd(F_Fed_B1out@meta.data$nFeature_RNA))
                               & mitoPercent < (mean(F_Fed_B1out@meta.data$mitoPercent) + 4*sd(F_Fed_B1out@meta.data$mitoPercent))
                               & riboPercent < (mean(F_Fed_B1out@meta.data$riboPercent) + 4*sd(F_Fed_B1out@meta.data$riboPercent))
                               & hemoglobin < (mean(F_Fed_B1out@meta.data$hemoglobin) + 4*sd(F_Fed_B1out@meta.data$hemoglobin)))


# plot new quality markers values
VlnPlot(F_Fed_B1out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fed_B1_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)


#scDblFinder

# prep object for scDbl
# set default assay as RNA
DefaultAssay(F_Fed_B1out_filtered) <- 'RNA'
# save seurat object with only RNA count info
F_Fed_B1out_filtered <- DietSeurat(F_Fed_B1out_filtered, assays = 'RNA')
# transform into a SC experiment object
F_Fed_B1out_filtered <- as.SingleCellExperiment(F_Fed_B1out_filtered)

# set doublet rate at number of cells / 50000
F_Fed_B1_DBR <- F_Fed_B1out_filtered@colData@nrows/50000
# run scDblFinder using clusters from before, 2000 features, 30 dimensions, a random proportion of 0.2 and 5 iterations
F_Fed_B1_scDbl <- scDblFinder(F_Fed_B1out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = F_Fed_B1_DBR, propRandom = 0.2, iter = 5)
# transform back into seurat object
F_Fed_B1_scDbl <- as.Seurat(F_Fed_B1_scDbl)

# normalize expression calculate variable features ignoring y chromosome, mitochondria genes, RBP genes and hemoglobin genes
F_Fed_B1_scDbl <- SCTransform(F_Fed_B1_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')

F_Fed_B1_scDbl@assays$SCT@var.features <- F_Fed_B1_scDbl@assays$SCT@var.features[(!F_Fed_B1_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features
F_Fed_B1_scDbl <- RunPCA(F_Fed_B1_scDbl, features = VariableFeatures(object = F_Fed_B1_scDbl))
# calculate UMAP projections based on first 30 PCs
F_Fed_B1_scDbl <- RunUMAP(F_Fed_B1_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs and clusters at 0.8 resolution
F_Fed_B1_scDbl <- FindNeighbors(F_Fed_B1_scDbl, dims = 1:30)
F_Fed_B1_scDbl <- FindClusters(F_Fed_B1_scDbl, resolution = 0.8)


# plot singlets and doublets in UMAP space
DimPlot(F_Fed_B1_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/F_Fed_B1_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)
#F_Fed_B1_scDbl <- DietSeurat(F_Fed_B1_scDbl, dimreducs = NULL)

# save only singlets
F_Fed_B1_scDbl <- subset(F_Fed_B1_scDbl, scDblFinder.class == 'singlet')

# normalize expression of new object
F_Fed_B1_scDbl <- SCTransform(F_Fed_B1_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')

F_Fed_B1_scDbl@assays$SCT@var.features <- F_Fed_B1_scDbl@assays$SCT@var.features[(!F_Fed_B1_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs of new object
F_Fed_B1_scDbl <- RunPCA(F_Fed_B1_scDbl, features = VariableFeatures(object = F_Fed_B1_scDbl))
# calculate UMAP of new object
F_Fed_B1_scDbl <- RunUMAP(F_Fed_B1_scDbl, reduction = "pca", dims = 1:30)


# find neighbors and clusters of new object
F_Fed_B1_scDbl <- FindNeighbors(F_Fed_B1_scDbl, dims = 1:30)
F_Fed_B1_scDbl <- FindClusters(F_Fed_B1_scDbl, resolution = 0.8)

#check purity of clusters
# plot canonical ARH neuron genes expression
FeaturePlot(F_Fed_B1_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B1_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)
# plot neuron and glia marker genes expression
FeaturePlot(F_Fed_B1_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B1_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# calculate marker genes for singlets object
F_Fed_B1_markers_post_scDblFnd <- FindAllMarkers(F_Fed_B1_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

#FindMarkers(F_Fed_B1_scDbl, ident.1 = c('4','8'), features = c('Agrp','Npy'), only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save F_Fed_B1_scDbl, markers
rm(F_Fed_B1_flt, F_Fed_B1_meta, F_Fed_B1_raw, F_Fed_B1_soup, F_Fed_B1_UMAP, 
   F_Fed_B1out, F_Fed_B1out_filtered, F_Fed_B1_singlets)


# repeat steps for next sample

#F_Fed_B2
#create basic clustering for sample
# read filtered sample in and create seurat object
F_Fed_B2_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fed_B2Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fed_B2_flt <- CreateSeuratObject(F_Fed_B2_flt)

# normalize expression, find variable features ignoring xy mt rbp and hemoglobin genes
F_Fed_B2_flt <- SCTransform(F_Fed_B2_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fed_B2_flt@assays$SCT@var.features <- F_Fed_B2_flt@assays$SCT@var.features[(!F_Fed_B2_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

# calculate PCs based on varibale features, calculate UMAP projections based on first 30 PCs
F_Fed_B2_flt <- RunPCA(F_Fed_B2_flt, features = VariableFeatures(object = F_Fed_B2_flt))
F_Fed_B2_flt <- RunUMAP(F_Fed_B2_flt, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fed_B2_flt <- FindNeighbors(F_Fed_B2_flt, dims = 1:30)
F_Fed_B2_flt <- FindClusters(F_Fed_B2_flt, resolution = 0.8)

# plot UMAP projection
DimPlot(F_Fed_B2_flt, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(F_Fed_B2_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B2_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(F_Fed_B2_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B2_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes min pct of 0.33, logfc threshold log2(1.33)
F_Fed_B2_markers_b4_soupx <- FindAllMarkers(F_Fed_B2_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# save UMAP projections and meta data for future use
F_Fed_B2_UMAP <- F_Fed_B2_flt@reductions$umap@cell.embeddings |> as.data.frame()
F_Fed_B2_meta <- F_Fed_B2_flt@meta.data[,c(2,3,7)] |> cbind(F_Fed_B2_UMAP)



#start SoupX in earnest 
# read in filtered and raw samples, create soupX object
F_Fed_B2_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fed_B2Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fed_B2_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fed_B2Solo.out/GeneFull_Ex50pAS/raw/')
F_Fed_B2_soup <-SoupChannel(F_Fed_B2_raw, F_Fed_B2_flt)

# set dimension reduction to saved UMAP projections, set clusters based on saved clusters in metadata
F_Fed_B2_soup = setDR(F_Fed_B2_soup, F_Fed_B2_meta[colnames(F_Fed_B2_flt), c("UMAP_1", "UMAP_2")])
F_Fed_B2_soup = setClusters(F_Fed_B2_soup, setNames(F_Fed_B2_meta$seurat_clusters, colnames(F_Fed_B2_flt)))

# plot UMAP projections
ggplot(F_Fed_B2_soup$metaData, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = clusters))

# estimate rho
F_Fed_B2_soup = autoEstCont(F_Fed_B2_soup)
# rho estimated at 0.08; however, 0.20 producing cleaner looking clusters based on Agrp expression

# set contamination fraction and adjust counts 
F_Fed_B2_soup <- setContaminationFraction(F_Fed_B2_soup, 0.2)
F_Fed_B2out <- adjustCounts(F_Fed_B2_soup, method = 'multinomial')

# SoupX finished, check results
# create seurat object
F_Fed_B2out <- CreateSeuratObject(F_Fed_B2out)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
F_Fed_B2out <- SCTransform(F_Fed_B2out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fed_B2out@assays$SCT@var.features <- F_Fed_B2out@assays$SCT@var.features[(!F_Fed_B2out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP projections based on first 30 PCs
F_Fed_B2out <- RunPCA(F_Fed_B2out, features = VariableFeatures(object = F_Fed_B2out))
F_Fed_B2out <- RunUMAP(F_Fed_B2out, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fed_B2out <- FindNeighbors(F_Fed_B2out, dims = 1:30)
F_Fed_B2out <- FindClusters(F_Fed_B2out, resolution = 0.8)

# plot UMAP projection
DimPlot(F_Fed_B2out, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(F_Fed_B2out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B2_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(F_Fed_B2out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B2_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all markers genes min pct 0.33 logfc threshold log2(1.33)
F_Fed_B2_markers_post_soupx <- FindAllMarkers(F_Fed_B2out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# calculate and save percent of mt gene expression
mitoPercent <- Matrix::colSums(F_Fed_B2out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(F_Fed_B2out@assays$RNA@counts)*100
F_Fed_B2out$mitoPercent <- mitoPercent

# calculate and save hemoglobin gene counts
hemoglobin <- Matrix::colSums(F_Fed_B2out@assays$RNA@counts[hemoglobin_genes, ])
F_Fed_B2out <- AddMetaData(object = F_Fed_B2out, metadata = hemoglobin, col.name = "hemoglobin")


# calculate and save pct of rbp gene expression
riboPercent <- Matrix::colSums(F_Fed_B2out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(F_Fed_B2out@assays$RNA@counts)*100
F_Fed_B2out$riboPercent <- riboPercent

# plot QC markers
VlnPlot(F_Fed_B2out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fed_B2_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

# filter extreme values that may be low quality or possible doublet, mean + 4sd
F_Fed_B2out_filtered <- subset(F_Fed_B2out, 
                               subset = nCount_RNA < (mean(F_Fed_B2out@meta.data$nCount_RNA) + 4*sd(F_Fed_B2out@meta.data$nCount_RNA))
                               & nFeature_RNA < (mean(F_Fed_B2out@meta.data$nFeature_RNA) + 4*sd(F_Fed_B2out@meta.data$nFeature_RNA))
                               & mitoPercent < (mean(F_Fed_B2out@meta.data$mitoPercent) + 4*sd(F_Fed_B2out@meta.data$mitoPercent))
                               & riboPercent < (mean(F_Fed_B2out@meta.data$riboPercent) + 4*sd(F_Fed_B2out@meta.data$riboPercent))
                               & hemoglobin < (mean(F_Fed_B2out@meta.data$hemoglobin) + 4*sd(F_Fed_B2out@meta.data$hemoglobin)))

# plot new QC values
VlnPlot(F_Fed_B2out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fed_B2_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
# set default assay to RNA, save only RNA counts, transform to a sc experiment
DefaultAssay(F_Fed_B2out_filtered) <- 'RNA'
F_Fed_B2out_filtered <- DietSeurat(F_Fed_B2out_filtered, assays = 'RNA')
F_Fed_B2out_filtered <- as.SingleCellExperiment(F_Fed_B2out_filtered)

# set doublet rate to number of cells / 50000
F_Fed_B2_DBR <- F_Fed_B2out_filtered@colData@nrows/50000
# run scdbl finder, using seurat clusters, 2000 features, 30 dims, random proportion 0.2 and 5 iterations
F_Fed_B2_scDbl <- scDblFinder(F_Fed_B2out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = F_Fed_B2_DBR, propRandom = 0.2, iter = 5)
# transform back to seurat object
F_Fed_B2_scDbl <- as.Seurat(F_Fed_B2_scDbl)

# normalize expression, find varible features excluding xy mt rbp hemoglobin genes
F_Fed_B2_scDbl <- SCTransform(F_Fed_B2_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fed_B2_scDbl@assays$SCT@var.features <- F_Fed_B2_scDbl@assays$SCT@var.features[(!F_Fed_B2_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on varible features, calculate UMAP based on first 30 PCs
F_Fed_B2_scDbl <- RunPCA(F_Fed_B2_scDbl, features = VariableFeatures(object = F_Fed_B2_scDbl))
F_Fed_B2_scDbl <- RunUMAP(F_Fed_B2_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fed_B2_scDbl <- FindNeighbors(F_Fed_B2_scDbl, dims = 1:30)
F_Fed_B2_scDbl <- FindClusters(F_Fed_B2_scDbl, resolution = 0.8)


# plot singlets and doublets in UMAP space
DimPlot(F_Fed_B2_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/F_Fed_B2_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)

# save only singlets
F_Fed_B2_scDbl <- subset(F_Fed_B2_scDbl, scDblFinder.class == 'singlet')

# normalize expression of new object, find variable features excluding xy mt rbp hemoglobin genes
F_Fed_B2_scDbl <- SCTransform(F_Fed_B2_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fed_B2_scDbl@assays$SCT@var.features <- F_Fed_B2_scDbl@assays$SCT@var.features[(!F_Fed_B2_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
F_Fed_B2_scDbl <- RunPCA(F_Fed_B2_scDbl, features = VariableFeatures(object = F_Fed_B2_scDbl))
F_Fed_B2_scDbl <- RunUMAP(F_Fed_B2_scDbl, reduction = "pca", dims = 1:30)


# find neighbors using first 30 PCs, find clusters at 0.8 resolution
F_Fed_B2_scDbl <- FindNeighbors(F_Fed_B2_scDbl, dims = 1:30)
F_Fed_B2_scDbl <- FindClusters(F_Fed_B2_scDbl, resolution = 0.8)

#check purity of clusters
# plot expression of canonical ARH neuron genes
FeaturePlot(F_Fed_B2_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B2_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot neuron and glia marker genes
FeaturePlot(F_Fed_B2_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B2_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot UMAP projections
DimPlot(F_Fed_B2_scDbl, label = TRUE) + NoLegend()
# find all marker genes min pct 0.33 logfc thershold log2(1.33)
F_Fed_B2_markers_post_scDblFnd <- FindAllMarkers(F_Fed_B2_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)



# remove unneeded objects, save F_Fed_B2_scDbl, markers
rm(F_Fed_B2_flt, F_Fed_B2_meta, F_Fed_B2_raw, F_Fed_B2_soup, F_Fed_B2_UMAP, 
   F_Fed_B2out, F_Fed_B2out_filtered, F_Fed_B2_singlets)


# repeat steps for next sample

#F_Fed_B3
#create basic clustering for sample
# read filtered sample in, create seurat object
F_Fed_B3_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fed_B3Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fed_B3_flt <- CreateSeuratObject(F_Fed_B3_flt)

# normalize expression, find varibale features excluding xy mt rbp  hemoglobin genes
F_Fed_B3_flt <- SCTransform(F_Fed_B3_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fed_B3_flt@assays$SCT@var.features <- F_Fed_B3_flt@assays$SCT@var.features[(!F_Fed_B3_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
F_Fed_B3_flt <- RunPCA(F_Fed_B3_flt, features = VariableFeatures(object = F_Fed_B3_flt))
F_Fed_B3_flt <- RunUMAP(F_Fed_B3_flt, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fed_B3_flt <- FindNeighbors(F_Fed_B3_flt, dims = 1:30)
F_Fed_B3_flt <- FindClusters(F_Fed_B3_flt, resolution = 0.8)

# plot UMAP projections
DimPlot(F_Fed_B3_flt, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(F_Fed_B3_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B3_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(F_Fed_B3_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B3_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes, min pct 0.33, logfc threshold log2(1.33)
F_Fed_B3_markers_b4_soupx <- FindAllMarkers(F_Fed_B3_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# save UMAP projections and metadata for future use
F_Fed_B3_UMAP <- F_Fed_B3_flt@reductions$umap@cell.embeddings |> as.data.frame()
F_Fed_B3_meta <- F_Fed_B3_flt@meta.data[,c(2,3,7)] |> cbind(F_Fed_B3_UMAP)



#start SoupX in earnest 
# read in filtered and raw samples, create soupx object
F_Fed_B3_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fed_B3Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fed_B3_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fed_B3Solo.out/GeneFull_Ex50pAS/raw/')
F_Fed_B3_soup <-SoupChannel(F_Fed_B3_raw, F_Fed_B3_flt)

# set dimension reduction to saved UMAP, set clusters to saved clusters
F_Fed_B3_soup = setDR(F_Fed_B3_soup, F_Fed_B3_meta[colnames(F_Fed_B3_flt), c("UMAP_1", "UMAP_2")])
F_Fed_B3_soup = setClusters(F_Fed_B3_soup, setNames(F_Fed_B3_meta$seurat_clusters, colnames(F_Fed_B3_flt)))

# estimate rho
F_Fed_B3_soup = autoEstCont(F_Fed_B3_soup)
# rho estimated at 0.06; however, 0.15 producing cleaner looking clusters based on Agrp expression

# set contamination fraction and adjust counts
F_Fed_B3_soup <- setContaminationFraction(F_Fed_B3_soup, 0.15)
F_Fed_B3out <- adjustCounts(F_Fed_B3_soup, method = 'multinomial')

# SoupX finished, check results
# create seurat object
F_Fed_B3out <- CreateSeuratObject(F_Fed_B3out)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
F_Fed_B3out <- SCTransform(F_Fed_B3out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fed_B3out@assays$SCT@var.features <- F_Fed_B3out@assays$SCT@var.features[(!F_Fed_B3out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
F_Fed_B3out <- RunPCA(F_Fed_B3out, features = VariableFeatures(object = F_Fed_B3out))
F_Fed_B3out <- RunUMAP(F_Fed_B3out, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fed_B3out <- FindNeighbors(F_Fed_B3out, dims = 1:30)
F_Fed_B3out <- FindClusters(F_Fed_B3out, resolution = 0.8)

# plot UMAP projections
DimPlot(F_Fed_B3out, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(F_Fed_B3out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B3_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(F_Fed_B3out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B3_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes, min pct 0.33, logfc threshold log2(1.33)
F_Fed_B3_markers_post_soupx <- FindAllMarkers(F_Fed_B3out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# calculate percent of mt gene expression and save
mitoPercent <- Matrix::colSums(F_Fed_B3out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(F_Fed_B3out@assays$RNA@counts)*100
F_Fed_B3out$mitoPercent <- mitoPercent

# calculate and save hemoglobin gene expression counts
hemoglobin <- Matrix::colSums(F_Fed_B3out@assays$RNA@counts[hemoglobin_genes, ])
F_Fed_B3out <- AddMetaData(object = F_Fed_B3out, metadata = hemoglobin, col.name = "hemoglobin")


# calculate and save pct rbp gene expression
riboPercent <- Matrix::colSums(F_Fed_B3out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(F_Fed_B3out@assays$RNA@counts)*100
F_Fed_B3out$riboPercent <- riboPercent

# plot QC features
VlnPlot(F_Fed_B3out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fed_B3_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

# filter out extreme values that may be low quality or possible doublet, mean + 4sd
F_Fed_B3out_filtered <- subset(F_Fed_B3out, 
                               subset = nCount_RNA < (mean(F_Fed_B3out@meta.data$nCount_RNA) + 4*sd(F_Fed_B3out@meta.data$nCount_RNA))
                               & nFeature_RNA < (mean(F_Fed_B3out@meta.data$nFeature_RNA) + 4*sd(F_Fed_B3out@meta.data$nFeature_RNA))
                               & mitoPercent < (mean(F_Fed_B3out@meta.data$mitoPercent) + 4*sd(F_Fed_B3out@meta.data$mitoPercent))
                               & riboPercent < (mean(F_Fed_B3out@meta.data$riboPercent) + 4*sd(F_Fed_B3out@meta.data$riboPercent))
                               & hemoglobin < (mean(F_Fed_B3out@meta.data$hemoglobin) + 4*sd(F_Fed_B3out@meta.data$hemoglobin)))

# plot new QC values
VlnPlot(F_Fed_B3out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fed_B3_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
# set default assay to RNA, save only RNA counts, transform into a sc experiment
DefaultAssay(F_Fed_B3out_filtered) <- 'RNA'
F_Fed_B3out_filtered <- DietSeurat(F_Fed_B3out_filtered, assays = 'RNA')
F_Fed_B3out_filtered <- as.SingleCellExperiment(F_Fed_B3out_filtered)

# set doublet rate to number of cells / 50000
F_Fed_B3_DBR <- F_Fed_B3out_filtered@colData@nrows/50000
# run scdbl finder, using saved clusters, 2000 features, 30 dims, random proportion 0.2 and 5 iterations
F_Fed_B3_scDbl <- scDblFinder(F_Fed_B3out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = F_Fed_B3_DBR, propRandom = 0.2, iter = 5)
# tranfrom back to seurat object
F_Fed_B3_scDbl <- as.Seurat(F_Fed_B3_scDbl)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
F_Fed_B3_scDbl <- SCTransform(F_Fed_B3_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fed_B3_scDbl@assays$SCT@var.features <- F_Fed_B3_scDbl@assays$SCT@var.features[(!F_Fed_B3_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
F_Fed_B3_scDbl <- RunPCA(F_Fed_B3_scDbl, features = VariableFeatures(object = F_Fed_B3_scDbl))
F_Fed_B3_scDbl <- RunUMAP(F_Fed_B3_scDbl, reduction = "pca", dims = 1:30)


# find neighbors using first 30 PCs, find clusters at 0.8 resolution
F_Fed_B3_scDbl <- FindNeighbors(F_Fed_B3_scDbl, dims = 1:30)
F_Fed_B3_scDbl <- FindClusters(F_Fed_B3_scDbl, resolution = 0.8)


# plot singlets and doublets in UMAP space
DimPlot(F_Fed_B3_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/F_Fed_B3_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)

# save only singlets
F_Fed_B3_scDbl <- subset(F_Fed_B3_scDbl, scDblFinder.class == 'singlet')

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
F_Fed_B3_scDbl <- SCTransform(F_Fed_B3_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fed_B3_scDbl@assays$SCT@var.features <- F_Fed_B3_scDbl@assays$SCT@var.features[(!F_Fed_B3_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
F_Fed_B3_scDbl <- RunPCA(F_Fed_B3_scDbl, features = VariableFeatures(object = F_Fed_B3_scDbl))
F_Fed_B3_scDbl <- RunUMAP(F_Fed_B3_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fed_B3_scDbl <- FindNeighbors(F_Fed_B3_scDbl, dims = 1:30)
F_Fed_B3_scDbl <- FindClusters(F_Fed_B3_scDbl, resolution = 0.8)

#check purity of clusters
# plot expression of canonical ARH neuron genes
FeaturePlot(F_Fed_B3_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fed_B3_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(F_Fed_B3_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fed_B3_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes, min pct 0.33, logfc threshold log2(1.33)
F_Fed_B3_markers_post_scDblFnd <- FindAllMarkers(F_Fed_B3_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save F_Fed_B3_scDbl, markers
rm(F_Fed_B3_flt, F_Fed_B3_meta, F_Fed_B3_raw, F_Fed_B3_soup, F_Fed_B3_UMAP, 
   F_Fed_B3out, F_Fed_B3out_filtered, F_Fed_B3_singlets)

# repeat steps for next sample

#F_Fast_B1
#create basic clustering for sample
# read filtered sample in, create seurat object
F_Fast_B1_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B1Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fast_B1_flt <- CreateSeuratObject(F_Fast_B1_flt)

# noramalize counts, find variable features excluding xy mt rbp hemoglobin genes
F_Fast_B1_flt <- SCTransform(F_Fast_B1_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fast_B1_flt@assays$SCT@var.features <- F_Fast_B1_flt@assays$SCT@var.features[(!F_Fast_B1_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

# calculate PCs based on variable genes, calculate UMAP based on first 30 PCs
F_Fast_B1_flt <- RunPCA(F_Fast_B1_flt, features = VariableFeatures(object = F_Fast_B1_flt))
F_Fast_B1_flt <- RunUMAP(F_Fast_B1_flt, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fast_B1_flt <- FindNeighbors(F_Fast_B1_flt, dims = 1:30)
F_Fast_B1_flt <- FindClusters(F_Fast_B1_flt, resolution = 0.8)

# plot UMAP projections
DimPlot(F_Fast_B1_flt, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(F_Fast_B1_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B1_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(F_Fast_B1_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B1_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes, min pct 0.33, logfc threshold log2(1.33)
F_Fast_B1_markers_b4_soupx <- FindAllMarkers(F_Fast_B1_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# save UMAP projections and metadata for later use
F_Fast_B1_UMAP <- F_Fast_B1_flt@reductions$umap@cell.embeddings |> as.data.frame()
F_Fast_B1_meta <- F_Fast_B1_flt@meta.data[,c(2,3,7)] |> cbind(F_Fast_B1_UMAP)



#start SoupX in earnest 
# read in filtered and raw samples in, create soupx object
F_Fast_B1_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B1Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fast_B1_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B1Solo.out/GeneFull_Ex50pAS/raw/')
F_Fast_B1_soup <-SoupChannel(F_Fast_B1_raw, F_Fast_B1_flt)

# set dimension reduction to saved UMAP, set clusters to saved clusters
F_Fast_B1_soup = setDR(F_Fast_B1_soup, F_Fast_B1_meta[colnames(F_Fast_B1_flt), c("UMAP_1", "UMAP_2")])
F_Fast_B1_soup = setClusters(F_Fast_B1_soup, setNames(F_Fast_B1_meta$seurat_clusters, colnames(F_Fast_B1_flt)))

# estimate rho
F_Fast_B1_soup = autoEstCont(F_Fast_B1_soup)
# rho estimated at 0.03; however, 0.06 producing cleaner looking clusters based on Agrp expression

# set contamination fraction, adjust counts
F_Fast_B1_soup <- setContaminationFraction(F_Fast_B1_soup, 0.06)
F_Fast_B1out <- adjustCounts(F_Fast_B1_soup, method = 'multinomial')

# SoupX finished, check results
# create seurat object
F_Fast_B1out <- CreateSeuratObject(F_Fast_B1out)

# normalize counts, find variable features excluding xy mt rbp hemoglobin genes
F_Fast_B1out <- SCTransform(F_Fast_B1out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fast_B1out@assays$SCT@var.features <- F_Fast_B1out@assays$SCT@var.features[(!F_Fast_B1out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
F_Fast_B1out <- RunPCA(F_Fast_B1out, features = VariableFeatures(object = F_Fast_B1out))
F_Fast_B1out <- RunUMAP(F_Fast_B1out, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fast_B1out <- FindNeighbors(F_Fast_B1out, dims = 1:30)
F_Fast_B1out <- FindClusters(F_Fast_B1out, resolution = 0.8)

# plot UMAP projections
DimPlot(F_Fast_B1out, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(F_Fast_B1out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B1_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(F_Fast_B1out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B1_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes, min pct 0.33, logfc threshold log2(1.33)
F_Fast_B1_markers_post_soupx <- FindAllMarkers(F_Fast_B1out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# calculate and save pct of mt genes expression
mitoPercent <- Matrix::colSums(F_Fast_B1out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(F_Fast_B1out@assays$RNA@counts)*100
F_Fast_B1out$mitoPercent <- mitoPercent

# calculate and save hemoglobin gene counts
hemoglobin <- Matrix::colSums(F_Fast_B1out@assays$RNA@counts[hemoglobin_genes, ])
F_Fast_B1out <- AddMetaData(object = F_Fast_B1out, metadata = hemoglobin, col.name = "hemoglobin")


# calculate and save pct of rbp gene expression
riboPercent <- Matrix::colSums(F_Fast_B1out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(F_Fast_B1out@assays$RNA@counts)*100
F_Fast_B1out$riboPercent <- riboPercent

# plot QC markers
VlnPlot(F_Fast_B1out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fast_B1_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

# filter extreme values that may be low quality or possible doublet
F_Fast_B1out_filtered <- subset(F_Fast_B1out, 
                                subset = nCount_RNA < (mean(F_Fast_B1out@meta.data$nCount_RNA) + 4*sd(F_Fast_B1out@meta.data$nCount_RNA))
                                & nFeature_RNA < (mean(F_Fast_B1out@meta.data$nFeature_RNA) + 4*sd(F_Fast_B1out@meta.data$nFeature_RNA))
                                & mitoPercent < (mean(F_Fast_B1out@meta.data$mitoPercent) + 4*sd(F_Fast_B1out@meta.data$mitoPercent))
                                & riboPercent < (mean(F_Fast_B1out@meta.data$riboPercent) + 4*sd(F_Fast_B1out@meta.data$riboPercent))
                                & hemoglobin < (mean(F_Fast_B1out@meta.data$hemoglobin) + 4*sd(F_Fast_B1out@meta.data$hemoglobin)))


# plot updated QC markers
VlnPlot(F_Fast_B1out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fast_B1_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
# set default assay to RNA, save only RNA counts, convert to sc experiment
DefaultAssay(F_Fast_B1out_filtered) <- 'RNA'
F_Fast_B1out_filtered <- DietSeurat(F_Fast_B1out_filtered, assays = 'RNA')
F_Fast_B1out_filtered <- as.SingleCellExperiment(F_Fast_B1out_filtered)

# set doublet rate to number of cells / 50000
F_Fast_B1_DBR <- F_Fast_B1out_filtered@colData@nrows/50000
# run scdbl finder, 2000 features, 30 dimensions, random porportion 0.2 and 5 iteration
F_Fast_B1_scDbl <- scDblFinder(F_Fast_B1out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = F_Fast_B1_DBR, propRandom = 0.2, iter = 5)
# convert back to seurat object
F_Fast_B1_scDbl <- as.Seurat(F_Fast_B1_scDbl)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
F_Fast_B1_scDbl <- SCTransform(F_Fast_B1_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fast_B1_scDbl@assays$SCT@var.features <- F_Fast_B1_scDbl@assays$SCT@var.features[(!F_Fast_B1_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
F_Fast_B1_scDbl <- RunPCA(F_Fast_B1_scDbl, features = VariableFeatures(object = F_Fast_B1_scDbl))
F_Fast_B1_scDbl <- RunUMAP(F_Fast_B1_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fast_B1_scDbl <- FindNeighbors(F_Fast_B1_scDbl, dims = 1:30)
F_Fast_B1_scDbl <- FindClusters(F_Fast_B1_scDbl, resolution = 0.8)


# plot singlets and doublets in UMAP space
DimPlot(F_Fast_B1_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/F_Fast_B1_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)

# save only singlets
F_Fast_B1_scDbl <- subset(F_Fast_B1_scDbl, scDblFinder.class == 'singlet')

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
F_Fast_B1_scDbl <- SCTransform(F_Fast_B1_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fast_B1_scDbl@assays$SCT@var.features <- F_Fast_B1_scDbl@assays$SCT@var.features[(!F_Fast_B1_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
F_Fast_B1_scDbl <- RunPCA(F_Fast_B1_scDbl, features = VariableFeatures(object = F_Fast_B1_scDbl))
F_Fast_B1_scDbl <- RunUMAP(F_Fast_B1_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fast_B1_scDbl <- FindNeighbors(F_Fast_B1_scDbl, dims = 1:30)
F_Fast_B1_scDbl <- FindClusters(F_Fast_B1_scDbl, resolution = 0.8)

#check purity of clusters
# plot expression of canonical ARH neuron genes
FeaturePlot(F_Fast_B1_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B1_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(F_Fast_B1_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B1_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes, min pct 0.33, logfc thershold log2(1.33)
F_Fast_B1_markers_post_scDblFnd <- FindAllMarkers(F_Fast_B1_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save F_Fast_B1_scDbl, markers
rm(F_Fast_B1_flt, F_Fast_B1_meta, F_Fast_B1_raw, F_Fast_B1_soup, F_Fast_B1_UMAP, 
   F_Fast_B1out, F_Fast_B1out_filtered, F_Fast_B1_singlets)

# repeat steps for next sample

#F_Fast_B2
#create basic clustering for sample
# read in filtered sample, create seurat object
F_Fast_B2_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B2Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fast_B2_flt <- CreateSeuratObject(F_Fast_B2_flt)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
F_Fast_B2_flt <- SCTransform(F_Fast_B2_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fast_B2_flt@assays$SCT@var.features <- F_Fast_B2_flt@assays$SCT@var.features[(!F_Fast_B2_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
F_Fast_B2_flt <- RunPCA(F_Fast_B2_flt, features = VariableFeatures(object = F_Fast_B2_flt))
F_Fast_B2_flt <- RunUMAP(F_Fast_B2_flt, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fast_B2_flt <- FindNeighbors(F_Fast_B2_flt, dims = 1:30)
F_Fast_B2_flt <- FindClusters(F_Fast_B2_flt, resolution = 0.8)

# plot UMAP projections
DimPlot(F_Fast_B2_flt, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(F_Fast_B2_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B2_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(F_Fast_B2_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B2_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes , min pct 0.33, logfc threshold log2(1.33)
F_Fast_B2_markers_b4_soupx <- FindAllMarkers(F_Fast_B2_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# save UMAP projections and metadata for future use
F_Fast_B2_UMAP <- F_Fast_B2_flt@reductions$umap@cell.embeddings |> as.data.frame()
F_Fast_B2_meta <- F_Fast_B2_flt@meta.data[,c(2,3,7)] |> cbind(F_Fast_B2_UMAP)



#start SoupX in earnest 
# read in filtered and raw samples, create soupx object
F_Fast_B2_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B2Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fast_B2_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B2Solo.out/GeneFull_Ex50pAS/raw/')
F_Fast_B2_soup <-SoupChannel(F_Fast_B2_raw, F_Fast_B2_flt)

# set dimension reduction to saved UMAP, set clusters to saved clusters
F_Fast_B2_soup = setDR(F_Fast_B2_soup, F_Fast_B2_meta[colnames(F_Fast_B2_flt), c("UMAP_1", "UMAP_2")])
F_Fast_B2_soup = setClusters(F_Fast_B2_soup, setNames(F_Fast_B2_meta$seurat_clusters, colnames(F_Fast_B2_flt)))

# estimate rho
F_Fast_B2_soup = autoEstCont(F_Fast_B2_soup)
# rho estimated at 0.06; however, 0.12 producing cleaner looking clusters based on Agrp expression

# set contamination fraction, adjust counts
F_Fast_B2_soup <- setContaminationFraction(F_Fast_B2_soup, 0.12)
F_Fast_B2out <- adjustCounts(F_Fast_B2_soup, method = 'multinomial')

# SoupX finished, check results
# create seurat object
F_Fast_B2out <- CreateSeuratObject(F_Fast_B2out)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
F_Fast_B2out <- SCTransform(F_Fast_B2out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fast_B2out@assays$SCT@var.features <- F_Fast_B2out@assays$SCT@var.features[(!F_Fast_B2out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculat UMAP based on first 30 PCs
F_Fast_B2out <- RunPCA(F_Fast_B2out, features = VariableFeatures(object = F_Fast_B2out))
F_Fast_B2out <- RunUMAP(F_Fast_B2out, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fast_B2out <- FindNeighbors(F_Fast_B2out, dims = 1:30)
F_Fast_B2out <- FindClusters(F_Fast_B2out, resolution = 0.8)

# plot UMAP projections
DimPlot(F_Fast_B2out, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(F_Fast_B2out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B2_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(F_Fast_B2out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B2_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes , min pct 0.33, logfc threshold log2(1.33)
F_Fast_B2_markers_post_soupx <- FindAllMarkers(F_Fast_B2out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# calculate and save pct of mt gene expression
mitoPercent <- Matrix::colSums(F_Fast_B2out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(F_Fast_B2out@assays$RNA@counts)*100
F_Fast_B2out$mitoPercent <- mitoPercent

# calculate and save hemoglobin gene expression counts
hemoglobin <- Matrix::colSums(F_Fast_B2out@assays$RNA@counts[hemoglobin_genes, ])
F_Fast_B2out <- AddMetaData(object = F_Fast_B2out, metadata = hemoglobin, col.name = "hemoglobin")


# calculate and save pct of rbp gene expression
riboPercent <- Matrix::colSums(F_Fast_B2out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(F_Fast_B2out@assays$RNA@counts)*100
F_Fast_B2out$riboPercent <- riboPercent

# plot QC markers
VlnPlot(F_Fast_B2out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fast_B2_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

# filter out extreme values that may be low quality or possible doublet
F_Fast_B2out_filtered <- subset(F_Fast_B2out, 
                                subset = nCount_RNA < (mean(F_Fast_B2out@meta.data$nCount_RNA) + 4*sd(F_Fast_B2out@meta.data$nCount_RNA))
                                & nFeature_RNA < (mean(F_Fast_B2out@meta.data$nFeature_RNA) + 4*sd(F_Fast_B2out@meta.data$nFeature_RNA))
                                & mitoPercent < (mean(F_Fast_B2out@meta.data$mitoPercent) + 4*sd(F_Fast_B2out@meta.data$mitoPercent))
                                & riboPercent < (mean(F_Fast_B2out@meta.data$riboPercent) + 4*sd(F_Fast_B2out@meta.data$riboPercent))
                                & hemoglobin < (mean(F_Fast_B2out@meta.data$hemoglobin) + 4*sd(F_Fast_B2out@meta.data$hemoglobin)))


# plot updated QC markers
VlnPlot(F_Fast_B2out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fast_B2_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
# set default assay to RNA, save only RNA counts, transform to sc experiment
DefaultAssay(F_Fast_B2out_filtered) <- 'RNA'
F_Fast_B2out_filtered <- DietSeurat(F_Fast_B2out_filtered, assays = 'RNA')
F_Fast_B2out_filtered <- as.SingleCellExperiment(F_Fast_B2out_filtered)

# set doublet rate to number of cells / 50000
F_Fast_B2_DBR <- F_Fast_B2out_filtered@colData@nrows/50000
# run scdbl finder, 2000 features, 30 dims, random porportion 0.2, 5 iterations
F_Fast_B2_scDbl <- scDblFinder(F_Fast_B2out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = F_Fast_B2_DBR, propRandom = 0.2, iter = 5)
# transform back to seurat object
F_Fast_B2_scDbl <- as.Seurat(F_Fast_B2_scDbl)

# normalize expression , find variable features excluding xy mt rbp hemoglobin genes
F_Fast_B2_scDbl <- SCTransform(F_Fast_B2_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fast_B2_scDbl@assays$SCT@var.features <- F_Fast_B2_scDbl@assays$SCT@var.features[(!F_Fast_B2_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
F_Fast_B2_scDbl <- RunPCA(F_Fast_B2_scDbl, features = VariableFeatures(object = F_Fast_B2_scDbl))
F_Fast_B2_scDbl <- RunUMAP(F_Fast_B2_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fast_B2_scDbl <- FindNeighbors(F_Fast_B2_scDbl, dims = 1:30)
F_Fast_B2_scDbl <- FindClusters(F_Fast_B2_scDbl, resolution = 0.8)


# plot singlet and doublet in UMAP space
DimPlot(F_Fast_B2_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/F_Fast_B2_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)

# save only singlets
F_Fast_B2_scDbl <- subset(F_Fast_B2_scDbl, scDblFinder.class == 'singlet')

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
F_Fast_B2_scDbl <- SCTransform(F_Fast_B2_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fast_B2_scDbl@assays$SCT@var.features <- F_Fast_B2_scDbl@assays$SCT@var.features[(!F_Fast_B2_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
F_Fast_B2_scDbl <- RunPCA(F_Fast_B2_scDbl, features = VariableFeatures(object = F_Fast_B2_scDbl))
F_Fast_B2_scDbl <- RunUMAP(F_Fast_B2_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fast_B2_scDbl <- FindNeighbors(F_Fast_B2_scDbl, dims = 1:30)
F_Fast_B2_scDbl <- FindClusters(F_Fast_B2_scDbl, resolution = 0.8)

#check purity of clusters
# plot expression of canonical ARH neuron marker genes
FeaturePlot(F_Fast_B2_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B2_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(F_Fast_B2_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B2_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes , min pct 0.33 , logfc threshold log2(1.33)
F_Fast_B2_markers_post_scDblFnd <- FindAllMarkers(F_Fast_B2_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save F_Fast_B2_scDbl, markers
rm(F_Fast_B2_flt, F_Fast_B2_meta, F_Fast_B2_raw, F_Fast_B2_soup, F_Fast_B2_UMAP, 
   F_Fast_B2out, F_Fast_B2out_filtered, F_Fast_B2_singlets)

# repeat steps for next sample

#F_Fast_B3
#create basic clustering for sample
# read in filtered sample, create seurat object
F_Fast_B3_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B3Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fast_B3_flt <- CreateSeuratObject(F_Fast_B3_flt)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
F_Fast_B3_flt <- SCTransform(F_Fast_B3_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fast_B3_flt@assays$SCT@var.features <- F_Fast_B3_flt@assays$SCT@var.features[(!F_Fast_B3_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
F_Fast_B3_flt <- RunPCA(F_Fast_B3_flt, features = VariableFeatures(object = F_Fast_B3_flt))
F_Fast_B3_flt <- RunUMAP(F_Fast_B3_flt, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fast_B3_flt <- FindNeighbors(F_Fast_B3_flt, dims = 1:30)
F_Fast_B3_flt <- FindClusters(F_Fast_B3_flt, resolution = 0.8)

# plot UMAP projections
DimPlot(F_Fast_B3_flt, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(F_Fast_B3_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B3_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(F_Fast_B3_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B3_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes , min pct 0.33 , logfc threshold log2(1.33)
F_Fast_B3_markers_b4_soupx <- FindAllMarkers(F_Fast_B3_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# save UMAP projections , and metadata for later use
F_Fast_B3_UMAP <- F_Fast_B3_flt@reductions$umap@cell.embeddings |> as.data.frame()
F_Fast_B3_meta <- F_Fast_B3_flt@meta.data[,c(2,3,7)] |> cbind(F_Fast_B3_UMAP)



#start SoupX in earnest
# read in filtered and raw samples and create soupx object
F_Fast_B3_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B3Solo.out/GeneFull_Ex50pAS/filtered/')
F_Fast_B3_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/F_Fast_B3Solo.out/GeneFull_Ex50pAS/raw/')
F_Fast_B3_soup <-SoupChannel(F_Fast_B3_raw, F_Fast_B3_flt)

# set dimension reduction to saved UMAP, set clusters to saved clusters
F_Fast_B3_soup = setDR(F_Fast_B3_soup, F_Fast_B3_meta[colnames(F_Fast_B3_flt), c("UMAP_1", "UMAP_2")])
F_Fast_B3_soup = setClusters(F_Fast_B3_soup, setNames(F_Fast_B3_meta$seurat_clusters, colnames(F_Fast_B3_flt)))

# estimate rho
F_Fast_B3_soup = autoEstCont(F_Fast_B3_soup)
# rho estimated at 0.04; however, 0.10 producing cleaner looking clusters based on Agrp expression

# set contamination fraction, and adjust counts
F_Fast_B3_soup <- setContaminationFraction(F_Fast_B3_soup, 0.10)
F_Fast_B3out <- adjustCounts(F_Fast_B3_soup, method = 'multinomial')

# SoupX finished, check results
# create seurat object
F_Fast_B3out <- CreateSeuratObject(F_Fast_B3out)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
F_Fast_B3out <- SCTransform(F_Fast_B3out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fast_B3out@assays$SCT@var.features <- F_Fast_B3out@assays$SCT@var.features[(!F_Fast_B3out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
F_Fast_B3out <- RunPCA(F_Fast_B3out, features = VariableFeatures(object = F_Fast_B3out))
F_Fast_B3out <- RunUMAP(F_Fast_B3out, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fast_B3out <- FindNeighbors(F_Fast_B3out, dims = 1:30)
F_Fast_B3out <- FindClusters(F_Fast_B3out, resolution = 0.8)

# plot UMAP
DimPlot(F_Fast_B3out, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(F_Fast_B3out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B3_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(F_Fast_B3out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B3_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes , min pct 0.33 , logfc threshold log2(1.33)
F_Fast_B3_markers_post_soupx <- FindAllMarkers(F_Fast_B3out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# calculate and save pct mt gene expression
mitoPercent <- Matrix::colSums(F_Fast_B3out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(F_Fast_B3out@assays$RNA@counts)*100
F_Fast_B3out$mitoPercent <- mitoPercent

# calculate and save hemoglobin gene counts
hemoglobin <- Matrix::colSums(F_Fast_B3out@assays$RNA@counts[hemoglobin_genes, ])
F_Fast_B3out <- AddMetaData(object = F_Fast_B3out, metadata = hemoglobin, col.name = "hemoglobin")


# calculate pct rbp gene expression
riboPercent <- Matrix::colSums(F_Fast_B3out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(F_Fast_B3out@assays$RNA@counts)*100
F_Fast_B3out$riboPercent <- riboPercent

# plot QC markers
VlnPlot(F_Fast_B3out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fast_B3_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

# filter out extreme values that may be low quality or possible doublet (mean + 4sd)
F_Fast_B3out_filtered <- subset(F_Fast_B3out, 
                                subset = nCount_RNA < (mean(F_Fast_B3out@meta.data$nCount_RNA) + 4*sd(F_Fast_B3out@meta.data$nCount_RNA))
                                & nFeature_RNA < (mean(F_Fast_B3out@meta.data$nFeature_RNA) + 4*sd(F_Fast_B3out@meta.data$nFeature_RNA))
                                & mitoPercent < (mean(F_Fast_B3out@meta.data$mitoPercent) + 4*sd(F_Fast_B3out@meta.data$mitoPercent))
                                & riboPercent < (mean(F_Fast_B3out@meta.data$riboPercent) + 4*sd(F_Fast_B3out@meta.data$riboPercent))
                                & hemoglobin < (mean(F_Fast_B3out@meta.data$hemoglobin) + 4*sd(F_Fast_B3out@meta.data$hemoglobin)))

# plot updated QC markers
VlnPlot(F_Fast_B3out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/F_Fast_B3_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
# set default assay to RNA, save only RNA counts, transform into sc experiment
DefaultAssay(F_Fast_B3out_filtered) <- 'RNA'
F_Fast_B3out_filtered <- DietSeurat(F_Fast_B3out_filtered, assays = 'RNA')
F_Fast_B3out_filtered <- as.SingleCellExperiment(F_Fast_B3out_filtered)

# set doublet rate at number of cells / 50000
F_Fast_B3_DBR <- F_Fast_B3out_filtered@colData@nrows/50000
# run scdbl finder with 2000 features , 30 dims , random porportion 0.2 and 5 iterations
F_Fast_B3_scDbl <- scDblFinder(F_Fast_B3out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = F_Fast_B3_DBR, propRandom = 0.2, iter = 5)
# transform back to seurat object
F_Fast_B3_scDbl <- as.Seurat(F_Fast_B3_scDbl)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
F_Fast_B3_scDbl <- SCTransform(F_Fast_B3_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fast_B3_scDbl@assays$SCT@var.features <- F_Fast_B3_scDbl@assays$SCT@var.features[(!F_Fast_B3_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
F_Fast_B3_scDbl <- RunPCA(F_Fast_B3_scDbl, features = VariableFeatures(object = F_Fast_B3_scDbl))
F_Fast_B3_scDbl <- RunUMAP(F_Fast_B3_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
F_Fast_B3_scDbl <- FindNeighbors(F_Fast_B3_scDbl, dims = 1:30)
F_Fast_B3_scDbl <- FindClusters(F_Fast_B3_scDbl, resolution = 0.8)


# plot singlets and doublets in UMAP space
DimPlot(F_Fast_B3_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/F_Fast_B3_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)

# save only singlets
F_Fast_B3_scDbl <- subset(F_Fast_B3_scDbl, scDblFinder.class == 'singlet')

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
F_Fast_B3_scDbl <- SCTransform(F_Fast_B3_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
F_Fast_B3_scDbl@assays$SCT@var.features <- F_Fast_B3_scDbl@assays$SCT@var.features[(!F_Fast_B3_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculat PCs based on variable features, calculate UMAP based on first 30 PCs
F_Fast_B3_scDbl <- RunPCA(F_Fast_B3_scDbl, features = VariableFeatures(object = F_Fast_B3_scDbl))
F_Fast_B3_scDbl <- RunUMAP(F_Fast_B3_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 1.5 resolution
F_Fast_B3_scDbl <- FindNeighbors(F_Fast_B3_scDbl, dims = 1:30)
F_Fast_B3_scDbl <- FindClusters(F_Fast_B3_scDbl, resolution = 1.5)

# plot UMAP projections
DimPlot(F_Fast_B3_scDbl, label = TRUE) + NoLegend()
#check purity of clusters
# plot expression of canonical ARH neuron genes
FeaturePlot(F_Fast_B3_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/F_Fast_B3_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(F_Fast_B3_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/F_Fast_B3_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all maeker genes , min pct 0.33 , logfc threshold log2(1.33)
F_Fast_B3_markers_post_scDblFnd <- FindAllMarkers(F_Fast_B3_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save F_Fast_B3_scDbl, markers
rm(F_Fast_B3_flt, F_Fast_B3_meta, F_Fast_B3_raw, F_Fast_B3_soup, F_Fast_B3_UMAP, 
   F_Fast_B3out, F_Fast_B3out_filtered, F_Fast_B3_singlets)

# repeat steps for next sample

#M_Fed_B1
#create basic clustering for sample
# read in filtered sample and create seurat object
M_Fed_B1_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B1Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fed_B1_flt <- CreateSeuratObject(M_Fed_B1_flt)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fed_B1_flt <- SCTransform(M_Fed_B1_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fed_B1_flt@assays$SCT@var.features <- M_Fed_B1_flt@assays$SCT@var.features[(!M_Fed_B1_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

# calculate PCs based on varibale features, calculate UMAP based on first 30 PCs
M_Fed_B1_flt <- RunPCA(M_Fed_B1_flt, features = VariableFeatures(object = M_Fed_B1_flt))
M_Fed_B1_flt <- RunUMAP(M_Fed_B1_flt, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
M_Fed_B1_flt <- FindNeighbors(M_Fed_B1_flt, dims = 1:30)
M_Fed_B1_flt <- FindClusters(M_Fed_B1_flt, resolution = 0.8)

# plot UMAP projections
DimPlot(M_Fed_B1_flt, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fed_B1_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B1_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fed_B1_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B1_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes , min pct 0.33 , logfc threshold log2(1.33)
M_Fed_B1_markers_b4_soupx <- FindAllMarkers(M_Fed_B1_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# save UMAP and metadata for future use
M_Fed_B1_UMAP <- M_Fed_B1_flt@reductions$umap@cell.embeddings |> as.data.frame()
M_Fed_B1_meta <- M_Fed_B1_flt@meta.data[,c(2,3,7)] |> cbind(M_Fed_B1_UMAP)



#start SoupX in earnest 
# read in filtered and raw samples, create soupx object
M_Fed_B1_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B1Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fed_B1_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B1Solo.out/GeneFull_Ex50pAS/raw/')
M_Fed_B1_soup <-SoupChannel(M_Fed_B1_raw, M_Fed_B1_flt)

# set dimension reduction to saved UMAP set clusters to saved clusters
M_Fed_B1_soup = setDR(M_Fed_B1_soup, M_Fed_B1_meta[colnames(M_Fed_B1_flt), c("UMAP_1", "UMAP_2")])
M_Fed_B1_soup = setClusters(M_Fed_B1_soup, setNames(M_Fed_B1_meta$seurat_clusters, colnames(M_Fed_B1_flt)))

# estimate rho
M_Fed_B1_soup = autoEstCont(M_Fed_B1_soup)
# rho estimated at 0.06; however, 0.15 producing cleaner looking clusters based on Agrp expression

# set contamination fraction and adjust counts
M_Fed_B1_soup <- setContaminationFraction(M_Fed_B1_soup, 0.15)
M_Fed_B1out <- adjustCounts(M_Fed_B1_soup, method = 'multinomial')

# SoupX finished, check results
# create seurat object
M_Fed_B1out <- CreateSeuratObject(M_Fed_B1out)

# normalize expression , find varibale features excluding xy mt rbp hemoglobin genes
M_Fed_B1out <- SCTransform(M_Fed_B1out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fed_B1out@assays$SCT@var.features <- M_Fed_B1out@assays$SCT@var.features[(!M_Fed_B1out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fed_B1out <- RunPCA(M_Fed_B1out, features = VariableFeatures(object = M_Fed_B1out))
M_Fed_B1out <- RunUMAP(M_Fed_B1out, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
M_Fed_B1out <- FindNeighbors(M_Fed_B1out, dims = 1:30)
M_Fed_B1out <- FindClusters(M_Fed_B1out, resolution = 0.8)

# plot UMAP projection
DimPlot(M_Fed_B1out, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fed_B1out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B1_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fed_B1out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B1_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes , min pct 0.33 , logfc threshold log2(1.33)
M_Fed_B1_markers_post_soupx <- FindAllMarkers(M_Fed_B1out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# calculate and save pct mt gene expression
mitoPercent <- Matrix::colSums(M_Fed_B1out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(M_Fed_B1out@assays$RNA@counts)*100
M_Fed_B1out$mitoPercent <- mitoPercent

# calculate and save hemoglobin gene counts
hemoglobin <- Matrix::colSums(M_Fed_B1out@assays$RNA@counts[hemoglobin_genes, ])
M_Fed_B1out <- AddMetaData(object = M_Fed_B1out, metadata = hemoglobin, col.name = "hemoglobin")


# calculate and save rbp pct gene expression
riboPercent <- Matrix::colSums(M_Fed_B1out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(M_Fed_B1out@assays$RNA@counts)*100
M_Fed_B1out$riboPercent <- riboPercent

# plot QC markers
VlnPlot(M_Fed_B1out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fed_B1_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

# filter extreme values that may be low quality or possible doublet
M_Fed_B1out_filtered <- subset(M_Fed_B1out, 
                               subset = nCount_RNA < (mean(M_Fed_B1out@meta.data$nCount_RNA) + 4*sd(M_Fed_B1out@meta.data$nCount_RNA))
                               & nFeature_RNA < (mean(M_Fed_B1out@meta.data$nFeature_RNA) + 4*sd(M_Fed_B1out@meta.data$nFeature_RNA))
                               & mitoPercent < (mean(M_Fed_B1out@meta.data$mitoPercent) + 4*sd(M_Fed_B1out@meta.data$mitoPercent))
                               & riboPercent < (mean(M_Fed_B1out@meta.data$riboPercent) + 4*sd(M_Fed_B1out@meta.data$riboPercent))
                               & hemoglobin < (mean(M_Fed_B1out@meta.data$hemoglobin) + 4*sd(M_Fed_B1out@meta.data$hemoglobin)))

# plot updated QC markers
VlnPlot(M_Fed_B1out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fed_B1_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
# set default assay to RNA, save only RNA counts, transform to sc experiment
DefaultAssay(M_Fed_B1out_filtered) <- 'RNA'
M_Fed_B1out_filtered <- DietSeurat(M_Fed_B1out_filtered, assays = 'RNA')
M_Fed_B1out_filtered <- as.SingleCellExperiment(M_Fed_B1out_filtered)

# set doublet rate at number of cells / 50000
M_Fed_B1_DBR <- M_Fed_B1out_filtered@colData@nrows/50000
# run scbbl finder , 2000 features , 30 dims, random porportion 0.2, 5 iterations
M_Fed_B1_scDbl <- scDblFinder(M_Fed_B1out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = M_Fed_B1_DBR, propRandom = 0.2, iter = 5)
M_Fed_B1_scDbl <- as.Seurat(M_Fed_B1_scDbl)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fed_B1_scDbl <- SCTransform(M_Fed_B1_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fed_B1_scDbl@assays$SCT@var.features <- M_Fed_B1_scDbl@assays$SCT@var.features[(!M_Fed_B1_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fed_B1_scDbl <- RunPCA(M_Fed_B1_scDbl, features = VariableFeatures(object = M_Fed_B1_scDbl))
M_Fed_B1_scDbl <- RunUMAP(M_Fed_B1_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
M_Fed_B1_scDbl <- FindNeighbors(M_Fed_B1_scDbl, dims = 1:30)
M_Fed_B1_scDbl <- FindClusters(M_Fed_B1_scDbl, resolution = 0.8)


# plot singlets and doublets in UMAP space
DimPlot(M_Fed_B1_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/M_Fed_B1_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)

# save only singlets
M_Fed_B1_scDbl <- subset(M_Fed_B1_scDbl, scDblFinder.class == 'singlet')

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fed_B1_scDbl <- SCTransform(M_Fed_B1_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fed_B1_scDbl@assays$SCT@var.features <- M_Fed_B1_scDbl@assays$SCT@var.features[(!M_Fed_B1_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fed_B1_scDbl <- RunPCA(M_Fed_B1_scDbl, features = VariableFeatures(object = M_Fed_B1_scDbl))
M_Fed_B1_scDbl <- RunUMAP(M_Fed_B1_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
M_Fed_B1_scDbl <- FindNeighbors(M_Fed_B1_scDbl, dims = 1:30)
M_Fed_B1_scDbl <- FindClusters(M_Fed_B1_scDbl, resolution = 0.8)

#check purity of clusters
# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fed_B1_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B1_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fed_B1_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B1_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes, min pct 0.33, logfc threshold log2(1.33)
M_Fed_B1_markers_post_scDblFnd <- FindAllMarkers(M_Fed_B1_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save M_Fed_B1_scDbl, markers
rm(M_Fed_B1_flt, M_Fed_B1_meta, M_Fed_B1_raw, M_Fed_B1_soup, M_Fed_B1_UMAP, 
   M_Fed_B1out, M_Fed_B1out_filtered, M_Fed_B1_singlets)

# repeat steps for next sample

#M_Fed_B2
#create basic clustering for sample
# read in filtered sample, create seurat object
M_Fed_B2_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B2Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fed_B2_flt <- CreateSeuratObject(M_Fed_B2_flt)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fed_B2_flt <- SCTransform(M_Fed_B2_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fed_B2_flt@assays$SCT@var.features <- M_Fed_B2_flt@assays$SCT@var.features[(!M_Fed_B2_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

# calculate PCs based on variable features , calculate UMAP based on first 30 PCs
M_Fed_B2_flt <- RunPCA(M_Fed_B2_flt, features = VariableFeatures(object = M_Fed_B2_flt))
M_Fed_B2_flt <- RunUMAP(M_Fed_B2_flt, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
M_Fed_B2_flt <- FindNeighbors(M_Fed_B2_flt, dims = 1:30)
M_Fed_B2_flt <- FindClusters(M_Fed_B2_flt, resolution = 0.8)

# plot UMAP projections 
DimPlot(M_Fed_B2_flt, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fed_B2_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B2_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fed_B2_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B2_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes, min pct 0.33 , logfc threshold log2(1.33)
M_Fed_B2_markers_b4_soupx <- FindAllMarkers(M_Fed_B2_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# save UMAP projections and metadata for future use
M_Fed_B2_UMAP <- M_Fed_B2_flt@reductions$umap@cell.embeddings |> as.data.frame()
M_Fed_B2_meta <- M_Fed_B2_flt@meta.data[,c(2,3,7)] |> cbind(M_Fed_B2_UMAP)



#start SoupX in earnest 
# read in filtered and raw samples , create soupx object
M_Fed_B2_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B2Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fed_B2_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B2Solo.out/GeneFull_Ex50pAS/raw/')
M_Fed_B2_soup <-SoupChannel(M_Fed_B2_raw, M_Fed_B2_flt)

# set dimension reduction to saved UMAP , set clusters to saved clusters
M_Fed_B2_soup = setDR(M_Fed_B2_soup, M_Fed_B2_meta[colnames(M_Fed_B2_flt), c("UMAP_1", "UMAP_2")])
M_Fed_B2_soup = setClusters(M_Fed_B2_soup, setNames(M_Fed_B2_meta$seurat_clusters, colnames(M_Fed_B2_flt)))

# estimate contamination fraction
M_Fed_B2_soup = autoEstCont(M_Fed_B2_soup)
# rho estimated at 0.03; however, 0.105 producing cleaner looking clusters

# set contamination fraction and adjust counts
M_Fed_B2_soup <- setContaminationFraction(M_Fed_B2_soup, 0.105)
M_Fed_B2out <- adjustCounts(M_Fed_B2_soup, method = 'multinomial')

# SoupX finished, check results
# create seurat object
M_Fed_B2out <- CreateSeuratObject(M_Fed_B2out)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fed_B2out <- SCTransform(M_Fed_B2out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fed_B2out@assays$SCT@var.features <- M_Fed_B2out@assays$SCT@var.features[(!M_Fed_B2out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fed_B2out <- RunPCA(M_Fed_B2out, features = VariableFeatures(object = M_Fed_B2out))
M_Fed_B2out <- RunUMAP(M_Fed_B2out, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 1.5 resolution
M_Fed_B2out <- FindNeighbors(M_Fed_B2out, dims = 1:30)
M_Fed_B2out <- FindClusters(M_Fed_B2out, resolution = 1.5)

# plot UMAP projections
DimPlot(M_Fed_B2out, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fed_B2out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B2_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fed_B2out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B2_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

M_Fed_B2_markers_post_soupx <- FindAllMarkers(M_Fed_B2out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# calculate and save mt pct gene expression
mitoPercent <- Matrix::colSums(M_Fed_B2out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(M_Fed_B2out@assays$RNA@counts)*100
M_Fed_B2out$mitoPercent <- mitoPercent

# calculate and save hemoglobin gene expression count
hemoglobin <- Matrix::colSums(M_Fed_B2out@assays$RNA@counts[hemoglobin_genes, ])
M_Fed_B2out <- AddMetaData(object = M_Fed_B2out, metadata = hemoglobin, col.name = "hemoglobin")


# calculate and save rbp pct gene expression
riboPercent <- Matrix::colSums(M_Fed_B2out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(M_Fed_B2out@assays$RNA@counts)*100
M_Fed_B2out$riboPercent <- riboPercent

# plot QC markers
VlnPlot(M_Fed_B2out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fed_B2_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

# filter out extreme values that may be low quality or possible doublet
M_Fed_B2out_filtered <- subset(M_Fed_B2out, 
                               subset = nCount_RNA < (mean(M_Fed_B2out@meta.data$nCount_RNA) + 4*sd(M_Fed_B2out@meta.data$nCount_RNA))
                               & nFeature_RNA < (mean(M_Fed_B2out@meta.data$nFeature_RNA) + 4*sd(M_Fed_B2out@meta.data$nFeature_RNA))
                               & mitoPercent < (mean(M_Fed_B2out@meta.data$mitoPercent) + 4*sd(M_Fed_B2out@meta.data$mitoPercent))
                               & riboPercent < (mean(M_Fed_B2out@meta.data$riboPercent) + 4*sd(M_Fed_B2out@meta.data$riboPercent))
                               & hemoglobin < (mean(M_Fed_B2out@meta.data$hemoglobin) + 4*sd(M_Fed_B2out@meta.data$hemoglobin)))

# plot updated QC markers
VlnPlot(M_Fed_B2out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fed_B2_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
# set default assay to RNA, save only RNA counts, transform to sc experiment
DefaultAssay(M_Fed_B2out_filtered) <- 'RNA'
M_Fed_B2out_filtered <- DietSeurat(M_Fed_B2out_filtered, assays = 'RNA')
M_Fed_B2out_filtered <- as.SingleCellExperiment(M_Fed_B2out_filtered)

# set doublet rate to number of cells / 50000
M_Fed_B2_DBR <- M_Fed_B2out_filtered@colData@nrows/50000
# run scdbl finder , 2000 features , 30 dims, random porportion 0.2 , 5 iteration
M_Fed_B2_scDbl <- scDblFinder(M_Fed_B2out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = M_Fed_B2_DBR, propRandom = 0.2, iter = 5)
# transform back to seurat object
M_Fed_B2_scDbl <- as.Seurat(M_Fed_B2_scDbl)

# normalize expression , find variable features excluding xy mt rbp hemoglobin genes
M_Fed_B2_scDbl <- SCTransform(M_Fed_B2_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fed_B2_scDbl@assays$SCT@var.features <- M_Fed_B2_scDbl@assays$SCT@var.features[(!M_Fed_B2_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fed_B2_scDbl <- RunPCA(M_Fed_B2_scDbl, features = VariableFeatures(object = M_Fed_B2_scDbl))
M_Fed_B2_scDbl <- RunUMAP(M_Fed_B2_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 0.8 resolution 
M_Fed_B2_scDbl <- FindNeighbors(M_Fed_B2_scDbl, dims = 1:30)
M_Fed_B2_scDbl <- FindClusters(M_Fed_B2_scDbl, resolution = 0.8)


# plot singlets and doublets in UMAP space
DimPlot(M_Fed_B2_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/M_Fed_B2_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)

# save singlets only
M_Fed_B2_scDbl <- subset(M_Fed_B2_scDbl, scDblFinder.class == 'singlet')

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fed_B2_scDbl <- SCTransform(M_Fed_B2_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fed_B2_scDbl@assays$SCT@var.features <- M_Fed_B2_scDbl@assays$SCT@var.features[(!M_Fed_B2_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fed_B2_scDbl <- RunPCA(M_Fed_B2_scDbl, features = VariableFeatures(object = M_Fed_B2_scDbl))
M_Fed_B2_scDbl <- RunUMAP(M_Fed_B2_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 1.2 resolution
M_Fed_B2_scDbl <- FindNeighbors(M_Fed_B2_scDbl, dims = 1:30)
M_Fed_B2_scDbl <- FindClusters(M_Fed_B2_scDbl, resolution = 1.2)
DimPlot(M_Fed_B2_scDbl, label = TRUE) + NoLegend()
#check purity of clusters
# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fed_B2_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B2_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fed_B2_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B2_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes , min pct 0.33, logfc threshold log2(1.33)
M_Fed_B2_markers_post_scDblFnd <- FindAllMarkers(M_Fed_B2_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save M_Fed_B2_scDbl, markers
rm(M_Fed_B2_flt, M_Fed_B2_meta, M_Fed_B2_raw, M_Fed_B2_soup, M_Fed_B2_UMAP, 
   M_Fed_B2out, M_Fed_B2out_filtered, M_Fed_B2_singlets)

# repeat steps for next sample

#M_Fed_B3
#create basic clustering for sample
# read in filtered sample, create seurat object
M_Fed_B3_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B3Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fed_B3_flt <- CreateSeuratObject(M_Fed_B3_flt)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fed_B3_flt <- SCTransform(M_Fed_B3_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fed_B3_flt@assays$SCT@var.features <- M_Fed_B3_flt@assays$SCT@var.features[(!M_Fed_B3_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fed_B3_flt <- RunPCA(M_Fed_B3_flt, features = VariableFeatures(object = M_Fed_B3_flt))
M_Fed_B3_flt <- RunUMAP(M_Fed_B3_flt, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
M_Fed_B3_flt <- FindNeighbors(M_Fed_B3_flt, dims = 1:30)
M_Fed_B3_flt <- FindClusters(M_Fed_B3_flt, resolution = 0.8)

# plot UMAP projections
DimPlot(M_Fed_B3_flt, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fed_B3_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B3_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fed_B3_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B3_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes , min pct 0.33, logfc threshold log2(1.33)
M_Fed_B3_markers_b4_soupx <- FindAllMarkers(M_Fed_B3_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# save UMAP and metadata for later use
M_Fed_B3_UMAP <- M_Fed_B3_flt@reductions$umap@cell.embeddings |> as.data.frame()
M_Fed_B3_meta <- M_Fed_B3_flt@meta.data[,c(2,3,7)] |> cbind(M_Fed_B3_UMAP)



#start SoupX in earnest 
# read in filtered and raw samples, create soupx object
M_Fed_B3_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B3Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fed_B3_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fed_B3Solo.out/GeneFull_Ex50pAS/raw/')
M_Fed_B3_soup <-SoupChannel(M_Fed_B3_raw, M_Fed_B3_flt)

# set dimension reduction to saved UMAP, set clusters to saved clusters
M_Fed_B3_soup = setDR(M_Fed_B3_soup, M_Fed_B3_meta[colnames(M_Fed_B3_flt), c("UMAP_1", "UMAP_2")])
M_Fed_B3_soup = setClusters(M_Fed_B3_soup, setNames(M_Fed_B3_meta$seurat_clusters, colnames(M_Fed_B3_flt)))

# estimate contamination fraction
M_Fed_B3_soup = autoEstCont(M_Fed_B3_soup)
# rho estimated at 0.07; however, 0.20 producing cleaner looking clusters

# set contamination fraction, adjust counts
M_Fed_B3_soup <- setContaminationFraction(M_Fed_B3_soup, 0.20)
M_Fed_B3out <- adjustCounts(M_Fed_B3_soup, method = 'multinomial')

# SoupX finished, check results
# create seurat object
M_Fed_B3out <- CreateSeuratObject(M_Fed_B3out)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fed_B3out <- SCTransform(M_Fed_B3out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fed_B3out@assays$SCT@var.features <- M_Fed_B3out@assays$SCT@var.features[(!M_Fed_B3out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fed_B3out <- RunPCA(M_Fed_B3out, features = VariableFeatures(object = M_Fed_B3out))
M_Fed_B3out <- RunUMAP(M_Fed_B3out, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 1.5 resolution
M_Fed_B3out <- FindNeighbors(M_Fed_B3out, dims = 1:30)
M_Fed_B3out <- FindClusters(M_Fed_B3out, resolution = 1.5)

# plot UMAP projections
DimPlot(M_Fed_B3out, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fed_B3out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B3_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fed_B3out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B3_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes , min pct 0.33, logfc threshold log2(1.33)
M_Fed_B3_markers_post_soupx <- FindAllMarkers(M_Fed_B3out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# calculate and save mt genes expression pct 
mitoPercent <- Matrix::colSums(M_Fed_B3out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(M_Fed_B3out@assays$RNA@counts)*100
M_Fed_B3out$mitoPercent <- mitoPercent

# calculate and save hemoglobin genes counts
hemoglobin <- Matrix::colSums(M_Fed_B3out@assays$RNA@counts[hemoglobin_genes, ])
M_Fed_B3out <- AddMetaData(object = M_Fed_B3out, metadata = hemoglobin, col.name = "hemoglobin")


# calculate and save rbp genes pct
riboPercent <- Matrix::colSums(M_Fed_B3out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(M_Fed_B3out@assays$RNA@counts)*100
M_Fed_B3out$riboPercent <- riboPercent

# plot QC markers
VlnPlot(M_Fed_B3out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fed_B3_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

# filtered out extreme values that may be low quality or possible doublet
M_Fed_B3out_filtered <- subset(M_Fed_B3out, 
                               subset = nCount_RNA < (mean(M_Fed_B3out@meta.data$nCount_RNA) + 4*sd(M_Fed_B3out@meta.data$nCount_RNA))
                               & nFeature_RNA < (mean(M_Fed_B3out@meta.data$nFeature_RNA) + 4*sd(M_Fed_B3out@meta.data$nFeature_RNA))
                               & mitoPercent < (mean(M_Fed_B3out@meta.data$mitoPercent) + 4*sd(M_Fed_B3out@meta.data$mitoPercent))
                               & riboPercent < (mean(M_Fed_B3out@meta.data$riboPercent) + 4*sd(M_Fed_B3out@meta.data$riboPercent))
                               & hemoglobin < (mean(M_Fed_B3out@meta.data$hemoglobin) + 4*sd(M_Fed_B3out@meta.data$hemoglobin)))

# plot updated QC markers
VlnPlot(M_Fed_B3out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fed_B3_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
# set default assay to RNA, save RNA counts only, transform to sc experiment
DefaultAssay(M_Fed_B3out_filtered) <- 'RNA'
M_Fed_B3out_filtered <- DietSeurat(M_Fed_B3out_filtered, assays = 'RNA')
M_Fed_B3out_filtered <- as.SingleCellExperiment(M_Fed_B3out_filtered)

# set doublet rate to number of cells / 50000
M_Fed_B3_DBR <- M_Fed_B3out_filtered@colData@nrows/50000
# run scdbl finder , 2000 features , 30 dims, random porportion 0.2 , 5 iterations
M_Fed_B3_scDbl <- scDblFinder(M_Fed_B3out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = M_Fed_B3_DBR, propRandom = 0.2, iter = 5)
# transform back to seurat object
M_Fed_B3_scDbl <- as.Seurat(M_Fed_B3_scDbl)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fed_B3_scDbl <- SCTransform(M_Fed_B3_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fed_B3_scDbl@assays$SCT@var.features <- M_Fed_B3_scDbl@assays$SCT@var.features[(!M_Fed_B3_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fed_B3_scDbl <- RunPCA(M_Fed_B3_scDbl, features = VariableFeatures(object = M_Fed_B3_scDbl))
M_Fed_B3_scDbl <- RunUMAP(M_Fed_B3_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
M_Fed_B3_scDbl <- FindNeighbors(M_Fed_B3_scDbl, dims = 1:30)
M_Fed_B3_scDbl <- FindClusters(M_Fed_B3_scDbl, resolution = 0.8)


# plot singlets and doublets in UMAP space
DimPlot(M_Fed_B3_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/M_Fed_B3_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)

# save singlets only
M_Fed_B3_scDbl <- subset(M_Fed_B3_scDbl, scDblFinder.class == 'singlet')

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fed_B3_scDbl <- SCTransform(M_Fed_B3_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fed_B3_scDbl@assays$SCT@var.features <- M_Fed_B3_scDbl@assays$SCT@var.features[(!M_Fed_B3_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fed_B3_scDbl <- RunPCA(M_Fed_B3_scDbl, features = VariableFeatures(object = M_Fed_B3_scDbl))
M_Fed_B3_scDbl <- RunUMAP(M_Fed_B3_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 1.2 resolution
M_Fed_B3_scDbl <- FindNeighbors(M_Fed_B3_scDbl, dims = 1:30)
M_Fed_B3_scDbl <- FindClusters(M_Fed_B3_scDbl, resolution = 1.2)
DimPlot(M_Fed_B3_scDbl, label = TRUE) + NoLegend()
#check purity of clusters
# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fed_B3_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fed_B3_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fed_B3_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fed_B3_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes , min pct 0.33, logfc threshold log2(1.33)
M_Fed_B3_markers_post_scDblFnd <- FindAllMarkers(M_Fed_B3_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save M_Fed_B3_scDbl, markers
rm(M_Fed_B3_flt, M_Fed_B3_meta, M_Fed_B3_raw, M_Fed_B3_soup, M_Fed_B3_UMAP, 
   M_Fed_B3out, M_Fed_B3out_filtered, M_Fed_B3_singlets)


# repeat steps for next sample

#M_Fast_B1
#create basic clustering for sample
# read in filtered sample, create seurat object
M_Fast_B1_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B1Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fast_B1_flt <- CreateSeuratObject(M_Fast_B1_flt)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fast_B1_flt <- SCTransform(M_Fast_B1_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fast_B1_flt@assays$SCT@var.features <- M_Fast_B1_flt@assays$SCT@var.features[(!M_Fast_B1_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fast_B1_flt <- RunPCA(M_Fast_B1_flt, features = VariableFeatures(object = M_Fast_B1_flt))
M_Fast_B1_flt <- RunUMAP(M_Fast_B1_flt, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
M_Fast_B1_flt <- FindNeighbors(M_Fast_B1_flt, dims = 1:30)
M_Fast_B1_flt <- FindClusters(M_Fast_B1_flt, resolution = 0.8)

# plot UMAP projections
DimPlot(M_Fast_B1_flt, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fast_B1_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B1_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fast_B1_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B1_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all markers , min pct 0.33, logfc threshold log2(1.33)
M_Fast_B1_markers_b4_soupx <- FindAllMarkers(M_Fast_B1_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# save UMAP, and metadata for future use
M_Fast_B1_UMAP <- M_Fast_B1_flt@reductions$umap@cell.embeddings |> as.data.frame()
M_Fast_B1_meta <- M_Fast_B1_flt@meta.data[,c(2,3,7)] |> cbind(M_Fast_B1_UMAP)



#start SoupX in earnest 
# read in filtered and raw samples, create soupx object
M_Fast_B1_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B1Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fast_B1_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B1Solo.out/GeneFull_Ex50pAS/raw/')
M_Fast_B1_soup <-SoupChannel(M_Fast_B1_raw, M_Fast_B1_flt)

# set dimension reduction to saved UMAP , set clusters to saved clusters
M_Fast_B1_soup = setDR(M_Fast_B1_soup, M_Fast_B1_meta[colnames(M_Fast_B1_flt), c("UMAP_1", "UMAP_2")])
M_Fast_B1_soup = setClusters(M_Fast_B1_soup, setNames(M_Fast_B1_meta$seurat_clusters, colnames(M_Fast_B1_flt)))

# estimate contamination fraction
M_Fast_B1_soup = autoEstCont(M_Fast_B1_soup)
# rho estimated at 0.07; however, 0.175 producing cleaner looking clusters based on Agrp expression 

# set contamination fraction, adjust counts
M_Fast_B1_soup <- setContaminationFraction(M_Fast_B1_soup, 0.175)
M_Fast_B1out <- adjustCounts(M_Fast_B1_soup, method = 'multinomial')

# SoupX finished, check results
# create seurat object
M_Fast_B1out <- CreateSeuratObject(M_Fast_B1out)

# normalize expression , find variable features excluding xy mt rbp hemoglobin genes
M_Fast_B1out <- SCTransform(M_Fast_B1out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fast_B1out@assays$SCT@var.features <- M_Fast_B1out@assays$SCT@var.features[(!M_Fast_B1out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fast_B1out <- RunPCA(M_Fast_B1out, features = VariableFeatures(object = M_Fast_B1out))
M_Fast_B1out <- RunUMAP(M_Fast_B1out, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 1.5 resolution
M_Fast_B1out <- FindNeighbors(M_Fast_B1out, dims = 1:30)
M_Fast_B1out <- FindClusters(M_Fast_B1out, resolution = 1.5)

# plot UMAP projections
DimPlot(M_Fast_B1out, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH genes
FeaturePlot(M_Fast_B1out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B1_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fast_B1out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B1_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes , min pct 0.33, logfc threshold log2(1.33)
M_Fast_B1_markers_post_soupx <- FindAllMarkers(M_Fast_B1out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# calculate and save mt genes expression pct
mitoPercent <- Matrix::colSums(M_Fast_B1out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(M_Fast_B1out@assays$RNA@counts)*100
M_Fast_B1out$mitoPercent <- mitoPercent

# calculate and save hemoglobin genes counts
hemoglobin <- Matrix::colSums(M_Fast_B1out@assays$RNA@counts[hemoglobin_genes, ])
M_Fast_B1out <- AddMetaData(object = M_Fast_B1out, metadata = hemoglobin, col.name = "hemoglobin")


# calculate and save rbp genes expression pct
riboPercent <- Matrix::colSums(M_Fast_B1out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(M_Fast_B1out@assays$RNA@counts)*100
M_Fast_B1out$riboPercent <- riboPercent

# plot QC markers
VlnPlot(M_Fast_B1out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fast_B1_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

# filter out extreme values that may be low quality or possible doublet
M_Fast_B1out_filtered <- subset(M_Fast_B1out, 
                                subset = nCount_RNA < (mean(M_Fast_B1out@meta.data$nCount_RNA) + 4*sd(M_Fast_B1out@meta.data$nCount_RNA))
                                & nFeature_RNA < (mean(M_Fast_B1out@meta.data$nFeature_RNA) + 4*sd(M_Fast_B1out@meta.data$nFeature_RNA))
                                & mitoPercent < (mean(M_Fast_B1out@meta.data$mitoPercent) + 4*sd(M_Fast_B1out@meta.data$mitoPercent))
                                & riboPercent < (mean(M_Fast_B1out@meta.data$riboPercent) + 4*sd(M_Fast_B1out@meta.data$riboPercent))
                                & hemoglobin < (mean(M_Fast_B1out@meta.data$hemoglobin) + 4*sd(M_Fast_B1out@meta.data$hemoglobin)))

# plot updated QC markers
VlnPlot(M_Fast_B1out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fast_B1_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
# set default assay to RNA, save only RNA counts, transform into sc experiment
DefaultAssay(M_Fast_B1out_filtered) <- 'RNA'
M_Fast_B1out_filtered <- DietSeurat(M_Fast_B1out_filtered, assays = 'RNA')
M_Fast_B1out_filtered <- as.SingleCellExperiment(M_Fast_B1out_filtered)

# set doublet rate to number of cells / 50000
M_Fast_B1_DBR <- M_Fast_B1out_filtered@colData@nrows/50000
# run scdbl finder , 2000 features, 30 dims, random porportion 0.2 , 5 iterations
M_Fast_B1_scDbl <- scDblFinder(M_Fast_B1out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = M_Fast_B1_DBR, propRandom = 0.2, iter = 5)
# transform back into seurat object
M_Fast_B1_scDbl <- as.Seurat(M_Fast_B1_scDbl)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fast_B1_scDbl <- SCTransform(M_Fast_B1_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fast_B1_scDbl@assays$SCT@var.features <- M_Fast_B1_scDbl@assays$SCT@var.features[(!M_Fast_B1_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculat PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fast_B1_scDbl <- RunPCA(M_Fast_B1_scDbl, features = VariableFeatures(object = M_Fast_B1_scDbl))
M_Fast_B1_scDbl <- RunUMAP(M_Fast_B1_scDbl, reduction = "pca", dims = 1:30)



#M_Fast_B1_scDbl <- FindNeighbors(M_Fast_B1_scDbl, dims = 1:30)
#M_Fast_B1_scDbl <- FindClusters(M_Fast_B1_scDbl, resolution = 0.8)


# plot singlets and doublets in UMAP space
DimPlot(M_Fast_B1_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/M_Fast_B1_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)

# save only singlets
M_Fast_B1_scDbl <- subset(M_Fast_B1_scDbl, scDblFinder.class == 'singlet')

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fast_B1_scDbl <- SCTransform(M_Fast_B1_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fast_B1_scDbl@assays$SCT@var.features <- M_Fast_B1_scDbl@assays$SCT@var.features[(!M_Fast_B1_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fast_B1_scDbl <- RunPCA(M_Fast_B1_scDbl, features = VariableFeatures(object = M_Fast_B1_scDbl))
M_Fast_B1_scDbl <- RunUMAP(M_Fast_B1_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 1.5 resolution, plot UMAP
M_Fast_B1_scDbl <- FindNeighbors(M_Fast_B1_scDbl, dims = 1:30)
M_Fast_B1_scDbl <- FindClusters(M_Fast_B1_scDbl, resolution = 1.5)
DimPlot(M_Fast_B1_scDbl, label = TRUE) + NoLegend()

#check purity of clusters
# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fast_B1_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B1_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fast_B1_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B1_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes, min pct 0.33, logfc threshold log2(1.33)
M_Fast_B1_markers_post_scDblFnd <- FindAllMarkers(M_Fast_B1_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save M_Fast_B1_scDbl, markers
rm(M_Fast_B1_flt, M_Fast_B1_meta, M_Fast_B1_raw, M_Fast_B1_soup, M_Fast_B1_UMAP, 
   M_Fast_B1out, M_Fast_B1out_filtered, M_Fast_B1_singlets)

# repeat steps for next sample

#M_Fast_B2
#create basic clustering for sample
# read in filtered sample , create seurat object
M_Fast_B2_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B2Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fast_B2_flt <- CreateSeuratObject(M_Fast_B2_flt)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fast_B2_flt <- SCTransform(M_Fast_B2_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fast_B2_flt@assays$SCT@var.features <- M_Fast_B2_flt@assays$SCT@var.features[(!M_Fast_B2_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fast_B2_flt <- RunPCA(M_Fast_B2_flt, features = VariableFeatures(object = M_Fast_B2_flt))
M_Fast_B2_flt <- RunUMAP(M_Fast_B2_flt, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
M_Fast_B2_flt <- FindNeighbors(M_Fast_B2_flt, dims = 1:30)
M_Fast_B2_flt <- FindClusters(M_Fast_B2_flt, resolution = 0.8)

# plot UMAP projections
DimPlot(M_Fast_B2_flt, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fast_B2_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B2_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fast_B2_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B2_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes ,  min pct 0.33 , logfc threshold  log2(1.33)
M_Fast_B2_markers_b4_soupx <- FindAllMarkers(M_Fast_B2_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# save UMAP and metadata for future use
M_Fast_B2_UMAP <- M_Fast_B2_flt@reductions$umap@cell.embeddings |> as.data.frame()
M_Fast_B2_meta <- M_Fast_B2_flt@meta.data[,c(2,3,7)] |> cbind(M_Fast_B2_UMAP)



#start SoupX in earnest 
# read in filtered and raw samples , create soupx object
M_Fast_B2_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B2Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fast_B2_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B2Solo.out/GeneFull_Ex50pAS/raw/')
M_Fast_B2_soup <-SoupChannel(M_Fast_B2_raw, M_Fast_B2_flt)

# set dimension reduction to saved UMAP, set clusters to saved clusters
M_Fast_B2_soup = setDR(M_Fast_B2_soup, M_Fast_B2_meta[colnames(M_Fast_B2_flt), c("UMAP_1", "UMAP_2")])
M_Fast_B2_soup = setClusters(M_Fast_B2_soup, setNames(M_Fast_B2_meta$seurat_clusters, colnames(M_Fast_B2_flt)))

# estimate contamination fraction
M_Fast_B2_soup = autoEstCont(M_Fast_B2_soup)
# rho estimated at 0.05; however, 0.175 producing cleaner looking clusters based on Agrp expression

# set contamination fraction, adjust counts
M_Fast_B2_soup <- setContaminationFraction(M_Fast_B2_soup, 0.175)
M_Fast_B2out <- adjustCounts(M_Fast_B2_soup, method = 'multinomial')

# SoupX finished, check results
# create seurat object
M_Fast_B2out <- CreateSeuratObject(M_Fast_B2out)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fast_B2out <- SCTransform(M_Fast_B2out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fast_B2out@assays$SCT@var.features <- M_Fast_B2out@assays$SCT@var.features[(!M_Fast_B2out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fast_B2out <- RunPCA(M_Fast_B2out, features = VariableFeatures(object = M_Fast_B2out))
M_Fast_B2out <- RunUMAP(M_Fast_B2out, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 1.5 resolution
M_Fast_B2out <- FindNeighbors(M_Fast_B2out, dims = 1:30)
M_Fast_B2out <- FindClusters(M_Fast_B2out, resolution = 1.5)

# plot UMAP proojectioions
DimPlot(M_Fast_B2out, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fast_B2out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B2_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fast_B2out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B2_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes , min pct 0.33, logfc threshold log2(1.33)
M_Fast_B2_markers_post_soupx <- FindAllMarkers(M_Fast_B2out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# calculate and save mt genes expression pct
mitoPercent <- Matrix::colSums(M_Fast_B2out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(M_Fast_B2out@assays$RNA@counts)*100
M_Fast_B2out$mitoPercent <- mitoPercent

# calculate and save hemoglobin genes expression counts
hemoglobin <- Matrix::colSums(M_Fast_B2out@assays$RNA@counts[hemoglobin_genes, ])
M_Fast_B2out <- AddMetaData(object = M_Fast_B2out, metadata = hemoglobin, col.name = "hemoglobin")


# calculate and save rbp genes expression pct
riboPercent <- Matrix::colSums(M_Fast_B2out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(M_Fast_B2out@assays$RNA@counts)*100
M_Fast_B2out$riboPercent <- riboPercent

# plot QC markers
VlnPlot(M_Fast_B2out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fast_B2_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

# filter out extreme values that may be low quality or possible doublet
M_Fast_B2out_filtered <- subset(M_Fast_B2out, 
                                subset = nCount_RNA < (mean(M_Fast_B2out@meta.data$nCount_RNA) + 4*sd(M_Fast_B2out@meta.data$nCount_RNA))
                                & nFeature_RNA < (mean(M_Fast_B2out@meta.data$nFeature_RNA) + 4*sd(M_Fast_B2out@meta.data$nFeature_RNA))
                                & mitoPercent < (mean(M_Fast_B2out@meta.data$mitoPercent) + 4*sd(M_Fast_B2out@meta.data$mitoPercent))
                                & riboPercent < (mean(M_Fast_B2out@meta.data$riboPercent) + 4*sd(M_Fast_B2out@meta.data$riboPercent))
                                & hemoglobin < (mean(M_Fast_B2out@meta.data$hemoglobin) + 4*sd(M_Fast_B2out@meta.data$hemoglobin)))

# plot updated QC markers
VlnPlot(M_Fast_B2out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fast_B2_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
# set default assay to RNA, save only RNA counts, transform to sc experiment
DefaultAssay(M_Fast_B2out_filtered) <- 'RNA'
M_Fast_B2out_filtered <- DietSeurat(M_Fast_B2out_filtered, assays = 'RNA')
M_Fast_B2out_filtered <- as.SingleCellExperiment(M_Fast_B2out_filtered)

# set doublet rate to number of cells / 50000
M_Fast_B2_DBR <- M_Fast_B2out_filtered@colData@nrows/50000
# run scdbl finder, 2000 features, 30 dims, 0.2 random proportion, 5 iterations
M_Fast_B2_scDbl <- scDblFinder(M_Fast_B2out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = M_Fast_B2_DBR, propRandom = 0.2, iter = 5)
# transform back to seurat object
M_Fast_B2_scDbl <- as.Seurat(M_Fast_B2_scDbl)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fast_B2_scDbl <- SCTransform(M_Fast_B2_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fast_B2_scDbl@assays$SCT@var.features <- M_Fast_B2_scDbl@assays$SCT@var.features[(!M_Fast_B2_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fast_B2_scDbl <- RunPCA(M_Fast_B2_scDbl, features = VariableFeatures(object = M_Fast_B2_scDbl))
M_Fast_B2_scDbl <- RunUMAP(M_Fast_B2_scDbl, reduction = "pca", dims = 1:30)



#M_Fast_B2_scDbl <- FindNeighbors(M_Fast_B2_scDbl, dims = 1:30)
#M_Fast_B2_scDbl <- FindClusters(M_Fast_B2_scDbl, resolution = 0.8)


# plot singlets and doublets in UMAP space
DimPlot(M_Fast_B2_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/M_Fast_B2_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)

# save only singlets
M_Fast_B2_scDbl <- subset(M_Fast_B2_scDbl, scDblFinder.class == 'singlet')

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fast_B2_scDbl <- SCTransform(M_Fast_B2_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fast_B2_scDbl@assays$SCT@var.features <- M_Fast_B2_scDbl@assays$SCT@var.features[(!M_Fast_B2_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fast_B2_scDbl <- RunPCA(M_Fast_B2_scDbl, features = VariableFeatures(object = M_Fast_B2_scDbl))
M_Fast_B2_scDbl <- RunUMAP(M_Fast_B2_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 1.5 resolution, plot UMAP
M_Fast_B2_scDbl <- FindNeighbors(M_Fast_B2_scDbl, dims = 1:30)
M_Fast_B2_scDbl <- FindClusters(M_Fast_B2_scDbl, resolution = 1.5)
DimPlot(M_Fast_B2_scDbl, label = TRUE) + NoLegend()
#check purity of clusters

# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fast_B2_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B2_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fast_B2_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B2_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes, min pct 0.33, logfc threshold log2(1.33)
M_Fast_B2_markers_post_scDblFnd <- FindAllMarkers(M_Fast_B2_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save M_Fast_B2_scDbl, markers
rm(M_Fast_B2_flt, M_Fast_B2_meta, M_Fast_B2_raw, M_Fast_B2_soup, M_Fast_B2_UMAP, 
   M_Fast_B2out, M_Fast_B2out_filtered, M_Fast_B2_singlets)

# repeat steps for next sample

#M_Fast_B3
#create basic clustering for sample
# read in filtered sample , create seurat object
M_Fast_B3_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B3Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fast_B3_flt <- CreateSeuratObject(M_Fast_B3_flt)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fast_B3_flt <- SCTransform(M_Fast_B3_flt, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fast_B3_flt@assays$SCT@var.features <- M_Fast_B3_flt@assays$SCT@var.features[(!M_Fast_B3_flt@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fast_B3_flt <- RunPCA(M_Fast_B3_flt, features = VariableFeatures(object = M_Fast_B3_flt))
M_Fast_B3_flt <- RunUMAP(M_Fast_B3_flt, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters at 0.8 resolution
M_Fast_B3_flt <- FindNeighbors(M_Fast_B3_flt, dims = 1:30)
M_Fast_B3_flt <- FindClusters(M_Fast_B3_flt, resolution = 0.8)

# plot UMAP projections
DimPlot(M_Fast_B3_flt, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fast_B3_flt, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B3_ARH_neurons_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fast_B3_flt, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B3_ARH_cells_b4_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes, min pct 0.33 , logfc threshold log2(1.33)
M_Fast_B3_markers_b4_soupx <- FindAllMarkers(M_Fast_B3_flt, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# save UMAP and metadata for later use
M_Fast_B3_UMAP <- M_Fast_B3_flt@reductions$umap@cell.embeddings |> as.data.frame()
M_Fast_B3_meta <- M_Fast_B3_flt@meta.data[,c(2,3,7)] |> cbind(M_Fast_B3_UMAP)



#start SoupX in earnest 
# read in filtered and raw samples , create soupx object
M_Fast_B3_flt <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B3Solo.out/GeneFull_Ex50pAS/filtered/')
M_Fast_B3_raw <- Read10X('../../ARH_Sex_by_Nutr/STARsolo/M_Fast_B3Solo.out/GeneFull_Ex50pAS/raw/')
M_Fast_B3_soup <-SoupChannel(M_Fast_B3_raw, M_Fast_B3_flt)

# set dimension reduction to saved UMAP, set clusters to saved clusters
M_Fast_B3_soup = setDR(M_Fast_B3_soup, M_Fast_B3_meta[colnames(M_Fast_B3_flt), c("UMAP_1", "UMAP_2")])
M_Fast_B3_soup = setClusters(M_Fast_B3_soup, setNames(M_Fast_B3_meta$seurat_clusters, colnames(M_Fast_B3_flt)))

# estimate contamination fraction
M_Fast_B3_soup = autoEstCont(M_Fast_B3_soup)
# rho estimated at 0.08; however, 0.16 producing cleaner looking clusters based of Agrp expression

# set contamination fraction, adjust counts
M_Fast_B3_soup <- setContaminationFraction(M_Fast_B3_soup, 0.16)
M_Fast_B3out <- adjustCounts(M_Fast_B3_soup, method = 'multinomial')

# SoupX finished, check results
# create seurat object
M_Fast_B3out <- CreateSeuratObject(M_Fast_B3out)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fast_B3out <- SCTransform(M_Fast_B3out, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fast_B3out@assays$SCT@var.features <- M_Fast_B3out@assays$SCT@var.features[(!M_Fast_B3out@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fast_B3out <- RunPCA(M_Fast_B3out, features = VariableFeatures(object = M_Fast_B3out))
M_Fast_B3out <- RunUMAP(M_Fast_B3out, reduction = "pca", dims = 1:30)

# find neighbors based on first 30 PCs, find clusters 0.8 resolution
M_Fast_B3out <- FindNeighbors(M_Fast_B3out, dims = 1:30)
M_Fast_B3out <- FindClusters(M_Fast_B3out, resolution = 0.8)

# plot UMAP projections
DimPlot(M_Fast_B3out, label = TRUE, label.size = 2) + NoLegend()

# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fast_B3out, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B3_ARH_neurons_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fast_B3out, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B3_ARH_cells_post_soupx.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all marker genes , min pct 0.33, logfc threshold log2(1.33)
M_Fast_B3_markers_post_soupx <- FindAllMarkers(M_Fast_B3out, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)

# calculate and save mt genes expression pct
mitoPercent <- Matrix::colSums(M_Fast_B3out@assays$RNA@counts[c(mito.genes), ])/Matrix::colSums(M_Fast_B3out@assays$RNA@counts)*100
M_Fast_B3out$mitoPercent <- mitoPercent

# calculate and save hemoglobin genes expression counts
hemoglobin <- Matrix::colSums(M_Fast_B3out@assays$RNA@counts[hemoglobin_genes, ])
M_Fast_B3out <- AddMetaData(object = M_Fast_B3out, metadata = hemoglobin, col.name = "hemoglobin")


# calculate and save rbp genes expression pct
riboPercent <- Matrix::colSums(M_Fast_B3out@assays$RNA@counts[c(ribo.genes), ])/Matrix::colSums(M_Fast_B3out@assays$RNA@counts)*100
M_Fast_B3out$riboPercent <- riboPercent

# plot QC markers
VlnPlot(M_Fast_B3out, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fast_B3_b4_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

# filter out extreme values that may be low quality or possible doublet
M_Fast_B3out_filtered <- subset(M_Fast_B3out, 
                                subset = nCount_RNA < (mean(M_Fast_B3out@meta.data$nCount_RNA) + 4*sd(M_Fast_B3out@meta.data$nCount_RNA))
                                & nFeature_RNA < (mean(M_Fast_B3out@meta.data$nFeature_RNA) + 4*sd(M_Fast_B3out@meta.data$nFeature_RNA))
                                & mitoPercent < (mean(M_Fast_B3out@meta.data$mitoPercent) + 4*sd(M_Fast_B3out@meta.data$mitoPercent))
                                & riboPercent < (mean(M_Fast_B3out@meta.data$riboPercent) + 4*sd(M_Fast_B3out@meta.data$riboPercent))
                                & hemoglobin < (mean(M_Fast_B3out@meta.data$hemoglobin) + 4*sd(M_Fast_B3out@meta.data$hemoglobin)))

# plot updated QC markers
VlnPlot(M_Fast_B3out_filtered, features = c('nCount_RNA', 'nFeature_RNA','riboPercent','hemoglobin','mitoPercent'), group.by = 'orig.ident')
#ggsave('soupx_figures/M_Fast_B3_post_filter_QC.tiff', units = 'in', width = 7, height = 7, dpi = 600)

#scDblFinder

#scDbl
# set default assay to RNA, save only RNA counts, transform into sc experiment
DefaultAssay(M_Fast_B3out_filtered) <- 'RNA'
M_Fast_B3out_filtered <- DietSeurat(M_Fast_B3out_filtered, assays = 'RNA')
M_Fast_B3out_filtered <- as.SingleCellExperiment(M_Fast_B3out_filtered)

# set doublet rate to number of cells / 50000
M_Fast_B3_DBR <- M_Fast_B3out_filtered@colData@nrows/50000
# run scdbl finder, 2000 features, 30 dims , 0.2 random proportion , 5 iterations
M_Fast_B3_scDbl <- scDblFinder(M_Fast_B3out_filtered, clusters = 'seurat_clusters', nfeatures = 2000, dims = 30, dbr = M_Fast_B3_DBR, propRandom = 0.2, iter = 5)
# transform back to seurat object
M_Fast_B3_scDbl <- as.Seurat(M_Fast_B3_scDbl)

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fast_B3_scDbl <- SCTransform(M_Fast_B3_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fast_B3_scDbl@assays$SCT@var.features <- M_Fast_B3_scDbl@assays$SCT@var.features[(!M_Fast_B3_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features, calculate UMAP based on first 30 PCs
M_Fast_B3_scDbl <- RunPCA(M_Fast_B3_scDbl, features = VariableFeatures(object = M_Fast_B3_scDbl))
M_Fast_B3_scDbl <- RunUMAP(M_Fast_B3_scDbl, reduction = "pca", dims = 1:30)



#M_Fast_B3_scDbl <- FindNeighbors(M_Fast_B3_scDbl, dims = 1:30)
#M_Fast_B3_scDbl <- FindClusters(M_Fast_B3_scDbl, resolution = 0.8)


# plot singlets and doublets in UMAP space
DimPlot(M_Fast_B3_scDbl, group.by = 'scDblFinder.class') + NoLegend()
#ggsave('soupx_figures/M_Fast_B3_scDblFinder.tiff', units = 'in', width = 3.5, height = 3.5, dpi = 600)

# save only singlets
M_Fast_B3_scDbl <- subset(M_Fast_B3_scDbl, scDblFinder.class == 'singlet')

# normalize expression, find variable features excluding xy mt rbp hemoglobin genes
M_Fast_B3_scDbl <- SCTransform(M_Fast_B3_scDbl, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')
M_Fast_B3_scDbl@assays$SCT@var.features <- M_Fast_B3_scDbl@assays$SCT@var.features[(!M_Fast_B3_scDbl@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]


# calculate PCs based on variable features , calculate UMAP based on first 30 PCs
M_Fast_B3_scDbl <- RunPCA(M_Fast_B3_scDbl, features = VariableFeatures(object = M_Fast_B3_scDbl))
M_Fast_B3_scDbl <- RunUMAP(M_Fast_B3_scDbl, reduction = "pca", dims = 1:30)


# find neighbors based on first 30 PCs, find clusters at 0.8 resolution, plot UMAP
M_Fast_B3_scDbl <- FindNeighbors(M_Fast_B3_scDbl, dims = 1:30)
M_Fast_B3_scDbl <- FindClusters(M_Fast_B3_scDbl, resolution = 0.8)
DimPlot(M_Fast_B3_scDbl, label = TRUE) + NoLegend()
#check purity of clusters

# plot expression of canonical ARH neuron genes
FeaturePlot(M_Fast_B3_scDbl, features = c('Agrp','Pomc','Tac2','Ghrh'))
#ggsave('soupx_figures/M_Fast_B3_ARH_neurons_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# plot expression of neuron and glia marker genes
FeaturePlot(M_Fast_B3_scDbl, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('soupx_figures/M_Fast_B3_ARH_cells_post_scDblFnd.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

# find all markers , min pct 0.33, logfc threshold log2(1.33)
M_Fast_B3_markers_post_scDblFnd <- FindAllMarkers(M_Fast_B3_scDbl, only.pos = TRUE, min.pct = 0.33, logfc.threshold = 0.415)


# remove unneeded objects, save M_Fast_B3_scDbl, markers
rm(M_Fast_B3_flt, M_Fast_B3_meta, M_Fast_B3_raw, M_Fast_B3_soup, M_Fast_B3_UMAP, 
   M_Fast_B3out, M_Fast_B3out_filtered, M_Fast_B3_singlets)




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


# prep samples

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

# Combine all samples into one object

ARH_Sex_by_Nutr <- merge(F_Fed_B1_scDbl, y= c(F_Fed_B2_scDbl, F_Fed_B3_scDbl, 
                                              F_Fast_B1_scDbl, F_Fast_B2_scDbl, F_Fast_B3_scDbl, 
                                              M_Fed_B1_scDbl, M_Fed_B2_scDbl, M_Fed_B3_scDbl, 
                                              M_Fast_B1_scDbl, M_Fast_B2_scDbl, M_Fast_B3_scDbl), 
                         add.cell.ids = c('ffdb1','ffdb2','ffdb3',
                                          'fftb1','fftb2','fftb3',
                                          'mfdb1','mfdb2','mfdb3',
                                          'mftb1','mftb2','mftb3'), 
                         project = 'ARH_Sex_by_Nutr')

# save object for future scripts
saveRDS(ARH_Sex_by_Nutr, file = 'data/ARH_Sex_by_Nutr.rds')

#remove all objects, libs
rm(cells2HL, cluster_celltype_scores_exon_only, dr10_eo_arh_cell_markers, dr10_eo_arh_cell_markers250pc, dr10_eo_arh_neurons_markers, F_Fast_B1_DBR, F_Fast_B1_markers,
   F_Fast_B1_scDbl, F_Fast_B2_DBR, F_Fast_B1_markers, F_Fast_B2_scDbl, F_Fast_B3_DBR, F_Fast_B3_markers, F_Fast_B3_scDbl, F_Fed_B1_DBR, F_Fed_B1_markers, 
   F_Fed_B1_scDbl, F_Fed_B1_scDbl2, F_Fed_B2_DBR, F_Fed_B2_markers, F_Fed_B2_scDbl, F_Fed_B3_DBR, F_Fed_B3_markers, F_Fed_B3_scDbl, M_Fast_B1_DBR, M_Fast_B1_markers,
   M_Fast_B1_scDbl, hemoglobin, arh_features, libs_min, M_Fast_B2_DBR, M_Fast_B2_markers, M_Fast_B2_scDbl, M_Fast_B3_DBR, M_Fast_B3_scDbl, M_Fed_B1_DBR,
   M_Fed_B1_markers, M_Fed_B1_scDbl, M_Fed_B1_markers, M_Fed_B1_scDbl, M_Fed_B2_DBR, M_Fed_B2_markers, M_Fed_B2_scDbl, M_Fed_B3_DBR, M_Fed_B3_markers, M_Fed_B3_scDbl,
   markers, mito.genes2, mitoPercent, riboPercent, select_arh_features, sn_cell_types, sn_conditions)
# remove any remaining objects



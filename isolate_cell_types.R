# the purpose of this script is to investigate a limited number of cell-types in more detail
# read in required libraries
libs <- c( 'gplots','stringi','reshape2','cowplot','RColorBrewer',
           'sctransform','stringr','org.Mm.eg.db','AnnotationDbi',
           'IRanges','S4Vectors','Biobase','BiocGenerics','clusterProfiler',
           'biomaRt','Matrix','DESeq2','RcppThread', 'extrafont', 'openxlsx',
           'Seurat','dplyr','tidyr','ggplot2','harmony','ggalluvial',
           'scDblFinder','SoupX','UpSetR','ComplexUpset','CellChat','gprofiler2','NMF','ggalluvial')

lapply(libs, require, character.only = TRUE)



# read in seurat object
ARH_Sex_by_Nutr <- readRDS('data/ARH_Sex_by_Nutr.rds')
#Agrp_Sex_by_Nutr <- readRDS('data/Agrp_Sex_by_Nutr.rds')
# subset Agrp cell-type into seurat object
Agrp_Sex_by_Nutr <- subset(ARH_Sex_by_Nutr, subset = cell_type3 == 'Agrp')
# normalize expression and find variable features
Agrp_Sex_by_Nutr <- SCTransform(Agrp_Sex_by_Nutr, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')

# exclude variable features on xy mt rbp and hemoglobin genes
Agrp_Sex_by_Nutr@assays$SCT@var.features <- Agrp_Sex_by_Nutr@assays$SCT@var.features[(!Agrp_Sex_by_Nutr@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]
# calculate PCs based on variable features, calculate TSNE based on first 30 PCs
Agrp_Sex_by_Nutr <- RunPCA(Agrp_Sex_by_Nutr, features = VariableFeatures(object = Agrp_Sex_by_Nutr))
Agrp_Sex_by_Nutr <- RunTSNE(Agrp_Sex_by_Nutr, reduction = "pca", dims = 1:30)
# find neighbors based on first 30 PCs, find clusters at 0.5 resolution
Agrp_Sex_by_Nutr <- FindNeighbors(Agrp_Sex_by_Nutr, dims = 1:30)
Agrp_Sex_by_Nutr <- FindClusters(Agrp_Sex_by_Nutr, resolution = 0.5)
# plot TSNE projections
DimPlot(Agrp_Sex_by_Nutr, label = TRUE, label.size = 8, repel = FALSE, reduction = 'tsne') + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())



# plot TSNE grouped by sexXnutr
DimPlot(Agrp_Sex_by_Nutr, label = F, label.size = 2,pt.size = 1, repel = TRUE, 
        shuffle = TRUE, group.by = 'sexXnutr', reduction = 'tsne',
        cols = c('#f9994f','#6a816a','#ff6e00','#19552b')) &
  theme(plot.title = element_blank(),
        legend.text = element_text(family = 'Arial', size = 8, color = 'black'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
ggsave('figures/agrp_sex_by_nutr_dimplot.tiff', device = 'tiff', units = 'in', width = 4.5,height = 3.5,dpi = 600)

# plot female fed condition
FeaturePlot(Agrp_Sex_by_Nutr, label = FALSE, 
            features = 'F_Fed', reduction = 'tsne',
            pt.size = 1, order = T,
            cols = c('lightgrey','#6a816a')) + 
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/agrp_f_fed_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)

# plot female fasted condition
FeaturePlot(Agrp_Sex_by_Nutr, label = FALSE, 
            features = 'F_Fast', reduction = 'tsne', 
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#f9994f')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/agrp_f_fast_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)
# plot male fed condition
FeaturePlot(Agrp_Sex_by_Nutr, label = FALSE, 
            features = 'M_Fed', reduction = 'tsne',
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#19552b')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/agrp_m_fed_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)

# plot male fasted condition
FeaturePlot(Agrp_Sex_by_Nutr, label = FALSE, 
            features = 'M_Fast', reduction = 'tsne',
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#ff6e00')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/agrp_m_fast_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)


# plot expression of Fos
FeaturePlot(Agrp_Sex_by_Nutr, label = FALSE, 
            features = 'Fos', reduction = 'tsne', max.cutoff = 1,
            pt.size = 1, order = TRUE, cols = c('grey','darkred')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))
ggsave('figures/agrp_fos_dimplot2.tiff', device = 'tiff', units = 'in', width = 3.25,height = 3.25,dpi = 600)


saveRDS(Agrp_Sex_by_Nutr, file = 'data/Agrp_Sex_by_Nutr.rds')

#Agrp.F.Fd.v.Fst$p_val_adj[1] = Agrp.F.Fd.v.Fst$p_val_adj[2]

# plot volcano plot with DE genes in color female fed vs female fasted
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Agrp', group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'F_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Agrp.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Agrp.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Agrp.F.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/agrp_female_fed_vs_fasted_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)

# plot volcano plot with DE genes in color male fed vs male fasted
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Agrp', group.by = 'sexXnutr',ident.1 = 'M_Fed', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0,test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Agrp.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Agrp.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Agrp.M.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 0.9), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/agrp_male_fed_vs_fasted_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)




####

# plot volcano plot with DE genes in color female fed vs male fed
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Agrp', group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'M_Fed', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Agrp.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Agrp.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Agrp.Fd.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/agrp_fed_female_vs_male_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot with DE genes in color female fasted vs male fasted
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Agrp', group.by = 'sexXnutr',ident.1 = 'F_Fast', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Agrp.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Agrp.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Agrp.Fst.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.7)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/agrp_fasted_female_vs_male_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)

library(gprofiler2)

# gene ontology molecular function 
Agrp.F.Fd.Fst.suppress.gost <- gost(query = (Agrp.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Agrp.F.Fd.Fst.induce.gost <- gost(query = (Agrp.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)

# save GO:MF to excel file
write.xlsx(Agrp.F.Fd.Fst.suppress.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Agrp/Agrp_F_Fd_Fst_suppress.xlsx")


# save select GO:MF fasting induced genes in females
write.xlsx(Agrp.F.Fd.Fst.induce.gost$result[c(1,4,8,9,12),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Agrp/Agrp_F_Fd_Fst_induce.xlsx")



# GO:MF for males
Agrp.M.Fd.Fst.suppress.gost <- gost(query = (Agrp.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Agrp.M.Fd.Fst.induce.gost <- gost(query = (Agrp.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)

# save select GO):MF terms for fasting suppressed genes in males
write.xlsx(Agrp.M.Fd.Fst.suppress.gost$result[c(4,5,11,14,15),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Agrp/Agrp_M_Fd_Fst_suppress.xlsx")

# save select GO:MF terms for fasting induced genes in  males to excel file
write.xlsx(Agrp.M.Fd.Fst.induce.gost$result[c(4,5,11,14,15),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Agrp/Agrp_M_Fd_Fst_induce.xlsx")


####
# GO:MF for female fed vs male fed
Agrp.Fd.female.higher.gost <- gost(query = (Agrp.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                   organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Agrp.Fd.male.higher.gost <- gost(query = (Agrp.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                 organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)

# save select GO:MF terms for fed females higher genes to excel file
write.xlsx(Agrp.Fd.female.higher.gost[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-01-06/GOMF_Tables/Agrp_Fd_Female_higher.xlsx")


# save select GO:MF terms for fed males higher to excel files
write.xlsx(Agrp.Fd.male.higher.gost$result[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-01-06/GOMF_Tables/Agrp_Fd_male_higher.xlsx")


#####
# GO:MF for female fasted vs male fasted
Agrp.Fst.female.higher.gost <- gost(query = (Agrp.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Agrp.Fst.male.higher.gost <- gost(query = (Agrp.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)

# save select GO:MF terms for fasted females higher to excel file
write.xlsx(Agrp.Fst.female.higher.gost$result[c(1),c(9,11,3)], 
           file = "../paper_figures/post_2025-01-06/GOMF_Tables/Agrp_Fst_Female_higher.xlsx")


# save select GO:MF terms for fasted males higher to excel file
write.xlsx(Agrp.Fst.male.higher.gost$result[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-01-06/GOMF_Tables/Agrp_Fst_male_higher.xlsx")





###barplot for fos in agrp####

Agrp_Sex_by_Nutr$SCT@data[c('Agrp','Fos'),] %>% 
  t() %>% 
  as.data.frame() %>% 
  cbind(Agrp_Sex_by_Nutr@meta.data[,c(1:5)]) %>% 
  ggplot(aes(sexXnutr, Fos)) +
  stat_summary(geom = 'bar', color = 'black', aes(fill = sexXnutr), width = 0.8) +
  stat_summary(geom = 'errorbar', color = 'black', width = 0.3) +
  scale_fill_manual(values = c('#f9994f','#6a816a','#ff6e00','#19552b')) +
  scale_x_discrete(limits = c('F_Fed', 'F_Fast','M_Fed', 'M_Fast'), labels = c('Fed','Fasted','Fed','Fasted')) +
  labs(y='Fos Expression', title = '', x = '') +
  theme_classic() +
  theme(plot.title = element_text(family = 'Arial', size = 10, color = 'black', hjust = 0.5),
        axis.title=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.y=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.x=element_text(family = 'Arial', size = 8, color = 'black'),
        legend.position = 'none')
ggsave('figures/Fos_in_Agrp_barplot', device = 'tiff', units = 'in', width = 4, height = 2, dpi = 600)


# statistical tests for Fos expression
FindMarkers(Agrp_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'F_Fed', 
            ident.2 = 'F_Fast', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)


FindMarkers(Agrp_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'M_Fed', 
            ident.2 = 'M_Fast', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)




FindMarkers(Agrp_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'F_Fed', 
            ident.2 = 'M_Fed', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)





####




##### KNDy
#KNDy_Sex_by_Nutr <- readRDS('data/KNDy_Sex_by_Nutr.rds')
# subset KNDy neurons into seurat object
KNDy_Sex_by_Nutr <- subset(ARH_Sex_by_Nutr, subset = cell_type3 == 'KNDy')
# normalize expression and find variable features
KNDy_Sex_by_Nutr <- SCTransform(KNDy_Sex_by_Nutr, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')

# exclude variable features on xy mt rbp or hemoglobin genes
KNDy_Sex_by_Nutr@assays$SCT@var.features <- KNDy_Sex_by_Nutr@assays$SCT@var.features[(!KNDy_Sex_by_Nutr@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]
# calculate PCs based on variable features, calculate TSNE based on first 30 PCs
KNDy_Sex_by_Nutr <- RunPCA(KNDy_Sex_by_Nutr, features = VariableFeatures(object = KNDy_Sex_by_Nutr))
KNDy_Sex_by_Nutr <- RunTSNE(KNDy_Sex_by_Nutr, reduction = "pca", dims = 1:30)
# find neighbors based on first 30 PCs, find clusters at 0.5 resolution
KNDy_Sex_by_Nutr <- FindNeighbors(KNDy_Sex_by_Nutr, dims = 1:30)
KNDy_Sex_by_Nutr <- FindClusters(KNDy_Sex_by_Nutr, resolution = 0.5)
# plot TSNE
DimPlot(KNDy_Sex_by_Nutr, label = TRUE, label.size = 8, repel = FALSE, reduction = 'tsne') + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())



# plot TSNE grouped by sexXnutr
DimPlot(KNDy_Sex_by_Nutr, label = F, label.size = 2,pt.size = 1, repel = TRUE, 
        shuffle = TRUE, group.by = 'sexXnutr', reduction = 'tsne',
        cols = c('#f9994f','#6a816a','#ff6e00','#19552b')) &
  theme(plot.title = element_blank(),
        legend.text = element_text(family = 'Arial', size = 8, color = 'black'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
ggsave('figures/KNDy_sex_by_nutr_dimplot.tiff', device = 'tiff', units = 'in', width = 4.5,height = 3.5,dpi = 600)

# plot TSNE for female fed
FeaturePlot(KNDy_Sex_by_Nutr, label = FALSE, 
            features = 'F_Fed', reduction = 'tsne',
            pt.size = 1, order = T,
            cols = c('lightgrey','#6a816a')) + 
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/KNDy_f_fed_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)

# plot TSNE for female fasted
FeaturePlot(KNDy_Sex_by_Nutr, label = FALSE, 
            features = 'F_Fast', reduction = 'tsne', 
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#f9994f')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/KNDy_f_fast_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)
# plot TSNE for male fed
FeaturePlot(KNDy_Sex_by_Nutr, label = FALSE, 
            features = 'M_Fed', reduction = 'tsne',
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#19552b')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/KNDy_m_fed_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)

# plot TSNE for male fasted
FeaturePlot(KNDy_Sex_by_Nutr, label = FALSE, 
            features = 'M_Fast', reduction = 'tsne',
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#ff6e00')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/KNDy_m_fast_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)

# plot Fos expression 
FeaturePlot(KNDy_Sex_by_Nutr, label = FALSE, 
            features = 'Fos', reduction = 'tsne', cols = c('grey', 'darkred'), max.cutoff = 1,
            pt.size = 1, order = TRUE) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/KNDy_fos_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)


###barplot for fos in KNDy####

KNDy_Sex_by_Nutr$SCT@data[c('Tac2','Fos'),] %>% 
  t() %>% 
  as.data.frame() %>% 
  cbind(KNDy_Sex_by_Nutr@meta.data[,c(1:5)]) %>% 
  ggplot(aes(sexXnutr, Fos)) +
  stat_summary(geom = 'bar', color = 'black', aes(fill = sexXnutr), width = 0.8) +
  stat_summary(geom = 'errorbar', color = 'black', width = 0.3) +
  scale_fill_manual(values = c('#f9994f','#6a816a','#ff6e00','#19552b')) +
  scale_x_discrete(limits = c('F_Fed', 'F_Fast','M_Fed', 'M_Fast'), labels = c('Fed','Fasted','Fed','Fasted')) +
  labs(y='Fos Expression', title = '', x = '') +
  theme_classic() +
  theme(plot.title = element_text(family = 'Arial', size = 10, color = 'black', hjust = 0.5),
        axis.title=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.y=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.x=element_text(family = 'Arial', size = 8, color = 'black'),
        legend.position = 'none')
ggsave('figures/Fos_in_KNDy_barplot', device = 'tiff', units = 'in', width = 4, height = 2, dpi = 600)

# statistical tests for Fos expression
FindMarkers(KNDy_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'F_Fed', 
            ident.2 = 'F_Fast', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)


FindMarkers(KNDy_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'M_Fed', 
            ident.2 = 'M_Fast', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)




FindMarkers(KNDy_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'F_Fed', 
            ident.2 = 'M_Fed', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)





####


# save seurat object
saveRDS(KNDy_Sex_by_Nutr, file = 'data/KNDy_Sex_by_Nutr.rds')

# plot volcano plot for female fed vs female fasted with DE genes in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'KNDy', group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'F_Fast',
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(KNDy.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(KNDy.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(KNDy.F.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/KNDy_female_fed_vs_fasted_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot for male fed vs male fasted with DE genes in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'KNDy', group.by = 'sexXnutr',ident.1 = 'M_Fed', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(KNDy.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(KNDy.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(KNDy.M.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > .5), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/KNDy_male_fed_vs_fasted_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)




####

# plot volcano plot for female fed vs male fed with DE genes in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'KNDy', group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'M_Fed', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(KNDy.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(KNDy.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(KNDy.Fd.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,3.01)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/KNDy_fed_female_vs_male_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot for female fasted vs male fasted with DE genes in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'KNDy', group.by = 'sexXnutr',ident.1 = 'F_Fast', ident.2 = 'M_Fast',
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(KNDy.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(KNDy.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(KNDy.Fst.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.8)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/KNDy_fasted_female_vs_male_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)

library(gprofiler2)

# GO:MF for female fed vs female fasted
KNDy.F.Fd.Fst.suppress.gost <- gost(query = (KNDy.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
KNDy.F.Fd.Fst.induce.gost <- gost(query = (KNDy.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)

# save select GO:MF for fasting suppressed genes in females to excel file
write.xlsx(KNDy.F.Fd.Fst.suppress.gost$result[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/KNDy_F_Fd_Fst_suppress.xlsx")




# save select GO:MF terms for fasting induced genes in females to excel file
write.xlsx(KNDy.F.Fd.Fst.induce.gost$result[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/KNDy_F_Fd_Fst_induce.xlsx")







# GO:MF for male fed vs male fasted
KNDy.M.Fd.Fst.suppress.gost <- gost(query = (KNDy.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
KNDy.M.Fd.Fst.induce.gost <- gost(query = (KNDy.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


# save select GO:MF terms for fasting suppressed genes in males to excel file 
write.xlsx(KNDy.M.Fd.Fst.suppress.gost$result[c(3,7,10,11,15),c(9,11,3)], 
           file = "../paper_figures/post_2025-01-06/GOMF_Tables/KNDy/KNDy_M_Fd_Fst_suppress.xlsx")






# save GO:MF terms for fasting induced genes in males to excel file
write.xlsx(KNDy.M.Fd.Fst.induce.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/KNDy_M_Fd_Fst_induce.xlsx")








#####
# GO:MF for female fed vs male fed
KNDy.Fd.female.higher.gost <- gost(query = (KNDy.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
KNDy.Fd.male.higher.gost <- gost(query = (KNDy.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


# save select GO:MF terms for genes higher in fed females to excel file
write.xlsx(KNDy.Fd.female.higher.gost$result[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/KNDy_Fd_Female_higher.xlsx")







# save GO:MF terms for genes higher in fed males to excel file
write.xlsx(KNDy.Fd.male.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/KNDy_Fd_male_higher.xlsx")




#####
# GO:MF for fasted females vs fasted males
KNDy.Fst.female.higher.gost <- gost(query = (KNDy.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                   organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
KNDy.Fst.male.higher.gost <- gost(query = (KNDy.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                 organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


# save GO:MF terms for genes higher in fasted females to excel file
write.xlsx(KNDy.Fst.female.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/KNDy_Fst_Female_higher.xlsx")




# save GO:MF terms for genes higher in fasted males
write.xlsx(KNDy.Fst.male.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/KNDy_Fst_male_higher.xlsx")






##### DA
#DA_Sex_by_Nutr <- readRDS('data/DA_Sex_by_Nutr.rds')
# subset DA neurons to seurat object
DA_Sex_by_Nutr <- subset(ARH_Sex_by_Nutr, subset = cell_type3 == 'DA')
# normalize expression and find variable features
DA_Sex_by_Nutr <- SCTransform(DA_Sex_by_Nutr, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')

# exclude variable features on xy mt rbp and hemoglobin genes
DA_Sex_by_Nutr@assays$SCT@var.features <- DA_Sex_by_Nutr@assays$SCT@var.features[(!DA_Sex_by_Nutr@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]
# calculate PCs based on variable features, calculate TSNE based on first 30 PCs
DA_Sex_by_Nutr <- RunPCA(DA_Sex_by_Nutr, features = VariableFeatures(object = DA_Sex_by_Nutr))
DA_Sex_by_Nutr <- RunTSNE(DA_Sex_by_Nutr, reduction = "pca", dims = 1:30)
# find neighbors based on first 30 PCs, find clusters at 0.5 resolution
DA_Sex_by_Nutr <- FindNeighbors(DA_Sex_by_Nutr, dims = 1:30)
DA_Sex_by_Nutr <- FindClusters(DA_Sex_by_Nutr, resolution = 0.5)
# plot TSNE
DimPlot(DA_Sex_by_Nutr, label = TRUE, label.size = 8, repel = FALSE, reduction = 'tsne') + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())



# plot TSNE grouped by sexXnutr
DimPlot(DA_Sex_by_Nutr, label = F, label.size = 2,pt.size = 1, repel = TRUE, 
        shuffle = TRUE, group.by = 'sexXnutr', reduction = 'tsne',
        cols = c('#f9994f','#6a816a','#ff6e00','#19552b')) &
  theme(plot.title = element_blank(),
        legend.text = element_text(family = 'Arial', size = 8, color = 'black'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
ggsave('figures/DA_sex_by_nutr_dimplot.tiff', device = 'tiff', units = 'in', width = 4.5,height = 3.5,dpi = 600)

# plot TSNE of female fed
FeaturePlot(DA_Sex_by_Nutr, label = FALSE, 
            features = 'F_Fed', reduction = 'tsne',
            pt.size = 1, order = T,
            cols = c('lightgrey','#6a816a')) + 
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/DA_f_fed_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)

# plot TSNE of female fasted
FeaturePlot(DA_Sex_by_Nutr, label = FALSE, 
            features = 'F_Fast', reduction = 'tsne', 
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#f9994f')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/DA_f_fast_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)
# plot TSNE of male fed
FeaturePlot(DA_Sex_by_Nutr, label = FALSE, 
            features = 'M_Fed', reduction = 'tsne',
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#19552b')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/DA_m_fed_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)

# plot TSNE of male fasted
FeaturePlot(DA_Sex_by_Nutr, label = FALSE, 
            features = 'M_Fast', reduction = 'tsne',
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#ff6e00')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/DA_m_fast_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)



# plot Fos expression
FeaturePlot(DA_Sex_by_Nutr, label = FALSE, 
            features = 'Fos', reduction = 'tsne', cols = c('grey','darkred'), max.cutoff = 1,
            pt.size = 1, order = TRUE) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/DA_fos_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)


###barplot for fos in DA ####

DA_Sex_by_Nutr$SCT@data[c('Th','Fos'),] %>% 
  t() %>% 
  as.data.frame() %>% 
  cbind(DA_Sex_by_Nutr@meta.data[,c(1:5)]) %>% 
  ggplot(aes(sexXnutr, Fos)) +
  stat_summary(geom = 'bar', color = 'black', aes(fill = sexXnutr), width = 0.8) +
  stat_summary(geom = 'errorbar', color = 'black', width = 0.3) +
  scale_fill_manual(values = c('#f9994f','#6a816a','#ff6e00','#19552b')) +
  scale_x_discrete(limits = c('F_Fed', 'F_Fast','M_Fed', 'M_Fast'), labels = c('Fed','Fasted','Fed','Fasted')) +
  labs(y='Fos Expression', title = '', x = '') +
  theme_classic() +
  theme(plot.title = element_text(family = 'Arial', size = 10, color = 'black', hjust = 0.5),
        axis.title=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.y=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.x=element_text(family = 'Arial', size = 8, color = 'black'),
        legend.position = 'none')
ggsave('figures/Fos_in_DA_barplot', device = 'tiff', units = 'in', width = 4, height = 2, dpi = 600)

# statistical tests for Fos expression
FindMarkers(DA_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'F_Fed', 
            ident.2 = 'F_Fast', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)


FindMarkers(DA_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'M_Fed', 
            ident.2 = 'M_Fast', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)




FindMarkers(DA_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'F_Fed', 
            ident.2 = 'M_Fed', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)





####





# save seurat object
saveRDS(DA_Sex_by_Nutr, file = 'data/DA_Sex_by_Nutr.rds')

# plot volcano plot for female fed vs female fasted with DE genes in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'DA', group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'F_Fast',
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(DA.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(DA.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(DA.F.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > .7), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/DA_female_fed_vs_fasted_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot for male fed vs male fasted with DE genes in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'DA', group.by = 'sexXnutr',ident.1 = 'M_Fed', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(DA.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(DA.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(DA.M.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/DA_male_fed_vs_fasted_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)




####

# plot volcano plot for female fed vs male fed with DE genes in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'DA', group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'M_Fed', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(DA.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(DA.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(DA.Fd.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.6)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/DA_fed_female_vs_male_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot for female fasted vs male fasted with DE genes in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'DA', group.by = 'sexXnutr',ident.1 = 'F_Fast', ident.2 = 'M_Fast',
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(DA.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(DA.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(DA.Fst.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.6)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/DA_fasted_female_vs_male_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)

library(gprofiler2)

# GO:MF for female fed vs female fasted
DA.F.Fd.Fst.suppress.gost <- gost(query = (DA.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
DA.F.Fd.Fst.induce.gost <- gost(query = (DA.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


# save GO:MF for fasting suppressed genes in females to excel file
write.xlsx(DA.F.Fd.Fst.suppress.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/DA_F_Fd_Fst_suppress.xlsx")





# save GO:MF for fasting induced genes in females
write.xlsx(DA.F.Fd.Fst.induce.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/DA_F_Fd_Fst_suppress.xlsx")






# GO:MF for male fed vs male fasted
DA.M.Fd.Fst.suppress.gost <- gost(query = (DA.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
DA.M.Fd.Fst.induce.gost <- gost(query = (DA.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


# no results to save or plot





#####
# GO:MF for female fed vs male fed
DA.Fd.female.higher.gost <- gost(query = (DA.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                   organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
DA.Fd.male.higher.gost <- gost(query = (DA.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                 organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


# save GO:MF for genes higher in females to excel file
write.xlsx(DA.Fd.female.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/DA_Fd_Female_higher.xlsx")





# no results to plot or save for higher in fed males


#####
# GO:MF for female fasted vs male fasted
DA.Fst.female.higher.gost <- gost(query = (DA.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
DA.Fst.male.higher.gost <- gost(query = (DA.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


# save GO:MF for genes higher in fasted females to excel file
write.xlsx(DA.Fst.female.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/DA_Fst_Female_higher.xlsx")




# save GO:MF for genes higher in fasted males to excel file
write.xlsx(DA.Fst.male.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-01-06/GOMF_Tables/KNDy/DA_Fst_Male_higher.xlsx")






####

###microglia####

ARH_Sex_by_Nutr <- readRDS('data/ARH_Sex_by_Nutr.rds')
#Microglia_Sex_by_Nutr <- readRDS('data/Microglia_Sex_by_Nutr.rds')
# subset microglia to seurat object
Microglia_Sex_by_Nutr <- subset(ARH_Sex_by_Nutr, subset = cell_type3 == 'Microglia')
# normalize expression and find variable features
Microglia_Sex_by_Nutr <- SCTransform(Microglia_Sex_by_Nutr, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')

# exclude variable features on xy mt rbp and hemoglobin genes
Microglia_Sex_by_Nutr@assays$SCT@var.features <- Microglia_Sex_by_Nutr@assays$SCT@var.features[(!Microglia_Sex_by_Nutr@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]
# calculate PCs based on variable features, calculate TSNE based on first 30 PCs
Microglia_Sex_by_Nutr <- RunPCA(Microglia_Sex_by_Nutr, features = VariableFeatures(object = Microglia_Sex_by_Nutr))
Microglia_Sex_by_Nutr <- RunTSNE(Microglia_Sex_by_Nutr, reduction = "pca", dims = 1:30)
# find neighbors based on first 30 PCs, find clusters at 0.5 resolution
Microglia_Sex_by_Nutr <- FindNeighbors(Microglia_Sex_by_Nutr, dims = 1:30)
Microglia_Sex_by_Nutr <- FindClusters(Microglia_Sex_by_Nutr, resolution = 0.5)
# plot TSNE
DimPlot(Microglia_Sex_by_Nutr, label = TRUE, label.size = 8, repel = FALSE, reduction = 'tsne') + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())



# plot TSNE grouped by sexXnutr
DimPlot(Microglia_Sex_by_Nutr, label = F, label.size = 2,pt.size = 1, repel = TRUE, 
        shuffle = TRUE, group.by = 'sexXnutr', reduction = 'tsne',
        cols = c('#f9994f','#6a816a','#ff6e00','#19552b')) &
  theme(plot.title = element_blank(),
        legend.text = element_text(family = 'Arial', size = 8, color = 'black'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
ggsave('figures/Microglia_sex_by_nutr_dimplot.tiff', device = 'tiff', units = 'in', width = 4.5,height = 3.5,dpi = 600)

# plot TSNE for female fed
FeaturePlot(Microglia_Sex_by_Nutr, label = FALSE, 
            features = 'F_Fed', reduction = 'tsne',
            pt.size = 1, order = T,
            cols = c('lightgrey','#6a816a')) + 
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/Microglia_f_fed_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)

# plot TSNE for female fasted
FeaturePlot(Microglia_Sex_by_Nutr, label = FALSE, 
            features = 'F_Fast', reduction = 'tsne', 
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#f9994f')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/Microglia_f_fast_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)
# plot TSNE for male fed
FeaturePlot(Microglia_Sex_by_Nutr, label = FALSE, 
            features = 'M_Fed', reduction = 'tsne',
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#19552b')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/Microglia_m_fed_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)

# plot TSNE for male fasted
FeaturePlot(Microglia_Sex_by_Nutr, label = FALSE, 
            features = 'M_Fast', reduction = 'tsne',
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#ff6e00')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/Microglia_m_fast_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)



# plot Fos expression
FeaturePlot(Microglia_Sex_by_Nutr, label = FALSE, 
            features = 'Fos', reduction = 'tsne', cols = c('grey','darkred'), max.cutoff = 1,
            pt.size = 1, order = TRUE) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))

ggsave('figures/Microglia_fos_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)



###barplot for fos in Microglia####

Microglia_Sex_by_Nutr$SCT@data[c('Th','Fos'),] %>% 
  t() %>% 
  as.data.frame() %>% 
  cbind(Microglia_Sex_by_Nutr@meta.data[,c(1:5)]) %>% 
  ggplot(aes(sexXnutr, Fos)) +
  stat_summary(geom = 'bar', color = 'black', aes(fill = sexXnutr), width = 0.8) +
  stat_summary(geom = 'errorbar', color = 'black', width = 0.3) +
  scale_fill_manual(values = c('#f9994f','#6a816a','#ff6e00','#19552b')) +
  scale_x_discrete(limits = c('F_Fed', 'F_Fast','M_Fed', 'M_Fast'), labels = c('Fed','Fasted','Fed','Fasted')) +
  labs(y='Fos Expression', title = '', x = '') +
  theme_classic() +
  theme(plot.title = element_text(family = 'Arial', size = 10, color = 'black', hjust = 0.5),
        axis.title=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.y=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.x=element_text(family = 'Arial', size = 8, color = 'black'),
        legend.position = 'none')
ggsave('figures/Fos_in_Microglia_barplot', device = 'tiff', units = 'in', width = 4, height = 2, dpi = 600)

# statistical tests for Fos expression
FindMarkers(Microglia_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'F_Fed', 
            ident.2 = 'F_Fast', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)


FindMarkers(Microglia_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'M_Fed', 
            ident.2 = 'M_Fast', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)




FindMarkers(Microglia_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'F_Fed', 
            ident.2 = 'M_Fed', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)





####




# save seurate object
saveRDS(Microglia_Sex_by_Nutr, file = 'data/Microglia_Sex_by_Nutr.rds')

# plot volcano plot for female fed vs female fasted with DE genes plotted in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Microglia', group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'F_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Microglia.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Microglia.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Microglia.F.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 0.7), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Microglia_female_fed_vs_fasted_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)



# plot volcano plot for male fed vs male fasted with DE genes plotted in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Microglia', group.by = 'sexXnutr',ident.1 = 'M_Fed', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Microglia.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Microglia.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Microglia.M.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 0.7), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Microglia_male_fed_vs_fasted_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)




####


# plot volcano plot for female fed vs male fed with DE genes plotted in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Microglia', group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'M_Fed', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Microglia.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Microglia.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Microglia.Fd.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,3.13)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Microglia_fed_female_vs_male_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)



# plot volcano plot for female fasted vs male fasted with DE genes plotted in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Microglia', group.by = 'sexXnutr',ident.1 = 'F_Fast', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Microglia.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Microglia.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Microglia.Fst.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,3.1)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Microglia_fasted_female_vs_male_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)

#library(gprofiler2)

# GO:MF for DE genes between female fed vs femlae fasted
Microglia.F.Fd.Fst.suppress.gost <- gost(query = (Microglia.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Microglia.F.Fd.Fst.induce.gost <- gost(query = (Microglia.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)

# no GO:MF results to save for fasting suppressed genes

# save GO:MF for fasting induced genes in females
write.xlsx(Microglia.F.Fd.Fst.induce.gost$result[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Microglia/Microglia_F_Fd_Fst_induce.xlsx")




# GO:MF for male fed vs male fasted 
Microglia.M.Fd.Fst.suppress.gost <- gost(query = (Microglia.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Microglia.M.Fd.Fst.induce.gost <- gost(query = (Microglia.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)



# no GO:MF results for fasting suppressed genes

# save GO:MF for fasting induced genes in males
write.xlsx(Microglia.M.Fd.Fst.induce.gost$result[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Microglia/Microglia_M_Fd_Fst_induce.xlsx")





####
# GO:MF for female fed vs male fed
Microglia.Fd.female.higher.gost <- gost(query = (Microglia.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                   organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Microglia.Fd.male.higher.gost <- gost(query = (Microglia.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                 organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)

# no results to save





#####
# GO:MF for female fasted vs male fasted
Microglia.Fst.female.higher.gost <- gost(query = (Microglia.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Microglia.Fst.male.higher.gost <- gost(query = (Microglia.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


# no GO:MF results for higher in fasting females


# save GO:MF results for higher in fasting males to excel file
write.xlsx(Microglia.Fst.male.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Microglia/Microglia_Fst_Male_higher.xlsx")








###Oligo####


ARH_Sex_by_Nutr <- readRDS('data/ARH_Sex_by_Nutr.rds')
#Oligo_Sex_by_Nutr <- readRDS('data/Oligo_Sex_by_Nutr.rds')
 # subset oligodendrocytes to seurat object
Oligo_Sex_by_Nutr <- subset(ARH_Sex_by_Nutr, subset = cell_type3 == 'Oligodendrocytes')
# normalize expression and find variable features
Oligo_Sex_by_Nutr <- SCTransform(Oligo_Sex_by_Nutr, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')

# exclude variable features on xy mt rbp and hemoglobin genes
Oligo_Sex_by_Nutr@assays$SCT@var.features <- Oligo_Sex_by_Nutr@assays$SCT@var.features[(!Oligo_Sex_by_Nutr@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]
# calculate PCs based on variable features, calculate TSNE based on first 30 PCs
Oligo_Sex_by_Nutr <- RunPCA(Oligo_Sex_by_Nutr, features = VariableFeatures(object = Oligo_Sex_by_Nutr))
Oligo_Sex_by_Nutr <- RunTSNE(Oligo_Sex_by_Nutr, reduction = "pca", dims = 1:30)
# find neighbors based on first 30 PCs, find clusters at 0.5 resolution
Oligo_Sex_by_Nutr <- FindNeighbors(Oligo_Sex_by_Nutr, dims = 1:30)
Oligo_Sex_by_Nutr <- FindClusters(Oligo_Sex_by_Nutr, resolution = 0.5)
# plot TSNE
DimPlot(Oligo_Sex_by_Nutr, label = TRUE, label.size = 8, repel = FALSE, reduction = 'tsne') + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())



# plot TSNE grouped by sexXnutr
DimPlot(Oligo_Sex_by_Nutr, label = F, label.size = 2,pt.size = 1, repel = TRUE, 
        shuffle = TRUE, group.by = 'sexXnutr', reduction = 'tsne',
        cols = c('#f9994f','#6a816a','#ff6e00','#19552b')) &
  theme(plot.title = element_blank(),
        legend.text = element_text(family = 'Arial', size = 8, color = 'black'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
ggsave('figures/Oligo_sex_by_nutr_dimplot.tiff', device = 'tiff', units = 'in', width = 4.5,height = 3.5,dpi = 600)

# plot TSNE for female fed
FeaturePlot(Oligo_Sex_by_Nutr, label = FALSE, 
            features = 'F_Fed', reduction = 'tsne',
            pt.size = 1, order = T,
            cols = c('lightgrey','#6a816a')) + 
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/Oligo_f_fed_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)

# plot TSNE for female fasted
FeaturePlot(Oligo_Sex_by_Nutr, label = FALSE, 
            features = 'F_Fast', reduction = 'tsne', 
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#f9994f')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/Oligo_f_fast_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)
# plot TSNE for male fed
FeaturePlot(Oligo_Sex_by_Nutr, label = FALSE, 
            features = 'M_Fed', reduction = 'tsne',
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#19552b')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/Oligo_m_fed_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)

# plot TSNE for male fasted
FeaturePlot(Oligo_Sex_by_Nutr, label = FALSE, 
            features = 'M_Fast', reduction = 'tsne',
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#ff6e00')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/Oligo_m_fast_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)


# plot Fos expression
FeaturePlot(Oligo_Sex_by_Nutr, label = FALSE, 
            features = 'Fos', reduction = 'tsne', cols = c('grey','darkred'), max.cutoff = 1,
            pt.size = 1, order = TRUE) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/Oligo_fos_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)



###barplot for fos in Microglia####

Oligo_Sex_by_Nutr$SCT@data[c('Th','Fos'),] %>% 
  t() %>% 
  as.data.frame() %>% 
  cbind(Oligo_Sex_by_Nutr@meta.data[,c(1:5)]) %>% 
  ggplot(aes(sexXnutr, Fos)) +
  stat_summary(geom = 'bar', color = 'black', aes(fill = sexXnutr), width = 0.8) +
  stat_summary(geom = 'errorbar', color = 'black', width = 0.3) +
  scale_fill_manual(values = c('#f9994f','#6a816a','#ff6e00','#19552b')) +
  scale_x_discrete(limits = c('F_Fed', 'F_Fast','M_Fed', 'M_Fast'), labels = c('Fed','Fasted','Fed','Fasted')) +
  labs(y='Fos Expression', title = '', x = '') +
  theme_classic() +
  theme(plot.title = element_text(family = 'Arial', size = 10, color = 'black', hjust = 0.5),
        axis.title=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.y=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.x=element_text(family = 'Arial', size = 8, color = 'black'),
        legend.position = 'none')
ggsave('figures/Fos_in_Oligo_barplot', device = 'tiff', units = 'in', width = 4, height = 2, dpi = 600)

# statistical tests for Fos expression
FindMarkers(Oligo_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'F_Fed', 
            ident.2 = 'F_Fast', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)


FindMarkers(Oligo_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'M_Fed', 
            ident.2 = 'M_Fast', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)




FindMarkers(Oligo_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'F_Fed', 
            ident.2 = 'M_Fed', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)





####


# save seurat object to RDS
saveRDS(Oligo_Sex_by_Nutr, file = 'data/Oligo_Sex_by_Nutr.rds')

# plot volcano plot for female fed vs females fasted with DE genes plotted in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Oligodendrocytes', group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'F_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Oligo.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Oligo.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Oligo.F.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Oligo_female_fed_vs_fasted_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot for male fed vs males fasted with DE genes plotted in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Oligodendrocytes', group.by = 'sexXnutr',ident.1 = 'M_Fed', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Oligo.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Oligo.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Oligo.M.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Oligo_male_fed_vs_fasted_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)




####

# plot volcano plot for female fed vs males fed with DE genes plotted in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Oligodendrocytes', group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'M_Fed', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Oligo.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Oligo.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Oligo.Fd.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Oligo_fed_female_vs_male_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot for female fasted vs males fasted with DE genes plotted in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Oligodendrocytes', group.by = 'sexXnutr',ident.1 = 'F_Fast', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Oligo.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Oligo.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Oligo.Fst.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Oligo_fasted_female_vs_male_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)

#library(gprofiler2)

# GO:MF for female fed vs female fasted
Oligo.F.Fd.Fst.suppress.gost <- gost(query = (Oligo.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                         organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Oligo.F.Fd.Fst.induce.gost <- gost(query = (Oligo.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                       organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)

# no results to save

#write.xlsx(Oligo.F.Fd.Fst.suppress.gost$result[,c(9,11,3)], 
 #          file = "../paper_figures/post_2025-01-06/GOMF_Tables/Oligo/Oligo_F_Fd_Fst_suppress.xlsx")


#write.xlsx(Oligo.F.Fd.Fst.induce.gost$result[c(1:3,5,6),c(9,11,3)], 
 #          file = "../paper_figures/post_2025-01-06/GOMF_Tables/Oligo/Oligo_F_Fd_Fst_induce.xlsx")




# GO:MF for male fed vs male fasted

Oligo.M.Fd.Fst.suppress.gost <- gost(query = (Oligo.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                         organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Oligo.M.Fd.Fst.induce.gost <- gost(query = (Oligo.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                       organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


# no results for fasting suppressed
#write.xlsx(Oligo.M.Fd.Fst.suppress.gost$result[c(2,3,6,7,8),c(9,11,3)], 
 #          file = "../paper_figures/post_2025-01-06/GOMF_Tables/Oligo/Oligo_M_Fd_Fst_suppress.xlsx")




# save GO:MF results for fasting induced genes in males to excel file
write.xlsx(Oligo.M.Fd.Fst.induce.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Oligo/Oligo_M_Fd_Fst_induce.xlsx")




####
# GO:MF for female fed vs male fed
Oligo.Fd.female.higher.gost <- gost(query = (Oligo.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                        organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)

Oligo.Fd.male.higher.gost <- gost(query = (Oligo.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                      organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)

# no results to save





#####
# GO:MF for female fasting vs male fasting
Oligo.Fst.female.higher.gost <- gost(query = (Oligo.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                         organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Oligo.Fst.male.higher.gost <- gost(query = (Oligo.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                       organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


# save GO:MF results for higher in fasting females to excel file
write.xlsx(Oligo.Fst.female.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Oligo/Oligo_Fst_Female_higher.xlsx")



# no results for higher in fasting males




####POMC####


# subset POMC to seurat object
Pomc_Sex_by_Nutr <- subset(ARH_Sex_by_Nutr, subset = cell_type3 == 'Pomc')
# normalize expressiong and find varibale features
Pomc_Sex_by_Nutr <- SCTransform(Pomc_Sex_by_Nutr, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')

# exclude variable features on xy mt rbp and hemoglobin genes
Pomc_Sex_by_Nutr@assays$SCT@var.features <- Pomc_Sex_by_Nutr@assays$SCT@var.features[(!Pomc_Sex_by_Nutr@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]
# calculate PCs based on variable features, calculate TSNE based on first 30 PCs
Pomc_Sex_by_Nutr <- RunPCA(Pomc_Sex_by_Nutr, features = VariableFeatures(object = Pomc_Sex_by_Nutr))
Pomc_Sex_by_Nutr <- RunTSNE(Pomc_Sex_by_Nutr, reduction = "pca", dims = 1:30)
# find neighbors based on first 30 PCs, find clusters at 0.5 resolution
Pomc_Sex_by_Nutr <- FindNeighbors(Pomc_Sex_by_Nutr, dims = 1:30)
Pomc_Sex_by_Nutr <- FindClusters(Pomc_Sex_by_Nutr, resolution = 0.5)
# plot TSNE
DimPlot(Pomc_Sex_by_Nutr, label = TRUE, label.size = 8, repel = FALSE, reduction = 'tsne') + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())



# plot TSNE grouped by sexXnutr
DimPlot(Pomc_Sex_by_Nutr, label = F, label.size = 2,pt.size = 1, repel = TRUE, 
        shuffle = TRUE, group.by = 'sexXnutr', reduction = 'tsne',
        cols = c('#f9994f','#6a816a','#ff6e00','#19552b')) &
  theme(plot.title = element_blank(),
        legend.text = element_text(family = 'Arial', size = 8, color = 'black'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
ggsave('figures/Pomc_sex_by_nutr_dimplot.tiff', device = 'tiff', units = 'in', width = 4.5,height = 3.5,dpi = 600)

# plot TSNE for female fed
FeaturePlot(Pomc_Sex_by_Nutr, label = FALSE, 
            features = 'F_Fed', reduction = 'tsne',
            pt.size = 1, order = T,
            cols = c('lightgrey','#6a816a')) + 
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/Pomc_f_fed_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)

# plot TSNE for female fasted
FeaturePlot(Pomc_Sex_by_Nutr, label = FALSE, 
            features = 'F_Fast', reduction = 'tsne', 
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#f9994f')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/Pomc_f_fast_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)
# plot TSNE for male fed
FeaturePlot(Pomc_Sex_by_Nutr, label = FALSE, 
            features = 'M_Fed', reduction = 'tsne',
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#19552b')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/Pomc_m_fed_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)

# plot TSNE for male fasted
FeaturePlot(Pomc_Sex_by_Nutr, label = FALSE, 
            features = 'M_Fast', reduction = 'tsne',
            pt.size = 1, order = TRUE,
            cols = c('lightgrey','#ff6e00')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


ggsave('figures/Pomc_m_fast_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)


# plot Fos expression
FeaturePlot(Pomc_Sex_by_Nutr, label = FALSE, 
            features = 'Fos', reduction = 'tsne', max.cutoff = 1,
            pt.size = 1, order = TRUE, cols = c('grey','darkred')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))
ggsave('figures/Pomc_fos_dimplot.tiff', device = 'tiff', units = 'in', width = 3.25,height = 3.25,dpi = 600)

# save seurat object to RDS
saveRDS(Pomc_Sex_by_Nutr, file = 'data/Pomc_Sex_by_Nutr.rds')



# plot volcano plot for female fed vs female fasted with DE genes plotted in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Pomc', group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'F_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Pomc.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Pomc.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Pomc.F.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 0.7), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Pomc_female_fed_vs_fasted_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)

# plot volcano plot for male fed vs male fasted with DE genes plotted in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Pomc', group.by = 'sexXnutr',ident.1 = 'M_Fed', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0,test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Pomc.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Pomc.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Pomc.M.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 0.7), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Pomc_male_fed_vs_fasted_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)




####

# plot volcano plot for female fed vs male fed with DE genes plotted in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Pomc', group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'M_Fed', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Pomc.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Pomc.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Pomc.Fd.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Pomc_fed_female_vs_male_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot for female fasted vs male fasted with DE genes plotted in color
FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Pomc', group.by = 'sexXnutr',ident.1 = 'F_Fast', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Pomc.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Pomc.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Pomc.Fst.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Pomc_fasted_female_vs_male_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)



###barplot for fos in POMC####

Pomc_Sex_by_Nutr$SCT@data[c('Pomc','Fos'),] %>% 
  t() %>% 
  as.data.frame() %>% 
  cbind(Pomc_Sex_by_Nutr@meta.data[,c(1:5)]) %>% 
  ggplot(aes(sexXnutr, Fos)) +
  stat_summary(geom = 'bar', color = 'black', aes(fill = sexXnutr), width = 0.8) +
  stat_summary(geom = 'errorbar', color = 'black', width = 0.3) +
  scale_fill_manual(values = c('#f9994f','#6a816a','#ff6e00','#19552b')) +
  scale_x_discrete(limits = c('F_Fed', 'F_Fast','M_Fed', 'M_Fast'), labels = c('Fed','Fasted','Fed','Fasted')) +
  labs(y='Fos Expression', title = '', x = '') +
  theme_classic() +
  theme(plot.title = element_text(family = 'Arial', size = 10, color = 'black', hjust = 0.5),
        axis.title=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.y=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.x=element_text(family = 'Arial', size = 8, color = 'black'),
        legend.position = 'none')
ggsave('figures/Fos_in_Pomc_barplot', device = 'tiff', units = 'in', width = 4, height = 2, dpi = 600)

# statistical tests for Fos expression
FindMarkers(Pomc_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'F_Fed', 
            ident.2 = 'F_Fast', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)


FindMarkers(Pomc_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'M_Fed', 
            ident.2 = 'M_Fast', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)




FindMarkers(Pomc_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'F_Fed', 
            ident.2 = 'M_Fed', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)





####

# GO:MF for female fed vs females fasted
Pomc.F.Fd.Fst.suppress.gost <- gost(query = (Pomc.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Pomc.F.Fd.Fst.induce.gost <- gost(query = (Pomc.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)

# save GO:MF results to excel files
write.xlsx(Pomc.F.Fd.Fst.suppress.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Pomc/Pomc_F_Fd_Fst_suppress.xlsx")
write.xlsx(Pomc.F.Fd.Fst.induce.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Pomc/Pomc_F_Fd_Fst_induce.xlsx")


# GO:MF for male fed vs male fasted
Pomc.M.Fd.Fst.suppress.gost <- gost(query = (Pomc.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Pomc.M.Fd.Fst.induce.gost <- gost(query = (Pomc.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)

# save GO:MF results to excel files
write.xlsx(Pomc.M.Fd.Fst.suppress.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Pomc/Pomc_M_Fd_Fst_suppress.xlsx")
write.xlsx(Pomc.M.Fd.Fst.induce.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Pomc/Pomc_M_Fd_Fst_induce.xlsx")

# GO:MF fir female fed vs male fed
Pomc.Fd.female.higher.gost <- gost(query = (Pomc.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                   organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Pomc.Fd.male.higher.gost <- gost(query = (Pomc.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                 organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)

# save GO:MF results to excel files
write.xlsx(Pomc.Fd.female.higher.gost$result[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Pomc/Pomc_Fd_Female_higher.xlsx")

write.xlsx(Pomc.Fd.male.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Pomc/Pomc_Fd_Male_higher.xlsx")

# GO:MF for female fasted vs male fasted
Pomc.Fst.female.higher.gost <- gost(query = (Pomc.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                   organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Pomc.Fst.male.higher.gost <- gost(query = (Pomc.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                 organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)

# save GO:MF results to excel files
write.xlsx(Pomc.Fst.female.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Pomc/Pomc_Fst_Female_higher.xlsx")

write.xlsx(Pomc.Fst.male.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Pomc/Pomc_Fst_Male_higher.xlsx")



# subclusters of Agrp and Pomc cell-types

# read in data
Agrp_Gm <- readRDS('data/Agrp_gm.rds')

# read in DE analyses from DE_analysis_for_subcluster.R
Agrp.Gm.F.Fd.v.Fst <- read.csv('data/Agrp_gm_F_Fd_v_Fst.csv')
Agrp.Gm.M.Fd.v.Fst <- read.csv('data/Agrp_gm_M_Fd_v_Fst.csv')
Agrp.Gm.Fd.F.v.M <- read.csv('data/Agrp_gm_Fd_F_v_M.csv')
Agrp.Gm.Fst.F.v.M <- read.csv('data/Agrp_gm_Fst_F_v_M.csv')

# plot volcano plot with DE genes in color female fed vs female fasted
FindMarkers(Agrp_Gm, group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'F_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Agrp.Gm.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Agrp.Gm.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Agrp.Gm.F.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-5.63,5.63)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/agrp_gm_female_fed_vs_fasted_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot with DE genes in color male fed vs male fasted
FindMarkers(Agrp_Gm, group.by = 'sexXnutr',ident.1 = 'M_Fed', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0,test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Agrp.Gm.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Agrp.Gm.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Agrp.Gm.M.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 0.9), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-6.11,6.11)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/agrp_gm_male_fed_vs_fasted_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot with DE genes in color female fed vs male fed
FindMarkers(Agrp_Gm, group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'M_Fed', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Agrp.Gm.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Agrp.Gm.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Agrp.Gm.Fd.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-12.34,12.34)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/agrp_gm_fed_female_vs_male_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot with DE genes in color female fasted vs male fasted
FindMarkers(Agrp_Gm, group.by = 'sexXnutr',ident.1 = 'F_Fast', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Agrp.Gm.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Agrp.Gm.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Agrp.Gm.Fst.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-10.39,10.39)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/agrp_gm_fasted_female_vs_male_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)



# read in data
Agrp_Sst <- readRDS('data/Agrp_sst.rds')

# read in DE analyses from DE_analysis_for_subcluster.R
Agrp.Sst.F.Fd.v.Fst <- read.csv('data/Agrp_Sst_F_Fd_v_Fst.csv')
Agrp.Sst.M.Fd.v.Fst <- read.csv('data/Agrp_Sst_M_Fd_v_Fst.csv')
Agrp.Sst.Fd.F.v.M <- read.csv('data/Agrp_Sst_Fd_F_v_M.csv')
Agrp.Sst.Fst.F.v.M <- read.csv('data/Agrp_Sst_Fst_F_v_M.csv')


# plot volcano plot with DE genes in color female fed vs female fasted
FindMarkers(Agrp_Sst, group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'F_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Agrp.Sst.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Agrp.Sst.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Agrp.Sst.F.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.73,2.73)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/agrp_sst_female_fed_vs_fasted_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot with DE genes in color male fed vs male fasted
FindMarkers(Agrp_Sst, group.by = 'sexXnutr',ident.1 = 'M_Fed', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0,test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Agrp.Sst.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Agrp.Sst.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Agrp.Sst.M.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 0.9), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-4.35,4.35)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/agrp_sst_male_fed_vs_fasted_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot with DE genes in color female fed vs male fed
FindMarkers(Agrp_Sst, group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'M_Fed', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Agrp.Sst.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Agrp.Sst.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Agrp.Sst.Fd.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-10.46,10.46)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/agrp_sst_fed_female_vs_male_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot with DE genes in color female fasted vs male fasted
FindMarkers(Agrp_Sst, group.by = 'sexXnutr',ident.1 = 'F_Fast', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Agrp.Sst.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Agrp.Sst.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Agrp.Sst.Fst.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-10.39,10.39)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/agrp_sst_fasted_female_vs_male_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)



# read in data
Pomc_Anxa2 <- readRDS('data/Pomc_Anxa2.rds')

# read in DE analyses from DE_analysis_for_subcluster.R
Pomc.Anxa2.F.Fd.v.Fst <- read.csv('data/Pomc_anxa2_F_Fd_v_Fst.csv')
Pomc.Anxa2.M.Fd.v.Fst <- read.csv('data/Pomc_anxa2_M_Fd_v_Fst.csv')
Pomc.Anxa2.Fd.F.v.M <- read.csv('data/Pomc_anxa2_Fd_F_v_M.csv')
Pomc.Anxa2.Fst.F.v.M <- read.csv('data/Pomc_anxa2_Fst_F_v_M.csv')


# plot volcano plot with DE genes in color female fed vs female fasted
FindMarkers(Pomc_Anxa2, group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'F_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Pomc.Anxa2.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Pomc.Anxa2.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Pomc.Anxa2.F.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-5.59,5.59)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/pomc_anxa2_female_fed_vs_fasted_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot with DE genes in color male fed vs male fasted
FindMarkers(Pomc_Anxa2, group.by = 'sexXnutr',ident.1 = 'M_Fed', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0,test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Pomc.Anxa2.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Pomc.Anxa2.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Pomc.Anxa2.M.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 0.9), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-3.05,3.05)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/pomc_anxa2_male_fed_vs_fasted_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot with DE genes in color female fed vs male fed
FindMarkers(Pomc_Anxa2, group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'M_Fed', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Pomc.Anxa2.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Pomc.Anxa2.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Pomc.Anxa2.Fd.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-11.31,11.31)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/pomc_anxa2_fed_female_vs_male_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot with DE genes in color female fasted vs male fasted
FindMarkers(Pomc_Anxa2, group.by = 'sexXnutr',ident.1 = 'F_Fast', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Pomc.Anxa2.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Pomc.Anxa2.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Pomc.Anxa2.Fst.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-9.35,9.35)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/pomc_anxa2_fasted_female_vs_male_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)



# read in data
Pomc_Glipr1 <- readRDS('data/Pomc_Glipr1.rds')

# read in DE analyses from DE_analysis_for_subcluster.R
Pomc.Glipr1.F.Fd.v.Fst <- read.csv('data/Pomc_glipr1_F_Fd_v_Fst.csv')
Pomc.Glipr1.M.Fd.v.Fst <- read.csv('data/Pomc_glipr1_M_Fd_v_Fst.csv')
Pomc.Glipr1.Fd.F.v.M <- read.csv('data/Pomc_glipr1_Fd_F_v_M.csv')
Pomc.Glipr1.Fst.F.v.M <- read.csv('data/Pomc_glipr1_Fst_F_v_M.csv')


# plot volcano plot with DE genes in color female fed vs female fasted
FindMarkers(Pomc_Glipr1, group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'F_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Pomc.Glipr1.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Pomc.Glipr1.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Pomc.Glipr1.F.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-3.3,3.3)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/pomc_glipr1_female_fed_vs_fasted_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot with DE genes in color male fed vs male fasted
FindMarkers(Pomc_Glipr1, group.by = 'sexXnutr',ident.1 = 'M_Fed', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0,test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Pomc.Glipr1.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Pomc.Glipr1.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Pomc.Glipr1.M.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 0.9), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-1.83,1.83)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/pomc_glipr1_male_fed_vs_fasted_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot with DE genes in color female fed vs male fed
FindMarkers(Pomc_Glipr1, group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'M_Fed', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Pomc.Glipr1.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Pomc.Glipr1.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Pomc.Glipr1.Fd.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-11.72,11.72)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/pomc_glipr1_fed_female_vs_male_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot with DE genes in color female fasted vs male fasted
FindMarkers(Pomc_Glipr1, group.by = 'sexXnutr',ident.1 = 'F_Fast', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Pomc.Glipr1.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Pomc.Glipr1.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Pomc.Glipr1.Fst.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-7.05,7.05)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/pomc_glipr1_fasted_female_vs_male_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)



# read in data
Pomc_Ttr <- readRDS('data/Pomc_Ttr.rds')

# read in DE analyses from DE_analysis_for_subcluster.R
Pomc.Ttr.F.Fd.v.Fst <- read.csv('data/Pomc_ttr_F_Fd_v_Fst.csv')
Pomc.Ttr.M.Fd.v.Fst <- read.csv('data/Pomc_ttr_M_Fd_v_Fst.csv')
Pomc.Ttr.Fd.F.v.M <- read.csv('data/Pomc_ttr_Fd_F_v_M.csv')
Pomc.Ttr.Fst.F.v.M <- read.csv('data/Pomc_ttr_Fst_F_v_M.csv')


# plot volcano plot with DE genes in color female fed vs female fasted
FindMarkers(Pomc_Ttr, group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'F_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Pomc.Ttr.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Pomc.Ttr.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Pomc.Ttr.F.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-3.07,3.07)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/pomc_ttr_female_fed_vs_fasted_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot with DE genes in color male fed vs male fasted
FindMarkers(Pomc_Ttr, group.by = 'sexXnutr',ident.1 = 'M_Fed', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0,test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Pomc.Ttr.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Pomc.Ttr.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Pomc.Ttr.M.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 0.9), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-1.23,1.23)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/pomc_ttr_male_fed_vs_fasted_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot with DE genes in color female fed vs male fed
FindMarkers(Pomc_Ttr, group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'M_Fed', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Pomc.Ttr.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Pomc.Ttr.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Pomc.Ttr.Fd.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-12.36,12.36)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/pomc_ttr_fed_female_vs_male_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)


# plot volcano plot with DE genes in color female fasted vs male fasted
FindMarkers(Pomc_Ttr, group.by = 'sexXnutr',ident.1 = 'F_Fast', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Pomc.Ttr.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Pomc.Ttr.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Pomc.Ttr.Fst.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-10.89,10.89)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures3/pomc_ttr_fasted_female_vs_male_volc.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)

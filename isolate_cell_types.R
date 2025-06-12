libs <- c( 'gplots','stringi','reshape2','cowplot','RColorBrewer',
           'sctransform','stringr','org.Mm.eg.db','AnnotationDbi',
           'IRanges','S4Vectors','Biobase','BiocGenerics','clusterProfiler',
           'biomaRt','Matrix','DESeq2','RcppThread', 'extrafont', 'openxlsx',
           'Seurat','dplyr','tidyr','ggplot2','harmony','ggalluvial',
           'scDblFinder','SoupX','UpSetR','ComplexUpset','CellChat','gprofiler2','NMF','ggalluvial')

lapply(libs, require, character.only = TRUE)



#create basic clustering for sample
ARH_Sex_by_Nutr <- readRDS('data/ARH_Sex_by_Nutr.rds')
Agrp_Sex_by_Nutr <- readRDS('data/Agrp_Sex_by_Nutr.rds')

Agrp_Sex_by_Nutr <- subset(ARH_Sex_by_Nutr, subset = cell_type3 == 'Agrp')

Agrp_Sex_by_Nutr <- SCTransform(Agrp_Sex_by_Nutr, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')


Agrp_Sex_by_Nutr@assays$SCT@var.features <- Agrp_Sex_by_Nutr@assays$SCT@var.features[(!Agrp_Sex_by_Nutr@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

Agrp_Sex_by_Nutr <- RunPCA(Agrp_Sex_by_Nutr, features = VariableFeatures(object = Agrp_Sex_by_Nutr))
Agrp_Sex_by_Nutr <- RunUMAP(Agrp_Sex_by_Nutr, reduction = "pca", dims = 1:30)

Agrp_Sex_by_Nutr <- FindNeighbors(Agrp_Sex_by_Nutr, dims = 1:30)
Agrp_Sex_by_Nutr <- FindClusters(Agrp_Sex_by_Nutr, resolution = 0.5)

DimPlot(Agrp_Sex_by_Nutr, label = TRUE, label.size = 8, repel = FALSE, reduction = 'tsne') + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())




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


Agrp.F.Fd.Fst.suppress.gost <- gost(query = (Agrp.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Agrp.F.Fd.Fst.induce.gost <- gost(query = (Agrp.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


write.xlsx(Agrp.F.Fd.Fst.suppress.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Agrp/Agrp_F_Fd_Fst_suppress.xlsx")

Agrp.F.Fd.Fst.suppress.gost$result[c(4,17,18,22, 23),] |> ggplot(aes(-log10(p_value), 
                                                                     factor(term_name, levels = rev(c('PDZ domain binding',
                                                                                                     'potassium channel activity',
                                                                                                     'voltage-gated calcium channel activity',
                                                                                                     'kinase binding',
                                                                                                     'transporter activity'))))) + 
  geom_col(fill = '#6a816a', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/agrp_fasting_suppressed_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)



write.xlsx(Agrp.F.Fd.Fst.induce.gost$result[c(1,4,8,9,12),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Agrp/Agrp_F_Fd_Fst_induce.xlsx")

Agrp.F.Fd.Fst.induce.gost$result[c(3,5,7,8,18),] |> ggplot(aes(-log10(p_value), 
                                                                     factor(term_name, levels =c('cAMP binding',
                                                                                                'protein kinase binding',
                                                                                                'GTPase regulator activity',
                                                                                                'transporter activity',
                                                                                                'enzyme binding')))) + 
  geom_col(fill = '#f9994f', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/agrp_fasting_induced_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)





Agrp.M.Fd.Fst.suppress.gost <- gost(query = (Agrp.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Agrp.M.Fd.Fst.induce.gost <- gost(query = (Agrp.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


write.xlsx(Agrp.M.Fd.Fst.suppress.gost$result[c(4,5,11,14,15),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Agrp/Agrp_M_Fd_Fst_suppress.xlsx")

Agrp.M.Fd.Fst.suppress.gost$result[c(3,7,10,11,15),] |> ggplot(aes(-log10(p_value), 
                                                                     factor(term_name, levels = c('protein kinase binding',
                                                                                                  'cell-cell adhesion mediator activity',
                                                                                                  'GPI-linked ephrin receptor activity',
                                                                                                  'voltage-gated channel activity',
                                                                                                  'protein domain specific binding')))) + 
  geom_col(fill = '#19552b', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/agrp_fasting_suppressed_in_males_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)



write.xlsx(Agrp.M.Fd.Fst.induce.gost$result[c(4,5,11,14,15),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Agrp/Agrp_M_Fd_Fst_induce.xlsx")

Agrp.M.Fd.Fst.induce.gost$result[c(5,6,9,10,14),] |> ggplot(aes(-log10(p_value), 
                                                               factor(term_name, levels = c('ion binding',
                                                                                            'cytoskeletal protein binding',
                                                                                            'cAMP binding',
                                                                                            'GTPase binding',
                                                                                            'kinase binding')))) + 
  geom_col(fill = '#ff6e00', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/agrp_fasting_induced_in_males_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)



####

Agrp.Fd.female.higher.gost <- gost(query = (Agrp.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                   organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Agrp.Fd.male.higher.gost <- gost(query = (Agrp.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                 organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


write.xlsx(Agrp.Fd.female.higher.gost[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-01-06/GOMF_Tables/Agrp_Fd_Female_higher.xlsx")

Agrp.Fd.female.higher.gost$result[c(1:5),] |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name))) + 
  geom_col(fill = '#6a816a', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Agrp_higher_in_fed_females_GOMF.tiff', device = 'tiff', units = 'in', width = 4, height = 1.5, dpi = 600)



write.xlsx(Agrp.Fd.male.higher.gost$result[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-01-06/GOMF_Tables/Agrp_Fd_male_higher.xlsx")

Agrp.Fd.male.higher.gost$result[c(1:5),] |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name, levels = c('intracellular ligand-gated monoatomic ion channel activity',
                                                                                             'ligand-gated channel activity',
                                                                                             'ligand-gated monoatomic ion channel activity',
                                                                                             'cadherin binding',
                                                                                             'gated channel activity')))) + 
  geom_col(fill = '#19552b', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Agrp_higher_in_fed_males_GOMF.tiff', device = 'tiff', units = 'in', width = 4, height = 1.5, dpi = 600)



#####

Agrp.Fst.female.higher.gost <- gost(query = (Agrp.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Agrp.Fst.male.higher.gost <- gost(query = (Agrp.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


write.xlsx(Agrp.Fst.female.higher.gost$result[c(1),c(9,11,3)], 
           file = "../paper_figures/post_2025-01-06/GOMF_Tables/Agrp_Fst_Female_higher.xlsx")

Agrp.Fst.female.higher.gost$result[c(1),] |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name))) + 
  geom_col(fill = '#f9994f', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Agrp_higher_in_fasted_females_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 0.75, dpi = 600)



write.xlsx(Agrp.Fst.male.higher.gost$result[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-01-06/GOMF_Tables/Agrp_Fst_male_higher.xlsx")

Agrp.Fst.male.higher.gost$result[c(1:5),] |> ggplot(aes(-log10(p_value), 
                                                               factor(term_name, levels = c('ligand-gated monoatomic ion channel activity involved in regulation of presynaptic membrane potential',
                                                                                            'protein binding',
                                                                                            'glutamate-gated receptor activity',
                                                                                            'protein domain specific binding',
                                                                                            'PDZ domain binding')))) + 
  geom_col(fill = '#ff6e00', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Agrp_higher_in_fasted_males_GOMF.tiff', device = 'tiff', units = 'in', width = 7, height = 1.5, dpi = 600)







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
KNDy_Sex_by_Nutr <- readRDS('data/KNDy_Sex_by_Nutr.rds')

KNDy_Sex_by_Nutr <- subset(ARH_Sex_by_Nutr, subset = cell_type3 == 'KNDy')

KNDy_Sex_by_Nutr <- SCTransform(KNDy_Sex_by_Nutr, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')


KNDy_Sex_by_Nutr@assays$SCT@var.features <- KNDy_Sex_by_Nutr@assays$SCT@var.features[(!KNDy_Sex_by_Nutr@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

KNDy_Sex_by_Nutr <- RunPCA(KNDy_Sex_by_Nutr, features = VariableFeatures(object = KNDy_Sex_by_Nutr))
KNDy_Sex_by_Nutr <- RunTSNE(KNDy_Sex_by_Nutr, reduction = "pca", dims = 1:30)

KNDy_Sex_by_Nutr <- FindNeighbors(KNDy_Sex_by_Nutr, dims = 1:30)
KNDy_Sex_by_Nutr <- FindClusters(KNDy_Sex_by_Nutr, resolution = 0.5)

DimPlot(KNDy_Sex_by_Nutr, label = TRUE, label.size = 8, repel = FALSE, reduction = 'tsne') + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())




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



saveRDS(KNDy_Sex_by_Nutr, file = 'data/KNDy_Sex_by_Nutr.rds')


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


KNDy.F.Fd.Fst.suppress.gost <- gost(query = (KNDy.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
KNDy.F.Fd.Fst.induce.gost <- gost(query = (KNDy.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


write.xlsx(KNDy.F.Fd.Fst.suppress.gost$result[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/KNDy_F_Fd_Fst_suppress.xlsx")



KNDy.F.Fd.Fst.suppress.gost$result[c(4,6,9,11,13),] |> ggplot(aes(-log10(p_value), 
                                                                     factor(term_name, levels =c('cAMP binding', 
                                                                                                 'GTPase regulator activity',
                                                                                                 '14-3-3 protein binding',
                                                                                                 'nuclear receptor activity',
                                                                                                 'voltage-gated channel activity') ))) + 
  geom_col(fill = '#6a816a', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/KNDy_fasting_suppressed_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)




write.xlsx(KNDy.F.Fd.Fst.induce.gost$result[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/KNDy_F_Fd_Fst_induce.xlsx")





KNDy.F.Fd.Fst.induce.gost$result[c(2,3,6,13,18),] |> ggplot(aes(-log10(p_value), 
                                                               factor(term_name, levels = c('potassium channel regulator activity',
                                                                                            'passive transmembrane transporter activity',
                                                                                            'ligand-gated channel activity',
                                                                                            'glutamate receptor activity',
                                                                                            'PDZ domain binding')))) + 
  geom_col(fill = '#f9994f', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/KNDy_fasting_induced_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)





KNDy.M.Fd.Fst.suppress.gost <- gost(query = (KNDy.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
KNDy.M.Fd.Fst.induce.gost <- gost(query = (KNDy.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)



write.xlsx(KNDy.M.Fd.Fst.suppress.gost$result[c(3,7,10,11,15),c(9,11,3)], 
           file = "../paper_figures/post_2025-01-06/GOMF_Tables/KNDy/KNDy_M_Fd_Fst_suppress.xlsx")





KNDy.M.Fd.Fst.suppress.gost$result[c(3,7,10,11,15),] |> ggplot(aes(-log10(p_value), 
                                                                   factor(term_name, levels = c('protein kinase binding',
                                                                                                'cell-cell adhesion mediator activity',
                                                                                                'GPI-linked ephrin receptor activity',
                                                                                                'voltage-gated channel activity',
                                                                                                'protein domain specific binding')))) + 
  geom_col(fill = '#19552b', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/KNDy_fasting_suppressed_in_males_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)




write.xlsx(KNDy.M.Fd.Fst.induce.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/KNDy_M_Fd_Fst_induce.xlsx")





KNDy.M.Fd.Fst.induce.gost$result[c(5,6,9,10,14),] |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name, levels = c('ion binding',
                                                                                             'cytoskeletal protein binding',
                                                                                             'cAMP binding',
                                                                                             'GTPase binding',
                                                                                             'kinase binding')))) + 
  geom_col(fill = '#ff6e00', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/KNDy_fasting_induced_in_males_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)





#####

KNDy.Fd.female.higher.gost <- gost(query = (KNDy.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
KNDy.Fd.male.higher.gost <- gost(query = (KNDy.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)



write.xlsx(KNDy.Fd.female.higher.gost$result[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/KNDy_Fd_Female_higher.xlsx")





KNDy.Fd.female.higher.gost$result[c(3,5,6,8,10),] |> ggplot(aes(-log10(p_value), 
                                                                   factor(term_name, levels = c('intracellular sodium-activated potassium channel activity',
                                                                                                'cell adhesion molecule binding',
                                                                                                'transcription coregulator binding',
                                                                                                'phosphoric diester hydrolase activity',
                                                                                                'nuclear receptor activity')))) + 
  geom_col(fill = '#6a816a', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/KNDy_higher_in_fed_females_GOMF.tiff', device = 'tiff', units = 'in', width = 4, height = 1.5, dpi = 600)




write.xlsx(KNDy.Fd.male.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/KNDy_Fd_male_higher.xlsx")



KNDy.Fd.male.higher.gost$result[c(3,5,10,12,14),] |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name, levels = c('PDZ domain binding',
                                                                                             'glutamate receptor activity',
                                                                                             'GTPase regulator activity',
                                                                                             'passive transmembrane transporter activity',
                                                                                             'monoatomic ion channel activity')))) + 
  geom_col(fill = '#19552b', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/KNDy_higher_in_fed_males_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)



#####

KNDy.Fst.female.higher.gost <- gost(query = (KNDy.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                   organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
KNDy.Fst.male.higher.gost <- gost(query = (KNDy.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                 organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)



write.xlsx(KNDy.Fst.female.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/KNDy_Fst_Female_higher.xlsx")



KNDy.Fst.female.higher.gost$result[c(1,4,5,7,8),] |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name, levels = c('AMPA glutamate receptor activity',
                                                                                             'pancreatic polypeptide receptor activity',
                                                                                             'PDZ domain binding',
                                                                                             'calmodulin binding',
                                                                                             'glutamate receptor activity')))) + 
  geom_col(fill = '#f9994f', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/KNDy_higher_in_fasted_females_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)



write.xlsx(KNDy.Fst.male.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/KNDy_Fst_male_higher.xlsx")




KNDy.Fst.male.higher.gost$result[c(1,2,4,5,10),] |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name, levels = c('voltage-gated calcium channel activity',
                                                                                             'cell adhesion molecule binding',
                                                                                             'transporter regulator activity',
                                                                                             'cadherin binding',
                                                                                             'actin binding')))) + 
  geom_col(fill = '#ff6e00', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/KNDy_higher_in_fasted_males_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)






##### DA
DA_Sex_by_Nutr <- readRDS('data/DA_Sex_by_Nutr.rds')

DA_Sex_by_Nutr <- subset(ARH_Sex_by_Nutr, subset = cell_type3 == 'DA')

DA_Sex_by_Nutr <- SCTransform(DA_Sex_by_Nutr, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')


DA_Sex_by_Nutr@assays$SCT@var.features <- DA_Sex_by_Nutr@assays$SCT@var.features[(!DA_Sex_by_Nutr@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

DA_Sex_by_Nutr <- RunPCA(DA_Sex_by_Nutr, features = VariableFeatures(object = DA_Sex_by_Nutr))
DA_Sex_by_Nutr <- RunTSNE(DA_Sex_by_Nutr, reduction = "pca", dims = 1:30)

DA_Sex_by_Nutr <- FindNeighbors(DA_Sex_by_Nutr, dims = 1:30)
DA_Sex_by_Nutr <- FindClusters(DA_Sex_by_Nutr, resolution = 0.5)

DimPlot(DA_Sex_by_Nutr, label = TRUE, label.size = 8, repel = FALSE, reduction = 'tsne') + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())




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






saveRDS(DA_Sex_by_Nutr, file = 'data/DA_Sex_by_Nutr.rds')


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


DA.F.Fd.Fst.suppress.gost <- gost(query = (DA.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
DA.F.Fd.Fst.induce.gost <- gost(query = (DA.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)



write.xlsx(DA.F.Fd.Fst.suppress.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/DA_F_Fd_Fst_suppress.xlsx")




DA.F.Fd.Fst.suppress.gost$result |> ggplot(aes(-log10(p_value), 
                                                                  factor(term_name, levels = c('signaling receptor binding',
                                                                                               'protein binding',
                                                                                               'PDZ domain binding') ))) + 
  geom_col(fill = '#6a816a', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/DA_fasting_suppressed_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)




write.xlsx(DA.F.Fd.Fst.induce.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/DA_F_Fd_Fst_suppress.xlsx")



DA.F.Fd.Fst.induce.gost$result |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name))) + 
  geom_col(fill = '#f9994f', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/DA_fasting_induced_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)





DA.M.Fd.Fst.suppress.gost <- gost(query = (DA.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
DA.M.Fd.Fst.induce.gost <- gost(query = (DA.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)








#####

DA.Fd.female.higher.gost <- gost(query = (DA.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                   organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
DA.Fd.male.higher.gost <- gost(query = (DA.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                 organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)



write.xlsx(DA.Fd.female.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/DA_Fd_Female_higher.xlsx")



DA.Fd.female.higher.gost$result |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name))) + 
  geom_col(fill = '#6a816a', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/DA_higher_in_fed_females_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = .75, dpi = 600)



#write.xlsx(DA.Fd.male.higher.gost$result[c(1,2,7,10,15),c(9,11,3)], 
#           file = "../paper_figures/post_2025-01-06/GOMF_Tables/KNDy/DA_Fd_Male_higher.xlsx")




DA.Fd.male.higher.gost$result[c(1,2,7,10,15),] |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name, levels =c('scaffold protein binding',
                                                                                            'potassium ion transmembrane transporter activity',
                                                                                            'transmitter-gated channel activity',
                                                                                            'cell adhesion molecule binding',
                                                                                            'glutamate receptor activity')))) + 
  geom_col(fill = '#19552b', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/DA_higher_in_fed_males_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)



#####

DA.Fst.female.higher.gost <- gost(query = (DA.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
DA.Fst.male.higher.gost <- gost(query = (DA.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)



write.xlsx(DA.Fst.female.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/KNDy/DA_Fst_Female_higher.xlsx")



DA.Fst.female.higher.gost$result |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name, levels = c('neurotrophin binding',
                                                                                             'brain-derived neurotrophic factor binding',
                                                                                             'netrin receptor activity')))) + 
  geom_col(fill = '#f9994f', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/DA_higher_in_fasted_females_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)



write.xlsx(DA.Fst.male.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-01-06/GOMF_Tables/KNDy/DA_Fst_Male_higher.xlsx")




DA.Fst.male.higher.gost$result |> ggplot(aes(-log10(p_value), 
                                                               factor(term_name))) + 
  geom_col(fill = '#ff6e00', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/DA_higher_in_fasted_males_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)





####

##### Ghrh.Chat
Ghrh.Chat_Sex_by_Nutr <- readRDS('data/Ghrh.Chat_Sex_by_Nutr.rds')

Ghrh.Chat_Sex_by_Nutr <- subset(ARH_Sex_by_Nutr, subset = cell_type3 == 'Ghrh/Chat')

Ghrh.Chat_Sex_by_Nutr <- SCTransform(Ghrh.Chat_Sex_by_Nutr, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')


Ghrh.Chat_Sex_by_Nutr@assays$SCT@var.features <- Ghrh.Chat_Sex_by_Nutr@assays$SCT@var.features[(!Ghrh.Chat_Sex_by_Nutr@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

Ghrh.Chat_Sex_by_Nutr <- RunPCA(Ghrh.Chat_Sex_by_Nutr, features = VariableFeatures(object = Ghrh.Chat_Sex_by_Nutr))
Ghrh.Chat_Sex_by_Nutr <- RunTSNE(Ghrh.Chat_Sex_by_Nutr, reduction = "pca", dims = 1:30)

Ghrh.Chat_Sex_by_Nutr <- FindNeighbors(Ghrh.Chat_Sex_by_Nutr, dims = 1:30)
Ghrh.Chat_Sex_by_Nutr <- FindClusters(Ghrh.Chat_Sex_by_Nutr, resolution = 0.5)

DimPlot(Ghrh.Chat_Sex_by_Nutr, label = TRUE, label.size = 8, repel = FALSE, reduction = 'tsne') + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())




DimPlot(Ghrh.Chat_Sex_by_Nutr, label = F, label.size = 2,pt.size = 1, repel = TRUE, 
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
ggsave('figures/Ghrh.Chat_sex_by_nutr_dimplot.tiff', device = 'tiff', units = 'in', width = 4.5,height = 3.5,dpi = 600)


FeaturePlot(Ghrh.Chat_Sex_by_Nutr, label = FALSE, 
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


ggsave('figures/Ghrh.Chat_f_fed_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)


FeaturePlot(Ghrh.Chat_Sex_by_Nutr, label = FALSE, 
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


ggsave('figures/Ghrh.Chat_f_fast_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)

FeaturePlot(Ghrh.Chat_Sex_by_Nutr, label = FALSE, 
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


ggsave('figures/Ghrh.Chat_m_fed_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)


FeaturePlot(Ghrh.Chat_Sex_by_Nutr, label = FALSE, 
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


ggsave('figures/Ghrh.Chat_m_fast_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)



FeaturePlot(Ghrh.Chat_Sex_by_Nutr, label = FALSE, 
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


ggsave('figures/Ghrh.Chat_fos_dimplot.tiff', device = 'tiff', units = 'in', width = 3.5,height = 3.5,dpi = 600)


###barplot for fos in Ghrh/Chat ####

Ghrh.Chat_Sex_by_Nutr$SCT@data[c('Th','Fos'),] %>% 
  t() %>% 
  as.data.frame() %>% 
  cbind(Ghrh.Chat_Sex_by_Nutr@meta.data[,c(1:5)]) %>% 
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
ggsave('figures/Fos_in_Ghrh.Chat_barplot', device = 'tiff', units = 'in', width = 4, height = 2, dpi = 600)

FindMarkers(Ghrh.Chat_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'F_Fed', 
            ident.2 = 'F_Fast', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)


FindMarkers(Ghrh.Chat_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'M_Fed', 
            ident.2 = 'M_Fast', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)




FindMarkers(Ghrh.Chat_Sex_by_Nutr, 
            features = 'Fos',
            group.by = 'sexXnutr', 
            ident.1 = 'F_Fed', 
            ident.2 = 'M_Fed', 
            logfc.threshold = 0, 
            min.pct = 0, 
            pseudocount.use = 0)





####






saveRDS(Ghrh.Chat_Sex_by_Nutr, file = 'data/Ghrh.Chat_Sex_by_Nutr.rds')


FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Ghrh/Chat', group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'F_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Ghrh.Chat.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Ghrh.Chat.F.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Ghrh.Chat.F.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Ghrh.Chat_female_fed_vs_fasted_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)



FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Ghrh/Chat', group.by = 'sexXnutr',ident.1 = 'M_Fed', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Ghrh.Chat.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Ghrh.Chat.M.Fd.v.Fst, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Ghrh.Chat.M.Fd.v.Fst, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Fed/Fasted)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Ghrh.Chat_male_fed_vs_fasted_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)




####


FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Ghrh/Chat', group.by = 'sexXnutr',ident.1 = 'F_Fed', ident.2 = 'M_Fed', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Ghrh.Chat.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#6a816a', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Ghrh.Chat.Fd.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#19552b', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Ghrh.Chat.Fd.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Ghrh.Chat_fed_female_vs_male_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)



FindMarkers(ARH_Sex_by_Nutr, subset.ident = 'Ghrh/Chat', group.by = 'sexXnutr',ident.1 = 'F_Fast', ident.2 = 'M_Fast', 
            logfc.threshold = 0, min.pct = 0, test.use = 'MAST', latent.vars = 'Sample_ID') |> 
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) + geom_point(shape = 21, color = 'grey', fill= 'grey', size =1) +
  geom_point(inherit.aes = F, data = filter(Ghrh.Chat.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC > 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#f9994f', color = 'black', shape = 21, size = 1.5) +
  geom_point(inherit.aes = F, data = filter(Ghrh.Chat.Fst.F.v.M, p_val_adj < 0.05, avg_log2FC < 0), 
             aes(avg_log2FC, -log10(p_val_adj)), fill = '#ff6e00', color = 'black', shape = 21, size = 1.5) +
  ggrepel::geom_text_repel(data = filter(Ghrh.Chat.Fst.F.v.M, p_val_adj < 0.05, abs(avg_log2FC) > 1), 
                           aes(avg_log2FC, -log10(p_val_adj), label = gene), size = 5/.pt) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  labs(x = 'log2(Female/Male)', y= '-log10(adj. p-value)') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Ghrh.Chat_fasted_female_vs_male_volc2.tiff', device = 'tiff', units = 'in', width = 3.5, height = 2, dpi = 600)

library(gprofiler2)
######20240516#######

Ghrh.Chat.F.Fd.Fst.suppress.gost <- gost(query = (Ghrh.Chat.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Ghrh.Chat.F.Fd.Fst.induce.gost <- gost(query = (Ghrh.Chat.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)



#write.xlsx(Ghrh.Chat.F.Fd.Fst.suppress.gost$result[c(2,5,6,8,10),c(9,11,3)], 
 #          file = "../paper_figures/post_2025-01-06/GOMF_Tables/Ghrh/Ghrh_F_Fd_Fst_suppress.xlsx")



Ghrh.Chat.F.Fd.Fst.suppress.gost$result[c(2,5,6,8,10),] |> ggplot(aes(-log10(p_value), 
                                                                  factor(term_name, levels = c('passive transmembrane transporter activity',
                                                                                               'potassium ion transmembrane transporter activity',
                                                                                               'GABA-gated chloride ion channel activity',
                                                                                               'voltage-gated monoatomic cation channel activity',
                                                                                               'monoatomic ion channel activity') ))) + 
  geom_col(fill = '#6a816a', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Ghrh.Chat_fasting_suppressed_in_female_GOMF.tiff', device = 'tiff', units = 'in', width = 3.75, height = 1.5, dpi = 600)




#write.xlsx(Ghrh.Chat.F.Fd.Fst.induce.gost$result[c(2,3,4,10,16),c(9,11,3)], 
 #          file = "../paper_figures/post_2025-01-06/GOMF_Tables/Ghrh/Ghrh_F_Fd_Fst_induce.xlsx")



Ghrh.Chat.F.Fd.Fst.induce.gost$result[c(2,3,4,10,16),] |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name, levels = c('PDZ domain binding',
                                                                                             'passive transmembrane transporter activity',
                                                                                             'gated channel activity',
                                                                                             'glutamate receptor activity',
                                                                                             'G protein-coupled glutamate receptor activity')))) + 
  geom_col(fill = '#f9994f', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Ghrh.Chat_fasting_induced_in_female_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)





Ghrh.Chat.M.Fd.Fst.suppress.gost <- gost(query = (Ghrh.Chat.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Ghrh.Chat.M.Fd.Fst.induce.gost <- gost(query = (Ghrh.Chat.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


####20240517#####


#write.xlsx(Ghrh.Chat.M.Fd.Fst.suppress.gost$result[,c(9,11,3)], 
 #          file = "../paper_figures/post_2025-01-06/GOMF_Tables/Ghrh/Ghrh_M_Fd_Fst_suppress.xlsx")



Ghrh.Chat.M.Fd.Fst.suppress.gost$result |> ggplot(aes(-log10(p_value), 
                                                                   factor(term_name))) + 
  geom_col(fill = '#19552b', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Ghrh.Chat_fasting_suppressed_in_males_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)



#no results

Ghrh.Chat.M.Fd.Fst.induce.gost$result |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name))) + 
  geom_col(fill = '#ff6e00', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
#ggsave(filename = 'figures/Ghrh.Chat_fasting_induced_in_males_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)





#####

Ghrh.Chat.Fd.female.higher.gost <- gost(query = (Ghrh.Chat.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                   organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Ghrh.Chat.Fd.male.higher.gost <- gost(query = (Ghrh.Chat.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                 organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)




#write.xlsx(Ghrh.Chat.Fd.female.higher.gost$result[,c(9,11,3)], 
 #          file = "../paper_figures/post_2025-01-06/GOMF_Tables/Ghrh/Ghrh_Fd_Female_higher.xlsx")




Ghrh.Chat.Fd.female.higher.gost$result |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name))) + 
  geom_col(fill = '#6a816a', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Ghrh.Chat_higher_in_fed_females_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)




#write.xlsx(Ghrh.Chat.Fd.male.higher.gost$result[,c(9,11,3)], 
 #          file = "../paper_figures/post_2025-01-06/GOMF_Tables/Ghrh/Ghrh_Fd_Male_higher.xlsx")



Ghrh.Chat.Fd.male.higher.gost$result |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name, levels = c(
                                                                  'PDZ domain binding',
                                                                  'hormone binding',
                                                                  'cell adhesion molecule binding'
                                                                )))) + 
  geom_col(fill = '#19552b', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Ghrh.Chat_higher_in_fed_males_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)



#####

Ghrh.Chat.Fst.female.higher.gost <- gost(query = (Ghrh.Chat.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Ghrh.Chat.Fst.male.higher.gost <- gost(query = (Ghrh.Chat.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)




#write.xlsx(Ghrh.Chat.Fst.female.higher.gost$result[,c(9,11,3)], 
 #          file = "../paper_figures/post_2025-01-06/GOMF_Tables/Ghrh/Ghrh_Fst_Female_higher.xlsx")




Ghrh.Chat.Fst.female.higher.gost$result |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name))) + 
  geom_col(fill = '#f9994f', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Ghrh.Chat_higher_in_fasted_females_GOMF.tiff', device = 'tiff', units = 'in', width = 5, height = 1.5, dpi = 600)




#write.xlsx(Ghrh.Chat.Fst.male.higher.gost$result[,c(9,11,3)], 
 #          file = "../paper_figures/post_2025-01-06/GOMF_Tables/Ghrh/Ghrh_Fst_Male_higher.xlsx")



Ghrh.Chat.Fst.male.higher.gost$result |> ggplot(aes(-log10(p_value), 
                                                               factor(term_name, levels = c(
                                                                 '[heparan sulfate]-glucosamine 3-sulfotransferase 1 activity',
                                                                 'netrin receptor activity',
                                                                 'heparan sulfate sulfotransferase activity'
                                                               )))) + 
  geom_col(fill = '#ff6e00', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Ghrh.Chat_higher_in_fasted_males_GOMF.tiff', device = 'tiff', units = 'in', width = 4, height = 1.5, dpi = 600)





###microglia####

#create basic clustering for sample
ARH_Sex_by_Nutr <- readRDS('data/ARH_Sex_by_Nutr.rds')
Microglia_Sex_by_Nutr <- readRDS('data/Microglia_Sex_by_Nutr.rds')

Microglia_Sex_by_Nutr <- subset(ARH_Sex_by_Nutr, subset = cell_type3 == 'Microglia')

Microglia_Sex_by_Nutr <- SCTransform(Microglia_Sex_by_Nutr, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')


Microglia_Sex_by_Nutr@assays$SCT@var.features <- Microglia_Sex_by_Nutr@assays$SCT@var.features[(!Microglia_Sex_by_Nutr@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

Microglia_Sex_by_Nutr <- RunPCA(Microglia_Sex_by_Nutr, features = VariableFeatures(object = Microglia_Sex_by_Nutr))
Microglia_Sex_by_Nutr <- RunTSNE(Microglia_Sex_by_Nutr, reduction = "pca", dims = 1:30)

Microglia_Sex_by_Nutr <- FindNeighbors(Microglia_Sex_by_Nutr, dims = 1:30)
Microglia_Sex_by_Nutr <- FindClusters(Microglia_Sex_by_Nutr, resolution = 0.5)

DimPlot(Microglia_Sex_by_Nutr, label = TRUE, label.size = 8, repel = FALSE, reduction = 'tsne') + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())




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



ggsave('figures/Microglia_sex_by_nutr_dimplot.tiff', device = 'tiff', units = 'in', width = 4.5,height = 3.5,dpi = 600)


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





saveRDS(Microglia_Sex_by_Nutr, file = 'data/Microglia_Sex_by_Nutr.rds')


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


Microglia.F.Fd.Fst.suppress.gost <- gost(query = (Microglia.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Microglia.F.Fd.Fst.induce.gost <- gost(query = (Microglia.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)



write.xlsx(Microglia.F.Fd.Fst.induce.gost$result[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Microglia/Microglia_F_Fd_Fst_induce.xlsx")



Microglia.F.Fd.Fst.induce.gost$result |> ggplot(aes(-log10(p_value), 
                                                               factor(term_name, levels = c('dystroglycan binding',
                                                                                            'transporter regulator activity',
                                                                                            'channel regulator activity',
                                                                                            'ion channel regulator activity')))) + 
  geom_col(fill = '#f9994f', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Microglia_fasting_induced_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)





Microglia.M.Fd.Fst.suppress.gost <- gost(query = (Microglia.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Microglia.M.Fd.Fst.induce.gost <- gost(query = (Microglia.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)






write.xlsx(Microglia.M.Fd.Fst.induce.gost$result[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Microglia/Microglia_M_Fd_Fst_induce.xlsx")




Microglia.M.Fd.Fst.induce.gost$result[c(1:5),] |> ggplot(aes(-log10(p_value), 
                                                                factor(term_name, levels = c('ligand-gated monoatomic ion channel activity',
                                                                                             'glutamate receptor activity',
                                                                                             'PDZ domain binding',
                                                                                             'ligand-gated sodium channel activity',
                                                                                             'gated channel activity'
                                                                                             )))) + 
  geom_col(fill = '#ff6e00', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Microglia_fasting_induced_in_males_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)



####

Microglia.Fd.female.higher.gost <- gost(query = (Microglia.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                   organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Microglia.Fd.male.higher.gost <- gost(query = (Microglia.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                 organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)






#####

Microglia.Fst.female.higher.gost <- gost(query = (Microglia.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Microglia.Fst.male.higher.gost <- gost(query = (Microglia.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)






write.xlsx(Microglia.Fst.male.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Microglia/Microglia_Fst_Male_higher.xlsx")




Microglia.Fst.male.higher.gost$result[c(1:5),] |> ggplot(aes(-log10(p_value), 
                                                        factor(term_name, levels  = c('structural constituent of synapse',
                                                                                      'structural constituent of postsynapse',
                                                                                      'structural constituent of postsynaptic specialization',
                                                                                      'PDZ domain binding',
                                                                                      'structural constituent of postsynaptic density')))) + 
  geom_col(fill = '#ff6e00', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Microglia_higher_in_fasted_males_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)






###Oligo####

#create basic clustering for sample
ARH_Sex_by_Nutr <- readRDS('data/ARH_Sex_by_Nutr.rds')
Oligo_Sex_by_Nutr <- readRDS('data/Oligo_Sex_by_Nutr.rds')

Oligo_Sex_by_Nutr <- subset(ARH_Sex_by_Nutr, subset = cell_type3 == 'Oligodendrocytes')

Oligo_Sex_by_Nutr <- SCTransform(Oligo_Sex_by_Nutr, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')


Oligo_Sex_by_Nutr@assays$SCT@var.features <- Oligo_Sex_by_Nutr@assays$SCT@var.features[(!Oligo_Sex_by_Nutr@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

Oligo_Sex_by_Nutr <- RunPCA(Oligo_Sex_by_Nutr, features = VariableFeatures(object = Oligo_Sex_by_Nutr))
Oligo_Sex_by_Nutr <- RunTSNE(Oligo_Sex_by_Nutr, reduction = "pca", dims = 1:30)

Oligo_Sex_by_Nutr <- FindNeighbors(Oligo_Sex_by_Nutr, dims = 1:30)
Oligo_Sex_by_Nutr <- FindClusters(Oligo_Sex_by_Nutr, resolution = 0.5)

DimPlot(Oligo_Sex_by_Nutr, label = TRUE, label.size = 8, repel = FALSE, reduction = 'tsne') + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())




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



saveRDS(Oligo_Sex_by_Nutr, file = 'data/Oligo_Sex_by_Nutr.rds')


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


Oligo.F.Fd.Fst.suppress.gost <- gost(query = (Oligo.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                         organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Oligo.F.Fd.Fst.induce.gost <- gost(query = (Oligo.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                       organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)



#write.xlsx(Oligo.F.Fd.Fst.suppress.gost$result[,c(9,11,3)], 
 #          file = "../paper_figures/post_2025-01-06/GOMF_Tables/Oligo/Oligo_F_Fd_Fst_suppress.xlsx")



Oligo.F.Fd.Fst.suppress.gost$result |> ggplot(aes(-log10(p_value), 
                                                                          factor(term_name, levels = c('protein binding',
                                                                                                       'cell adhesion molecule binding')))) + 
  geom_col(fill = '#6a816a', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Oligo_fasting_suppressed_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)




#write.xlsx(Oligo.F.Fd.Fst.induce.gost$result[c(1:3,5,6),c(9,11,3)], 
 #          file = "../paper_figures/post_2025-01-06/GOMF_Tables/Oligo/Oligo_F_Fd_Fst_induce.xlsx")



Oligo.F.Fd.Fst.induce.gost$result[c(1:3,5,6),] |> ggplot(aes(-log10(p_value), 
                                                                    factor(term_name, levels = c('potassium channel regulator activity',
                                                                                                 'protein domain specific binding',
                                                                                                 'transporter regulator activity',
                                                                                                 'channel regulator activity',
                                                                                                 'ion channel regulator activity')))) + 
  geom_col(fill = '#f9994f', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Oligo_fasting_induced_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)





Oligo.M.Fd.Fst.suppress.gost <- gost(query = (Oligo.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                         organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Oligo.M.Fd.Fst.induce.gost <- gost(query = (Oligo.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                       organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)



#write.xlsx(Oligo.M.Fd.Fst.suppress.gost$result[c(2,3,6,7,8),c(9,11,3)], 
 #          file = "../paper_figures/post_2025-01-06/GOMF_Tables/Oligo/Oligo_M_Fd_Fst_suppress.xlsx")



Oligo.M.Fd.Fst.suppress.gost$result[c(2,3,6,7,8),] |> ggplot(aes(-log10(p_value), 
                                                                        factor(term_name, levels = c('Rho GDP-dissociation inhibitor binding',
                                                                                                     'dATP binding',
                                                                                                     'UTP binding',
                                                                                                     'cell adhesion molecule binding',
                                                                                                     'CTP binding')))) + 
  geom_col(fill = '#19552b', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Oligo_fasting_suppressed_in_males_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)




write.xlsx(Oligo.M.Fd.Fst.induce.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Oligo/Oligo_M_Fd_Fst_induce.xlsx")



Oligo.M.Fd.Fst.induce.gost$result |> ggplot(aes(-log10(p_value), 
                                                                     factor(term_name, levels = c('voltage-gated monoatomic ion channel activity involved in regulation of postsynaptic membrane potential',
                                                                                                  'potassium channel regulator activity',
                                                                                                  'transporter regulator activity',
                                                                                                  'channel regulator activity',
                                                                                                  'ion channel regulator activity')))) + 
  geom_col(fill = '#ff6e00', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Oligo_fasting_induced_in_males_GOMF.tiff', device = 'tiff', units = 'in', width = 7, height = 1.5, dpi = 600)



####

Oligo.Fd.female.higher.gost <- gost(query = (Oligo.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                        organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)

Oligo.Fd.male.higher.gost <- gost(query = (Oligo.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                      organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)







#####

Oligo.Fst.female.higher.gost <- gost(query = (Oligo.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                         organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Oligo.Fst.male.higher.gost <- gost(query = (Oligo.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                       organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)



write.xlsx(Oligo.Fst.female.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Oligo/Oligo_Fst_Female_higher.xlsx")



Oligo.Fst.female.higher.gost$result |> ggplot(aes(-log10(p_value), 
                                                             factor(term_name, levels = c('structural constituent of synapse',
                                                                                          'protein domain specific binding',
                                                                                          'protein binding')))) + 
  geom_col(fill = '#f9994f', color = 'black') +
  theme_classic() +
  labs(y = '') +
  theme(axis.text = element_text(family = 'Arial', color = 'black',size = 8),
        axis.title = element_text(family = 'Arial', color = 'black',size = 8))
ggsave(filename = 'figures/Oligo_higher_in_fasted_females_GOMF.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.5, dpi = 600)






####POMC####



Pomc_Sex_by_Nutr <- subset(ARH_Sex_by_Nutr, subset = cell_type3 == 'Pomc')

Pomc_Sex_by_Nutr <- SCTransform(Pomc_Sex_by_Nutr, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')


Pomc_Sex_by_Nutr@assays$SCT@var.features <- Pomc_Sex_by_Nutr@assays$SCT@var.features[(!Pomc_Sex_by_Nutr@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

Pomc_Sex_by_Nutr <- RunPCA(Pomc_Sex_by_Nutr, features = VariableFeatures(object = Pomc_Sex_by_Nutr))
Pomc_Sex_by_Nutr <- RunTSNE(Pomc_Sex_by_Nutr, reduction = "pca", dims = 1:30)

Pomc_Sex_by_Nutr <- FindNeighbors(Pomc_Sex_by_Nutr, dims = 1:30)
Pomc_Sex_by_Nutr <- FindClusters(Pomc_Sex_by_Nutr, resolution = 0.5)

DimPlot(Pomc_Sex_by_Nutr, label = TRUE, label.size = 8, repel = FALSE, reduction = 'tsne') + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())




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


saveRDS(Pomc_Sex_by_Nutr, file = 'data/Pomc_Sex_by_Nutr.rds')

#Pomc.F.Fd.v.Fst$p_val_adj[1] = Pomc.F.Fd.v.Fst$p_val_adj[2]


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


Pomc.F.Fd.Fst.suppress.gost <- gost(query = (Pomc.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Pomc.F.Fd.Fst.induce.gost <- gost(query = (Pomc.F.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


write.xlsx(Pomc.F.Fd.Fst.suppress.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Pomc/Pomc_F_Fd_Fst_suppress.xlsx")
write.xlsx(Pomc.F.Fd.Fst.induce.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Pomc/Pomc_F_Fd_Fst_induce.xlsx")



Pomc.M.Fd.Fst.suppress.gost <- gost(query = (Pomc.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                    organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Pomc.M.Fd.Fst.induce.gost <- gost(query = (Pomc.M.Fd.v.Fst |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                  organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


write.xlsx(Pomc.M.Fd.Fst.suppress.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Pomc/Pomc_M_Fd_Fst_suppress.xlsx")
write.xlsx(Pomc.M.Fd.Fst.induce.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Pomc/Pomc_M_Fd_Fst_induce.xlsx")


Pomc.Fd.female.higher.gost <- gost(query = (Pomc.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                   organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Pomc.Fd.male.higher.gost <- gost(query = (Pomc.Fd.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                 organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


write.xlsx(Pomc.Fd.female.higher.gost$result[c(1:5),c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Pomc/Pomc_Fd_Female_higher.xlsx")

write.xlsx(Pomc.Fd.male.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Pomc/Pomc_Fd_Male_higher.xlsx")


Pomc.Fst.female.higher.gost <- gost(query = (Pomc.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC > 0) |> select(gene) |> as.list()) ,
                                   organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)
Pomc.Fst.male.higher.gost <- gost(query = (Pomc.Fst.F.v.M |> filter(p_val_adj < 0.05, avg_log2FC < 0) |> select(gene) |> as.list()) ,
                                 organism = 'mmusculus', sources = 'GO:MF', exclude_iea = T)


write.xlsx(Pomc.Fst.female.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Pomc/Pomc_Fst_Female_higher.xlsx")

write.xlsx(Pomc.Fst.male.higher.gost$result[,c(9,11,3)], 
           file = "../paper_figures/post_2025-04-07/GOMF_Tables/Pomc/Pomc_Fst_Male_higher.xlsx")

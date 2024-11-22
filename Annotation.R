libs <- c( 'gplots','stringi','reshape2','cowplot','RColorBrewer',
           'sctransform','stringr','org.Mm.eg.db','AnnotationDbi',
           'IRanges','S4Vectors','Biobase','BiocGenerics','clusterProfiler',
           'biomaRt','Matrix','DESeq2','RcppThread', 'extrafont', 'openxlsx',
           'Seurat','dplyr','tidyr','ggplot2','harmony','ggalluvial',
           'scDblFinder','SoupX','UpSetR','ComplexUpset','glmGamPoi')
BiocManager::install('glmGamPoi')
devtools::install_github('immunogenomics/presto')
lapply(libs, require, character.only = TRUE)


#create basic clustering for sample
ARH_Sex_by_Nutr <- readRDS('data/ARH_Sex_by_Nutr.rds')

ARH_Sex_by_Nutr <- SCTransform(ARH_Sex_by_Nutr, verbose = TRUE, do.scale = FALSE, do.center = FALSE, assay = 'RNA')

yx_chrom_genes <- c('Eif2s3y','Sry','Zfy','Rps4y1','Amely','Tbl1y','Pcdh11y','Tgif2ly',
                    'Tspy1','Tspy2','Azfa','Usp9y','Ddx3y','Uty','Tb4y','Azfb',
                    'Cyorf15','Rps4y2','Eif1ay','Kdm5d','Xkry','Hsfy1','Hsfy2',
                    'Pry','Pry2','Rbmy1a1','Azfc','Daz1','Daz2','Daz3','Daz4',
                    'Cdy1','Cdy2','Vcy1','Vcy2','Xist','Tsix'
)

mito.genes <- grep(pattern = "^mt-", x = rownames(x = ARH_Sex_by_Nutr@assays$RNA@data), value = TRUE)
hemoglobin_genes <- c('Hbq1a','Hbb-y','Hbb-bt','Hba-a2','Hba-a1')
ribo.genes <- grep(pattern = '^Rp[sl]', x = rownames(x = ARH_Sex_by_Nutr@assays$RNA@data), value = TRUE)


ARH_Sex_by_Nutr@assays$SCT@var.features <- ARH_Sex_by_Nutr@assays$SCT@var.features[(!ARH_Sex_by_Nutr@assays$SCT@var.features %in% c(yx_chrom_genes,mito.genes,hemoglobin_genes,ribo.genes))]

ARH_Sex_by_Nutr <- RunPCA(ARH_Sex_by_Nutr, features = VariableFeatures(object = ARH_Sex_by_Nutr))
ARH_Sex_by_Nutr <- RunTSNE(ARH_Sex_by_Nutr, reduction = "pca", dims = 1:30)

ARH_Sex_by_Nutr <- FindNeighbors(ARH_Sex_by_Nutr, dims = 1:30)
ARH_Sex_by_Nutr <- FindClusters(ARH_Sex_by_Nutr, resolution = 1.5)

DimPlot(ARH_Sex_by_Nutr, label = TRUE, label.size = 8, repel = FALSE, reduction = 'tsne') + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
#ggsave('post_2023-05-02_figures/dimplot_by_number.tiff', device = 'tiff', units = 'in', width = 9, height = 9, dpi = 600)
#ggsave('post_2023-05-02_figures/dimplot_by_number.svg', device = 'svg', width = 9, height = 9)

DimPlot(ARH_Sex_by_Nutr, label = FALSE, split.by = 'Batch', reduction = 'tsne') & 
  NoLegend() &
  theme(strip.text.x = element_text(hjust = 0.5,family = "Arial", face = 'plain', size = 10, color = 'black'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
#ggsave('all_samples_figures/dimplot_split_by_batch', device = 'tiff', units = 'in', width = 7, height = 3, dpi = 600)
#ggsave('post_2023-05-02_figures/dimplot_by_batch.svg', device = 'svg', width = 9, height = 3)





DimPlot(ARH_Sex_by_Nutr, label = TRUE, label.size = 2,pt.size = 0.01 , 
        shuffle = TRUE, group.by = 'Batch') + 
  NoLegend() +
  theme(title = element_text(family = 'Arial', size = 8, color = 'black'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
#ggsave('all_samples_figures/dimplot_group_by_batch', device = 'tiff', units = 'in', width = 3, height = 3, dpi = 600)


DimPlot(ARH_Sex_by_Nutr, label = TRUE, label.size = 2, split.by = 'Sex') + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
#ggsave('all_samples_figures/dimplot_split_by_sex', device = 'tiff', units = 'in', width = 5, height = 3, dpi = 600)


DimPlot(ARH_Sex_by_Nutr, label = TRUE, label.size = 2,pt.size = 0.01, repel = TRUE, 
        shuffle = TRUE, group.by = 'Sex') + 
  NoLegend() +
  theme(title = element_text(family = 'Arial', size = 8, color = 'black'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
#ggsave('all_samples_figures/dimplot_group_by_sex', device = 'tiff', units = 'in', width = 3, height = 3, dpi = 600)


DimPlot(ARH_Sex_by_Nutr, label = TRUE, label.size = 2, split.by = 'Nutr_State') + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
#ggsave('all_samples_figures/dimplot_split_by_nutr', device = 'tiff', units = 'in', width = 5, height = 3, dpi = 600)


DimPlot(ARH_Sex_by_Nutr, label = TRUE, label.size = 2,pt.size = 0.01, repel = TRUE, 
        shuffle = TRUE, group.by = 'Nutr_State') + 
  NoLegend() +
  theme(title = element_text(family = 'Arial', size = 8, color = 'black'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
#ggsave('all_samples_figures/dimplot_group_by_nutr', device = 'tiff', units = 'in', width = 3, height = 3, dpi = 600)


DimPlot(ARH_Sex_by_Nutr,split.by = 'sexXnutr', ncol = 2) + 
  NoLegend() +
  theme(strip.text.x = element_text(hjust = 0.5,family = "Arial", face = 'plain', size = 12, color = 'black'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
#ggsave('all_samples_figures/dimplot_split_by_sexXnutr', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)
#ggsave('post_2023-05-02_figures/dimplot_by_sexXnutr.svg', device = 'svg', width = 6, height = 6)







DimPlot(ARH_Sex_by_Nutr, label = TRUE, label.size = 2,pt.size = 0.01, repel = TRUE, 
        shuffle = TRUE, group.by = 'sexXnutr') + 
  theme(title = element_text(family = 'Arial', size = 8, color = 'black'),
        legend.text = element_text(family = 'Arial', size = 8, color = 'black'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
#ggsave('all_samples_figures/dimplot_group_by_sexXnutr', device = 'tiff', units = 'in', width = 4, height = 3, dpi = 600)



FeaturePlot(ARH_Sex_by_Nutr, features = c('Agrp','Pomc','Tac2','Ghrh'), reduction = 'tsne')
#ggsave('all_samples_figures/featureplot_arh_neurons', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

FeaturePlot(ARH_Sex_by_Nutr, features = c('Snap25','Slc1a3','Ptprc','Mbp'))
#ggsave('all_samples_figures/featureplot_arh_cells', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)


FeaturePlot(ARH_Sex_by_Nutr, features = c('Slc17a6','Gad2','Th','Hdc'))
#ggsave('all_samples_figures/featureplot_neurotransmitters', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)

all_sample_markers <- FindAllMarkers(ARH_Sex_by_Nutr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log2(1.25), pseudocount.use = 0)


#saveRDS(ARH_Sex_by_Nutr, file = 'data/ARH_Sex_by_Nutr.rds')
####2024-04-28####

######Cell Types######

neuro.genes <- c('Snap25','Rbfox3','Syp','Dlg4','Slc17a6','Gad1','Gad2','Slc32a1','Th','Slc6a3','Slc18a2','Hdc','Chat','Chga')

neuroPercent <- Matrix::colSums(ARH_Sex_by_Nutr@assays$RNA@counts[c(neuro.genes), ])/
  Matrix::colSums(ARH_Sex_by_Nutr@assays$RNA@counts)*100


ARH_Sex_by_Nutr$neuroPercent <- neuroPercent

FeaturePlot(ARH_Sex_by_Nutr, features ='neuroPercent', cols = c('lightgrey','red','purple'), max.cutoff = 0.3)
DimPlot(ARH_Sex_by_Nutr, label = TRUE) + NoLegend()


astro.genes <- c('Gfap','S100b','Aldh1l1','Aldoc','Slc1a2','Slc1a3')

astroPercent <- Matrix::colSums(ARH_Sex_by_Nutr@assays$RNA@counts[c(astro.genes), ])/
  Matrix::colSums(ARH_Sex_by_Nutr@assays$RNA@counts)*100


ARH_Sex_by_Nutr$astroPercent <- astroPercent


oligo.genes <- c('Mbp','Mal','Olig1','Olig2','Mog','Cnp')

oligoPercent <- Matrix::colSums(ARH_Sex_by_Nutr@assays$RNA@counts[c(oligo.genes), ])/
  Matrix::colSums(ARH_Sex_by_Nutr@assays$RNA@counts)*100


ARH_Sex_by_Nutr$oligoPercent <- oligoPercent




micro.genes <- c('Cx3cr1','Tmem119','Aif1','Itgam','Cd68','Ptprc')

microPercent <- Matrix::colSums(ARH_Sex_by_Nutr@assays$RNA@counts[c(micro.genes), ])/
  Matrix::colSums(ARH_Sex_by_Nutr@assays$RNA@counts)*100


ARH_Sex_by_Nutr$microPercent <- microPercent


epend.genes <- c('Foxj1','Meig1','Fam183b','Tmem107','Tppp3',"Pltp",'Mt3')

ependPercent <- Matrix::colSums(ARH_Sex_by_Nutr@assays$RNA@counts[c(epend.genes), ])/
  Matrix::colSums(ARH_Sex_by_Nutr@assays$RNA@counts)*100


ARH_Sex_by_Nutr$ependPercent <- ependPercent


FeaturePlot(ARH_Sex_by_Nutr, features =c('neuroPercent','astroPercent','oligoPercent','microPercent','ependPercent'), cols = c('lightgrey','red','purple'), max.cutoff = 0.2)
#ggsave('all_samples_figures/cell_type_percent.tiff', device = 'tiff',units = 'in', width = 7, height = 7, dpi = 600)



#Tanycytes
FeaturePlot(ARH_Sex_by_Nutr, features = c('Ppp1r1b','Vim','Dio2','Rax','Slc16a2'))

#tany alpha
FeaturePlot(ARH_Sex_by_Nutr, features = c('Cd59a','Slc17a8','Crym','Vcan'))
#tany beta
FeaturePlot(ARH_Sex_by_Nutr, features = c('Col25a1','Cacna2d2'))

#OPC
FeaturePlot(ARH_Sex_by_Nutr, features = c('Pdgfra','Cspg4'))
#Endothe
FeaturePlot(ARH_Sex_by_Nutr, features = c('Cldn5','Adgrf5'))
#Epend
FeaturePlot(ARH_Sex_by_Nutr, features = c('Ccdc153','Pltp','Meig1','Foxj1'))
FeaturePlot(ARH_Sex_by_Nutr, features = epend.genes)
#vascular and leptomeningeal cells (VLMC)
FeaturePlot(ARH_Sex_by_Nutr, features = c('Slc6a13','Kcnj8','Dcn'))
#####



sn_cell_types <- ARH_Sex_by_Nutr@active.ident |> as.data.frame()
sn_cell_types <- rename(sn_cell_types, clusters = 'ARH_Sex_by_Nutr@active.ident')


sn_cell_types <- sn_cell_types |> 
  mutate(cell_type = if_else((clusters == 2 | clusters == 4 | clusters == 5 | clusters == 6 |
                                clusters == 7 |  clusters == 8 | clusters == 9 |  
                                clusters == 10 | clusters == 11 | clusters == 13 |clusters == 14 | 
                                        clusters == 15 |clusters == 16 |clusters == 17 | clusters == 18 | 
                                        clusters == 19 | clusters == 21 |
                                        clusters == 22 | clusters == 23 | clusters == 24 | clusters == 25 | 
                                 clusters == 27 | clusters == 32), 'Neurons',
                                                             if_else((clusters == 29), 'Pars Tuberalis',
                                                                     if_else((clusters == 0),'Astrocytes',
                                                                             if_else((clusters == 1), 'α-Tanycytes',
                                                                                     if_else((clusters == 3), 'β-Tanycytes',
                                                                                             if_else((clusters == 12),'Oligodendrocytes',
                                                                                                     if_else((clusters == 33), 'TOP',
                                                                                                          if_else((clusters == 28), 'Microglia',
                                                                                                             if_else((clusters == 26), 'Ependymal',
                                                                                                                     if_else((clusters == 20), 'OPC',
                                                                                                                             if_else((clusters == 30), 'VLMC',
                                                                                                                                     if_else((clusters == 31), 'Endothelial', 'Unidentified'
                                                                                                                                             )
                                                                                                                                     )
                                                                                                                             )
                                                                                                                     )
                                                                                                             )
                                                                                                     )
                                                                                             )
                                                                                     )
                                                                             )
                                                                     )
                                                             )))


ARH_Sex_by_Nutr$cell_type <- sn_cell_types$cell_type



######


sn_conditions <- ARH_Sex_by_Nutr@meta.data[,c(1:4)] |> as.data.frame()
sn_conditions <- sn_conditions |> 
  mutate(F_Fed = (if_else(sexXnutr == "F_Fed", 1, 0)))

sn_conditions <- sn_conditions |> 
  mutate(F_Fast = (if_else(sexXnutr == "F_Fast", 1, 0)))

sn_conditions <- sn_conditions |> 
  mutate(M_Fed = (if_else(sexXnutr == "M_Fed", 1, 0)))

sn_conditions <- sn_conditions |> 
  mutate(M_Fast = (if_else(sexXnutr == "M_Fast", 1, 0)))


ARH_Sex_by_Nutr$F_Fed <- sn_conditions$F_Fed
ARH_Sex_by_Nutr$F_Fast <- sn_conditions$F_Fast
ARH_Sex_by_Nutr$M_Fed <- sn_conditions$M_Fed
ARH_Sex_by_Nutr$M_Fast <- sn_conditions$M_Fast




FeaturePlot(ARH_Sex_by_Nutr, label = FALSE, 
            features = 'F_Fed', reduction = 'tsne',
            pt.size = 0.01, order = T,
        cols = c('lightgrey','#6a816a')) + 
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


#ggsave('figures/f_fed_dimplot.tiff', device = 'tiff', units = 'in', width = 5,height = 5,dpi = 600)


FeaturePlot(ARH_Sex_by_Nutr, label = FALSE, 
            features = 'F_Fast', reduction = 'tsne', 
            pt.size = 0.01, order = TRUE,
        cols = c('lightgrey','#f9994f')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


#ggsave('figures/f_fast_dimplot.tiff', device = 'tiff', units = 'in', width = 5,height = 5,dpi = 600)

FeaturePlot(ARH_Sex_by_Nutr, label = FALSE, 
            features = 'M_Fed', reduction = 'tsne',
            pt.size = 0.01, order = TRUE,
        cols = c('lightgrey','#19552b')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


#ggsave('figures/m_fed_dimplot.tiff', device = 'tiff', units = 'in', width = 5,height = 5,dpi = 600)


FeaturePlot(ARH_Sex_by_Nutr, label = FALSE, 
            features = 'M_Fast', reduction = 'tsne',
            pt.size = 0.01, order = TRUE,
        cols = c('lightgrey','#ff6e00')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


#ggsave('figures/m_fast_dimplot.tiff', device = 'tiff', units = 'in', width = 5,height = 5,dpi = 600)




sn_batches <- ARH_Sex_by_Nutr@meta.data[,c(1,5)] |> as.data.frame()
sn_batches <- sn_batches |> 
  mutate(B1 = (if_else(Batch == "B1", 1, 0)))

sn_batches <- sn_batches |> 
  mutate(B2 = (if_else(Batch == "B2", 1, 0)))

sn_batches <- sn_batches |> 
  mutate(B3 = (if_else(Batch == "B3", 1, 0)))



ARH_Sex_by_Nutr$B1 <- sn_batches$B1
ARH_Sex_by_Nutr$B2 <- sn_batches$B2
ARH_Sex_by_Nutr$B3 <- sn_batches$B3



FeaturePlot(ARH_Sex_by_Nutr, label = FALSE, 
            features = 'B1', reduction = 'tsne', 
            pt.size = 0.01, order = TRUE,
            cols = c('lightgrey','#FF5733')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


#ggsave('figures/B1_dimplot.tiff', device = 'tiff', units = 'in', width = 5,height = 5,dpi = 600)



FeaturePlot(ARH_Sex_by_Nutr, label = FALSE, 
            features = 'B2', reduction = 'tsne',
            pt.size = 0.01, order = TRUE,
            cols = c('lightgrey','#C70039')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


#ggsave('figures/B2_dimplot.tiff', device = 'tiff', units = 'in', width = 5,height = 5,dpi = 600)


FeaturePlot(ARH_Sex_by_Nutr, label = FALSE, 
            features = 'B3', reduction = 'tsne',
            pt.size = 0.01, order = TRUE,
            cols = c('lightgrey','#900C3F')) +
  NoLegend() + 
  labs(title = "", x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10))


#ggsave('figures/B3_dimplot.tiff', device = 'tiff', units = 'in', width = 5,height = 5,dpi = 600)

#####





DotPlot(ARH_Sex_by_Nutr, features = (c(neuro.genes[c(1,2,5,11)],
                                                                
                                                                'Vipr2','Pomc','Agrp','Tac2',neuro.genes[c(8,10)],'Ghrh',neuro.genes[9],'Sst',
                                                                'Lepr',
                                                                'Asb4','Irs4','Ecel1','Tbx3','Nr5a2',
                                                                'Nr5a1','Fezf1','Nr2f2','Gda','Slit3',
                                                                'Tac1','Nos1','Lmx1b','Hdc',
                                                                'Postn','Adarb2','Tanc1',
                                                                'Rgs6',
                                                                'Adcy8','Avp','Oxt',
                                                                astro.genes[c(1,4,5)],
                                                                'Vim','Dio2','Rax','Slc16a2','Crym','Vcan','Col25a1',
                                                                oligo.genes[c(1,2,5,3)],'Cspg4','Pdgfra',
                                                                micro.genes[c(1,4,6)],
                                                                
                                                                'Ccdc153','Pltp','Foxj1',
                                                                'Tbx18','Slc6a13','Dcn','Foxc1','Cldn5','Adgrf5','Chga','Timeless','Tshb','Cck')), 
        scale = FALSE, scale.by = 'radius') &
  #coord_flip() &
  scale_y_discrete(limits = factor(c(
    32,
    22,
    10,
    18,14,4,
    23,19,7,
    27,25,24,23,21,17,16,15,14,
    13,11,10,9,8,7,6,5,2))) &
  scale_radius(limits = c(5,100)) &
  scale_color_continuous(limits = c(0,NA), type = "viridis") &
  theme_bw() &
  #labels = c(31,28,33,20,4,1,29,10,35,0,30,22,32,17,21,2,3,8,11,12,26,25,9,16,23,19,34,27,18,15,6,24,14,13,7,5)) &
  labs(x = '', y = '') &
  theme(title = element_text(family = 'Arial', size = 12, color = 'black'),
        axis.title=element_text(family = 'Arial', size = 12, color = 'black'),
        axis.text.y=element_text(family = 'Arial', size = 12, color = 'black'),
        axis.text.x=element_text(family = 'Arial', size = 12, color = 'black',angle = 90, hjust = 1, vjust =0.5)) &
  
  NoLegend()
ggsave('post_2023-06-02_figures/dotplot_celltypes_markers3.svg', device = 'svg', units = 'in', width = 8.48, height = 4, dpi = 600)






sn_cell_types <- ARH_Sex_by_Nutr@active.ident |> as.data.frame()
sn_cell_types <- rename(sn_cell_types, clusters = 'ARH_Sex_by_Nutr@active.ident')




sn_cell_types2 <- sn_cell_types |> 
  mutate(cell_type2 = if_else((clusters == 2 | clusters == 5 | clusters == 6 | clusters == 8 |
                                clusters == 9 |  clusters == 11 | clusters == 13 |  
                                clusters == 15 | clusters == 16 | clusters == 17 | clusters == 21 | 
                                clusters == 24 | clusters == 25 | clusters == 27), 'ARH Neurons',
                             if_else((clusters == 23 | clusters == 19 |clusters == 7), 'PVp Neurons',
                                     if_else((clusters == 10),'VMH Neurons',
                                             if_else((clusters == 22), 'SCN Neurons',
                                                     if_else((clusters == 32), 'TU Neurons',
                                                             if_else((clusters == 29), 'Pars Tuberalis',
                                                                     if_else((clusters == 0),'Astrocytes',
                                                                             if_else((clusters == 1), 'α-Tanycytes',
                                                                                     if_else((clusters == 3), 'β-Tanycytes',
                                                                                             if_else((clusters == 12),'Oligodendrocytes',
                                                                                                     if_else((clusters == 28), 'Microglia',
                                                                                                             if_else((clusters == 26), 'Ependymal',
                                                                                                                     if_else((clusters == 20), 'OPC',
                                                                                                                             if_else((clusters == 30), 'VLMC',
                                                                                                                                     if_else((clusters == 31), 'Endothelial', 
                                                                                                                                             if_else((clusters == 4 | clusters == 14 | clusters == 18), 'MM Neurons', 
                                                                                                                                                     if_else((clusters == 33),'TOP' ,'unidentified'
                                                                                                                                             )
                                                                                                                                     )
                                                                                                                             )
                                                                                                                     )
                                                                                                             )
                                                                                                     )
                                                                                             )
                                                                                     )
                                                                             )
                                                                     )
                                                             ))))))))


ARH_Sex_by_Nutr$cell_type2 <- sn_cell_types2$cell_type2


DimPlot(ARH_Sex_by_Nutr, label = T) & NoLegend()
DimPlot(ARH_Sex_by_Nutr, group.by = 'cell_type2', label = T, label.size = 3, repel = FALSE) + 
  NoLegend() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
#ggsave('post_2023-05-02_figures/dimplot_by_number.tiff', device = 'tiff', units = 'in', width = 9, height = 9, dpi = 600)
#ggsave('post_2023-05-02_figures/dimplot_by_number.svg', device = 'svg', width = 9, height = 9)




#saveRDS(ARH_Sex_by_Nutr, file = 'data/ARH_Sex_by_Nutr.rds')





#####2024-04-29######



campbell_clust_markers <- read.csv('../../scARC_Lowell/seurat/campbell_clust_markers2.csv', row.names = 1)
bean_clust_markers <- all_sample_markers |> filter(cluster != 3, cluster != 4, cluster != 5, cluster != 10, cluster != 12, 
                                                   cluster != 13, cluster != 21, 
                                           cluster != 26, cluster != 31, cluster != 36, cluster != 38,
                                           cluster != 40, cluster != 41, cluster != 43)

campbell_by_bean <- tibble(Bean_clust = 'a', Campbell_clust = 'b', 
                           Overlap = 0, Bean_count = 0, Campbell_count = 0, 
                           Bean_port = 0.00, Campbell_port = 0.00, .rows = 1054)
rownum = 0

for(i in unique(bean_clust_markers$cluster)){
  for(j in unique(campbell_clust_markers$cluster)){
    
    rownum = rownum + 1
    
    campbell_by_bean$Bean_clust[rownum] = i
    campbell_by_bean$Campbell_clust[rownum] = j
    
    campbell_by_bean$Overlap[rownum] = inner_join(
      bean_clust_markers[(bean_clust_markers$cluster == i & 
                            bean_clust_markers$pct.1 > 0.25 & 
                            bean_clust_markers$avg_log2FC > log2(1.25) & 
                            bean_clust_markers$p_val_adj < 0.05),], 
      campbell_clust_markers[(campbell_clust_markers$cluster == j & 
                                campbell_clust_markers$pct.1 > 0.25 & 
                                campbell_clust_markers$avg_log2FC > log2(1.25) & 
                                campbell_clust_markers$p_val_adj < 0.05),], 
      by = 'gene') |> 
      count()
    
    campbell_by_bean$Bean_count[rownum] = bean_clust_markers[(bean_clust_markers$cluster == i & 
                                                                bean_clust_markers$pct.1 > 0.25 & 
                                                                bean_clust_markers$avg_log2FC > log2(1.25) & 
                                                                bean_clust_markers$p_val_adj < 0.05),] |> count()
    campbell_by_bean$Campbell_count[rownum] = campbell_clust_markers[(campbell_clust_markers$cluster == j & 
                                                                        campbell_clust_markers$pct.1 > 0.25 & 
                                                                        campbell_clust_markers$avg_log2FC > log2(1.25) & 
                                                                        campbell_clust_markers$p_val_adj < 0.05),] |> count()
    
    campbell_by_bean$Bean_port[rownum] = as.numeric(campbell_by_bean$Overlap[rownum]) / as.numeric(campbell_by_bean$Bean_count[rownum])
    campbell_by_bean$Campbell_port[rownum] = as.numeric(campbell_by_bean$Overlap[rownum]) / as.numeric(campbell_by_bean$Campbell_count[rownum])
    
    print(rownum/1054*100)
    
  }
}



####
#scale_y_discrete(limits = c('Nfia/Nfib','Tafa4','Vgll3','ND','Lhx6','Gm33791','Zeb2','Ntng2','Trhr','Pdzrn3','Tbx19',
 #                           'Htr3b/Vgll3','Tac1/Reln','Lef1','Col25a1','Calcrl','Sst','Lamp5','Sgpp2','Tac2','Pomc','Agrp','Itpr2','Ghrh/Chat',
  #                          'Slc18a2/Dach2','Th/Slc6a3','Slc18a2/Satb2','Dach2','Tac1/Npsr1')) +
 # scale_x_discrete(limits = c('n02','n04','n07','n08','n09','n03','n19','n05','n10','n11','n12','n13','n14','n15',
   #                           'n20','n21','n22','n23','n24','n25','n26','n27','n32','n33','n34','n28'),
    #               labels = c('n02.Gm8773/Tac1','n04.Sst/Nts','n07.Arx/Nr5a2',
     #                         'n08.Th/Slc6a3','n09.Th/Nfib','n03.Th/Sst','n19.Gpr50','n05.Nfix/Nfib','n10.Ghrh','n11.Trh/Cxcl12','n12.Agrp/Sst',
      #                        'n13.Agrp/Gm8873','n14.Pomc/Ttr','n15.Pomc/Anxa2','n20.Kiss1/Tac2',
       #                       'n21.Pomc/Glipr1','n22.Tmem215','n23.Sst/Unc13c','n24.Sst/Pthlh',
        #                      'n25.Trh/Lef1','n26.Htr3b','n27.Tbx19','n32.Slc17a6/Trhr',
  #                            'n33.unassigned','n34.unassifned','n28.Qrfp')) +
####


campbell_by_bean <- campbell_by_bean |> as.data.frame()

  

ggplot(campbell_by_bean, aes(Campbell_clust, as.factor(Bean_clust) )) + 
  
  scale_y_discrete(limits = as.factor(c(44,6,
                                        39,33,30,29,27,
                                        24,18,17,14,
                                        1,7,32,25,42,
                                        23,11,16,19,20,
                                        9,8,15,2,0,
                                        28,22,34,35,37)),
                   labels = c(
                     'Tbx15','Gad2/Htr2c',
                     'Tac1/Reln','Ros1/Alk','Ebf3/Htr2c','Klhl1/Ebf3','MM.03',
                     'PVp.04','PVp.02','MM.02','PVp.01',
                     'MM.01','VMH.01','VMH.02','SCN','Tu',
                     'PVp.03','Slc17a6/Alk','Tbx19','Htr3b','Lef1',
                     'Sst/Unc13c','Lamp5/Npy5r','KNDy','Pomc','Agrp',
                     'Erg/Lepr','Ghrh/Chat','DA','Satb2/Slc18a2','Coch/Slc18a2'
                   )) +
  scale_x_discrete(limits = c('n03','n07','n08','n09','n19','n10','n11','n12','n13','n14','n15',
                              'n20','n21','n22','n04','n23','n24','n25','n26','n27','n32',
                              'n01','n02','n31','n06','n16','n17','n18','n29','n30','n33','n34','n05','n28'),
                   labels = c('n03.Th/Sst','n07.Arx/Nr5a2',
                              'n08.Th/Slc6a3','n09.Th/Nfib','n19.Gpr50','n10.Ghrh','n11.Trh/Cxcl12','n12.Agrp/Sst',
                              'n13.Agrp/Gm8873','n14.Pomc/Ttr','n15.Pomc/Anxa2','n20.Kiss1/Tac2',
                              'n21.Pomc/Glipr1','n22.Tmem215','n04.Sst/Nts','n23.Sst/Unc13c','n24.Sst/Pthlh',
                              'n25.Trh/Lef1','n26.Htr3b','n27.Tbx19','n32.Slc17a6/Trhr',
                              
                              'n01.Hdc','n02.Gm8773/Tac1','n31.Slc17a6/Fam19a2','n06.Oxt','n16.Rgs16/Vip','n17.Rgs16/Dlx1',
                              'n18.Rgs16/Nmu','n29.Nr5a1/Bdnf','n30.Nr5a1/Nfib',
                              'n33.unassigned','n34.unassifned','n05.Nfix/Nfib','n28.Qrfp')) +
  
  geom_bin_2d(aes(fill = pmin(pmax(Campbell_port, .15), .35))) + 
  scale_fill_continuous(type = "viridis") +

  theme_classic() +
  labs(x= 'Campbell et al. Clusters', y = 'Current Study Clusters') +
  theme(axis.title = element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.x =  element_text(family = 'Arial', size = 8, color = 'black', angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(family = 'Arial', size = 8, color = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(family = 'Arial', size = 8, color = 'black'))

#ggsave(filename = 'figures/correlation_campbell_by_bean_clust.svg', device = 'svg', width = 7.24, height = 7, dpi = 600)

####2024-04-30#####

##########'Asb4','Irs4','Ecel1','Tbx3','Nr5a2',
DotPlot(ARH_Sex_by_Nutr, features = rev(c(neuro.genes[c(1,2,5,8,10,11)],
                                       
                                       'Agrp','Pomc','Glipr1','Tac2','Pdyn','Kiss1','Ghrh','Th','Satb2','Coch','Sst',
                                       'Htr2c','Lamp5','Alk','Tbx19','Htr3b','Lef1','Erg',
                                       'Lepr','Klhl1','Ebf3','Ros1','Reln','Tbx15',
                                       'Vipr2',
                                       'Nr5a1','Fezf1','Nr2f2','Gda','Slit3',
                                       'Tac1','Nos1','Lmx1b','Hdc',
                                       'Postn','Adarb2','Tanc1',
                                       'Rgs6',
                                       'Adcy8','Avp','Oxt',
                                       astro.genes[c(1,4,5)],
                                       'Vim','Dio2','Rax','Slc16a2','Crym','Vcan','Col25a1',
                                       oligo.genes[c(1,2,5,3)],'Cspg4','Pdgfra',
                                       micro.genes[c(1,4,6)],
                                       
                                       'Ccdc153','Pltp','Foxj1',
                                       'Tbx18','Slc6a13','Dcn','Foxc1','Cldn5','Adgrf5','Chga','Timeless','Tshb','Cck')), 
        scale = FALSE, scale.by = 'radius') &
  coord_flip() &
  scale_y_discrete(limits = factor(c(0,2,15,22,34,35,37, 9,
                                     6,
                                     8,
                                     11,
                                     16,
                                     19,
                                     20,
                                     28,
                                     29,
                                     30,
                                     33,
                                     39,
                                     44,
                                     25,
                                     7,
                                     32,
                                     14,
                                     18,
                                     23,
                                     24,
                                     1,
                                     17,
                                     27,
                                     42,
                                     4,
                                     5,
                                     12,
                                     13,
                                     3,
                                     26,
                                     31,
                                     10,
                                     43,
                                     21,
                                     36,
                                     40,
                                     41,
                                     38)),
                   labels = c(
                     'Agrp',
                     'Pomc',
                     'KNDy',
                     'Ghrh/Chat',
                     'DA',
                     'Satb2/Slc18a2',
                     'Coch/Slc18a2',
                     'Sst/Unc13c',
                     'Gad2/Htr2c',
                     'Lamp5/Npy5r',
                     'Slc17a6/Alk',
                     'Tbx19',
                     'Htr3b',
                     'Lef1',
                     'Erg/Lepr',
                     'Klhl1/Ebf3',
                     'Ebf3/Htr2c',
                     'Ros1/Alk',
                     'Tac1/Reln',
                     'Tbx15',
                     'SCN',
                     'VMH.01',
                     'VMH.02',
                     'PVP.01',
                     'PVP.02',
                     'PVP.03',
                     'PVP.04',
                     'MM.01',
                     'MM.02',
                     'MM.03',
                     'Tu',
                     'Astrocytes',
                     'Astrocytes',
                     'Alpha-Tanycytes',
                     'Alpha-Tanycytes',
                     'Beta-Tanycytes',
                     'Beta-Tanycytes',
                     'Ependymal',
                     'Oligodendrocytes',
                     'TOP',
                     'OPC',
                     'Microglia',
                     'VLMC',
                     'Endothelial',
                     'Pars Tuberalis'
                     
                   )) &
  scale_radius(limits = c(5,100)) &
  scale_color_continuous(limits = c(0,NA), type = "viridis") &
  #labels = c(31,28,33,20,4,1,29,10,35,0,30,22,32,17,21,2,3,8,11,12,26,25,9,16,23,19,34,27,18,15,6,24,14,13,7,5)) &
  labs(x = 'Marker Genes', y = 'Cell-Type') &
  theme(title = element_text(family = 'Arial', size = 8, color = 'black'),
        axis.title=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.y=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.x=element_text(family = 'Arial', size = 8, color = 'black',angle = 90, hjust = 1, vjust =0.5)) 
#ggsave('figures/dotplot_celltypes_markers.svg', device = 'svg', units = 'in', width = 7.24, height = 10, dpi = 600)







sn_cell_types <- ARH_Sex_by_Nutr@active.ident |> as.data.frame()
sn_cell_types <- rename(sn_cell_types, cluster = 'ARH_Sex_by_Nutr@active.ident')




sn_cell_types3 <- sn_cell_types |> 
  mutate(cell_type3 = if_else(cluster == 0, 'Agrp',
                              if_else(cluster == 2, 'Pomc',
                                      if_else(cluster == 15, 'KNDy',
                                              if_else(cluster == 22, 'Ghrh/Chat',
                                                      if_else(cluster == 34, 'DA',
                                                              if_else(cluster == 35, 'Satb2/Slc18a2',
                                                                      if_else(cluster==37, 'Coch/Slc18a2',
                                                                              if_else(cluster == 9, 'Sst/Unc13c',
                                                                                      if_else(cluster == 6, 'Gad2/Htr2c',
                                                                                              if_else(cluster == 8, 'Lamp5/Npy5r',
                                                                                                      if_else(cluster == 11, 'Slc17a6/Alk',
                                                                                                              if_else(cluster == 16, 'Tbx19',
                                                                                                                      if_else(cluster==19, 'Htr3b',
                                                                                                                              if_else(cluster == 20, 'Lef1',
                                                                                                                                      if_else(cluster == 28, 'Erg/Lepr',
                                                                                                                                              if_else(cluster == 29, 'Klhl1/Ebf3',
                                                                                                                                                      if_else(cluster == 30, 'Ebf3/Htr2c',
                                                                                                                                                              if_else(cluster == 33, 'Ros1/Alk',
                                                                                                                                                                      if_else(cluster == 39, 'Tac1/Reln',
                                                                                                                                                                              if_else(cluster == 44, 'Tbx15',
                                                                                                                                                                                      if_else(cluster == 25, 'SCN',
                                                                                                                                                                                              if_else(cluster == 7, 'VMH.01',
                                                                                                                                                                                                      if_else(cluster == 32, 'VMH.02',
                                                                                                                                                                                                              if_else(cluster == 14, 'PVp.01',
                                                                                                                                                                                                                      if_else(cluster== 18, 'PVp.02',
                                                                                                                                                                                                                              if_else(cluster == 23, 'PVp.03',
                                                                                                                                                                                                                                      if_else(cluster == 24, 'PVp.04',
                                                                                                                                                                                                                                              if_else(cluster == 1, 'MM.01',
                                                                                                                                                                                                                                                      if_else(cluster == 17, 'MM.02',
                                                                                                                                                                                                                                                              if_else(cluster== 27, 'MM.03',
                                                                                                                                                                                                                                                                      if_else(cluster == 42, 'Tu',
                                                                                                                                                                                                                                                                              if_else(cluster == 4 | cluster == 5, 'Astrocytes',
                                                                                                                                                                                                                                                                                      if_else(cluster == 12 | cluster == 13, 'α-Tanycytes',
                                                                                                                                                                                                                                                                                              if_else(cluster == 3 | cluster == 26, 'β-Tanycytes',
                                                                                                                                                                                                                                                                                                      if_else(cluster == 31, 'Ependymal',
                                                                                                                                                                                                                                                                                                              if_else(cluster == 10 , 'Oligodendrocytes',
                                                                                                                                                                                                                                                                                                                      if_else(cluster == 43, 'TOP',
                                                                                                                                                                                                                                                                                                                              if_else(cluster == 21, 'OPC',
                                                                                                                                                                                                                                                                                                                                      if_else(cluster == 36, 'Microglia',
                                                                                                                                                                                                                                                                                                                                              if_else(cluster == 40, 'VLMC',
                                                                                                                                                                                                                                                                                                                                                      if_else(cluster == 41, 'Endothelial',
                                                                                                                                                                                                                                                                                                                                                              if_else(cluster == 38, 'Pars Tuberalis', 'Missed'))))))))))))))))))))))))))
                                                                                                                                                      )))))))))))))))))








ARH_Sex_by_Nutr$cell_type3 <- sn_cell_types3$cell_type3







#saveRDS(ARH_Sex_by_Nutr, file = 'data/ARH_Sex_by_Nutr.rds')

######2024-05-01######




DimPlot(ARH_Sex_by_Nutr, label = T, group.by = 'cell_type3', 
        repel = T, shuffle = F,
        label.size = 12/.pt
        ) &
  
  labs(title = "", x = '', y = '') &
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10)) & NoLegend()


ggsave('figures/cell_type_dimplot3.tiff', device = 'tiff', units = 'in', width = 7,height = 7,dpi = 600)
ggsave('figures/cell_type_dimplot.svg', device = 'svg', units = 'in', width = 7,height = 7,dpi = 600)


#####2024-05-01######



DotPlot(ARH_Sex_by_Nutr,group.by = 'cell_type3', features = rev(c(neuro.genes[c(1,2,5,8,10,11)], 'Ghr','Esr1','Prlr','Pgr',
                                          
                                          'Agrp','Pomc','Glipr1','Tac2','Pdyn','Kiss1','Ghrh','Th','Satb2','Coch','Sst',
                                          'Htr2c','Lamp5','Alk','Tbx19','Htr3b','Lef1','Erg',
                                          'Lepr','Klhl1','Ebf3','Ros1','Reln','Tbx15',
                                          'Vipr2',
                                          'Nr5a1','Fezf1','Nr2f2','Gda','Slit3',
                                          'Tac1','Nos1','Lmx1b','Hdc',
                                          'Postn','Adarb2','Tanc1',
                                          'Rgs6',
                                          'Adcy8','Avp','Oxt',
                                          astro.genes[c(1,4,5)],
                                          'Vim','Dio2','Rax','Slc16a2','Crym','Vcan','Col25a1',
                                          oligo.genes[c(1,2,5,3)],'Cspg4','Pdgfra',
                                          micro.genes[c(1,4,6)],
                                          
                                          'Ccdc153','Pltp','Foxj1',
                                          'Tbx18','Slc6a13','Dcn','Foxc1','Cldn5','Adgrf5','Chga','Timeless','Tshb','Cck')), 
        scale = FALSE, scale.by = 'radius') &
  coord_flip() &
  scale_y_discrete(limits = c(
                     'Agrp',
                     'Pomc',
                     'KNDy',
                     'Ghrh/Chat',
                     'DA',
                     'Satb2/Slc18a2',
                     'Coch/Slc18a2',
                     'Sst/Unc13c',
                     'Gad2/Htr2c',
                     'Lamp5/Npy5r',
                     'Slc17a6/Alk',
                     'Tbx19',
                     'Htr3b',
                     'Lef1',
                     'Erg/Lepr',
                     'Klhl1/Ebf3',
                     'Ebf3/Htr2c',
                     'Ros1/Alk',
                     'Tac1/Reln',
                     'Tbx15',
                     'SCN',
                     'VMH.01',
                     'VMH.02',
                     'PVp.01',
                     'PVp.02',
                     'PVp.03',
                     'PVp.04',
                     'MM.01',
                     'MM.02',
                     'MM.03',
                     'Tu',
                     'Astrocytes',
                     'α-Tanycytes',
                     'β-Tanycytes',
                     'Ependymal',
                     'Oligodendrocytes',
                     'TOP',
                     'OPC',
                     'Microglia',
                     'VLMC',
                     'Endothelial',
                     'Pars Tuberalis'
                     
                   )) &
  scale_radius(limits = c(5,100)) &
  scale_color_continuous(limits = c(0,NA), type = "viridis") &
  #labels = c(31,28,33,20,4,1,29,10,35,0,30,22,32,17,21,2,3,8,11,12,26,25,9,16,23,19,34,27,18,15,6,24,14,13,7,5)) &
  labs(x = 'Marker Genes', y = 'Cell-Type') &
  theme_classic() &
  theme(title = element_text(family = 'Arial', size = 8, color = 'black'),
        axis.title=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.y=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.x=element_text(family = 'Arial', size = 8, color = 'black',angle = 90, hjust = 1, vjust =0.5)) 
ggsave('figures/dotplot_celltypes_markers2.svg', device = 'svg', units = 'in', width = 7.24, height = 10, dpi = 600)







DotPlot(ARH_Sex_by_Nutr,group.by = 'cell_type3', features = rev(c(neuro.genes[c(1,2,5,8,10,11)], 'Lepr','Insr','Ghr','Esr1','Prlr','Pgr',
                                                                  
                                                                  'Agrp','Pomc','Glipr1','Tac2','Pdyn','Kiss1','Ghrh','Th','Satb2','Coch','Sst',
                                                                  'Htr2c','Lamp5','Alk','Tbx19','Htr3b','Lef1','Erg',
                                                                  'Klhl1','Ebf3','Ros1','Reln'
                                                                  )), 
        scale = FALSE, scale.by = 'radius') &
  coord_flip() &
  scale_y_discrete(limits = c(
    'Agrp',
    'Pomc',
    'KNDy',
    'Ghrh/Chat',
    'DA',
    'Satb2/Slc18a2',
    'Coch/Slc18a2',
    'Sst/Unc13c',
    'Gad2/Htr2c',
    'Lamp5/Npy5r',
    'Slc17a6/Alk',
    'Tbx19',
    'Htr3b',
    'Lef1',
    'Erg/Lepr',
    'Klhl1/Ebf3',
    'Ebf3/Htr2c',
    'Ros1/Alk',
    'Tac1/Reln'
  )) &
  scale_radius(limits = c(5,100)) &
  scale_color_continuous(limits = c(0,NA), type = "viridis") &
  #labels = c(31,28,33,20,4,1,29,10,35,0,30,22,32,17,21,2,3,8,11,12,26,25,9,16,23,19,34,27,18,15,6,24,14,13,7,5)) &
  labs(x = 'Marker Genes', y = 'Cell-Type') &
  theme_classic() &
  theme(title = element_text(family = 'Arial', size = 8, color = 'black'),
        axis.title=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.y=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.x=element_text(family = 'Arial', size = 8, color = 'black',angle = 90, hjust = 1, vjust =0.5)) 
ggsave('figures/dotplot_celltypes_markersA.svg', device = 'svg', units = 'in', width = 7, height = 7, dpi = 600)
ggsave('figures/dotplot_celltypes_markersA.tiff', device = 'tiff', units = 'in', width = 7, height = 7, dpi = 600)






DotPlot(ARH_Sex_by_Nutr,group.by = 'cell_type3', features = rev(c(neuro.genes[c(1,2,5,8)],
                                                                  
                                                                 
                                                                  'Vipr2',
                                                                  'Nr5a1','Fezf1','Nr2f2','Gda','Slit3',
                                                                  'Tac1','Nos1','Lmx1b','Hdc',
                                                                  'Postn','Adarb2','Tanc1',
                                                                  'Rgs6',
                                                                  'Adcy8','Avp','Oxt'
                                                                  )), 
        scale = FALSE, scale.by = 'radius') &
  coord_flip() &
  scale_y_discrete(limits = c(
    
    'SCN',
    'VMH.01',
    'VMH.02',
    'PVp.01',
    'PVp.02',
    'PVp.03',
    'PVp.04',
    'MM.01',
    'MM.02',
    'MM.03',
    'Tu'
    
  )) &
  scale_radius(limits = c(5,100)) &
  scale_color_continuous(limits = c(0,NA), type = "viridis") &
  #labels = c(31,28,33,20,4,1,29,10,35,0,30,22,32,17,21,2,3,8,11,12,26,25,9,16,23,19,34,27,18,15,6,24,14,13,7,5)) &
  labs(x = 'Marker Genes', y = 'Cell-Type') &
  theme_classic() &
  theme(title = element_text(family = 'Arial', size = 8, color = 'black'),
        axis.title=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.y=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.x=element_text(family = 'Arial', size = 8, color = 'black',angle = 90, hjust = 1, vjust =0.5)) 
ggsave('figures/dotplot_celltypes_markersB.svg', device = 'svg', units = 'in', width = 7, height = 4.5, dpi = 600)
ggsave('figures/dotplot_celltypes_markersB.tiff', device = 'tiff', units = 'in', width = 7, height = 4.5, dpi = 600)






DotPlot(ARH_Sex_by_Nutr,group.by = 'cell_type3', features = rev(c(
                                                                  astro.genes[c(1,4,5)],
                                                                  'Vim','Dio2','Rax','Slc16a2','Crym','Vcan','Col25a1',
                                                                  oligo.genes[c(1,2,5,3)],'Cspg4','Pdgfra',
                                                                  micro.genes[c(1,4,6)],
                                                                  
                                                                  'Ccdc153','Pltp','Foxj1',
                                                                  'Tbx18','Slc6a13','Dcn','Foxc1','Cldn5','Adgrf5','Chga','Timeless','Tshb','Cck')), 
        scale = FALSE, scale.by = 'radius') &
  coord_flip() &
  scale_y_discrete(limits = c(
    
    'Astrocytes',
    'α-Tanycytes',
    'β-Tanycytes',
    'Ependymal',
    'Oligodendrocytes',
    'TOP',
    'OPC',
    'Microglia',
    'VLMC',
    'Endothelial',
    'Pars Tuberalis'
    
  )) &
  scale_radius(limits = c(5,100)) &
  scale_color_continuous(limits = c(0,NA), type = "viridis") &
  #labels = c(31,28,33,20,4,1,29,10,35,0,30,22,32,17,21,2,3,8,11,12,26,25,9,16,23,19,34,27,18,15,6,24,14,13,7,5)) &
  labs(x = 'Marker Genes', y = 'Cell-Type') &
  theme_classic() &
  theme(title = element_text(family = 'Arial', size = 8, color = 'black'),
        axis.title=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.y=element_text(family = 'Arial', size = 8, color = 'black'),
        axis.text.x=element_text(family = 'Arial', size = 8, color = 'black',angle = 90, hjust = 1, vjust =0.5)) 
ggsave('figures/dotplot_celltypes_markersC.svg', device = 'svg', units = 'in', width = 7, height = 5, dpi = 600)
ggsave('figures/dotplot_celltypes_markersC.tiff', device = 'tiff', units = 'in', width = 7, height = 5, dpi = 600)








DimPlot(ARH_Sex_by_Nutr, label = F, group.by = 'cell_type3',
        repel = F, shuffle = F,
        label.size = 4,
        cols = c('#228B22',
                 '#228B22',
                 '#228B22',
                 '#55FF55',
                 '#55FF55',
                 '#008080',
                 '#005455',
                 '#008080',
                 '#556B2F',
                 '#556B2F',
                 '#556B2F',
                 '#55FF55',
                 '#005455',
                 '#005455',
                 '#55FF55',
                 '#228B22',
                 
                 '#C0392B',
                 
                 '#2980B9',
                 '#34495E',
                 
                 '#8E44AD',
                 
                 '#4169E1',
                 '#1E90FF',
                 '#87CEFA',
                 
                 '#F5B041',
                 '#FFD700',
                 '#CA6F1E',
                 
                 '#BC8F8F',
                 
                 '#9370DB',
                 '#BA55D3',
                 '#DA70D6',
                 
                 '#CD853F',
                 
                 '#556B2F',
                 
                 '#BDC3C7',
                 '#F4A460',
                 '#EC7063',
                 '#AF7AC5',
                 'grey',
                 'grey',
                 'grey',
                 'grey',
                 'grey',
                 'grey'
                 )) &
  
  labs(title = "", x = '', y = '') &
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10)) &
  NoLegend()


#ggsave('post_2023-07-18_figures/cell_type4_dimplotv3.tiff', device = 'tiff', units = 'in', width = 9,height = 9,dpi = 600)
#ggsave('post_2023-07-18_figures/cell_type4_dimplotv2.svg', device = 'svg', units = 'in', width = 9,height = 9,dpi = 600)


ARH_Sex_by_Nutr <- SetIdent(ARH_Sex_by_Nutr, value = 'cell_type3')

saveRDS(ARH_Sex_by_Nutr, file = 'data/ARH_Sex_by_Nutr.rds')


cell_type_markers <- FindAllMarkers(ARH_Sex_by_Nutr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log2(1.25), pseudocount.use = 0)

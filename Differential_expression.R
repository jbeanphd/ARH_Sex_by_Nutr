libs <- c( 'gplots','stringi','reshape2','cowplot','RColorBrewer',
           'sctransform','stringr','org.Mm.eg.db','AnnotationDbi',
           'IRanges','S4Vectors','Biobase','BiocGenerics','clusterProfiler',
           'biomaRt','Matrix','DESeq2','RcppThread', 'extrafont', 'openxlsx',
           'Seurat','dplyr','tidyr','ggplot2','harmony','ggalluvial',
           'scDblFinder','SoupX','UpSetR','ComplexUpset','glmGamPoi')
BiocManager::install('glmGamPoi')
devtools::install_github('immunogenomics/presto')
lapply(libs, require, character.only = TRUE)


#Differentially expressed genes 

ARH_Sex_by_Nutr <- readRDS('data/ARH_Sex_by_Nutr.rds')

missing_gene_chromosomes <- readODS::read_ods('../missing_gene_chromosomes.ods')
missing_gene_chromosomes <- missing_gene_chromosomes |> rename(gene = external_gene_name)

#Agrp

Agrp.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                   subset.ident = 'Agrp', 
                                   group.by = 'sexXnutr', 
                                   ident.1 = 'F_Fed', 
                                   ident.2 = 'F_Fast', 
                                   logfc.threshold = log2(1.25), 
                                   min.pct = 0.25, 
                                   #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')


Agrp.F.Fd.v.Fst$gene <- Agrp.F.Fd.v.Fst |> row.names()
ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Agrp.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Agrp.F.Fd.v.Fst <- left_join(Agrp.F.Fd.v.Fst, temp01, by = "gene")

write.xlsx(Agrp.F.Fd.v.Fst, file = '../paper_figures/post_2025-01-06/DE_genes/Agrp/Agrp_F_Fd_v_Fst.xlsx')


Agrp.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                   subset.ident = 'Agrp', 
                                   group.by = 'sexXnutr', 
                                   ident.1 = 'M_Fed', 
                                   ident.2 = 'M_Fast', 
                                   logfc.threshold = log2(1.25), 
                                   min.pct = 0.25, 
                                   #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Agrp.M.Fd.v.Fst$gene <- Agrp.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Agrp.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Agrp.M.Fd.v.Fst <- left_join(Agrp.M.Fd.v.Fst, temp01, by = "gene")

write.xlsx(Agrp.M.Fd.v.Fst, file = '../paper_figures/post_2025-01-06/DE_genes/Agrp/Agrp_M_Fd_v_Fst.xlsx')

Agrp.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Agrp', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'M_Fed', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Agrp.Fd.F.v.M$gene <- Agrp.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Agrp.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Agrp.Fd.F.v.M <- left_join(Agrp.Fd.F.v.M, temp01, by = "gene")

write.xlsx(Agrp.Fd.F.v.M, file = '../paper_figures/post_2025-01-06/DE_genes/Agrp/Agrp_Fd_F_v_M.xlsx')

Agrp.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Agrp', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fast', 
                             ident.2 = 'M_Fast', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Agrp.Fst.F.v.M$gene <- Agrp.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Agrp.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Agrp.Fst.F.v.M <- left_join(Agrp.Fst.F.v.M, temp01, by = "gene")

write.xlsx(Agrp.Fst.F.v.M, file = '../paper_figures/post_2025-01-06/DE_genes/Agrp/Agrp_Fst_F_v_M.xlsx')

#Pomc

Pomc.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Pomc', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Pomc.F.Fd.v.Fst$gene <- Pomc.F.Fd.v.Fst |> row.names()
ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Pomc.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Pomc.F.Fd.v.Fst <- left_join(Pomc.F.Fd.v.Fst, temp01, by = "gene")

write.xlsx(Pomc.F.Fd.v.Fst, file = '../paper_figures/post_2025-04-07/DE_genes/Pomc/Pomc_F_Fd_v_Fst.xlsx')


Pomc.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Pomc', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Pomc.M.Fd.v.Fst$gene <- Pomc.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Pomc.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Pomc.M.Fd.v.Fst <- left_join(Pomc.M.Fd.v.Fst, temp01, by = "gene")

write.xlsx(Pomc.M.Fd.v.Fst, file = '../paper_figures/post_2025-04-07/DE_genes/Pomc/Pomc_M_Fd_v_Fst.xlsx')


Pomc.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Pomc', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Pomc.Fd.F.v.M$gene <- Pomc.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Pomc.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Pomc.Fd.F.v.M <- left_join(Pomc.Fd.F.v.M, temp01, by = "gene")

write.xlsx(Pomc.Fd.F.v.M, file = '../paper_figures/post_2025-04-07/DE_genes/Pomc/Pomc_Fd_F_v_M.xlsx')


Pomc.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'Pomc', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

Pomc.Fst.F.v.M$gene <- Pomc.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Pomc.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Pomc.Fst.F.v.M <- left_join(Pomc.Fst.F.v.M, temp01, by = "gene")

write.xlsx(Pomc.Fst.F.v.M, file = '../paper_figures/post_2025-04-07/DE_genes/Pomc/Pomc_Fst_F_v_M.xlsx')


#KNDy


KNDy.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'KNDy', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

KNDy.F.Fd.v.Fst$gene <- KNDy.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(KNDy.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
KNDy.F.Fd.v.Fst <- left_join(KNDy.F.Fd.v.Fst, temp01, by = "gene")

write.xlsx(KNDy.F.Fd.v.Fst, file = '../paper_figures/post_2025-01-06/DE_genes/KNDy/KNDy_F_Fd_v_Fst.xlsx')

KNDy.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'KNDy', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

KNDy.M.Fd.v.Fst$gene <- KNDy.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(KNDy.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
KNDy.M.Fd.v.Fst <- left_join(KNDy.M.Fd.v.Fst, temp01, by = "gene")


write.xlsx(KNDy.M.Fd.v.Fst, file = '../paper_figures/post_2025-01-06/DE_genes/KNDy/KNDy_M_Fd_v_Fst.xlsx')

KNDy.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'KNDy', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

KNDy.Fd.F.v.M$gene <- KNDy.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(KNDy.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
KNDy.Fd.F.v.M <- left_join(KNDy.Fd.F.v.M, temp01, by = "gene")


write.xlsx(KNDy.Fd.F.v.M, file = '../paper_figures/post_2025-01-06/DE_genes/KNDy/KNDy_Fd_F_v_M.xlsx')

KNDy.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'KNDy', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

KNDy.Fst.F.v.M$gene <- KNDy.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(KNDy.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
KNDy.Fst.F.v.M <- left_join(KNDy.Fst.F.v.M, temp01, by = "gene")


write.xlsx(KNDy.Fst.F.v.M, file = '../paper_figures/post_2025-01-06/DE_genes/KNDy/KNDy_Fst_F_v_M.xlsx')

#DA

DA.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'DA', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

DA.F.Fd.v.Fst$gene <- DA.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(DA.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
DA.F.Fd.v.Fst <- left_join(DA.F.Fd.v.Fst, temp01, by = "gene")

write.xlsx(DA.F.Fd.v.Fst, file = '../paper_figures/post_2025-01-06/DE_genes/DA/DA_F_Fd_v_Fst.xlsx')

DA.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'DA', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

DA.M.Fd.v.Fst$gene <- DA.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(DA.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
DA.M.Fd.v.Fst <- left_join(DA.M.Fd.v.Fst, temp01, by = "gene")


write.xlsx(DA.M.Fd.v.Fst, file = '../paper_figures/post_2025-01-06/DE_genes/DA/DA_M_Fd_v_Fst.xlsx')

DA.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'DA', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                           test.use = 'MAST',
                           latent.vars = 'Sample_ID')

DA.Fd.F.v.M$gene <- DA.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(DA.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
DA.Fd.F.v.M <- left_join(DA.Fd.F.v.M, temp01, by = "gene")


write.xlsx(DA.Fd.F.v.M, file = '../paper_figures/post_2025-01-06/DE_genes/DA/DA_Fd_F_v_M.xlsx')

DA.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'DA', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                            test.use = 'MAST',
                            latent.vars = 'Sample_ID')

DA.Fst.F.v.M$gene <- DA.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(DA.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
DA.Fst.F.v.M <- left_join(DA.Fst.F.v.M, temp01, by = "gene")


write.xlsx(DA.Fst.F.v.M, file = '../paper_figures/post_2025-01-06/DE_genes/DA/DA_Fst_F_v_M.xlsx')

#Ghrh/Chat



Ghrh.Chat.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Ghrh/Chat', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'F_Fast', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Ghrh.Chat.F.Fd.v.Fst$gene <- Ghrh.Chat.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Ghrh.Chat.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Ghrh.Chat.F.Fd.v.Fst <- left_join(Ghrh.Chat.F.Fd.v.Fst, temp01, by = "gene")


write.xlsx(Ghrh.Chat.F.Fd.v.Fst, file = '../paper_figures/post_2025-01-06/DE_genes/Ghrh/Ghrh_F_Fd_v_Fst.xlsx')

Ghrh.Chat.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Ghrh/Chat', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'M_Fed', 
                             ident.2 = 'M_Fast', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Ghrh.Chat.M.Fd.v.Fst$gene <- Ghrh.Chat.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Ghrh.Chat.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Ghrh.Chat.M.Fd.v.Fst <- left_join(Ghrh.Chat.M.Fd.v.Fst, temp01, by = "gene")


write.xlsx(Ghrh.Chat.M.Fd.v.Fst, file = '../paper_figures/post_2025-01-06/DE_genes/Ghrh/Ghrh_M_Fd_v_Fst.xlsx')

Ghrh.Chat.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                           subset.ident = 'Ghrh/Chat', 
                           group.by = 'sexXnutr', 
                           ident.1 = 'F_Fed', 
                           ident.2 = 'M_Fed', 
                           logfc.threshold = log2(1.25), 
                           min.pct = 0.25, 
                           #pseudocount.use = 0,
                           test.use = 'MAST',
                           latent.vars = 'Sample_ID')

Ghrh.Chat.Fd.F.v.M$gene <- Ghrh.Chat.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Ghrh.Chat.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Ghrh.Chat.Fd.F.v.M <- left_join(Ghrh.Chat.Fd.F.v.M, temp01, by = "gene")


write.xlsx(Ghrh.Chat.Fd.F.v.M, file = '../paper_figures/post_2025-01-06/DE_genes/Ghrh/Ghrh_Fd_F_v_M.xlsx')

Ghrh.Chat.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                            subset.ident = 'Ghrh/Chat', 
                            group.by = 'sexXnutr', 
                            ident.1 = 'F_Fast', 
                            ident.2 = 'M_Fast', 
                            logfc.threshold = log2(1.25), 
                            min.pct = 0.25, 
                            #pseudocount.use = 0,
                            test.use = 'MAST',
                            latent.vars = 'Sample_ID')

Ghrh.Chat.Fst.F.v.M$gene <- Ghrh.Chat.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Ghrh.Chat.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Ghrh.Chat.Fst.F.v.M <- left_join(Ghrh.Chat.Fst.F.v.M, temp01, by = "gene")


write.xlsx(Ghrh.Chat.Fst.F.v.M, file = '../paper_figures/post_2025-01-06/DE_genes/Ghrh/Ghrh_Fst_F_v_M.xlsx')


#Sst/Unc13c

Sst.Unc13c.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                    subset.ident = 'Sst/Unc13c', 
                                    group.by = 'sexXnutr', 
                                    ident.1 = 'F_Fed', 
                                    ident.2 = 'F_Fast', 
                                    logfc.threshold = log2(1.25), 
                                    min.pct = 0.25, 
                                    #pseudocount.use = 0,
                                    test.use = 'MAST',
                                    latent.vars = 'Sample_ID')

Sst.Unc13c.F.Fd.v.Fst$gene <- Sst.Unc13c.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Sst.Unc13c.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Sst.Unc13c.F.Fd.v.Fst <- left_join(Sst.Unc13c.F.Fd.v.Fst, temp01, by = "gene")



Sst.Unc13c.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                    subset.ident = 'Sst/Unc13c', 
                                    group.by = 'sexXnutr', 
                                    ident.1 = 'M_Fed', 
                                    ident.2 = 'M_Fast', 
                                    logfc.threshold = log2(1.25), 
                                    min.pct = 0.25, 
                                    #pseudocount.use = 0,
                                    test.use = 'MAST',
                                    latent.vars = 'Sample_ID')

Sst.Unc13c.M.Fd.v.Fst$gene <- Sst.Unc13c.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Sst.Unc13c.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Sst.Unc13c.M.Fd.v.Fst <- left_join(Sst.Unc13c.M.Fd.v.Fst, temp01, by = "gene")



Sst.Unc13c.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                  subset.ident = 'Sst/Unc13c', 
                                  group.by = 'sexXnutr', 
                                  ident.1 = 'F_Fed', 
                                  ident.2 = 'M_Fed', 
                                  logfc.threshold = log2(1.25), 
                                  min.pct = 0.25, 
                                  #pseudocount.use = 0,
                                  test.use = 'MAST',
                                  latent.vars = 'Sample_ID')

Sst.Unc13c.Fd.F.v.M$gene <- Sst.Unc13c.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Sst.Unc13c.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Sst.Unc13c.Fd.F.v.M <- left_join(Sst.Unc13c.Fd.F.v.M, temp01, by = "gene")



Sst.Unc13c.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                   subset.ident = 'Sst/Unc13c', 
                                   group.by = 'sexXnutr', 
                                   ident.1 = 'F_Fast', 
                                   ident.2 = 'M_Fast', 
                                   logfc.threshold = log2(1.25), 
                                   min.pct = 0.25, 
                                   #pseudocount.use = 0,
                                   test.use = 'MAST',
                                   latent.vars = 'Sample_ID')

Sst.Unc13c.Fst.F.v.M$gene <- Sst.Unc13c.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Sst.Unc13c.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Sst.Unc13c.Fst.F.v.M <- left_join(Sst.Unc13c.Fst.F.v.M, temp01, by = "gene")



#Lamp5/Npy5r

Lamp5.Npy5r.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                    subset.ident = 'Lamp5/Npy5r', 
                                    group.by = 'sexXnutr', 
                                    ident.1 = 'F_Fed', 
                                    ident.2 = 'F_Fast', 
                                    logfc.threshold = log2(1.25), 
                                    min.pct = 0.25, 
                                    #pseudocount.use = 0,
                                    test.use = 'MAST',
                                    latent.vars = 'Sample_ID')

Lamp5.Npy5r.F.Fd.v.Fst$gene <- Lamp5.Npy5r.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Lamp5.Npy5r.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Lamp5.Npy5r.F.Fd.v.Fst <- left_join(Lamp5.Npy5r.F.Fd.v.Fst, temp01, by = "gene")



Lamp5.Npy5r.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                    subset.ident = 'Lamp5/Npy5r', 
                                    group.by = 'sexXnutr', 
                                    ident.1 = 'M_Fed', 
                                    ident.2 = 'M_Fast', 
                                    logfc.threshold = log2(1.25), 
                                    min.pct = 0.25, 
                                    #pseudocount.use = 0,
                                    test.use = 'MAST',
                                    latent.vars = 'Sample_ID')

Lamp5.Npy5r.M.Fd.v.Fst$gene <- Lamp5.Npy5r.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Lamp5.Npy5r.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Lamp5.Npy5r.M.Fd.v.Fst <- left_join(Lamp5.Npy5r.M.Fd.v.Fst, temp01, by = "gene")



Lamp5.Npy5r.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                  subset.ident = 'Lamp5/Npy5r', 
                                  group.by = 'sexXnutr', 
                                  ident.1 = 'F_Fed', 
                                  ident.2 = 'M_Fed', 
                                  logfc.threshold = log2(1.25), 
                                  min.pct = 0.25, 
                                  #pseudocount.use = 0,
                                  test.use = 'MAST',
                                  latent.vars = 'Sample_ID')

Lamp5.Npy5r.Fd.F.v.M$gene <- Lamp5.Npy5r.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Lamp5.Npy5r.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Lamp5.Npy5r.Fd.F.v.M <- left_join(Lamp5.Npy5r.Fd.F.v.M, temp01, by = "gene")



Lamp5.Npy5r.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                   subset.ident = 'Lamp5/Npy5r', 
                                   group.by = 'sexXnutr', 
                                   ident.1 = 'F_Fast', 
                                   ident.2 = 'M_Fast', 
                                   logfc.threshold = log2(1.25), 
                                   min.pct = 0.25, 
                                   #pseudocount.use = 0,
                                   test.use = 'MAST',
                                   latent.vars = 'Sample_ID')

Lamp5.Npy5r.Fst.F.v.M$gene <- Lamp5.Npy5r.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Lamp5.Npy5r.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Lamp5.Npy5r.Fst.F.v.M <- left_join(Lamp5.Npy5r.Fst.F.v.M, temp01, by = "gene")


#Lef1

Lef1.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                    subset.ident = 'Lef1', 
                                    group.by = 'sexXnutr', 
                                    ident.1 = 'F_Fed', 
                                    ident.2 = 'F_Fast', 
                                    logfc.threshold = log2(1.25), 
                                    min.pct = 0.25, 
                                    #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Lef1.F.Fd.v.Fst$gene <- Lef1.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Lef1.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Lef1.F.Fd.v.Fst <- left_join(Lef1.F.Fd.v.Fst, temp01, by = "gene")



Lef1.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                    subset.ident = 'Lef1', 
                                    group.by = 'sexXnutr', 
                                    ident.1 = 'M_Fed', 
                                    ident.2 = 'M_Fast', 
                                    logfc.threshold = log2(1.25), 
                                    min.pct = 0.25, 
                                    #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Lef1.M.Fd.v.Fst$gene <- Lef1.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Lef1.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Lef1.M.Fd.v.Fst <- left_join(Lef1.M.Fd.v.Fst, temp01, by = "gene")



Lef1.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                  subset.ident = 'Lef1', 
                                  group.by = 'sexXnutr', 
                                  ident.1 = 'F_Fed', 
                                  ident.2 = 'M_Fed', 
                                  logfc.threshold = log2(1.25), 
                                  min.pct = 0.25, 
                                  #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Lef1.Fd.F.v.M$gene <- Lef1.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Lef1.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Lef1.Fd.F.v.M <- left_join(Lef1.Fd.F.v.M, temp01, by = "gene")



Lef1.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                   subset.ident = 'Lef1', 
                                   group.by = 'sexXnutr', 
                                   ident.1 = 'F_Fast', 
                                   ident.2 = 'M_Fast', 
                                   logfc.threshold = log2(1.25), 
                                   min.pct = 0.25, 
                                   #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

Lef1.Fst.F.v.M$gene <- Lef1.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Lef1.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Lef1.Fst.F.v.M <- left_join(Lef1.Fst.F.v.M, temp01, by = "gene")



#Htr3b

Htr3b.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Htr3b', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Htr3b.F.Fd.v.Fst$gene <- Htr3b.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Htr3b.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Htr3b.F.Fd.v.Fst <- left_join(Htr3b.F.Fd.v.Fst, temp01, by = "gene")



Htr3b.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Htr3b', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Htr3b.M.Fd.v.Fst$gene <- Htr3b.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Htr3b.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Htr3b.M.Fd.v.Fst <- left_join(Htr3b.M.Fd.v.Fst, temp01, by = "gene")



Htr3b.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Htr3b', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Htr3b.Fd.F.v.M$gene <- Htr3b.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Htr3b.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Htr3b.Fd.F.v.M <- left_join(Htr3b.Fd.F.v.M, temp01, by = "gene")



Htr3b.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'Htr3b', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

Htr3b.Fst.F.v.M$gene <- Htr3b.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Htr3b.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Htr3b.Fst.F.v.M <- left_join(Htr3b.Fst.F.v.M, temp01, by = "gene")


#Gad2/Htr2c

Gad2.Htr2c.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Gad2/Htr2c', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Gad2.Htr2c.F.Fd.v.Fst$gene <- Gad2.Htr2c.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Gad2.Htr2c.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Gad2.Htr2c.F.Fd.v.Fst <- left_join(Gad2.Htr2c.F.Fd.v.Fst, temp01, by = "gene")



Gad2.Htr2c.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Gad2/Htr2c', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Gad2.Htr2c.M.Fd.v.Fst$gene <- Gad2.Htr2c.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Gad2.Htr2c.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Gad2.Htr2c.M.Fd.v.Fst <- left_join(Gad2.Htr2c.M.Fd.v.Fst, temp01, by = "gene")



Gad2.Htr2c.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Gad2/Htr2c', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Gad2.Htr2c.Fd.F.v.M$gene <- Gad2.Htr2c.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Gad2.Htr2c.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Gad2.Htr2c.Fd.F.v.M <- left_join(Gad2.Htr2c.Fd.F.v.M, temp01, by = "gene")



Gad2.Htr2c.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'Gad2/Htr2c', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

Gad2.Htr2c.Fst.F.v.M$gene <- Gad2.Htr2c.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Gad2.Htr2c.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Gad2.Htr2c.Fst.F.v.M <- left_join(Gad2.Htr2c.Fst.F.v.M, temp01, by = "gene")



#Ros1/Alk

Ros1.Alk.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Ros1/Alk', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Ros1.Alk.F.Fd.v.Fst$gene <- Ros1.Alk.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Ros1.Alk.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Ros1.Alk.F.Fd.v.Fst <- left_join(Ros1.Alk.F.Fd.v.Fst, temp01, by = "gene")



Ros1.Alk.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Ros1/Alk', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Ros1.Alk.M.Fd.v.Fst$gene <- Ros1.Alk.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Ros1.Alk.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Ros1.Alk.M.Fd.v.Fst <- left_join(Ros1.Alk.M.Fd.v.Fst, temp01, by = "gene")



Ros1.Alk.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Ros1/Alk', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Ros1.Alk.Fd.F.v.M$gene <- Ros1.Alk.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Ros1.Alk.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Ros1.Alk.Fd.F.v.M <- left_join(Ros1.Alk.Fd.F.v.M, temp01, by = "gene")



Ros1.Alk.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'Ros1/Alk', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

Ros1.Alk.Fst.F.v.M$gene <- Ros1.Alk.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Ros1.Alk.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Ros1.Alk.Fst.F.v.M <- left_join(Ros1.Alk.Fst.F.v.M, temp01, by = "gene")



#Satb2/Slc18a2


Satb2.Slc18a2.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Satb2/Slc18a2', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Satb2.Slc18a2.F.Fd.v.Fst$gene <- Satb2.Slc18a2.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Satb2.Slc18a2.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Satb2.Slc18a2.F.Fd.v.Fst <- left_join(Satb2.Slc18a2.F.Fd.v.Fst, temp01, by = "gene")



Satb2.Slc18a2.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Satb2/Slc18a2', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Satb2.Slc18a2.M.Fd.v.Fst$gene <- Satb2.Slc18a2.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Satb2.Slc18a2.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Satb2.Slc18a2.M.Fd.v.Fst <- left_join(Satb2.Slc18a2.M.Fd.v.Fst, temp01, by = "gene")



Satb2.Slc18a2.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Satb2/Slc18a2', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Satb2.Slc18a2.Fd.F.v.M$gene <- Satb2.Slc18a2.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Satb2.Slc18a2.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Satb2.Slc18a2.Fd.F.v.M <- left_join(Satb2.Slc18a2.Fd.F.v.M, temp01, by = "gene")



Satb2.Slc18a2.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'Satb2/Slc18a2', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

Satb2.Slc18a2.Fst.F.v.M$gene <- Satb2.Slc18a2.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Satb2.Slc18a2.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Satb2.Slc18a2.Fst.F.v.M <- left_join(Satb2.Slc18a2.Fst.F.v.M, temp01, by = "gene")


#Coch/Slc18a2

Coch.Slc18a2.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Coch/Slc18a2', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Coch.Slc18a2.F.Fd.v.Fst$gene <- Coch.Slc18a2.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Coch.Slc18a2.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Coch.Slc18a2.F.Fd.v.Fst <- left_join(Coch.Slc18a2.F.Fd.v.Fst, temp01, by = "gene")



Coch.Slc18a2.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Coch/Slc18a2', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Coch.Slc18a2.M.Fd.v.Fst$gene <- Coch.Slc18a2.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Coch.Slc18a2.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Coch.Slc18a2.M.Fd.v.Fst <- left_join(Coch.Slc18a2.M.Fd.v.Fst, temp01, by = "gene")



Coch.Slc18a2.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Coch/Slc18a2', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Coch.Slc18a2.Fd.F.v.M$gene <- Coch.Slc18a2.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Coch.Slc18a2.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Coch.Slc18a2.Fd.F.v.M <- left_join(Coch.Slc18a2.Fd.F.v.M, temp01, by = "gene")



Coch.Slc18a2.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'Coch/Slc18a2', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

Coch.Slc18a2.Fst.F.v.M$gene <- Coch.Slc18a2.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Coch.Slc18a2.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Coch.Slc18a2.Fst.F.v.M <- left_join(Coch.Slc18a2.Fst.F.v.M, temp01, by = "gene")



#Tbx19

Tbx19.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Tbx19', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Tbx19.F.Fd.v.Fst$gene <- Tbx19.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Tbx19.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Tbx19.F.Fd.v.Fst <- left_join(Tbx19.F.Fd.v.Fst, temp01, by = "gene")



Tbx19.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Tbx19', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Tbx19.M.Fd.v.Fst$gene <- Tbx19.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Tbx19.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Tbx19.M.Fd.v.Fst <- left_join(Tbx19.M.Fd.v.Fst, temp01, by = "gene")



Tbx19.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Tbx19', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Tbx19.Fd.F.v.M$gene <- Tbx19.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Tbx19.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Tbx19.Fd.F.v.M <- left_join(Tbx19.Fd.F.v.M, temp01, by = "gene")



Tbx19.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'Tbx19', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

Tbx19.Fst.F.v.M$gene <- Tbx19.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Tbx19.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Tbx19.Fst.F.v.M <- left_join(Tbx19.Fst.F.v.M, temp01, by = "gene")



#Tbx15

Tbx15.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Tbx15', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Tbx15.F.Fd.v.Fst$gene <- Tbx15.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Tbx15.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Tbx15.F.Fd.v.Fst <- left_join(Tbx15.F.Fd.v.Fst, temp01, by = "gene")



Tbx15.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Tbx15', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Tbx15.M.Fd.v.Fst$gene <- Tbx15.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Tbx15.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Tbx15.M.Fd.v.Fst <- left_join(Tbx15.M.Fd.v.Fst, temp01, by = "gene")



Tbx15.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Tbx15', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Tbx15.Fd.F.v.M$gene <- Tbx15.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Tbx15.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Tbx15.Fd.F.v.M <- left_join(Tbx15.Fd.F.v.M, temp01, by = "gene")



Tbx15.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'Tbx15', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

Tbx15.Fst.F.v.M$gene <- Tbx15.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Tbx15.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Tbx15.Fst.F.v.M <- left_join(Tbx15.Fst.F.v.M, temp01, by = "gene")


#Ebf3/Htr2c

Ebf3.Htr2c.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Ebf3/Htr2c', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Ebf3.Htr2c.F.Fd.v.Fst$gene <- Ebf3.Htr2c.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Ebf3.Htr2c.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Ebf3.Htr2c.F.Fd.v.Fst <- left_join(Ebf3.Htr2c.F.Fd.v.Fst, temp01, by = "gene")



Ebf3.Htr2c.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Ebf3/Htr2c', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Ebf3.Htr2c.M.Fd.v.Fst$gene <- Ebf3.Htr2c.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Ebf3.Htr2c.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Ebf3.Htr2c.M.Fd.v.Fst <- left_join(Ebf3.Htr2c.M.Fd.v.Fst, temp01, by = "gene")



Ebf3.Htr2c.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Ebf3/Htr2c', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Ebf3.Htr2c.Fd.F.v.M$gene <- Ebf3.Htr2c.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Ebf3.Htr2c.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Ebf3.Htr2c.Fd.F.v.M <- left_join(Ebf3.Htr2c.Fd.F.v.M, temp01, by = "gene")



Ebf3.Htr2c.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'Ebf3/Htr2c', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

Ebf3.Htr2c.Fst.F.v.M$gene <- Ebf3.Htr2c.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Ebf3.Htr2c.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Ebf3.Htr2c.Fst.F.v.M <- left_join(Ebf3.Htr2c.Fst.F.v.M, temp01, by = "gene")



#Slc17a6/Alk


Slc17a6.Alk.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Slc17a6/Alk', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Slc17a6.Alk.F.Fd.v.Fst$gene <- Slc17a6.Alk.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Slc17a6.Alk.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Slc17a6.Alk.F.Fd.v.Fst <- left_join(Slc17a6.Alk.F.Fd.v.Fst, temp01, by = "gene")



Slc17a6.Alk.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Slc17a6/Alk', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Slc17a6.Alk.M.Fd.v.Fst$gene <- Slc17a6.Alk.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Slc17a6.Alk.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Slc17a6.Alk.M.Fd.v.Fst <- left_join(Slc17a6.Alk.M.Fd.v.Fst, temp01, by = "gene")



Slc17a6.Alk.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Slc17a6/Alk', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Slc17a6.Alk.Fd.F.v.M$gene <- Slc17a6.Alk.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Slc17a6.Alk.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Slc17a6.Alk.Fd.F.v.M <- left_join(Slc17a6.Alk.Fd.F.v.M, temp01, by = "gene")



Slc17a6.Alk.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'Slc17a6/Alk', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

Slc17a6.Alk.Fst.F.v.M$gene <- Slc17a6.Alk.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Slc17a6.Alk.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Slc17a6.Alk.Fst.F.v.M <- left_join(Slc17a6.Alk.Fst.F.v.M, temp01, by = "gene")



#Erg/Lepr


Erg.Lepr.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Erg/Lepr', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Erg.Lepr.F.Fd.v.Fst$gene <- Erg.Lepr.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Erg.Lepr.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Erg.Lepr.F.Fd.v.Fst <- left_join(Erg.Lepr.F.Fd.v.Fst, temp01, by = "gene")



Erg.Lepr.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Erg/Lepr', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Erg.Lepr.M.Fd.v.Fst$gene <- Erg.Lepr.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Erg.Lepr.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Erg.Lepr.M.Fd.v.Fst <- left_join(Erg.Lepr.M.Fd.v.Fst, temp01, by = "gene")



Erg.Lepr.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Erg/Lepr', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Erg.Lepr.Fd.F.v.M$gene <- Erg.Lepr.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Erg.Lepr.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Erg.Lepr.Fd.F.v.M <- left_join(Erg.Lepr.Fd.F.v.M, temp01, by = "gene")



Erg.Lepr.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'Erg/Lepr', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

Erg.Lepr.Fst.F.v.M$gene <- Erg.Lepr.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Erg.Lepr.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Erg.Lepr.Fst.F.v.M <- left_join(Erg.Lepr.Fst.F.v.M, temp01, by = "gene")


#Tac1/Reln

Tac1.Reln.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Tac1/Reln', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Tac1.Reln.F.Fd.v.Fst$gene <- Tac1.Reln.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Tac1.Reln.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Tac1.Reln.F.Fd.v.Fst <- left_join(Tac1.Reln.F.Fd.v.Fst, temp01, by = "gene")



Tac1.Reln.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Tac1/Reln', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Tac1.Reln.M.Fd.v.Fst$gene <- Tac1.Reln.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Tac1.Reln.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Tac1.Reln.M.Fd.v.Fst <- left_join(Tac1.Reln.M.Fd.v.Fst, temp01, by = "gene")



Tac1.Reln.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Tac1/Reln', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Tac1.Reln.Fd.F.v.M$gene <- Tac1.Reln.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Tac1.Reln.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Tac1.Reln.Fd.F.v.M <- left_join(Tac1.Reln.Fd.F.v.M, temp01, by = "gene")



Tac1.Reln.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'Tac1/Reln', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

Tac1.Reln.Fst.F.v.M$gene <- Tac1.Reln.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Tac1.Reln.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Tac1.Reln.Fst.F.v.M <- left_join(Tac1.Reln.Fst.F.v.M, temp01, by = "gene")



#Klhl1/Ebf3


Klhl1.Ebf3.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Klhl1/Ebf3', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Klhl1.Ebf3.F.Fd.v.Fst$gene <- Klhl1.Ebf3.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Klhl1.Ebf3.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Klhl1.Ebf3.F.Fd.v.Fst <- left_join(Klhl1.Ebf3.F.Fd.v.Fst, temp01, by = "gene")



Klhl1.Ebf3.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Klhl1/Ebf3', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Klhl1.Ebf3.M.Fd.v.Fst$gene <- Klhl1.Ebf3.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Klhl1.Ebf3.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Klhl1.Ebf3.M.Fd.v.Fst <- left_join(Klhl1.Ebf3.M.Fd.v.Fst, temp01, by = "gene")



Klhl1.Ebf3.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Klhl1/Ebf3', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Klhl1.Ebf3.Fd.F.v.M$gene <- Klhl1.Ebf3.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Klhl1.Ebf3.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Klhl1.Ebf3.Fd.F.v.M <- left_join(Klhl1.Ebf3.Fd.F.v.M, temp01, by = "gene")



Klhl1.Ebf3.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'Klhl1/Ebf3', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

Klhl1.Ebf3.Fst.F.v.M$gene <- Klhl1.Ebf3.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Klhl1.Ebf3.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Klhl1.Ebf3.Fst.F.v.M <- left_join(Klhl1.Ebf3.Fst.F.v.M, temp01, by = "gene")




#SCN


SCN.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'SCN', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

SCN.F.Fd.v.Fst$gene <- SCN.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(SCN.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
SCN.F.Fd.v.Fst <- left_join(SCN.F.Fd.v.Fst, temp01, by = "gene")



SCN.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'SCN', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

SCN.M.Fd.v.Fst$gene <- SCN.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(SCN.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
SCN.M.Fd.v.Fst <- left_join(SCN.M.Fd.v.Fst, temp01, by = "gene")



SCN.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'SCN', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                            test.use = 'MAST',
                            latent.vars = 'Sample_ID')

SCN.Fd.F.v.M$gene <- SCN.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(SCN.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
SCN.Fd.F.v.M <- left_join(SCN.Fd.F.v.M, temp01, by = "gene")



SCN.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'SCN', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

SCN.Fst.F.v.M$gene <- SCN.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(SCN.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
SCN.Fst.F.v.M <- left_join(SCN.Fst.F.v.M, temp01, by = "gene")


#VMH.01


VMH.01.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'VMH.01', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

VMH.01.F.Fd.v.Fst$gene <- VMH.01.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(VMH.01.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
VMH.01.F.Fd.v.Fst <- left_join(VMH.01.F.Fd.v.Fst, temp01, by = "gene")



VMH.01.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'VMH.01', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

VMH.01.M.Fd.v.Fst$gene <- VMH.01.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(VMH.01.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
VMH.01.M.Fd.v.Fst <- left_join(VMH.01.M.Fd.v.Fst, temp01, by = "gene")



VMH.01.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'VMH.01', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             pseudocount.use = 0)

VMH.01.Fd.F.v.M$gene <- VMH.01.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(VMH.01.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
VMH.01.Fd.F.v.M <- left_join(VMH.01.Fd.F.v.M, temp01, by = "gene")



VMH.01.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'VMH.01', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

VMH.01.Fst.F.v.M$gene <- VMH.01.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(VMH.01.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
VMH.01.Fst.F.v.M <- left_join(VMH.01.Fst.F.v.M, temp01, by = "gene")



#VMH.02


VMH.02.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'VMH.02', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

VMH.02.F.Fd.v.Fst$gene <- VMH.02.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(VMH.02.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
VMH.02.F.Fd.v.Fst <- left_join(VMH.02.F.Fd.v.Fst, temp01, by = "gene")



VMH.02.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'VMH.02', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               pseudocount.use = 0)

VMH.02.M.Fd.v.Fst$gene <- VMH.02.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(VMH.02.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
VMH.02.M.Fd.v.Fst <- left_join(VMH.02.M.Fd.v.Fst, temp01, by = "gene")



VMH.02.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'VMH.02', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

VMH.02.Fd.F.v.M$gene <- VMH.02.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(VMH.02.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
VMH.02.Fd.F.v.M <- left_join(VMH.02.Fd.F.v.M, temp01, by = "gene")



VMH.02.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'VMH.02', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

VMH.02.Fst.F.v.M$gene <- VMH.02.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(VMH.02.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
VMH.02.Fst.F.v.M <- left_join(VMH.02.Fst.F.v.M, temp01, by = "gene")


#PVp.01

PVp.01.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'PVp.01', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

PVp.01.F.Fd.v.Fst$gene <- PVp.01.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(PVp.01.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
PVp.01.F.Fd.v.Fst <- left_join(PVp.01.F.Fd.v.Fst, temp01, by = "gene")



PVp.01.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'PVp.01', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               pseudocount.use = 0)

PVp.01.M.Fd.v.Fst$gene <- PVp.01.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(PVp.01.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
PVp.01.M.Fd.v.Fst <- left_join(PVp.01.M.Fd.v.Fst, temp01, by = "gene")



PVp.01.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'PVp.01', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

PVp.01.Fd.F.v.M$gene <- PVp.01.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(PVp.01.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
PVp.01.Fd.F.v.M <- left_join(PVp.01.Fd.F.v.M, temp01, by = "gene")



PVp.01.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'PVp.01', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

PVp.01.Fst.F.v.M$gene <- PVp.01.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(PVp.01.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
PVp.01.Fst.F.v.M <- left_join(PVp.01.Fst.F.v.M, temp01, by = "gene")



#PVp.02

PVp.02.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'PVp.02', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

PVp.02.F.Fd.v.Fst$gene <- PVp.02.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(PVp.02.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
PVp.02.F.Fd.v.Fst <- left_join(PVp.02.F.Fd.v.Fst, temp01, by = "gene")



PVp.02.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'PVp.02', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

PVp.02.M.Fd.v.Fst$gene <- PVp.02.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(PVp.02.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
PVp.02.M.Fd.v.Fst <- left_join(PVp.02.M.Fd.v.Fst, temp01, by = "gene")



PVp.02.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'PVp.02', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

PVp.02.Fd.F.v.M$gene <- PVp.02.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(PVp.02.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
PVp.02.Fd.F.v.M <- left_join(PVp.02.Fd.F.v.M, temp01, by = "gene")



PVp.02.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'PVp.02', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

PVp.02.Fst.F.v.M$gene <- PVp.02.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(PVp.02.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
PVp.02.Fst.F.v.M <- left_join(PVp.02.Fst.F.v.M, temp01, by = "gene")


#PVp.03

PVp.03.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'PVp.03', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

PVp.03.F.Fd.v.Fst$gene <- PVp.03.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(PVp.03.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
PVp.03.F.Fd.v.Fst <- left_join(PVp.03.F.Fd.v.Fst, temp01, by = "gene")



PVp.03.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'PVp.03', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

PVp.03.M.Fd.v.Fst$gene <- PVp.03.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(PVp.03.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
PVp.03.M.Fd.v.Fst <- left_join(PVp.03.M.Fd.v.Fst, temp01, by = "gene")



PVp.03.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'PVp.03', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

PVp.03.Fd.F.v.M$gene <- PVp.03.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(PVp.03.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
PVp.03.Fd.F.v.M <- left_join(PVp.03.Fd.F.v.M, temp01, by = "gene")



PVp.03.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'PVp.03', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

PVp.03.Fst.F.v.M$gene <- PVp.03.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(PVp.03.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
PVp.03.Fst.F.v.M <- left_join(PVp.03.Fst.F.v.M, temp01, by = "gene")



#PVp.04


PVp.04.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'PVp.04', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

PVp.04.F.Fd.v.Fst$gene <- PVp.04.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(PVp.04.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
PVp.04.F.Fd.v.Fst <- left_join(PVp.04.F.Fd.v.Fst, temp01, by = "gene")



PVp.04.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'PVp.04', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

PVp.04.M.Fd.v.Fst$gene <- PVp.04.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(PVp.04.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
PVp.04.M.Fd.v.Fst <- left_join(PVp.04.M.Fd.v.Fst, temp01, by = "gene")



PVp.04.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'PVp.04', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

PVp.04.Fd.F.v.M$gene <- PVp.04.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(PVp.04.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
PVp.04.Fd.F.v.M <- left_join(PVp.04.Fd.F.v.M, temp01, by = "gene")



PVp.04.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'PVp.04', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

PVp.04.Fst.F.v.M$gene <- PVp.04.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(PVp.04.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
PVp.04.Fst.F.v.M <- left_join(PVp.04.Fst.F.v.M, temp01, by = "gene")



#MM.01

MM.01.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'MM.01', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

MM.01.F.Fd.v.Fst$gene <- MM.01.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(MM.01.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
MM.01.F.Fd.v.Fst <- left_join(MM.01.F.Fd.v.Fst, temp01, by = "gene")



MM.01.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'MM.01', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

MM.01.M.Fd.v.Fst$gene <- MM.01.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(MM.01.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
MM.01.M.Fd.v.Fst <- left_join(MM.01.M.Fd.v.Fst, temp01, by = "gene")



MM.01.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'MM.01', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

MM.01.Fd.F.v.M$gene <- MM.01.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(MM.01.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
MM.01.Fd.F.v.M <- left_join(MM.01.Fd.F.v.M, temp01, by = "gene")



MM.01.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'MM.01', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

MM.01.Fst.F.v.M$gene <- MM.01.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(MM.01.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
MM.01.Fst.F.v.M <- left_join(MM.01.Fst.F.v.M, temp01, by = "gene")



#MM.02


MM.02.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'MM.02', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

MM.02.F.Fd.v.Fst$gene <- MM.02.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(MM.02.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
MM.02.F.Fd.v.Fst <- left_join(MM.02.F.Fd.v.Fst, temp01, by = "gene")



MM.02.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'MM.02', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

MM.02.M.Fd.v.Fst$gene <- MM.02.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(MM.02.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
MM.02.M.Fd.v.Fst <- left_join(MM.02.M.Fd.v.Fst, temp01, by = "gene")



MM.02.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'MM.02', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

MM.02.Fd.F.v.M$gene <- MM.02.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(MM.02.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
MM.02.Fd.F.v.M <- left_join(MM.02.Fd.F.v.M, temp01, by = "gene")



MM.02.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'MM.02', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

MM.02.Fst.F.v.M$gene <- MM.02.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(MM.02.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
MM.02.Fst.F.v.M <- left_join(MM.02.Fst.F.v.M, temp01, by = "gene")



#MM.03


MM.03.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'MM.03', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

MM.03.F.Fd.v.Fst$gene <- MM.03.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(MM.03.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
MM.03.F.Fd.v.Fst <- left_join(MM.03.F.Fd.v.Fst, temp01, by = "gene")



MM.03.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'MM.03', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

MM.03.M.Fd.v.Fst$gene <- MM.03.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(MM.03.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
MM.03.M.Fd.v.Fst <- left_join(MM.03.M.Fd.v.Fst, temp01, by = "gene")



MM.03.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'MM.03', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

MM.03.Fd.F.v.M$gene <- MM.03.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(MM.03.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
MM.03.Fd.F.v.M <- left_join(MM.03.Fd.F.v.M, temp01, by = "gene")



MM.03.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'MM.03', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

MM.03.Fst.F.v.M$gene <- MM.03.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(MM.03.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
MM.03.Fst.F.v.M <- left_join(MM.03.Fst.F.v.M, temp01, by = "gene")


#Tu

Tu.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Tu', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Tu.F.Fd.v.Fst$gene <- Tu.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Tu.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Tu.F.Fd.v.Fst <- left_join(Tu.F.Fd.v.Fst, temp01, by = "gene")



Tu.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Tu', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Tu.M.Fd.v.Fst$gene <- Tu.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Tu.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Tu.M.Fd.v.Fst <- left_join(Tu.M.Fd.v.Fst, temp01, by = "gene")



Tu.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Tu', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                           test.use = 'MAST',
                           latent.vars = 'Sample_ID')

Tu.Fd.F.v.M$gene <- Tu.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Tu.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Tu.Fd.F.v.M <- left_join(Tu.Fd.F.v.M, temp01, by = "gene")



Tu.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'Tu', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                            test.use = 'MAST',
                            latent.vars = 'Sample_ID')

Tu.Fst.F.v.M$gene <- Tu.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Tu.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Tu.Fst.F.v.M <- left_join(Tu.Fst.F.v.M, temp01, by = "gene")



#Pars Tuberalis

ParsTub.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Pars Tuberalis', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

ParsTub.F.Fd.v.Fst$gene <- ParsTub.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(ParsTub.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
ParsTub.F.Fd.v.Fst <- left_join(ParsTub.F.Fd.v.Fst, temp01, by = "gene")



ParsTub.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Pars Tuberalis', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

ParsTub.M.Fd.v.Fst$gene <- ParsTub.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(ParsTub.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
ParsTub.M.Fd.v.Fst <- left_join(ParsTub.M.Fd.v.Fst, temp01, by = "gene")



ParsTub.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Pars Tuberalis', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

ParsTub.Fd.F.v.M$gene <- ParsTub.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(ParsTub.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
ParsTub.Fd.F.v.M <- left_join(ParsTub.Fd.F.v.M, temp01, by = "gene")



ParsTub.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'Pars Tuberalis', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

ParsTub.Fst.F.v.M$gene <- ParsTub.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(ParsTub.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
ParsTub.Fst.F.v.M <- left_join(ParsTub.Fst.F.v.M, temp01, by = "gene")



#Astrocytes

Astrocytes.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Astrocytes', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Astrocytes.F.Fd.v.Fst$gene <- Astrocytes.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Astrocytes.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Astrocytes.F.Fd.v.Fst <- left_join(Astrocytes.F.Fd.v.Fst, temp01, by = "gene")



Astrocytes.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'Astrocytes', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Astrocytes.M.Fd.v.Fst$gene <- Astrocytes.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Astrocytes.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Astrocytes.M.Fd.v.Fst <- left_join(Astrocytes.M.Fd.v.Fst, temp01, by = "gene")



Astrocytes.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'Astrocytes', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

Astrocytes.Fd.F.v.M$gene <- Astrocytes.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Astrocytes.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Astrocytes.Fd.F.v.M <- left_join(Astrocytes.Fd.F.v.M, temp01, by = "gene")



Astrocytes.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'Astrocytes', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

Astrocytes.Fst.F.v.M$gene <- Astrocytes.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Astrocytes.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Astrocytes.Fst.F.v.M <- left_join(Astrocytes.Fst.F.v.M, temp01, by = "gene")



# alpha-Tanycytes

α.Tanycytes.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'α-Tanycytes', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'F_Fed', 
                               ident.2 = 'F_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

α.Tanycytes.F.Fd.v.Fst$gene <- α.Tanycytes.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(α.Tanycytes.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
α.Tanycytes.F.Fd.v.Fst <- left_join(α.Tanycytes.F.Fd.v.Fst, temp01, by = "gene")



α.Tanycytes.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                               subset.ident = 'α-Tanycytes', 
                               group.by = 'sexXnutr', 
                               ident.1 = 'M_Fed', 
                               ident.2 = 'M_Fast', 
                               logfc.threshold = log2(1.25), 
                               min.pct = 0.25, 
                               #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

α.Tanycytes.M.Fd.v.Fst$gene <- α.Tanycytes.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(α.Tanycytes.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
α.Tanycytes.M.Fd.v.Fst <- left_join(α.Tanycytes.M.Fd.v.Fst, temp01, by = "gene")



α.Tanycytes.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                             subset.ident = 'α-Tanycytes', 
                             group.by = 'sexXnutr', 
                             ident.1 = 'F_Fed', 
                             ident.2 = 'M_Fed', 
                             logfc.threshold = log2(1.25), 
                             min.pct = 0.25, 
                             #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

α.Tanycytes.Fd.F.v.M$gene <- α.Tanycytes.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(α.Tanycytes.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
α.Tanycytes.Fd.F.v.M <- left_join(α.Tanycytes.Fd.F.v.M, temp01, by = "gene")



α.Tanycytes.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                              subset.ident = 'α-Tanycytes', 
                              group.by = 'sexXnutr', 
                              ident.1 = 'F_Fast', 
                              ident.2 = 'M_Fast', 
                              logfc.threshold = log2(1.25), 
                              min.pct = 0.25, 
                              #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

α.Tanycytes.Fst.F.v.M$gene <- α.Tanycytes.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(α.Tanycytes.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
α.Tanycytes.Fst.F.v.M <- left_join(α.Tanycytes.Fst.F.v.M, temp01, by = "gene")



#beta-Tanycytes


β.Tanycytes.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                      subset.ident = 'β-Tanycytes', 
                                      group.by = 'sexXnutr', 
                                      ident.1 = 'F_Fed', 
                                      ident.2 = 'F_Fast', 
                                      logfc.threshold = log2(1.25), 
                                      min.pct = 0.25, 
                                      #pseudocount.use = 0,
                                      test.use = 'MAST',
                                      latent.vars = 'Sample_ID')

β.Tanycytes.F.Fd.v.Fst$gene <- β.Tanycytes.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(β.Tanycytes.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
β.Tanycytes.F.Fd.v.Fst <- left_join(β.Tanycytes.F.Fd.v.Fst, temp01, by = "gene")



β.Tanycytes.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                      subset.ident = 'β-Tanycytes', 
                                      group.by = 'sexXnutr', 
                                      ident.1 = 'M_Fed', 
                                      ident.2 = 'M_Fast', 
                                      logfc.threshold = log2(1.25), 
                                      min.pct = 0.25, 
                                      #pseudocount.use = 0,
                                      test.use = 'MAST',
                                      latent.vars = 'Sample_ID')

β.Tanycytes.M.Fd.v.Fst$gene <- β.Tanycytes.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(β.Tanycytes.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
β.Tanycytes.M.Fd.v.Fst <- left_join(β.Tanycytes.M.Fd.v.Fst, temp01, by = "gene")



β.Tanycytes.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                    subset.ident = 'β-Tanycytes', 
                                    group.by = 'sexXnutr', 
                                    ident.1 = 'F_Fed', 
                                    ident.2 = 'M_Fed', 
                                    logfc.threshold = log2(1.25), 
                                    min.pct = 0.25, 
                                    #pseudocount.use = 0,
                                    test.use = 'MAST',
                                    latent.vars = 'Sample_ID')

β.Tanycytes.Fd.F.v.M$gene <- β.Tanycytes.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(β.Tanycytes.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
β.Tanycytes.Fd.F.v.M <- left_join(β.Tanycytes.Fd.F.v.M, temp01, by = "gene")



β.Tanycytes.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                     subset.ident = 'β-Tanycytes', 
                                     group.by = 'sexXnutr', 
                                     ident.1 = 'F_Fast', 
                                     ident.2 = 'M_Fast', 
                                     logfc.threshold = log2(1.25), 
                                     min.pct = 0.25, 
                                     #pseudocount.use = 0,
                                     test.use = 'MAST',
                                     latent.vars = 'Sample_ID')

β.Tanycytes.Fst.F.v.M$gene <- β.Tanycytes.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(β.Tanycytes.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
β.Tanycytes.Fst.F.v.M <- left_join(β.Tanycytes.Fst.F.v.M, temp01, by = "gene")



#Ependymal


Ependymal.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                     subset.ident = 'Ependymal', 
                                     group.by = 'sexXnutr', 
                                     ident.1 = 'F_Fed', 
                                     ident.2 = 'F_Fast', 
                                     logfc.threshold = log2(1.25), 
                                     min.pct = 0.25, 
                                     #pseudocount.use = 0,
                                    test.use = 'MAST',
                                    latent.vars = 'Sample_ID')

Ependymal.F.Fd.v.Fst$gene <- Ependymal.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Ependymal.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Ependymal.F.Fd.v.Fst <- left_join(Ependymal.F.Fd.v.Fst, temp01, by = "gene")



Ependymal.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                     subset.ident = 'Ependymal', 
                                     group.by = 'sexXnutr', 
                                     ident.1 = 'M_Fed', 
                                     ident.2 = 'M_Fast', 
                                     logfc.threshold = log2(1.25), 
                                     min.pct = 0.25, 
                                     #pseudocount.use = 0,
                                    test.use = 'MAST',
                                    latent.vars = 'Sample_ID')

Ependymal.M.Fd.v.Fst$gene <- Ependymal.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Ependymal.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Ependymal.M.Fd.v.Fst <- left_join(Ependymal.M.Fd.v.Fst, temp01, by = "gene")



Ependymal.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                   subset.ident = 'Ependymal', 
                                   group.by = 'sexXnutr', 
                                   ident.1 = 'F_Fed', 
                                   ident.2 = 'M_Fed', 
                                   logfc.threshold = log2(1.25), 
                                   min.pct = 0.25, 
                                   #pseudocount.use = 0,
                                  test.use = 'MAST',
                                  latent.vars = 'Sample_ID')

Ependymal.Fd.F.v.M$gene <- Ependymal.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Ependymal.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Ependymal.Fd.F.v.M <- left_join(Ependymal.Fd.F.v.M, temp01, by = "gene")



Ependymal.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                    subset.ident = 'Ependymal', 
                                    group.by = 'sexXnutr', 
                                    ident.1 = 'F_Fast', 
                                    ident.2 = 'M_Fast', 
                                    logfc.threshold = log2(1.25), 
                                    min.pct = 0.25, 
                                    #pseudocount.use = 0,
                                   test.use = 'MAST',
                                   latent.vars = 'Sample_ID')

Ependymal.Fst.F.v.M$gene <- Ependymal.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Ependymal.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Ependymal.Fst.F.v.M <- left_join(Ependymal.Fst.F.v.M, temp01, by = "gene")


#Oligodendrocytes


Oligo.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                     subset.ident = 'Oligodendrocytes', 
                                     group.by = 'sexXnutr', 
                                     ident.1 = 'F_Fed', 
                                     ident.2 = 'F_Fast', 
                                     logfc.threshold = log2(1.25), 
                                     min.pct = 0.25, 
                                     #pseudocount.use = 0,
                                test.use = 'MAST',
                                latent.vars = 'Sample_ID')

Oligo.F.Fd.v.Fst$gene <- Oligo.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Oligo.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Oligo.F.Fd.v.Fst <- left_join(Oligo.F.Fd.v.Fst, temp01, by = "gene")

write.xlsx(Oligo.F.Fd.v.Fst, file = '../paper_figures/post_2025-01-06/DE_genes/Oligo/Oligo_F_Fd_v_Fst.xlsx')

Oligo.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                     subset.ident = 'Oligodendrocytes', 
                                     group.by = 'sexXnutr', 
                                     ident.1 = 'M_Fed', 
                                     ident.2 = 'M_Fast', 
                                     logfc.threshold = log2(1.25), 
                                     min.pct = 0.25, 
                                     #pseudocount.use = 0,
                                test.use = 'MAST',
                                latent.vars = 'Sample_ID')

Oligo.M.Fd.v.Fst$gene <- Oligo.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Oligo.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Oligo.M.Fd.v.Fst <- left_join(Oligo.M.Fd.v.Fst, temp01, by = "gene")


write.xlsx(Oligo.M.Fd.v.Fst, file = '../paper_figures/post_2025-01-06/DE_genes/Oligo/Oligo_M_Fd_v_Fst.xlsx')

Oligo.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                   subset.ident = 'Oligodendrocytes', 
                                   group.by = 'sexXnutr', 
                                   ident.1 = 'F_Fed', 
                                   ident.2 = 'M_Fed', 
                                   logfc.threshold = log2(1.25), 
                                   min.pct = 0.25, 
                                   #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

Oligo.Fd.F.v.M$gene <- Oligo.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Oligo.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Oligo.Fd.F.v.M <- left_join(Oligo.Fd.F.v.M, temp01, by = "gene")


write.xlsx(Oligo.Fd.F.v.M, file = '../paper_figures/post_2025-01-06/DE_genes/Oligo/Oligo_Fd_F_v_M.xlsx')

Oligo.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                    subset.ident = 'Oligodendrocytes', 
                                    group.by = 'sexXnutr', 
                                    ident.1 = 'F_Fast', 
                                    ident.2 = 'M_Fast', 
                                    logfc.threshold = log2(1.25), 
                                    min.pct = 0.25, 
                                    #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

Oligo.Fst.F.v.M$gene <- Oligo.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Oligo.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Oligo.Fst.F.v.M <- left_join(Oligo.Fst.F.v.M, temp01, by = "gene")

write.xlsx(Oligo.Fst.F.v.M, file = '../paper_figures/post_2025-01-06/DE_genes/Oligo/Oligo_Fst_F_v_M.xlsx')

#TOP


TOP.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                     subset.ident = 'TOP', 
                                     group.by = 'sexXnutr', 
                                     ident.1 = 'F_Fed', 
                                     ident.2 = 'F_Fast', 
                                     logfc.threshold = log2(1.25), 
                                     min.pct = 0.25, 
                                     #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

TOP.F.Fd.v.Fst$gene <- TOP.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(TOP.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)

temp01 <- temp01 |> rename('gene' = "external_gene_name")
TOP.F.Fd.v.Fst <- left_join(TOP.F.Fd.v.Fst, temp01, by = "gene")



TOP.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                     subset.ident = 'TOP', 
                                     group.by = 'sexXnutr', 
                                     ident.1 = 'M_Fed', 
                                     ident.2 = 'M_Fast', 
                                     logfc.threshold = log2(1.25), 
                                     min.pct = 0.25, 
                                     #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

TOP.M.Fd.v.Fst$gene <- TOP.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(TOP.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
TOP.M.Fd.v.Fst <- left_join(TOP.M.Fd.v.Fst, temp01, by = "gene")



TOP.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                   subset.ident = 'TOP', 
                                   group.by = 'sexXnutr', 
                                   ident.1 = 'F_Fed', 
                                   ident.2 = 'M_Fed', 
                                   logfc.threshold = log2(1.25), 
                                   min.pct = 0.25, 
                                   #pseudocount.use = 0,
                            test.use = 'MAST',
                            latent.vars = 'Sample_ID')

TOP.Fd.F.v.M$gene <- TOP.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(TOP.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
TOP.Fd.F.v.M <- left_join(TOP.Fd.F.v.M, temp01, by = "gene")



TOP.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                    subset.ident = 'TOP', 
                                    group.by = 'sexXnutr', 
                                    ident.1 = 'F_Fast', 
                                    ident.2 = 'M_Fast', 
                                    logfc.threshold = log2(1.25), 
                                    min.pct = 0.25, 
                                    #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

TOP.Fst.F.v.M$gene <- TOP.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(TOP.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
TOP.Fst.F.v.M <- left_join(TOP.Fst.F.v.M, temp01, by = "gene")



#OPC


OPC.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                     subset.ident = 'OPC', 
                                     group.by = 'sexXnutr', 
                                     ident.1 = 'F_Fed', 
                                     ident.2 = 'F_Fast', 
                                     logfc.threshold = log2(1.25), 
                                     min.pct = 0.25, 
                                     #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

OPC.F.Fd.v.Fst$gene <- OPC.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(OPC.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
OPC.F.Fd.v.Fst <- left_join(OPC.F.Fd.v.Fst, temp01, by = "gene")



OPC.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                     subset.ident = 'OPC', 
                                     group.by = 'sexXnutr', 
                                     ident.1 = 'M_Fed', 
                                     ident.2 = 'M_Fast', 
                                     logfc.threshold = log2(1.25), 
                                     min.pct = 0.25, 
                                     #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

OPC.M.Fd.v.Fst$gene <- OPC.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(OPC.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
OPC.M.Fd.v.Fst <- left_join(OPC.M.Fd.v.Fst, temp01, by = "gene")



OPC.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                   subset.ident = 'OPC', 
                                   group.by = 'sexXnutr', 
                                   ident.1 = 'F_Fed', 
                                   ident.2 = 'M_Fed', 
                                   logfc.threshold = log2(1.25), 
                                   min.pct = 0.25, 
                                   #pseudocount.use = 0,
                            test.use = 'MAST',
                            latent.vars = 'Sample_ID')

OPC.Fd.F.v.M$gene <- OPC.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(OPC.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
OPC.Fd.F.v.M <- left_join(OPC.Fd.F.v.M, temp01, by = "gene")



OPC.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                    subset.ident = 'OPC', 
                                    group.by = 'sexXnutr', 
                                    ident.1 = 'F_Fast', 
                                    ident.2 = 'M_Fast', 
                                    logfc.threshold = log2(1.25), 
                                    min.pct = 0.25, 
                                    #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

OPC.Fst.F.v.M$gene <- OPC.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(OPC.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
OPC.Fst.F.v.M <- left_join(OPC.Fst.F.v.M, temp01, by = "gene")


#Microglia


Microglia.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                     subset.ident = 'Microglia', 
                                     group.by = 'sexXnutr', 
                                     ident.1 = 'F_Fed', 
                                     ident.2 = 'F_Fast', 
                                     logfc.threshold = log2(1.25), 
                                     min.pct = 0.25, 
                                     #pseudocount.use = 0,
                                    test.use = 'MAST',
                                    latent.vars = 'Sample_ID')

Microglia.F.Fd.v.Fst$gene <- Microglia.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Microglia.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Microglia.F.Fd.v.Fst <- left_join(Microglia.F.Fd.v.Fst, temp01, by = "gene")


write.xlsx(Microglia.F.Fd.v.Fst, file = '../paper_figures/post_2025-01-06/DE_genes/Microglia/Microglia_F_Fd_v_Fst.xlsx')

Microglia.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                     subset.ident = 'Microglia', 
                                     group.by = 'sexXnutr', 
                                     ident.1 = 'M_Fed', 
                                     ident.2 = 'M_Fast', 
                                     logfc.threshold = log2(1.25), 
                                     min.pct = 0.25, 
                                     #pseudocount.use = 0,
                                    test.use = 'MAST',
                                    latent.vars = 'Sample_ID')

Microglia.M.Fd.v.Fst$gene <- Microglia.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Microglia.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Microglia.M.Fd.v.Fst <- left_join(Microglia.M.Fd.v.Fst, temp01, by = "gene")


write.xlsx(Microglia.M.Fd.v.Fst, file = '../paper_figures/post_2025-01-06/DE_genes/Microglia/Microglia_M_Fd_v_Fst.xlsx')

Microglia.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                   subset.ident = 'Microglia', 
                                   group.by = 'sexXnutr', 
                                   ident.1 = 'F_Fed', 
                                   ident.2 = 'M_Fed', 
                                   logfc.threshold = log2(1.25), 
                                   min.pct = 0.25, 
                                   #pseudocount.use = 0,
                                  test.use = 'MAST',
                                  latent.vars = 'Sample_ID')

Microglia.Fd.F.v.M$gene <- Microglia.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Microglia.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Microglia.Fd.F.v.M <- left_join(Microglia.Fd.F.v.M, temp01, by = "gene")


write.xlsx(Microglia.Fd.F.v.M, file = '../paper_figures/post_2025-01-06/DE_genes/Microglia/Microglia_Fd_F_v_M.xlsx')

Microglia.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                    subset.ident = 'Microglia', 
                                    group.by = 'sexXnutr', 
                                    ident.1 = 'F_Fast', 
                                    ident.2 = 'M_Fast', 
                                    logfc.threshold = log2(1.25), 
                                    min.pct = 0.25, 
                                    #pseudocount.use = 0,
                                   test.use = 'MAST',
                                   latent.vars = 'Sample_ID')

Microglia.Fst.F.v.M$gene <- Microglia.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Microglia.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Microglia.Fst.F.v.M <- left_join(Microglia.Fst.F.v.M, temp01, by = "gene")


write.xlsx(Microglia.Fst.F.v.M, file = '../paper_figures/post_2025-01-06/DE_genes/Microglia/Microglia_Fst_F_v_M.xlsx')

#Endothelial 


Endothelial.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                     subset.ident = 'Endothelial', 
                                     group.by = 'sexXnutr', 
                                     ident.1 = 'F_Fed', 
                                     ident.2 = 'F_Fast', 
                                     logfc.threshold = log2(1.25), 
                                     min.pct = 0.25, 
                                     #pseudocount.use = 0,
                                     test.use = 'MAST',
                                     latent.vars = 'Sample_ID')

Endothelial.F.Fd.v.Fst$gene <- Endothelial.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Endothelial.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Endothelial.F.Fd.v.Fst <- left_join(Endothelial.F.Fd.v.Fst, temp01, by = "gene")



Endothelial.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                     subset.ident = 'Endothelial', 
                                     group.by = 'sexXnutr', 
                                     ident.1 = 'M_Fed', 
                                     ident.2 = 'M_Fast', 
                                     logfc.threshold = log2(1.25), 
                                     min.pct = 0.25, 
                                     #pseudocount.use = 0,
                                     test.use = 'MAST',
                                     latent.vars = 'Sample_ID')

Endothelial.M.Fd.v.Fst$gene <- Endothelial.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Endothelial.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Endothelial.M.Fd.v.Fst <- left_join(Endothelial.M.Fd.v.Fst, temp01, by = "gene")



Endothelial.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                   subset.ident = 'Endothelial', 
                                   group.by = 'sexXnutr', 
                                   ident.1 = 'F_Fed', 
                                   ident.2 = 'M_Fed', 
                                   logfc.threshold = log2(1.25), 
                                   min.pct = 0.25, 
                                   #pseudocount.use = 0,
                                   test.use = 'MAST',
                                   latent.vars = 'Sample_ID')

Endothelial.Fd.F.v.M$gene <- Endothelial.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Endothelial.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Endothelial.Fd.F.v.M <- left_join(Endothelial.Fd.F.v.M, temp01, by = "gene")



Endothelial.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                    subset.ident = 'Endothelial', 
                                    group.by = 'sexXnutr', 
                                    ident.1 = 'F_Fast', 
                                    ident.2 = 'M_Fast', 
                                    logfc.threshold = log2(1.25), 
                                    min.pct = 0.25, 
                                    #pseudocount.use = 0,
                                    test.use = 'MAST',
                                    latent.vars = 'Sample_ID')

Endothelial.Fst.F.v.M$gene <- Endothelial.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(Endothelial.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
Endothelial.Fst.F.v.M <- left_join(Endothelial.Fst.F.v.M, temp01, by = "gene")




#VLMC

VLMC.F.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                     subset.ident = 'VLMC', 
                                     group.by = 'sexXnutr', 
                                     ident.1 = 'F_Fed', 
                                     ident.2 = 'F_Fast', 
                                     logfc.threshold = log2(1.25), 
                                     min.pct = 0.25, 
                                     #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

VLMC.F.Fd.v.Fst$gene <- VLMC.F.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(VLMC.F.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
VLMC.F.Fd.v.Fst <- left_join(VLMC.F.Fd.v.Fst, temp01, by = "gene")



VLMC.M.Fd.v.Fst <- FindMarkers(ARH_Sex_by_Nutr, 
                                     subset.ident = 'VLMC', 
                                     group.by = 'sexXnutr', 
                                     ident.1 = 'M_Fed', 
                                     ident.2 = 'M_Fast', 
                                     logfc.threshold = log2(1.25), 
                                     min.pct = 0.25, 
                                     #pseudocount.use = 0,
                               test.use = 'MAST',
                               latent.vars = 'Sample_ID')

VLMC.M.Fd.v.Fst$gene <- VLMC.M.Fd.v.Fst |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(VLMC.M.Fd.v.Fst), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
VLMC.M.Fd.v.Fst <- left_join(VLMC.M.Fd.v.Fst, temp01, by = "gene")



VLMC.Fd.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                   subset.ident = 'VLMC', 
                                   group.by = 'sexXnutr', 
                                   ident.1 = 'F_Fed', 
                                   ident.2 = 'M_Fed', 
                                   logfc.threshold = log2(1.25), 
                                   min.pct = 0.25, 
                                   #pseudocount.use = 0,
                             test.use = 'MAST',
                             latent.vars = 'Sample_ID')

VLMC.Fd.F.v.M$gene <- VLMC.Fd.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(VLMC.Fd.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
VLMC.Fd.F.v.M <- left_join(VLMC.Fd.F.v.M, temp01, by = "gene")



VLMC.Fst.F.v.M <- FindMarkers(ARH_Sex_by_Nutr, 
                                    subset.ident = 'VLMC', 
                                    group.by = 'sexXnutr', 
                                    ident.1 = 'F_Fast', 
                                    ident.2 = 'M_Fast', 
                                    logfc.threshold = log2(1.25), 
                                    min.pct = 0.25, 
                                    #pseudocount.use = 0,
                              test.use = 'MAST',
                              latent.vars = 'Sample_ID')

VLMC.Fst.F.v.M$gene <- VLMC.Fst.F.v.M |> row.names()
temp01 <- getBM(attributes = c("external_gene_name",'chromosome_name','description') , 
                filters = "external_gene_name", values = row.names(VLMC.Fst.F.v.M), mart = ensembl,
                checkFilters = TRUE, verbose = TRUE, uniqueRows = TRUE, bmHeader = FALSE)
temp01 <- temp01 |> rename('gene' = "external_gene_name")
VLMC.Fst.F.v.M <- left_join(VLMC.Fst.F.v.M, temp01, by = "gene")



cell_types <- c(
  

'Agrp',
'Pomc',
'KNDy',
'DA',
'Ghrh.Chat',
'Sst.Unc13c',
'Lamp5.Npy5r',
'Lef1',
'Htr3b',
'Gad2.Htr2c',
'Ros1.Alk',
'Satb2.Slc18a2',
'Coch.Slc18a2',
'Tbx19',
'Tbx15',
'Ebf3.Htr2c',
'Slc17a6.Alk',
'Erg.Lepr',
'Tac1.Reln',
'Klhl1.Ebf3',
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
'ParsTub',
'Astrocytes',
'α.Tanycytes',
'β.Tanycytes',
'Ependymal',
'Oligodendrocytes',
'TOP',
'OPC',
'Microglia',
'Endothelial',
'VLMC'

)

#####2024-05-02#####

Sex_by_Nutr_DE_nums <- tibble(cell_type = 'a', comparison = 'a', condition = 'a', DE = 0, .rows = 252)
row.number <- 0
{
 {
    #Agrp
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Agrp' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Agrp.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Agrp.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
  
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Agrp' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Agrp.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Agrp.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Agrp' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Agrp.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                     (Agrp.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Agrp' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Agrp.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Agrp.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Agrp' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Agrp.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Agrp.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Agrp' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Agrp.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Agrp.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #Pomc
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Pomc' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Pomc.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Pomc.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Pomc' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Pomc.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Pomc.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Pomc' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Pomc.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Pomc.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Pomc' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Pomc.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Pomc.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Pomc' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Pomc.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Pomc.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Pomc' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Pomc.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Pomc.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #KNDy
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'KNDy' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((KNDy.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (KNDy.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'KNDy' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((KNDy.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (KNDy.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'KNDy' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((KNDy.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (KNDy.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'KNDy' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((KNDy.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (KNDy.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'KNDy' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((KNDy.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (KNDy.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'KNDy' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((KNDy.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (KNDy.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    
    #DA
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'DA' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((DA.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (DA.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'DA' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((DA.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (DA.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'DA' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((DA.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (DA.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'DA' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((DA.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (DA.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'DA' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((DA.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (DA.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'DA' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((DA.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (DA.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #Ghrh/Chat
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ghrh/Chat' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Ghrh.Chat.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Ghrh.Chat.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ghrh/Chat' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Ghrh.Chat.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Ghrh.Chat.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ghrh/Chat' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Ghrh.Chat.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Ghrh.Chat.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ghrh/Chat' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Ghrh.Chat.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Ghrh.Chat.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ghrh/Chat' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Ghrh.Chat.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Ghrh.Chat.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ghrh/Chat' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Ghrh.Chat.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Ghrh.Chat.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    
    #Sst/Unc13c
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Sst/Unc13c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Sst.Unc13c.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Sst.Unc13c.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Sst/Unc13c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Sst.Unc13c.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Sst.Unc13c.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Sst/Unc13c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Sst.Unc13c.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Sst.Unc13c.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Sst/Unc13c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Sst.Unc13c.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Sst.Unc13c.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Sst/Unc13c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Sst.Unc13c.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Sst.Unc13c.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Sst/Unc13c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Sst.Unc13c.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Sst.Unc13c.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #Lamp5/Npy5r
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Lamp5/Npy5r' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Lamp5.Npy5r.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Lamp5.Npy5r.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Lamp5/Npy5r' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Lamp5.Npy5r.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Lamp5.Npy5r.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Lamp5/Npy5r' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Lamp5.Npy5r.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Lamp5.Npy5r.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Lamp5/Npy5r' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Lamp5.Npy5r.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Lamp5.Npy5r.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Lamp5/Npy5r' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Lamp5.Npy5r.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Lamp5.Npy5r.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Lamp5/Npy5r' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Lamp5.Npy5r.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Lamp5.Npy5r.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #Lef1
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Lef1' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Lef1.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Lef1.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Lef1' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Lef1.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Lef1.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Lef1' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Lef1.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Lef1.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Lef1' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Lef1.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Lef1.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Lef1' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Lef1.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Lef1.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Lef1' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Lef1.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Lef1.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #Htr3b
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Htr3b' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Htr3b.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Htr3b.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Htr3b' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Htr3b.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Htr3b.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Htr3b' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Htr3b.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Htr3b.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Htr3b' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Htr3b.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Htr3b.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Htr3b' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Htr3b.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Htr3b.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Htr3b' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Htr3b.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Htr3b.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    #Gad2/Htr2c
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Gad2/Htr2c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Gad2.Htr2c.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Gad2.Htr2c.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Gad2/Htr2c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Gad2.Htr2c.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Gad2.Htr2c.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Gad2/Htr2c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Gad2.Htr2c.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Gad2.Htr2c.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Gad2/Htr2c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Gad2.Htr2c.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Gad2.Htr2c.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Gad2/Htr2c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Gad2.Htr2c.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Gad2.Htr2c.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Gad2/Htr2c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Gad2.Htr2c.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Gad2.Htr2c.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    #Ros1/Alk
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ros1/Alk' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Ros1.Alk.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Ros1.Alk.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ros1/Alk' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Ros1.Alk.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Ros1.Alk.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ros1/Alk' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Ros1.Alk.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Ros1.Alk.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ros1/Alk' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Ros1.Alk.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Ros1.Alk.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ros1/Alk' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Ros1.Alk.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Ros1.Alk.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ros1/Alk' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Ros1.Alk.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Ros1.Alk.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    #Satb2/Slc18a2
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Satb2/Slc18a2' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Satb2.Slc18a2.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Satb2.Slc18a2.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Satb2/Slc18a2' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Satb2.Slc18a2.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Satb2.Slc18a2.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Satb2/Slc18a2' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Satb2.Slc18a2.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Satb2.Slc18a2.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Satb2/Slc18a2' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Satb2.Slc18a2.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Satb2.Slc18a2.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Satb2/Slc18a2' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Satb2.Slc18a2.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Satb2.Slc18a2.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Satb2/Slc18a2' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Satb2.Slc18a2.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Satb2.Slc18a2.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    
    #Coch/Slc18a2
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Coch/Slc18a2' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Coch.Slc18a2.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Coch.Slc18a2.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Coch/Slc18a2' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Coch.Slc18a2.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Coch.Slc18a2.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Coch/Slc18a2' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Coch.Slc18a2.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Coch.Slc18a2.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Coch/Slc18a2' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Coch.Slc18a2.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Coch.Slc18a2.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Coch/Slc18a2' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Coch.Slc18a2.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Coch.Slc18a2.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Coch/Slc18a2' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Coch.Slc18a2.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Coch.Slc18a2.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #Tbx19
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tbx19' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Tbx19.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Tbx19.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tbx19' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Tbx19.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Tbx19.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tbx19' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Tbx19.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Tbx19.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tbx19' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Tbx19.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Tbx19.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tbx19' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Tbx19.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Tbx19.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tbx19' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Tbx19.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Tbx19.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #Tbx15
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tbx15' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Tbx15.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Tbx15.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tbx15' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Tbx15.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Tbx15.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tbx15' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Tbx15.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Tbx15.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tbx15' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Tbx15.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Tbx15.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tbx15' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Tbx15.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Tbx15.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tbx15' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Tbx15.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Tbx15.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #Ebf3/Htr2c
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ebf3/Htr2c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Ebf3.Htr2c.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Ebf3.Htr2c.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ebf3/Htr2c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Ebf3.Htr2c.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Ebf3.Htr2c.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ebf3/Htr2c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Ebf3.Htr2c.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Ebf3.Htr2c.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ebf3/Htr2c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Ebf3.Htr2c.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Ebf3.Htr2c.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ebf3/Htr2c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Ebf3.Htr2c.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Ebf3.Htr2c.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ebf3/Htr2c' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Ebf3.Htr2c.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Ebf3.Htr2c.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #Slc17a6/Alk
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Slc17a6/Alk' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Slc17a6.Alk.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Slc17a6.Alk.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Slc17a6/Alk' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Slc17a6.Alk.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Slc17a6.Alk.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Slc17a6/Alk' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Slc17a6.Alk.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Slc17a6.Alk.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Slc17a6/Alk' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Slc17a6.Alk.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Slc17a6.Alk.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Slc17a6/Alk' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Slc17a6.Alk.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Slc17a6.Alk.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Slc17a6/Alk' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Slc17a6.Alk.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Slc17a6.Alk.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #Erg/Lepr
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Erg/Lepr' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Erg.Lepr.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Erg.Lepr.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Erg/Lepr' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Erg.Lepr.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Erg.Lepr.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Erg/Lepr' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Erg.Lepr.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Erg.Lepr.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Erg/Lepr' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Erg.Lepr.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Erg.Lepr.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Erg/Lepr' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Erg.Lepr.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Erg.Lepr.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Erg/Lepr' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Erg.Lepr.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Erg.Lepr.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #Tac1/Reln
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tac1/Reln' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Tac1.Reln.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Tac1.Reln.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tac1/Reln' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Tac1.Reln.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Tac1.Reln.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tac1/Reln' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Tac1.Reln.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Tac1.Reln.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tac1/Reln' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Tac1.Reln.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Tac1.Reln.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tac1/Reln' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Tac1.Reln.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Tac1.Reln.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tac1/Reln' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Tac1.Reln.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Tac1.Reln.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #Klhl1/Ebf3
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Klhl1/Ebf3' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Klhl1.Ebf3.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Klhl1.Ebf3.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Klhl1/Ebf3' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Klhl1.Ebf3.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Klhl1.Ebf3.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Klhl1/Ebf3' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Klhl1.Ebf3.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Klhl1.Ebf3.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Klhl1/Ebf3' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Klhl1.Ebf3.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Klhl1.Ebf3.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Klhl1/Ebf3' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Klhl1.Ebf3.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Klhl1.Ebf3.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Klhl1/Ebf3' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Klhl1.Ebf3.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Klhl1.Ebf3.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #SCN
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'SCN' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((SCN.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (SCN.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'SCN' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((SCN.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (SCN.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'SCN' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((SCN.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (SCN.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'SCN' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((SCN.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (SCN.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'SCN' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((SCN.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (SCN.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'SCN' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((SCN.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (SCN.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
  #VMH.01
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VMH.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((VMH.01.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (VMH.01.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VMH.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((VMH.01.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (VMH.01.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VMH.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((VMH.01.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (VMH.01.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VMH.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((VMH.01.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (VMH.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VMH.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((VMH.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (VMH.01.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VMH.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((VMH.01.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (VMH.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    #VMH.02
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VMH.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((VMH.02.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (VMH.02.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VMH.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((VMH.02.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (VMH.02.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VMH.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((VMH.02.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (VMH.02.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VMH.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((VMH.02.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (VMH.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VMH.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((VMH.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (VMH.02.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VMH.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((VMH.02.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (VMH.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #PVp.01
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((PVp.01.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (PVp.01.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((PVp.01.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (PVp.01.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((PVp.01.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (PVp.01.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((PVp.01.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (PVp.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((PVp.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (PVp.01.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((PVp.01.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (PVp.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #PVp.02
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((PVp.02.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (PVp.02.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((PVp.02.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (PVp.02.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((PVp.02.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (PVp.02.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((PVp.02.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (PVp.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((PVp.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (PVp.02.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((PVp.02.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (PVp.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
   #PVp.03
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.03' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((PVp.03.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (PVp.03.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.03' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((PVp.03.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (PVp.03.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.03' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((PVp.03.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (PVp.03.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.03' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((PVp.03.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (PVp.03.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.03' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((PVp.03.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (PVp.03.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.03' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((PVp.03.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (PVp.03.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #PVp.04
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.04' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((PVp.04.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (PVp.04.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.04' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((PVp.04.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (PVp.04.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.04' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((PVp.04.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (PVp.04.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.04' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((PVp.04.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (PVp.04.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.04' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((PVp.04.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (PVp.04.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'PVp.04' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((PVp.04.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (PVp.04.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #MM.01
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((MM.01.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (MM.01.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((MM.01.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (MM.01.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((MM.01.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (MM.01.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((MM.01.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (MM.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((MM.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (MM.01.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.01' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((MM.01.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (MM.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #MM.02
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((MM.02.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (MM.02.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((MM.02.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (MM.02.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((MM.02.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (MM.02.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((MM.02.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (MM.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((MM.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (MM.02.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.02' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((MM.02.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (MM.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #MM.03
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.03' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((MM.03.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (MM.03.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.03' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((MM.03.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (MM.03.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.03' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((MM.03.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (MM.03.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.03' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((MM.03.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (MM.03.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.03' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((MM.03.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (MM.03.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'MM.03' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((MM.03.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (MM.03.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    #Tu
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tu' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Tu.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Tu.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tu' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Tu.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Tu.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tu' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Tu.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Tu.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tu' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Tu.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Tu.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tu' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Tu.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Tu.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Tu' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Tu.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Tu.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
  #Pars Tuberalis
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Pars Tuberalis' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((ParsTub.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (ParsTub.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Pars Tuberalis' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((ParsTub.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (ParsTub.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Pars Tuberalis' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((ParsTub.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (ParsTub.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Pars Tuberalis' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((ParsTub.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (ParsTub.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Pars Tuberalis' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((ParsTub.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (ParsTub.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Pars Tuberalis' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((ParsTub.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (ParsTub.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #Astrocytes
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Astrocytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Astrocytes.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Astrocytes.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Astrocytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Astrocytes.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Astrocytes.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Astrocytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Astrocytes.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Astrocytes.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Astrocytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Astrocytes.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Astrocytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Astrocytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Astrocytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Astrocytes.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Astrocytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Astrocytes.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Astrocytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    
   #alpha-Tanycytes
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'α-Tanycytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((α.Tanycytes.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (α.Tanycytes.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'α-Tanycytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((α.Tanycytes.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (α.Tanycytes.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'α-Tanycytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((α.Tanycytes.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (α.Tanycytes.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'α-Tanycytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((α.Tanycytes.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (α.Tanycytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'α-Tanycytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((α.Tanycytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (α.Tanycytes.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'α-Tanycytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((α.Tanycytes.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (α.Tanycytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    #Beta-Tanycytes
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'β-Tanycytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((β.Tanycytes.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (β.Tanycytes.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'β-Tanycytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((β.Tanycytes.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (β.Tanycytes.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'β-Tanycytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((β.Tanycytes.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (β.Tanycytes.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'β-Tanycytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((β.Tanycytes.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (β.Tanycytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'β-Tanycytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((β.Tanycytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (β.Tanycytes.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'β-Tanycytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((β.Tanycytes.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (β.Tanycytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    #Ependymal
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ependymal' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Ependymal.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Ependymal.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ependymal' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Ependymal.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Ependymal.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ependymal' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Ependymal.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Ependymal.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ependymal' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Ependymal.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Ependymal.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ependymal' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Ependymal.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Ependymal.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Ependymal' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Ependymal.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Ependymal.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    #oligodendrocytes
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Oligodendrocytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Oligo.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Oligo.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Oligodendrocytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Oligo.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Oligo.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Oligodendrocytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Oligo.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Oligo.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Oligodendrocytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Oligo.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Oligo.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Oligodendrocytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Oligo.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Oligo.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Oligodendrocytes' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Oligo.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Oligo.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #TOP
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'TOP' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((TOP.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (TOP.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'TOP' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((TOP.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (TOP.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'TOP' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((TOP.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (TOP.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'TOP' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((TOP.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (TOP.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'TOP' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((TOP.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (TOP.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'TOP' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((TOP.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (TOP.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
  
    
    #OPC
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'OPC' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((OPC.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (OPC.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'OPC' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((OPC.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (OPC.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'OPC' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((OPC.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (OPC.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'OPC' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((OPC.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (OPC.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'OPC' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((OPC.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (OPC.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'OPC' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((OPC.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (OPC.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
  
    
    #Microglia
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Microglia' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Microglia.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Microglia.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Microglia' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Microglia.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Microglia.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Microglia' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Microglia.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Microglia.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Microglia' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Microglia.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Microglia.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Microglia' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Microglia.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Microglia.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Microglia' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Microglia.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Microglia.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    #Endothelial
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Endothelial' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Endothelial.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Endothelial.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Endothelial' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Endothelial.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (Endothelial.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Endothelial' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Endothelial.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (Endothelial.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Endothelial' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Endothelial.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Endothelial.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Endothelial' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((Endothelial.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Endothelial.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'Endothelial' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((Endothelial.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Endothelial.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    #VLMC
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VLMC' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Female'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((VLMC.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (VLMC.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VLMC' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Male'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((VLMC.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                   (VLMC.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VLMC' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Nurt'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'F.M.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((VLMC.F.Fd.v.Fst |> filter(p_val_adj < 0.05)), 
                                                    (VLMC.M.Fd.v.Fst |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VLMC' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((VLMC.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (VLMC.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VLMC' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_nums$DE[row.number] = anti_join((VLMC.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (VLMC.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_nums$cell_type[row.number] = 'VLMC' 
    Sex_by_Nutr_DE_nums$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_nums$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_nums$DE[row.number] = inner_join((VLMC.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (VLMC.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count()
    
    
    
    }
  
  
}

Sex_by_Nutr_DE_nums$DE <- Sex_by_Nutr_DE_nums$DE |> as.numeric()

Sex_by_Nutr_DE_nums <- Sex_by_Nutr_DE_nums %>% group_by(cell_type) %>% mutate(total = sum(DE))


Sex_by_Nutr_DE_nums |> filter(comparison == 'Nurt') |> 
  ggplot(aes(x = cell_type, y = DE)) + 
  geom_col(aes(fill = factor(condition, levels = c('Female','F.M.Ovlp','Male'))), color = 'black', size = 0.25) + 
  scale_fill_manual(values = c('#f8c471','#D35400','#E67E22'), name = '') +
  scale_x_discrete(limits = c('KNDy',
                              'Agrp',
                              'Lamp5/Npy5r',
                              'β-Tanycytes',
                              'Ghrh/Chat',
                              'Astrocytes',
                              'Pomc',
                              'MM.01',
                              'Oligodendrocytes',
                              'α-Tanycytes',
                              'DA',
                              'VMH.01',
                              'Sst/Unc13c',
                              'Gad2/Htr2c',
                              'Microglia',
                              'Coch/Slc18a2',
                              'Htr3b',
                              'Lef1',
                              'Ependymal',
                              'PVp.01',
                              'Tbx19',
                              'Erg/Lepr',
                              'Pars Tuberalis',
                              'Slc17a6/Alk',
                              'PVp.02',
                              'OPC',
                              'VLMC',
                              'PVp.03',
                              'MM.02',
                              'Satb2/Slc18a2',
                              'MM.03',
                              'Klhl1/Ebf3',
                              'Ebf3/Htr2c',
                              'VMH.02',
                              'Ros1/Alk',
                              'PVp.04',
                              'SCN',
                              'Tu',
                              'Tac1/Reln',
                              'Endothelial',
                              'TOP',
                              'Tbx15'
  ),
                   labels = function(x) {
                     max_length <- max(nchar(x))
                     lapply(x, function(label) {
                       diff_length <- max_length - nchar(label)
                       if (diff_length > 0) {
                         num_periods <- diff_length %/% 2 + 4
                         left_padding <- paste(rep(".", num_periods), collapse = "")
                         right_padding <- paste(rep(".", diff_length - num_periods + 7), collapse = "")
                         return(paste0(left_padding, label, right_padding))
                       } else {
                         return(label)
                       }
                     })
                   }) +
  labs(x = '', y= 'Nutritionally Regulated DE Genes', fill = '') +
  coord_cartesian(ylim = c(0,325)) +
  theme_classic() +
  theme(text = element_text(family = "Arial", size = 6, color = 'black'),
        axis.text.x= element_text(family = 'Arial', color = 'black', size = 6,angle = 90, hjust = 1, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y= element_text(family = 'Arial', color = 'black', size = 6),
        legend.key.height = unit(.1, 'in'),
        legend.direction = 'horizontal',
        legend.position = c(0.5,0.75))
ggsave(filename = 'figures/celltype_number_of_Nutr_DE3.tiff', device = 'tiff', units = 'in', width = 7.24, height = 2.25, dpi = 600)


######
Sex_by_Nutr_DE_nums |> filter(comparison == 'Sex') |> 
  ggplot(aes(x = cell_type, y = DE)) + 
  geom_col(aes(fill = factor(condition, levels = c('Fed','Fd.Fst.Ovlp','Fasted'))), color = 'black', size = 0.25) + 
  scale_fill_manual(values = c('#d7bde2','#8E44AD','#A569BD'), name = '') +
  scale_x_discrete(limits = c('KNDy',
                              'Agrp',
                              'Lamp5/Npy5r',
                              'β-Tanycytes',
                              'Ghrh/Chat',
                              'Astrocytes',
                              'Pomc',
                              'MM.01',
                              'Oligodendrocytes',
                              'α-Tanycytes',
                              'DA',
                              'VMH.01',
                              'Sst/Unc13c',
                              'Gad2/Htr2c',
                              'Microglia',
                              'Coch/Slc18a2',
                              'Htr3b',
                              'Lef1',
                              'Ependymal',
                              'PVp.01',
                              'Tbx19',
                              'Erg/Lepr',
                              'Pars Tuberalis',
                              'Slc17a6/Alk',
                              'PVp.02',
                              'OPC',
                              'VLMC',
                              'PVp.03',
                              'MM.02',
                              'Satb2/Slc18a2',
                              'MM.03',
                              'Klhl1/Ebf3',
                              'Ebf3/Htr2c',
                              'VMH.02',
                              'Ros1/Alk',
                              'PVp.04',
                              'SCN',
                              'Tu',
                              'Tac1/Reln',
                              'Endothelial',
                              'TOP',
                              'Tbx15'),
                   labels = function(x) {
                     max_length <- max(nchar(x))
                     lapply(x, function(label) {
                       diff_length <- max_length - nchar(label)
                       if (diff_length > 0) {
                         num_periods <- diff_length %/% 2 + 4
                         left_padding <- paste(rep(".", num_periods), collapse = "")
                         right_padding <- paste(rep(".", diff_length - num_periods + 7), collapse = "")
                         return(paste0(left_padding, label, right_padding))
                       } else {
                         return(label)
                       }
                     })
                   }
                   , position="top") +
  labs(x = '', y= 'Sexually Regulated DE Genes', fill = '') +
  scale_y_reverse() +
  coord_cartesian(ylim = c(325,0)) +
  theme_classic() +
  theme(text = element_text(family = "Arial", size = 6, color = 'black'),
        axis.text.x= element_text(family = 'Arial', color = 'black', size = 6,angle = 90, hjust = 1, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y= element_text(family = 'Arial', color = 'black', size = 6),
        legend.key.height = unit(.1, 'in'),
        legend.direction = 'horizontal',
        legend.position = c(0.5,0.25))
ggsave(filename = 'figures/celltype_number_of_Sex_DE3.tiff', device = 'tiff', units = 'in', width = 7.24, height = 2.25, dpi = 600)







#####2024-05-03######


# hdWGCNA analysis
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)

theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 8)


##load obj
seurat_obj <- readRDS('GSE282955_ARH_Sex_by_Nutr.rds')

#setup 
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", 
  fraction = 0.05, 
  wgcna_name = "ARH" 
)

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type3", "Batch"), 
  reduction = 'pca', 
  k = 25, 
  ident.group = 'cell_type3' 
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

metacell_obj <- GetMetacellObject(seurat_obj)


##set up the expression matrix
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "Agrp", 
  group.by='cell_type3', 
  assay = 'SCT', 
  layer = 'data' 
)

#Select soft-power threshold
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' 
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)


pdf('Agrp_all_softpower.pdf',height = 7, width = 10)
PlotSoftPowers(seurat_obj)
plot_list

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
dev.off()

power_table <- GetPowerTable(seurat_obj)
head(power_table)

write.csv(power_table,'Agrp_all_power_table.csv',row.names = F) 

#Construct co-expression network
# Reset module names and colors
library(MetBrewer)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module)
mod_colors <- dplyr::select(modules, c(module, color)) %>%
  distinct %>% arrange(module) %>% .$color
n_colors <- length(mod_colors) -1

new_colors <- paste0(met.brewer("Signac", n=n_colors, type='discrete'))
new_colors <- sample(new_colors)
seurat_obj <- ResetModuleColors(seurat_obj, new_colors)

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=10,
  setDatExpr=FALSE,
  tom_name = 'Agrp' 
)

pdf("Agrp_all_hdWGCNA_Dendrogram.pdf", width = 12, height = 6)
PlotDendrogram(seurat_obj, main='Agrp hdWGCNA Dendrogram')
p
dev.off()

png("Agrp_all_hdWGCNA_Dendrogram.png", width =4600, height = 2400, res = 300)
PlotDendrogram(seurat_obj, main='Agrp hdWGCNA Dendrogram')
p
dev.off()

#Compute harmonized module eigengenes


# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="Batch"
)

saveRDS(seurat_obj,'all_Agrp_WGCNA.rds')

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

#Compute module connectivity
# compute eigengene-based connectivity (kME):
seurat_obj<- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type3', group_name = 'Agrp'
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Agrp"
)

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=4)

p

ggsave("Agrp_PlotKMEs.pdf", p, width = 9, height = 5,dpi = 300)
ggsave("Agrp_PlotKMEs.png", p, width = 9, height = 5,dpi = 300)

# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)


# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

#Module Feature Plots
# make a featureplot of hMEs for each module
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)  # Or whatever dims are appropriate
pdf('Agrp_all_ModuleFeaturePlot_hub_hMEs.pdf',height = 7, width = 5.5)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs',# plot the hMEs
  reduction = 'tsne',
  order=TRUE, # order so the points with highest hMEs are on top
  point_size =0.1)

# stitch together with patchwork
wrap_plots(plot_list, ncol=4) + plot_annotation(title = "hME Module Expression")

dev.off()

png('Agrp_all_ModuleFeaturePlot_hub_hMEs.png',height = 700, width = 1000)
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs',# plot the hMEs
  reduction = 'tsne',
  order=TRUE, # order so the points with highest hMEs are on top
  point_size =0.1)

# stitch together with patchwork
wrap_plots(plot_list, ncol=4) + plot_annotation(title = "hME Module Expression")

dev.off()

# make a featureplot of hub scores for each module
pdf('Agrp_all_ModuleFeaturePlot_hub_score.pdf',height = 7, width = 10)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle',
  reduction = 'tsne',# order so cells are shuffled
  ucell = TRUE,# depending on Seurat vs UCell for gene scoring
  point_size =0.1)
wrap_plots(plot_list, ncol=4)+ plot_annotation(title = "Hub genes signature scores")

dev.off()

png('Agrp_all_ModuleFeaturePlot_hub_score.png',height = 700, width = 1000)
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle',
  reduction = 'tsne',# order so cells are shuffled
  ucell = TRUE,# depending on Seurat vs UCell for gene scoring
  point_size =0.1)
wrap_plots(plot_list, ncol=4)+ plot_annotation(title = "Hub genes signature scores")

dev.off()

# plot module correlagram
pdf('Agrp_all_module_correlagram.pdf',height = 7, width = 7)
ModuleCorrelogram(seurat_obj)
dev.off()

png('Agrp_all_module_correlagram.png',height = 700, width = 700)
ModuleCorrelogram(seurat_obj)
dev.off()

# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
#plot with Seurat DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'cell_type3')

p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
ggsave('Agrp_all_allcelltype_dotplot_module.pdf',p,height = 8,width = 8,dpi = 300)
ggsave('Agrp_all_allcelltype_dotplot_module.png',p,height = 8,width = 8,dpi = 300)



# NetWork_Visualization

theme_set(theme_cowplot())
set.seed(12345)

ModuleNetworkPlot(
  seurat_obj,
  outdir = 'Agrp_all_ModuleNetworks'
)

ModuleNetworkPlot(
  seurat_obj, 
  outdir='Agrp_all_ModuleNetworks2', # new folder name
  n_inner = 20, # number of genes in inner ring
  n_outer = 30, # number of genes in outer ring
  n_conns = Inf, # show all of the connections
  plot_size=c(10,10), # larger plotting area
  vertex.label.cex=1 # font size
)

options(future.globals.maxSize = 50 * 1024^3)  # 50 GiB

png("Agrp_All_HubGeneNetworkPlot_all.png", width = 2400, height = 2400, res = 300)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other = 5,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()

# Run HubGeneNetworkPlot again, this time only selecting 5 specific modules:
# get the list of modules:
modules <- GetModules(fast)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# hubgene network
png("HubGeneNetworkPlot_5module.png", width = 2400, height = 2400, res = 300)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 10, n_other=20,
  edge_prop = 0.75,
  mods = mods[1:5] # only select 5 modules
)
dev.off()

##Applying UMAP to co-expression networks
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, 
  n_neighbors=15, 
  min_dist=0.3 
)

#Next, make a simple visualization of the UMAP using ggplot2:
# get the hub gene UMAP table from the seurat object
graphics.off()
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
png("Agrp_All_ModuleUMAPPlot.png", width = 2400, height = 2400, res = 300)
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, 
    size=umap_df$kME*2 
  ) +
  umap_theme()

ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.15, 
  label_hubs=2 ,
  keep_grey_edges=FALSE
  
)
dev.off()

ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()
pdf("Agrp_All_ModuleUMAPPlot.pdf", width = 12, height = 12)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=3 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)
dev.off()

png("Agrp_All_ModuleUMAPPlot.png", width = 4800, height = 5000,res=500)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=3 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)
dev.off()

#DME analysis comparing two groups

theme_set(theme_cowplot())

##Fast
group1 <- seurat_obj@meta.data %>% subset(sexXnutr == 'F_Fast' & cell_type3 == 'Agrp') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(sexXnutr == 'M_Fast' & cell_type3 == 'Agrp') %>% rownames

head(group1)

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='ARH'
)

head(DMEs)


pdf("Agrp_All_DMEsLollipop_Fast_F_M.pdf", width = 6, height = 7)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name='ARH', 
  pvalue = "p_val_adj"
)
dev.off()

png("Agrp_All_DMEsLollipop_Fast_F_M.png", width = 1200, height = 1400,res = 300)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name='ARH', 
  pvalue = "p_val_adj"
)
dev.off()


##Fed
group1 <- seurat_obj@meta.data %>% subset(sexXnutr == 'F_Fed' & cell_type3 == 'Agrp') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(sexXnutr == 'M_Fed' & cell_type3 == 'Agrp') %>% rownames

head(group1)
head(group2)

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='ARH'
)

head(DMEs)

pdf("Agrp_All_DMEsLollipop_Fed_F_M.pdf", width = 6, height = 7)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name='ARH', 
  pvalue = "p_val_adj"
)
dev.off()

png("Agrp_All_DMEsLollipop_Fed_F_M.png", width = 1200, height = 1400,res = 300)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name='ARH', 
  pvalue = "p_val_adj"
)
dev.off()

##Female
group1 <- seurat_obj@meta.data %>% subset(sexXnutr == 'F_Fast' & cell_type3 == 'Agrp') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(sexXnutr == 'F_Fed' & cell_type3 == 'Agrp') %>% rownames

head(group1)
head(group2)

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='ARH'
)

head(DMEs)

pdf("Agrp_All_DMEsLollipop_Female_Fast_Fed.pdf", width = 6, height = 7)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name='ARH', 
  pvalue = "p_val_adj"
)
dev.off()

png("Agrp_All_DMEsLollipop_Female_Fast_Fed.png", width = 1200, height = 1400,res = 300)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name='ARH', 
  pvalue = "p_val_adj"
)
dev.off()

##Male
group1 <- seurat_obj@meta.data %>% subset(sexXnutr == 'M_Fast' & cell_type3 == 'Agrp') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(sexXnutr == 'M_Fed' & cell_type3 == 'Agrp') %>% rownames

head(group1)
head(group2)

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='ARH'
)

head(DMEs)

pdf("Agrp_All_DMEsLollipop_Male_Fast_Fed.pdf", width = 6, height = 7)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name='ARH', 
  pvalue = "p_val_adj"
)
dev.off()

png("Agrp_All_DMEsLollipop_Male_Fast_Fed.png", width = 1200, height = 1400,res = 300)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name='ARH', 
  pvalue = "p_val_adj",
)
dev.off()



# Pathway enrichment analysis

library(enrichR)

dbs <-c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021', 'WikiPathway_2021_Mouse', 'KEGG_2021_Mouse')

# compute GO terms:
seurat_obj <- RunEnrichr(seurat_obj, dbs=dbs)

enrichr_df <- GetEnrichrTable(seurat_obj) %>% subset(P.value < 0.05)
write.table(enrichr_df, quote=FALSE, sep='\t', row.names=FALSE, file = 'mouse_Agrp_enrichr.tsv')

# make GO term plots:
EnrichrBarPlot(
  seurat_obj,
  outdir = 
  n_terms = 10, 
  plot_size = c(5,7), 
  logscale=TRUE 
)


# GO_BP
pdf("Agrp_top3_GO_BP_enrichr_dotplot.pdf", width = 7, height = 8)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules
  database = "GO_Biological_Process_2021", 
  n_terms=3, # number of terms per module
  term_size=7, 
  p_adj = FALSE 
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

dev.off()


png("Agrp_top3_GO_BP_enrichr_dotplot.png", width = 2000, height = 2400,res = 300)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules 
  database = "GO_Biological_Process_2021", # 
  n_terms=3, # number of terms per module
  term_size=7, 
  p_adj = FALSE 
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

dev.off()

##GO_MF
pdf("Agrp_top3_GO_MF_enrichr_dotplot.pdf", width = 7, height = 8)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (default)
  database = "GO_Molecular_Function_2021", 
  n_terms=3, # number of terms per module
  term_size=7,
  p_adj = FALSE 
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

dev.off()


png("Agrp_top3_GO_MF_enrichr_dotplot.png", width = 2000, height = 2400,res = 300)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (default)
  database = "GO_Molecular_Function_2021",
  n_terms=3, # number of terms per module
  term_size=7, 
  p_adj = FALSE 
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

dev.off()

##KEGG
pdf("Agrp_top1_KEGG_enrichr_dotplot.pdf", width = 6, height = 5)
EnrichrDotPlot(
  seurat_obj,
  mods = 'all', # use all modules
  database = "KEGG_2021_Mouse", 
  n_terms=1, # number of terms per module
  term_size=7, 
  p_adj = FALSE, 
  break_ties=TRUE
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

dev.off()

png("Agrp_top1_KEGG_enrichr_dotplot.png", width = 2400, height = 2000,res = 300)
EnrichrDotPlot(
  seurat_obj,
  mods = 'all', # use all modules 
  database = "KEGG_2021_Mouse", 
  n_terms=1, # number of terms per module
  term_size=7, 
  p_adj = FALSE, 
  break_ties=TRUE
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

dev.off()

##GO_CC
pdf("Agrp_top1_GO_CC_enrichr_dotplot.pdf", width = 6, height = 5)
EnrichrDotPlot(
  seurat_obj,
  mods = 'all', # use all modules
  database = "GO_Cellular_Component_2021",
  n_terms=1, # number of terms per module
  term_size=7, 
  p_adj = FALSE, 
  break_ties=TRUE
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

dev.off()

png("Agrp_top1_GO_CC_dotplot.png", width = 2400, height = 2000,res = 300)
EnrichrDotPlot(
  seurat_obj,
  mods = 'all', # use all modules
  database = "GO_Cellular_Component_2021", 
  n_terms=1, # number of terms per module
  term_size=7,
  p_adj = FALSE,
  break_ties=TRUE
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

dev.off()

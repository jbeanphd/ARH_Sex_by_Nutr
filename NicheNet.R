##NicheNet analysis

library(nichenetr)
library(Seurat)
library(SeuratObject)
library(tidyverse)

#Agrp neurons

ligand_activities <- predict_ligand_activities(
  geneset = target_genes,
  background_expressed_genes = expressed_genes_receiver,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = ligand
)
ligand_activities

ligand_activities_long <- ligand_activities %>%
  select(test_ligand, pearson, aupr_corrected) %>%
  pivot_longer(cols = c(pearson, aupr_corrected),
               names_to = "score_type", values_to = "score")

##plot bar plot of ligends person score
score <- ggplot(ligand_activities_long, aes(y = reorder(test_ligand, score), x = score, fill = score_type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.7) +
  coord_flip() +
  labs(title = 'Agrp neurons',
       y = "Potential ligand",
       x = "Activity score") +
  scale_fill_manual(values = c("pearson" = "lightblue", "aupr_corrected" = "salmon"),
                    name = "Score type",
                    labels = c("Pearson", "AUPR (corrected)")) +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.spacing.y = unit(0.1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.ticks.length = unit(0.15, "cm"),                            
    axis.line = element_line(color = "black", size = 0.6)             
  ) +
  xlim(0, max(ligand_activities_long$score) + 0.05)
ggsave("Agrp_ligand_activated_score_barplot.pdf", score, width = 7, height = 4,dpi = 300)

best_upstream_ligands <- ligand_activities %>% top_n(10, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()

ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

max_val <- max(vis_ligand_aupr, na.rm = TRUE)

p_ligand_aupr <- make_heatmap_ggplot(
  vis_ligand_aupr,
  "Potential ligands", "Ligand activity", 
  legend_title = "AUPR", 
  color = "darkorange"
) + 
  scale_fill_gradient(low = "white", high = "darkorange", limits = c(0, max_val)) +
  theme(
    axis.text.x.top = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )
p_ligand_aupr
ggsave("Agrp_ligand_activated_score_barplot2.pdf", p_ligand_aupr, width = 2.5, height = 3,dpi = 300)

##potential receptors
receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors <- intersect(receptors,expressed_genes_receiver)

ligands <- lr_network %>% pull(from) %>% unique()
expressed_ligands <- intersect(ligands,expressed_genes_sender)

potential_ligands <-  lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
  pull(from) %>% unique()
ligand <-c('Nrg1','Nrg3','Efna5','Lrrc4c')
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  ligand, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

custom_ligand_order <- c( "Nrg1","Lrrc4c","Efna5","Nrg3" ) 
vis_ligand_receptor_network_KNDy <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  ligand,
  order_hclust = "both")
vis_ligand_receptor_network_ordered<- vis_ligand_receptor_network[, custom_ligand_order]

p <- make_heatmap_ggplot(
  t(vis_ligand_receptor_network_ordered),
  y_name = "Potential ligands",
  x_name = "Receptors expressed by Agrp neurons",
  color = "mediumvioletred",
  legend_title = "Prior interaction potential"
) +
  theme(
    axis.ticks = element_blank()
  )
p


ggsave("Agrp_ligand_target_receptors.pdf", p, width = 9, height = 4,dpi = 300)
ggsave("Agrp_ligand_target_receptors.png", p, width = 9, height = 4,dpi = 300)

sender_celltype <- c('Htr3b', 'Klhl1/Ebf3', 'Sst/Unc13c', 'Agrp', 'Tbx15', 'Oligodendrocytes')

# Dotplot of sender-focused approach(Agrp is receiver)
p_dotplot <- DotPlot(
  subset(female, cell_type3 %in% sender_celltype & sexXnutr == 'F_Fast'),
  features = rev(best_upstream_ligands),
  cols = "RdYlBu"
) + 
  coord_flip() +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.0, vjust = 0.5))

p_dotplot
ggsave("dotplot_ligand_for_Agrp.pdf", p, width = 4.5, height = 5,dpi = 300)


###########LFC heatmap
celltype_order <- levels(Idents(female)) 

# Use this if cell type labels are the identities of your Seurat object
# if not: indicate the celltype_col properly
DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltype],
  get_lfc_celltype, 
  seurat_obj = female,
  condition_colname = "sexXnutr",
  condition_oi = 'F_Fast',
  condition_reference = 'F_Fed',
  celltype_col = "cell_type3",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands 
) 

DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% 
  column_to_rownames("gene") 

vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])

p_lfc <- make_threecolor_heatmap_ggplot(
  vis_ligand_lfc,
  "Prioritized ligands", "LFC in Sender",
  low_color = "midnightblue", mid_color = "white",
  mid = median(vis_ligand_lfc), high_color = "red",
  legend_title = "LFC"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),axis.ticks = element_blank()
  )
p_lfc

##dotplot for all ligands and sender celltype in ARH
fasting_df <- read.csv("net_F_Fast.csv")
ligands_fasting <- unique(fasting_df$ligand)
##subset female data 
cells_to_remove <- rownames(female@meta.data)[female@meta.data$cell_type2 %in% c('VMH Neurons', 'MM Neurons', 'PVp Neurons', 'SCN Neurons', 'TU Neurons')]
female <- subset(female, cells = setdiff(colnames(female), cells_to_remove))
cells_to_remove <- rownames(female@meta.data)[female@meta.data$cell_type3 %in% c('PVp.02', 'PVp.03', 'PVp.04', 'VMH.01', 'VMH.02','Tu','Klhl1/Ebf3','SCN')]
female <- subset(female, cells = setdiff(colnames(female), cells_to_remove))

p_dotplot <- DotPlot(
  subset(female, sexXnutr == 'F_Fast'),
  features = rev(ligands_fasting),
  cols = "RdYlBu"
) + 
  coord_flip() +
  scale_y_discrete(position = "right") +
  labs(x = "Potential ligands", y = "Cell type") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.0, vjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

p_dotplot

ggsave("Fasting_F_dotplot_all_ligand_for_all_celltype.pdf", p_dotplot, width = 15, height = 5.3,dpi = 300)
ggsave("Fasting_F_dotplot_all_ligand_for_all_celltype.tif", p_dotplot, width = 15, height = 5.3,dpi = 300)


##LFC heatmap for all celltype in ARH
ligands_fasting <- c("Nrg1", "Nrg3", "Fgf1", "Sema3d", "Slit2", "Cadm1", "Efna5", "Lrrc4c",
                     "Nrxn1", "Tenm1", "Tenm3", "Ptn",
                     "Cdh4", "Mag", "Negr1", "Nrxn3", "Sema6a", "Tenm4", "Lrfn5")
Idents(female) <- "cell_type3"
celltype_order <- levels(Idents(female)) 
DE_table_top_ligands <- lapply(
  celltype_order,
  get_lfc_celltype, 
  seurat_obj = female,
  condition_colname = "sexXnutr",
  condition_oi = 'F_Fast',
  condition_reference = 'F_Fed',
  celltype_col = "cell_type3",
  min.pct = 0, logfc.threshold = 0,
  features = ligands_fasting
) 

DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% 
  column_to_rownames("gene") 

vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(ligands_fasting), , drop = FALSE])

p_lfc <- make_threecolor_heatmap_ggplot(
  vis_ligand_lfc,
  "Potential ligands", "LFC in Sender",
  low_color = "midnightblue", mid_color = "white",
  mid = median(vis_ligand_lfc), high_color = "red",
  legend_title = "LFC"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "Black", size = 12, face = "plain"),
    axis.text.y = element_text(color = "black", size = 12),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

p_lfc
ggsave("Fasting_F_LFC_heatmap_all_ligand_for_all_celltype.pdf", p_lfc, width = 9, height = 4.8,dpi = 300)
ggsave("Fasting_F_LFC_heatmap_all_ligand_for_all_celltype.tif", p_lfc, width = 9, height = 4.8,dpi = 300)


##reciptor
##LFC heatmap for all celltype in ARH
Idents(female) <- "cell_type3"
celltype_order <- levels(Idents(female)) 
receptors_fasting <- c("Erbb4", "Fgfr2", "Robo1", "Robo2", "Cadm1", "Epha5", "Ntng1", "Nlgn1",
                       "Adgrl2", "Adgrl3", "Lrrtm4", "Sdc2", "Cdh4", "Mag", "Negr1", "Plxna4", "Ptprd")
DE_table_top_receptors <- lapply(
  celltype_order,
  get_lfc_celltype, 
  seurat_obj = female,
  condition_colname = "sexXnutr",
  condition_oi = 'F_Fast',
  condition_reference = 'F_Fed',
  celltype_col = "cell_type3",
  min.pct = 0, logfc.threshold = 0,
  features = receptors_fasting
) 

DE_table_top_receptors<- DE_table_top_receptors %>%  reduce(., full_join) %>% 
  column_to_rownames("gene") 

vis_receptor_lfc <- as.matrix(DE_table_top_receptors[rev(receptors_fasting), , drop = FALSE])

p_lfc <- make_threecolor_heatmap_ggplot(
  vis_receptor_lfc,
  "Potential Recivers", "LFC in Reciver",
  low_color = "midnightblue", mid_color = "white",
  mid = median(vis_ligand_lfc), high_color = "red",
  legend_title = "LFC"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "Black", size = 12, face = "plain"),
    axis.text.y = element_text(color = "black", size = 12),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

p_lfc
ggsave("Fasting_F_LFC_heatmap_all_receptor_for_all_celltype.pdf", p_lfc, width = 9, height = 5.8,dpi = 300)
ggsave("Fasting_F_LFC_heatmap_all_receptor_for_all_celltype.tif", p_lfc, width = 9, height = 5.8,dpi = 300)

##KNDy neurons
##recepotr
expressed_genes_receiver_KNDy <- get_expressed_genes('KNDy', female, pct = 0.05)
receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors_KNDy <- intersect(receptors,expressed_genes_receiver_KNDy)
ligand_F_KNDy <-c('Nrg1','Nrg3','Slit2')
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  ligand_F_KNDy, expressed_receptors_KNDy,
  lr_network, weighted_networks$lr_sig) 

custom_ligand_order_KNDy <- c( 'Slit2','Nrg1','Nrg3') 
vis_ligand_receptor_network_KNDy <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  ligand,
  order_hclust = "both")
vis_ligand_receptor_network_ordered_KNDy<- vis_ligand_receptor_network_KNDy[, custom_ligand_order_KNDy]

p <- make_heatmap_ggplot(
  t(vis_ligand_receptor_network_ordered_KNDy),
  y_name = "Potential ligands",
  x_name = "Receptors expressed by KNDy neurons",
  color = "mediumvioletred",
  legend_title = "Prior interaction potential"
) +
  theme(
    axis.ticks = element_blank()
  )

ggsave("KNDy_ligand_target_receptors.pdf", p, width = 4, height = 4,dpi = 300)
ggsave("KNDy_ligand_target_receptors.png", p, width = 4, height = 4,dpi = 300)

##ligand_activates for KNDy
KNDy_cells <- WhichCells(seurat_obj, expression = cell_type3 == 'KNDy' & sexXnutr %in% c(group_fasting, group_fed))
KNDy_subset <- subset(seurat_obj, cells = KNDy_cells)
Idents(KNDy_subset) <- "sexXnutr"
deg_KNDy <- FindMarkers(KNDy_subset, ident.1 = group_fasting, ident.2 = group_fed, logfc.threshold = 0.10)
write.csv(deg_KNDy,'KNDy_DEG.csv')
target_genes <- rownames(deg_KNDy)[deg_KNDy$p_val_adj < 0.05]

ligand_activities <- predict_ligand_activities(
  geneset = target_genes,
  background_expressed_genes = expressed_genes_receiver_KNDy,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = ligand_F_KNDy
)
ligand_activities

ligand_activities_long <- ligand_activities %>%
  select(test_ligand, pearson, aupr_corrected) %>%
  pivot_longer(cols = c(pearson, aupr_corrected),
               names_to = "score_type", values_to = "score")

##plot bar plot of ligends person score
score <- ggplot(ligand_activities_long, aes(y = reorder(test_ligand, score), x = score, fill = score_type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.7) +
  coord_flip() +
  labs(title = 'KNDy neurons',
       y = "Potential ligand",
       x = "Activity score") +
  scale_fill_manual(values = c("pearson" = "lightblue", "aupr_corrected" = "salmon"),
                    name = "Score type",
                    labels = c("Pearson", "AUPR (corrected)")) +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.spacing.y = unit(0.1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.ticks.length = unit(0.15, "cm"),                            
    axis.line = element_line(color = "black", size = 0.6)             
  ) +
  xlim(0, max(ligand_activities_long$score) + 0.05)
ggsave("KNDy_ligand_activated_score_barplot.pdf", score, width = 7, height = 4,dpi = 300)
ggsave("KNDy_ligand_activated_score_barplot.png", score, width = 7, height = 4,dpi = 300)

##KNDy_ligand_activated_score_barplot2
best_upstream_ligands <- ligand_activities %>% top_n(10, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()

ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% ligand_F_KNDy) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

max_val <- max(vis_ligand_aupr, na.rm = TRUE)

p_ligand_aupr <- make_heatmap_ggplot(
  vis_ligand_aupr,
  "Potential ligands", "Ligand activity", 
  legend_title = "AUPR", 
  color = "darkorange"
) + 
  scale_fill_gradient(low = "white", high = "darkorange", limits = c(0, max_val)) +
  theme(
    axis.text.x.top = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )
p_ligand_aupr
ggsave("KNDy_ligand_activated_score_barplot2.pdf", p_ligand_aupr, width = 2.5, height = 3,dpi = 300)
ggsave("KNDy_ligand_activated_score_barplot2.png", p_ligand_aupr, width = 2.5, height = 3,dpi = 300)

# Dotplot of sender-focused approach(KNDy is receiver)
sender_celltype_KNDy <- c("Htr3b", "Sst/Unc13c", "Agrp", "Endothelial")
p_dotplot <- DotPlot(
  subset(female, cell_type3 %in% sender_celltype_KNDy & sexXnutr == 'F_Fast'),
  features = rev(best_upstream_ligands),
  cols = "RdYlBu"
) + 
  coord_flip() +
  scale_y_discrete(position = "right") +
  labs(x = "Potential ligands", y = "Cell type") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.0, vjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

p_dotplot
ggsave("dotplot_ligand_for_KNDy.pdf", p_dotplot, width = 4.5, height = 5,dpi = 300)
ggsave("dotplot_ligand_for_KNDy.png", p_dotplot, width = 4.5, height = 5,dpi = 300)

# LFC heatmap
celltype_order <- levels(Idents(female)) 

DE_table_top_ligands_KNDy <- lapply(
  celltype_order[celltype_order %in% sender_celltype_KNDy],
  get_lfc_celltype, 
  seurat_obj = female,
  condition_colname = "sexXnutr",
  condition_oi = 'F_Fast',
  condition_reference = 'F_Fed',
  celltype_col = "cell_type3",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands 
) 

DE_table_top_ligands_KNDy <- DE_table_top_ligands_KNDy %>%  reduce(., full_join) %>% 
  column_to_rownames("gene") 

vis_ligand_lfc <- as.matrix(DE_table_top_ligands_KNDy[rev(best_upstream_ligands), , drop = FALSE])

p_lfc <- make_threecolor_heatmap_ggplot(
  vis_ligand_lfc,
  "Prioritized ligands", "LFC in Sender",
  low_color = "midnightblue", mid_color = "white",
  mid = median(vis_ligand_lfc), high_color = "red",
  legend_title = "LFC"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),axis.ticks = element_blank(),legend.position = "right"
  )
p_lfc
ggsave("Fasting_F_LFC_heatmap_all_receptor_for_KNDy.pdf", p_lfc, width = 3.5, height = 4.5,dpi = 300)
ggsave("Fasting_F_LFC_heatmap_all_receptor_for_KNDy.tif", p_lfc, width = 3.5, height = 4.5,dpi = 300)

##Agrp neurons

ligand_activities <- predict_ligand_activities(
  geneset = target_genes,
  background_expressed_genes = expressed_genes_receiver,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = ligand
)
ligand_activities

ligand_activities_long <- ligand_activities %>%
  select(test_ligand, pearson, aupr_corrected) %>%
  pivot_longer(cols = c(pearson, aupr_corrected),
               names_to = "score_type", values_to = "score")

# plot bar plot of ligends person score
score <- ggplot(ligand_activities_long, aes(y = reorder(test_ligand, score), x = score, fill = score_type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.7) +
  coord_flip() +
  labs(title = 'Agrp neurons',
       y = "Potential ligand",
       x = "Activity score") +
  scale_fill_manual(values = c("pearson" = "lightblue", "aupr_corrected" = "salmon"),
                    name = "Score type",
                    labels = c("Pearson", "AUPR (corrected)")) +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.spacing.y = unit(0.1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.ticks.length = unit(0.15, "cm"),                            
    axis.line = element_line(color = "black", size = 0.6)             
  ) +
  xlim(0, max(ligand_activities_long$score) + 0.05)
ggsave("Agrp_ligand_activated_score_barplot.pdf", score, width = 7, height = 4,dpi = 300)

best_upstream_ligands <- ligand_activities %>% top_n(10, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()

ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

max_val <- max(vis_ligand_aupr, na.rm = TRUE)

p_ligand_aupr <- make_heatmap_ggplot(
  vis_ligand_aupr,
  "Potential ligands", "Ligand activity", 
  legend_title = "AUPR", 
  color = "darkorange"
) + 
  scale_fill_gradient(low = "white", high = "darkorange", limits = c(0, max_val)) +
  theme(
    axis.text.x.top = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )
p_ligand_aupr
ggsave("Agrp_ligand_activated_score_barplot2.pdf", p_ligand_aupr, width = 2.5, height = 3,dpi = 300)

# potential receptors
receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors <- intersect(receptors,expressed_genes_receiver)

ligands <- lr_network %>% pull(from) %>% unique()
expressed_ligands <- intersect(ligands,expressed_genes_sender)

potential_ligands <-  lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
  pull(from) %>% unique()
ligand <-c('Nrg1','Nrg3','Efna5','Lrrc4c')
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  ligand, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

custom_ligand_order <- c( "Nrg1","Lrrc4c","Efna5","Nrg3" ) 
vis_ligand_receptor_network_KNDy <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  ligand,
  order_hclust = "both")
vis_ligand_receptor_network_ordered<- vis_ligand_receptor_network[, custom_ligand_order]

p <- make_heatmap_ggplot(
  t(vis_ligand_receptor_network_ordered),
  y_name = "Potential ligands",
  x_name = "Receptors expressed by Agrp neurons",
  color = "mediumvioletred",
  legend_title = "Prior interaction potential"
) +
  theme(
    axis.ticks = element_blank()
  )
p


ggsave("Agrp_ligand_target_receptors.pdf", p, width = 9, height = 4,dpi = 300)
ggsave("Agrp_ligand_target_receptors.png", p, width = 9, height = 4,dpi = 300)

sender_celltype <- c('Htr3b', 'Klhl1/Ebf3', 'Sst/Unc13c', 'Agrp', 'Tbx15', 'Oligodendrocytes')

# Dotplot of sender-focused approach(Agrp is receiver)
p_dotplot <- DotPlot(
  subset(female, cell_type3 %in% sender_celltype & sexXnutr == 'F_Fast'),
  features = rev(best_upstream_ligands),
  cols = "RdYlBu"
) + 
  coord_flip() +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.0, vjust = 0.5))

p_dotplot
ggsave("dotplot_ligand_for_Agrp.pdf", p, width = 4.5, height = 5,dpi = 300)


# LFC heatmap
celltype_order <- levels(Idents(female)) 

DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltype],
  get_lfc_celltype, 
  seurat_obj = female,
  condition_colname = "sexXnutr",
  condition_oi = 'F_Fast',
  condition_reference = 'F_Fed',
  celltype_col = "cell_type3",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands 
) 

DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% 
  column_to_rownames("gene") 

vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])

p_lfc <- make_threecolor_heatmap_ggplot(
  vis_ligand_lfc,
  "Prioritized ligands", "LFC in Sender",
  low_color = "midnightblue", mid_color = "white",
  mid = median(vis_ligand_lfc), high_color = "red",
  legend_title = "LFC"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),axis.ticks = element_blank()
  )
p_lfc

##dotplot for all ligands and sender celltype in ARH
fasting_df <- read.csv("net_F_Fast.csv")
ligands_fasting <- unique(fasting_df$ligand)

##subset female data 
cells_to_remove <- rownames(female@meta.data)[female@meta.data$cell_type2 %in% c('VMH Neurons', 'MM Neurons', 'PVp Neurons', 'SCN Neurons', 'TU Neurons')]
female <- subset(female, cells = setdiff(colnames(female), cells_to_remove))
cells_to_remove <- rownames(female@meta.data)[female@meta.data$cell_type3 %in% c('PVp.02', 'PVp.03', 'PVp.04', 'VMH.01', 'VMH.02','Tu','Klhl1/Ebf3','SCN')]
female <- subset(female, cells = setdiff(colnames(female), cells_to_remove))

p_dotplot <- DotPlot(
  subset(female, sexXnutr == 'F_Fast'),
  features = rev(ligands_fasting),
  cols = "RdYlBu"
) + 
  coord_flip() +
  scale_y_discrete(position = "right") +
  labs(x = "Potential ligands", y = "Cell type") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.0, vjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

p_dotplot

ggsave("Fasting_F_dotplot_all_ligand_for_all_celltype.pdf", p_dotplot, width = 15, height = 5.3,dpi = 300)
ggsave("Fasting_F_dotplot_all_ligand_for_all_celltype.tif", p_dotplot, width = 15, height = 5.3,dpi = 300)


##LFC heatmap for all celltype in ARH
ligands_fasting <- c("Nrg1", "Nrg3", "Fgf1", "Sema3d", "Slit2", "Cadm1", "Efna5", "Lrrc4c",
                     "Nrxn1", "Tenm1", "Tenm3", "Ptn",
                     "Cdh4", "Mag", "Negr1", "Nrxn3", "Sema6a", "Tenm4", "Lrfn5")
Idents(female) <- "cell_type3"
celltype_order <- levels(Idents(female)) 
DE_table_top_ligands <- lapply(
  celltype_order,
  get_lfc_celltype, 
  seurat_obj = female,
  condition_colname = "sexXnutr",
  condition_oi = 'F_Fast',
  condition_reference = 'F_Fed',
  celltype_col = "cell_type3",
  min.pct = 0, logfc.threshold = 0,
  features = ligands_fasting
) 

DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% 
  column_to_rownames("gene") 

vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(ligands_fasting), , drop = FALSE])

p_lfc <- make_threecolor_heatmap_ggplot(
  vis_ligand_lfc,
  "Potential ligands", "LFC in Sender",
  low_color = "midnightblue", mid_color = "white",
  mid = median(vis_ligand_lfc), high_color = "red",
  legend_title = "LFC"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "Black", size = 12, face = "plain"),
    axis.text.y = element_text(color = "black", size = 12),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

p_lfc
ggsave("Fasting_F_LFC_heatmap_all_ligand_for_all_celltype.pdf", p_lfc, width = 9, height = 4.8,dpi = 300)
ggsave("Fasting_F_LFC_heatmap_all_ligand_for_all_celltype.tif", p_lfc, width = 9, height = 4.8,dpi = 300)


##reciptor
##LFC heatmap for all celltype in ARH
Idents(female) <- "cell_type3"
celltype_order <- levels(Idents(female)) 
receptors_fasting <- c("Erbb4", "Fgfr2", "Robo1", "Robo2", "Cadm1", "Epha5", "Ntng1", "Nlgn1",
                       "Adgrl2", "Adgrl3", "Lrrtm4", "Sdc2", "Cdh4", "Mag", "Negr1", "Plxna4", "Ptprd")
DE_table_top_receptors <- lapply(
  celltype_order,
  get_lfc_celltype, 
  seurat_obj = female,
  condition_colname = "sexXnutr",
  condition_oi = 'F_Fast',
  condition_reference = 'F_Fed',
  celltype_col = "cell_type3",
  min.pct = 0, logfc.threshold = 0,
  features = receptors_fasting
) 

DE_table_top_receptors<- DE_table_top_receptors %>%  reduce(., full_join) %>% 
  column_to_rownames("gene") 

vis_receptor_lfc <- as.matrix(DE_table_top_receptors[rev(receptors_fasting), , drop = FALSE])

p_lfc <- make_threecolor_heatmap_ggplot(
  vis_receptor_lfc,
  "Potential Recivers", "LFC in Reciver",
  low_color = "midnightblue", mid_color = "white",
  mid = median(vis_ligand_lfc), high_color = "red",
  legend_title = "LFC"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "Black", size = 12, face = "plain"),
    axis.text.y = element_text(color = "black", size = 12),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

p_lfc
ggsave("Fasting_F_LFC_heatmap_all_receptor_for_all_celltype.pdf", p_lfc, width = 9, height = 5.8,dpi = 300)
ggsave("Fasting_F_LFC_heatmap_all_receptor_for_all_celltype.tif", p_lfc, width = 9, height = 5.8,dpi = 300)

##KNDy neurons

expressed_genes_receiver_KNDy <- get_expressed_genes('KNDy', female, pct = 0.05)
receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors_KNDy <- intersect(receptors,expressed_genes_receiver_KNDy)
ligand_F_KNDy <-c('Nrg1','Nrg3','Slit2')
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  ligand_F_KNDy, expressed_receptors_KNDy,
  lr_network, weighted_networks$lr_sig) 

custom_ligand_order_KNDy <- c( 'Slit2','Nrg1','Nrg3') 
vis_ligand_receptor_network_KNDy <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  ligand,
  order_hclust = "both")
vis_ligand_receptor_network_ordered_KNDy<- vis_ligand_receptor_network_KNDy[, custom_ligand_order_KNDy]

p <- make_heatmap_ggplot(
  t(vis_ligand_receptor_network_ordered_KNDy),
  y_name = "Potential ligands",
  x_name = "Receptors expressed by KNDy neurons",
  color = "mediumvioletred",
  legend_title = "Prior interaction potential"
) +
  theme(
    axis.ticks = element_blank()
  )
p


ggsave("KNDy_ligand_target_receptors.pdf", p, width = 4, height = 4,dpi = 300)
ggsave("KNDy_ligand_target_receptors.png", p, width = 4, height = 4,dpi = 300)

# ligand_activates for KNDy
KNDy_cells <- WhichCells(seurat_obj, expression = cell_type3 == 'KNDy' & sexXnutr %in% c(group_fasting, group_fed))
KNDy_subset <- subset(seurat_obj, cells = KNDy_cells)
Idents(KNDy_subset) <- "sexXnutr"
deg_KNDy <- FindMarkers(KNDy_subset, ident.1 = group_fasting, ident.2 = group_fed, logfc.threshold = 0.10)
write.csv(deg_KNDy,'KNDy_DEG.csv')
target_genes <- rownames(deg_KNDy)[deg_KNDy$p_val_adj < 0.05]

ligand_activities <- predict_ligand_activities(
  geneset = target_genes,
  background_expressed_genes = expressed_genes_receiver_KNDy,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = ligand_F_KNDy
)
ligand_activities

ligand_activities_long <- ligand_activities %>%
  select(test_ligand, pearson, aupr_corrected) %>%
  pivot_longer(cols = c(pearson, aupr_corrected),
               names_to = "score_type", values_to = "score")

##plot bar plot of ligends person score
score <- ggplot(ligand_activities_long, aes(y = reorder(test_ligand, score), x = score, fill = score_type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.7) +
  coord_flip() +
  labs(title = 'KNDy neurons',
       y = "Potential ligand",
       x = "Activity score") +
  scale_fill_manual(values = c("pearson" = "lightblue", "aupr_corrected" = "salmon"),
                    name = "Score type",
                    labels = c("Pearson", "AUPR (corrected)")) +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.spacing.y = unit(0.1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.ticks.length = unit(0.15, "cm"),                            
    axis.line = element_line(color = "black", size = 0.6)             
  ) +
  xlim(0, max(ligand_activities_long$score) + 0.05)
ggsave("KNDy_ligand_activated_score_barplot.pdf", score, width = 7, height = 4,dpi = 300)
ggsave("KNDy_ligand_activated_score_barplot.png", score, width = 7, height = 4,dpi = 300)

##KNDy_ligand_activated_score_barplot2
best_upstream_ligands <- ligand_activities %>% top_n(10, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()

ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% ligand_F_KNDy) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

max_val <- max(vis_ligand_aupr, na.rm = TRUE)

p_ligand_aupr <- make_heatmap_ggplot(
  vis_ligand_aupr,
  "Potential ligands", "Ligand activity", 
  legend_title = "AUPR", 
  color = "darkorange"
) + 
  scale_fill_gradient(low = "white", high = "darkorange", limits = c(0, max_val)) +
  theme(
    axis.text.x.top = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )
p_ligand_aupr
ggsave("KNDy_ligand_activated_score_barplot2.pdf", p_ligand_aupr, width = 2.5, height = 3,dpi = 300)
ggsave("KNDy_ligand_activated_score_barplot2.png", p_ligand_aupr, width = 2.5, height = 3,dpi = 300)

# Dotplot of sender-focused approach
sender_celltype_KNDy <- c("Htr3b", "Sst/Unc13c", "Agrp", "Endothelial")
p_dotplot <- DotPlot(
  subset(female, cell_type3 %in% sender_celltype_KNDy & sexXnutr == 'F_Fast'),
  features = rev(best_upstream_ligands),
  cols = "RdYlBu"
) + 
  coord_flip() +
  scale_y_discrete(position = "right") +
  labs(x = "Potential ligands", y = "Cell type") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.0, vjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

p_dotplot
ggsave("dotplot_ligand_for_KNDy.pdf", p_dotplot, width = 4.5, height = 5,dpi = 300)
ggsave("dotplot_ligand_for_KNDy.png", p_dotplot, width = 4.5, height = 5,dpi = 300)

# LFC heatmap
celltype_order <- levels(Idents(female)) 

DE_table_top_ligands_KNDy <- lapply(
  celltype_order[celltype_order %in% sender_celltype_KNDy],
  get_lfc_celltype, 
  seurat_obj = female,
  condition_colname = "sexXnutr",
  condition_oi = 'F_Fast',
  condition_reference = 'F_Fed',
  celltype_col = "cell_type3",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands 
) 

DE_table_top_ligands_KNDy <- DE_table_top_ligands_KNDy %>%  reduce(., full_join) %>% 
  column_to_rownames("gene") 

vis_ligand_lfc <- as.matrix(DE_table_top_ligands_KNDy[rev(best_upstream_ligands), , drop = FALSE])

p_lfc <- make_threecolor_heatmap_ggplot(
  vis_ligand_lfc,
  "Prioritized ligands", "LFC in Sender",
  low_color = "midnightblue", mid_color = "white",
  mid = median(vis_ligand_lfc), high_color = "red",
  legend_title = "LFC"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),axis.ticks = element_blank(),legend.position = "right"
  )
p_lfc
ggsave("Fasting_F_LFC_heatmap_all_receptor_for_KNDy.pdf", p_lfc, width = 3.5, height = 4.5,dpi = 300)
ggsave("Fasting_F_LFC_heatmap_all_receptor_for_KNDy.tif", p_lfc, width = 3.5, height = 4.5,dpi = 300)

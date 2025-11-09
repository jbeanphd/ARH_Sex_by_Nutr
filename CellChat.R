# read in required libraries
devtools::install_github("jinworks/CellChat")
libs <- c( 'gplots','stringi','reshape2','cowplot','RColorBrewer',
           'sctransform','stringr','org.Mm.eg.db','AnnotationDbi',
           'IRanges','S4Vectors','Biobase','BiocGenerics','clusterProfiler',
           'biomaRt','Matrix','DESeq2','RcppThread', 'extrafont', 'openxlsx',
           'Seurat','dplyr','tidyr','ggplot2','harmony','ggalluvial',
           'scDblFinder','SoupX','UpSetR','ComplexUpset','CellChat','NeuronChat','NMF','ggalluvial')

lapply(libs, require, character.only = TRUE)

#ARH_CellChat <- readRDS('data3/ARH_CellChat.rds')




# increase future globals
options(future.globals.maxSize = 8000 * 1024^2)


# read in seurat object
ARH_Sex_by_Nutr <- readRDS('data/ARH_Sex_by_Nutr.rds')



####
# subset female fed condition to seurat object
F_Fed_ARH <- subset(ARH_Sex_by_Nutr, subset = sexXnutr == 'F_Fed')

# set database to mouse
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# create cellchat object
F_Fed_ARH_CC <- createCellChat(F_Fed_ARH, group.by = 'cell_type3', assay = 'SCT')

# set database for cellchat object
F_Fed_ARH_CC@DB <- CellChatDB
F_Fed_ARH_CC <- subsetData(F_Fed_ARH_CC) # This step is necessary even if using the whole database

# indentify over epressed genes and interactions
F_Fed_ARH_CC <- identifyOverExpressedGenes(F_Fed_ARH_CC)
F_Fed_ARH_CC <- identifyOverExpressedInteractions(F_Fed_ARH_CC)


# calculate probabilities 
F_Fed_ARH_CC <- computeCommunProb(F_Fed_ARH_CC)

F_Fed_ARH_CC <- computeCommunProbPathway(F_Fed_ARH_CC)
F_Fed_ARH_CC <- aggregateNet(F_Fed_ARH_CC)

# calculate centrality
F_Fed_ARH_CC <- netAnalysis_computeCentrality(F_Fed_ARH_CC)

# plot signaling role scatter plot and network plot
netAnalysis_signalingRole_scatter(F_Fed_ARH_CC)
netAnalysis_signalingRole_network(F_Fed_ARH_CC, width = 10, height = 3, font.size = 8)





####


####
# subset female fasted condition to seurat object
F_Fast_ARH <- subset(ARH_Sex_by_Nutr, subset = sexXnutr == 'F_Fast')


# create cellchat object
F_Fast_ARH_CC <- createCellChat(F_Fast_ARH, group.by = 'cell_type3', assay = 'SCT')

# set database for cellchat object
F_Fast_ARH_CC@DB <- CellChatDB
F_Fast_ARH_CC <- subsetData(F_Fast_ARH_CC) # This step is necessary even if using the whole database

# identify over expressed genes and interactions
F_Fast_ARH_CC <- identifyOverExpressedGenes(F_Fast_ARH_CC)
F_Fast_ARH_CC <- identifyOverExpressedInteractions(F_Fast_ARH_CC)

# calculate probabilities 
F_Fast_ARH_CC <- computeCommunProb(F_Fast_ARH_CC)

F_Fast_ARH_CC <- computeCommunProbPathway(F_Fast_ARH_CC)
F_Fast_ARH_CC <- aggregateNet(F_Fast_ARH_CC)


# calculate centrality 
F_Fast_ARH_CC <- netAnalysis_computeCentrality(F_Fast_ARH_CC)

# plot results in scatter plot and network plot
netAnalysis_signalingRole_scatter(F_Fast_ARH_CC)
netAnalysis_signalingRole_network(F_Fast_ARH_CC, width = 10, height = 3, font.size = 8)



####
####
# subset male fed condition to seurat object
M_Fed_ARH <- subset(ARH_Sex_by_Nutr, subset = sexXnutr == 'M_Fed')


# create cellchat object
M_Fed_ARH_CC <- createCellChat(M_Fed_ARH, group.by = 'cell_type3', assay = 'SCT')

# set database for cellchat object
M_Fed_ARH_CC@DB <- CellChatDB
M_Fed_ARH_CC <- subsetData(M_Fed_ARH_CC) # This step is necessary even if using the whole database

# identify over expressed genes and interactions
M_Fed_ARH_CC <- identifyOverExpressedGenes(M_Fed_ARH_CC)
M_Fed_ARH_CC <- identifyOverExpressedInteractions(M_Fed_ARH_CC)

# calculate probabilities
M_Fed_ARH_CC <- computeCommunProb(M_Fed_ARH_CC)


M_Fed_ARH_CC <- computeCommunProbPathway(M_Fed_ARH_CC)
M_Fed_ARH_CC <- aggregateNet(M_Fed_ARH_CC)

# calculate centrality
M_Fed_ARH_CC <- netAnalysis_computeCentrality(M_Fed_ARH_CC)

# plot results in scatter plot and network plot
netAnalysis_signalingRole_scatter(M_Fed_ARH_CC)
netAnalysis_signalingRole_network(M_Fed_ARH_CC, width = 10, height = 3, font.size = 8)




####
# subset male fasted condition to seurat object
M_Fast_ARH <- subset(ARH_Sex_by_Nutr, subset = sexXnutr == 'M_Fast')


# create cellchat object
M_Fast_ARH_CC <- createCellChat(M_Fast_ARH, group.by = 'cell_type3', assay = 'SCT')

# set daatbase for cellchat object
M_Fast_ARH_CC@DB <- CellChatDB
M_Fast_ARH_CC <- subsetData(M_Fast_ARH_CC) # This step is necessary even if using the whole database

# indentify over expressed genes and interactions
M_Fast_ARH_CC <- identifyOverExpressedGenes(M_Fast_ARH_CC)
M_Fast_ARH_CC <- identifyOverExpressedInteractions(M_Fast_ARH_CC)

# calculate probabilities
M_Fast_ARH_CC <- computeCommunProb(M_Fast_ARH_CC)


M_Fast_ARH_CC <- computeCommunProbPathway(M_Fast_ARH_CC)
M_Fast_ARH_CC <- aggregateNet(M_Fast_ARH_CC)

# calculate centrality 
M_Fast_ARH_CC <- netAnalysis_computeCentrality(M_Fast_ARH_CC)

# plot results in scatter plot and network plot
netAnalysis_signalingRole_scatter(M_Fast_ARH_CC)
netAnalysis_signalingRole_network(M_Fast_ARH_CC, width = 10, height = 3, font.size = 8)



#####
# combine all conditions into one object
ARH.F.M.Fd.Fst.CC <- mergeCellChat(list(F_Fed = F_Fed_ARH_CC, 
                                        F_Fast = F_Fast_ARH_CC, 
                                        M_Fed = M_Fed_ARH_CC, 
                                        M_Fast = M_Fast_ARH_CC), 
                                   add.names = names(list(F_Fed = F_Fed_ARH_CC, 
                                                          F_Fast = F_Fast_ARH_CC, 
                                                          M_Fed = M_Fed_ARH_CC, 
                                                          M_Fast = M_Fast_ARH_CC)))

# combare interactions between conditions
compareInteractions(ARH.F.M.Fd.Fst.CC, show.legend = F, group = c(1,2,3,4), color.use = c('#6a816a','#f9994f','#19552b','#ff6e00'))
ggsave(filename = 'figures/num_interactionsARH_cells_CC.tiff', width = 3, height = 3, dpi = 600)
compareInteractions(ARH.F.M.Fd.Fst.CC, show.legend = F, group = c(1,2,3,4), measure = "weight", color.use = c('#6a816a','#f9994f','#19552b','#ff6e00'))
ggsave(filename = 'figures/weight_interactionsARH_cells_CC.tiff', width = 3, height = 3, dpi = 600)

# plot differential interactions
netVisual_diffInteraction(ARH.F.M.Fd.Fst.CC, weight.scale = T)
netVisual_diffInteraction(ARH.F.M.Fd.Fst.CC, weight.scale = T, measure = "weight")

# plot heatmap
netVisual_heatmap(ARH.F.M.Fd.Fst.CC)
netVisual_heatmap(ARH.F.M.Fd.Fst.CC, measure = 'weight')




# save individual cellchat objects to RDS
saveRDS(F_Fed_ARH_CC, file = 'data/F_Fed_ARH_CC.rds')
saveRDS(F_Fast_ARH_CC, file = 'data/F_Fast_ARH_CC.rds')
saveRDS(M_Fed_ARH_CC, file = 'data/M_Fed_ARH_CC.rds')
saveRDS(M_Fast_ARH_CC, file = 'data/M_Fast_ARH_CC.rds')
saveRDS(ARH.F.M.Fd.Fst.CC, file = 'data/ARH.F.M.Fd.Fst.CC.rds')


# read in cellchat objects
#ARH.F.M.Fd.Fst.CC <- readRDS('data3/ARH.F.M.Fd.Fst.CC.rds')
F_Fed_ARH_CC <- readRDS('data/F_Fed_ARH_CC.rds')
F_Fast_ARH_CC <- readRDS('data/F_Fast_ARH_CC.rds')
M_Fed_ARH_CC <- readRDS('data/M_Fed_ARH_CC.rds')
M_Fast_ARH_CC <- readRDS('data/M_Fast_ARH_CC.rds')

# create object list with all conditions
object.list = list(F_Fed = F_Fed_ARH_CC, 
                   F_Fast = F_Fast_ARH_CC, 
                   M_Fed = M_Fed_ARH_CC, 
                   M_Fast = M_Fast_ARH_CC)



# create a unified pathway object
pathway.union <- union(object.list[[1]]@netP$pathways, c(object.list[[2]]@netP$pathways,object.list[[3]]@netP$pathways, object.list[[4]]@netP$pathways ))

# plot outgoing signaling pathways for female fed
netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "outgoing", signaling = pathway.union, 
                                  title = names(object.list)[1], width = 7, height = 12, font.size = 6, font.size.title = 8)
# plot outgoing signaling pathways for male fed
netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "outgoing", signaling = pathway.union, 
                                  title = names(object.list)[3], width = 7, height = 12, font.size = 6, font.size.title = 8)

# plot outgoing signaling pathways for female fasted
netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing", signaling = pathway.union, 
                                  title = names(object.list)[2], width = 7, height = 12, font.size = 6, font.size.title = 8)
# plot outgoing signaling pathways for male fasted
netAnalysis_signalingRole_heatmap(object.list[[4]], pattern = "outgoing", signaling = pathway.union, 
                                  title = names(object.list)[4], width = 7, height = 12, font.size = 6, font.size.title = 8)


# plot incoming signaling pathways for female fed
netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = pathway.union, 
                                  title = names(object.list)[1], width = 7, height = 12, font.size = 6, font.size.title = 8, color.heatmap = "GnBu")
# plot incoming signaling pathways for male fed
netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "incoming", signaling = pathway.union, 
                                  title = names(object.list)[3], width = 7, height = 12, font.size = 6, font.size.title = 8, color.heatmap = "GnBu")

# plot incoming signaling pathways for female fasted
netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = pathway.union, 
                                  title = names(object.list)[2], width = 7, height = 12, font.size = 6, font.size.title = 8, color.heatmap = "GnBu")
# plot incoming signaling pathways for male fasted
netAnalysis_signalingRole_heatmap(object.list[[4]], pattern = "incoming", signaling = pathway.union, 
                                  title = names(object.list)[4], width = 7, height = 12, font.size = 6, font.size.title = 8, color.heatmap = "GnBu")





# merge female fed and female fasted conditions to one cellchat object
ARH.F.Fd.v.Fst.CC <- mergeCellChat(list(F_Fed = F_Fed_ARH_CC, 
                                        F_Fast = F_Fast_ARH_CC), 
                                   add.names = names(list(F_Fed = F_Fed_ARH_CC, 
                                                          F_Fast = F_Fast_ARH_CC)))
# identify over expressed genes
ARH.F.Fd.v.Fst.CC <- identifyOverExpressedGenes(ARH.F.Fd.v.Fst.CC, group.dataset = 'datasets', 
                                                pos.dataset = 'F_Fast', only.pos = FALSE, thresh.pc = 0.25, thresh.fc = log2(1.25), thresh.p = 0.05)


# save DEG mapping
F_net <- netMappingDEG(ARH.F.Fd.v.Fst.CC,features.name = 'features')
# extract the ligand-receptor pairs with upregulated ligands in LS
net.F_Fast <- subsetCommunication(ARH.F.Fd.v.Fst.CC, net = F_net, datasets = "F_Fast",ligand.logFC = log2(1.25), receptor.logFC = NULL, receptor.pct.1 = 0.001)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.F_Fed <- subsetCommunication(ARH.F.Fd.v.Fst.CC, net = F_net, datasets = "F_Fed",ligand.logFC = -log2(1.25), receptor.logFC = NULL, receptor.pct.1 = 0.001)

gene.F_Fast <- extractGeneSubsetFromPair(net.F_Fast, ARH.F.Fd.v.Fst.CC)
gene.F_Fed <- extractGeneSubsetFromPair(net.F_Fed, ARH.F.Fd.v.Fst.CC)

# Chord diagram for female fed and female fasted for SOURCES 
#par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[1]], slot.name = 'net', 
                     net = net.F_Fed, lab.cex = 1, small.gap = 3.5,
                     sources.use = c('Agrp','KNDy','DA','Pomc','Oligodendrocytes','Microglia'),
                     title.name = 'Fasting Suppressed in Females')


netVisual_chord_gene(object.list[[2]], slot.name = 'net', 
                     net = net.F_Fast, lab.cex = 1, small.gap = 3.5,
                     sources.use = c('Agrp','KNDy','DA','Pomc','Oligodendrocytes','Microglia'),
                     title.name = 'Fasting Induced in Females')




netR.F_Fast <- subsetCommunication(ARH.F.Fd.v.Fst.CC, net = F_net, datasets = "F_Fast",ligand.logFC = NULL, receptor.logFC = log2(1.25), receptor.pct.1 = 25, ligand.pct.1 = 5)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
netR.F_Fed <- subsetCommunication(ARH.F.Fd.v.Fst.CC, net = F_net, datasets = "F_Fed",ligand.logFC = NULL, receptor.logFC = -log2(1.25), receptor.pct.2 = 25, ligand.pct.2 = 5)

# Chord diagram for female fed and female fasted for TARGETS 

netVisual_chord_gene(object.list[[1]], slot.name = 'net', 
                     net = netR.F_Fed, lab.cex = 1, small.gap = 3.5,
                     targets.use = c('Agrp','KNDy','DA','Pomc','Oligodendrocytes','Microglia'),
                     title.name = 'Fasting Suppressed in Females')


netVisual_chord_gene(object.list[[2]], slot.name = 'net', 
                     net = netR.F_Fast, lab.cex = 1, small.gap = 3.5,
                     targets.use = c('Agrp','KNDy','DA','Pomc','Oligodendrocytes','Microglia'),
                     title.name = 'Fasting Induced in Females')



######




# merge male fed and male fasted cellchat objects
ARH.M.Fd.v.Fst.CC <- mergeCellChat(list(M_Fed = M_Fed_ARH_CC, 
                                        M_Fast = M_Fast_ARH_CC), 
                                   add.names = names(list(M_Fed = M_Fed_ARH_CC, 
                                                          M_Fast = M_Fast_ARH_CC)))
# indentify over expressed genes
ARH.M.Fd.v.Fst.CC <- identifyOverExpressedGenes(ARH.M.Fd.v.Fst.CC, group.dataset = 'datasets', 
                                                     pos.dataset = 'M_Fast', only.pos = FALSE, thresh.pc = 0.25, thresh.fc = log2(1.25), thresh.p = 0.05)


# save DEG mapping
M_net <- netMappingDEG(ARH.M.Fd.v.Fst.CC,features.name = 'features')
# extract the ligand-receptor pairs with upregulated ligands in LS
net.M_Fast <- subsetCommunication(ARH.M.Fd.v.Fst.CC, net = M_net, datasets = "M_Fast",ligand.logFC = log2(1.25), receptor.logFC = NULL, receptor.pct.1 = 0.001)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.M_Fed <- subsetCommunication(ARH.M.Fd.v.Fst.CC, net = M_net, datasets = "M_Fed",ligand.logFC = -log2(1.25), receptor.logFC = NULL, receptor.pct.1 = 0.001)

gene.M_Fast <- extractGeneSubsetFromPair(net.M_Fast, ARH.M.Fd.v.Fst.CC)
gene.M_Fed <- extractGeneSubsetFromPair(net.M_Fed, ARH.M.Fd.v.Fst.CC)



# chord plots for male fed and male fasted for SOURCES
netVisual_chord_gene(object.list[[3]], slot.name = 'net', 
                     net = net.M_Fed, lab.cex = 1, small.gap = 3.5,
                     sources.use = c('Agrp','KNDy','DA','Pomc','Oligodendrocytes','Microglia'),
                     title.name = 'Fasting Suppressed in Males')


netVisual_chord_gene(object.list[[4]], slot.name = 'net', 
                     net = net.M_Fast, lab.cex = 1, small.gap = 3.5,
                     sources.use = c('Agrp','KNDy','DA','Pomc','Oligodendrocytes','Microglia'),
                     title.name = 'Fasting Induced in Males')








netR.M_Fast <- subsetCommunication(ARH.M.Fd.v.Fst.CC, net = M_net, datasets = "M_Fast",ligand.logFC = NULL, receptor.logFC = log2(1.25), receptor.pct.1 = 25, ligand.pct.1 = 5)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
netR.M_Fed <- subsetCommunication(ARH.M.Fd.v.Fst.CC, net = M_net, datasets = "M_Fed",ligand.logFC = NULL, receptor.logFC = -log2(1.25), receptor.pct.2 = 25, ligand.pct.2 = 5)

# chord plots for male fed and male fasted for TARGETS

netVisual_chord_gene(object.list[[1]], slot.name = 'net', 
                     net = netR.M_Fed, lab.cex = 1, small.gap = 3.5,
                     targets.use = c('Agrp','KNDy','DA','Pomc','Oligodendrocytes','Microglia'),
                     title.name = 'Fasting Suppressed in Males')


netVisual_chord_gene(object.list[[2]], slot.name = 'net', 
                     net = netR.M_Fast, lab.cex = 1, small.gap = 3.5,
                     targets.use = c('Agrp','KNDy','DA','Pomc','Oligodendrocytes','Microglia'),
                     title.name = 'Fasting Induced in Males')






# save objects to RDS
saveRDS(F_Fed_ARH_CC, file = 'data3/F_Fed_ARH_CC.rds')
saveRDS(F_Fast_ARH_CC, file = 'data3/F_Fast_ARH_CC.rds')
saveRDS(ARH.F.Fd.v.Fst.CC, file = 'data/ARH.F.Fd.v.Fst.CC.rds')

saveRDS(M_Fed_ARH_CC, file = 'data3/M_Fed_ARH_CC.rds')
saveRDS(M_Fast_ARH_CC, file = 'data3/M_Fast_ARH_CC.rds')
saveRDS(ARH.M.Fd.v.Fst.CC, file = 'data/ARH.M.Fd.v.Fst.CC.rds')

saveRDS(ARH.Fd.F.v.M.CC, file = 'data/ARH.Fd.F.v.M.CC.rds')
saveRDS(ARH.Fst.F.v.M.CC, file = 'data/ARH.Fst.F.v.M.CC.rds')


#### 20231117 ####
# read in saved objects
F_Fed_ARH_CC <- readRDS(file = 'data/F_Fed_ARH_CC.rds')
F_Fast_ARH_CC <- readRDS(file = 'data/F_Fast_ARH_CC.rds')
ARH.F.Fd.v.Fst.CC <- readRDS(file = 'data/ARH.F.Fd.v.Fst.CC.rds')

M_Fed_ARH_CC <- readRDS(file = 'data/M_Fed_ARH_CC.rds')
M_Fast_ARH_CC <- readRDS(file = 'data/M_Fast_ARH_CC.rds')
ARH.M.Fd.v.Fst.CC <- readRDS(file = 'data/ARH.M.Fd.v.Fst.CC.rds')

#ARH.Fd.F.v.M.CC <- readRDS('data/ARH.Fd.F.v.M.CC.rds')
#ARH.Fst.F.v.M.CC <- readRDS('data/ARH.Fst.F.v.M.CC.rds')

# save conditions to object list
object.list = list(F_Fed = F_Fed_ARH_CC, 
                   F_Fast = F_Fast_ARH_CC, 
                   M_Fed = M_Fed_ARH_CC, 
                   M_Fast = M_Fast_ARH_CC)


# merge femle fed and male fed into one cellchat object
ARH.Fd.F.v.M.CC <- mergeCellChat(list(F_Fed = F_Fed_ARH_CC, 
                                        M_Fed = M_Fed_ARH_CC), 
                                   add.names = names(list(F_Fed = F_Fed_ARH_CC, 
                                                          M_Fed = M_Fed_ARH_CC)))
# identify over expressed genes
ARH.Fd.F.v.M.CC <- identifyOverExpressedGenes(ARH.Fd.F.v.M.CC, group.dataset = 'datasets', 
                                                pos.dataset = 'M_Fed', only.pos = FALSE, thresh.pc = 0.25, thresh.fc = log2(1.25), thresh.p = 0.05)


# save DEG mapping
Fd_net <- netMappingDEG(ARH.Fd.F.v.M.CC,features.name = 'features')
# extract the ligand-receptor pairs with upregulated ligands in LS
net.F_Fed_sx <- subsetCommunication(ARH.Fd.F.v.M.CC, net = Fd_net, datasets = "F_Fed",ligand.logFC = -log2(1.25), receptor.logFC = NULL, receptor.pct.1 = 0.001)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.M_Fed_sx <- subsetCommunication(ARH.Fd.F.v.M.CC, net = Fd_net, datasets = "M_Fed",ligand.logFC = log2(1.25), receptor.logFC = NULL, receptor.pct.1 = 0.001)

gene.F_Fed_sx <- extractGeneSubsetFromPair(net.F_Fed_sx, ARH.Fd.F.v.M.CC)
gene.M_Fed_sx <- extractGeneSubsetFromPair(net.M_Fed_sx, ARH.Fd.F.v.M.CC)

# chord plot for female fed for SOURCES
netVisual_chord_gene(object.list[[1]], slot.name = 'net',
                     net = net.F_Fed_sx, lab.cex = 1, small.gap = 3.5,
                     sources.use = c('Agrp','KNDy','DA','Pomc','Oligodendrocytes','Microglia'),
                     title.name = 'Greater in Fed Females',
                     link.visible = T)


# chord plot for male fed condition for SOURCES
netVisual_chord_gene(object.list[[3]], slot.name = 'net', 
                     net = net.M_Fed_sx, lab.cex = 1, small.gap = 3.5,
                     sources.use = c('Agrp','KNDy','DA','Pomc','Oligodendrocytes','Microglia'),
                     title.name = 'Greater in Fed Males')




netR.F_Fed_sx <- subsetCommunication(ARH.Fd.F.v.M.CC, net = Fd_net, datasets = "F_Fed",ligand.logFC = NULL, receptor.logFC = -log2(1.25), receptor.pct.1 = 25, ligand.pct.1 = 5)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
netR.M_Fed_sx <- subsetCommunication(ARH.Fd.F.v.M.CC, net = Fd_net, datasets = "M_Fed",ligand.logFC = NULL, receptor.logFC = log2(1.25), receptor.pct.1 = 25, ligand.pct.1 = 5)


# chord plot for female fed condition for TARGETS
netVisual_chord_gene(object.list[[1]], slot.name = 'net',
                     net = netR.F_Fed_sx, lab.cex = 1, small.gap = 3.5, big.gap = 10,
                     targets.use = c('Agrp','KNDy','DA','Pomc','Oligodendrocytes','Microglia'),
                     title.name = 'Greater in Fed Females',
                     link.visible = T)


# chord plot for male fed condition for TARGETS
netVisual_chord_gene(object.list[[3]], slot.name = 'net', 
                     net = netR.M_Fed_sx, lab.cex = 1, small.gap = 3.5,
                     targets.use = c('Agrp','KNDy','DA','Pomc','Oligodendrocytes','Microglia'),
                     title.name = 'Greater in Fed Males')



# merge female fasted and male fasted conditions into on cellchat object
ARH.Fst.F.v.M.CC <- mergeCellChat(list(F_Fast = F_Fast_ARH_CC, 
                                      M_Fast = M_Fast_ARH_CC), 
                                 add.names = names(list(F_Fast = F_Fast_ARH_CC, 
                                                        M_Fast = M_Fast_ARH_CC)))
# identify over expressed genes
ARH.Fst.F.v.M.CC <- identifyOverExpressedGenes(ARH.Fst.F.v.M.CC, group.dataset = 'datasets', 
                                              pos.dataset = 'M_Fast', only.pos = FALSE, thresh.pc = 0.25, thresh.fc = log2(1.25), thresh.p = 0.05)


# save DEG mapping 
Fst_net <- netMappingDEG(ARH.Fst.F.v.M.CC,features.name = 'features')
# extract the ligand-receptor pairs with upregulated ligands in LS
net.F_Fast_sx <- subsetCommunication(ARH.Fst.F.v.M.CC, net = Fst_net, datasets = "F_Fast",ligand.logFC = -log2(1.25), receptor.logFC = NULL, receptor.pct.1 = 0.001)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.M_Fast_sx <- subsetCommunication(ARH.Fst.F.v.M.CC, net = Fst_net, datasets = "M_Fast",ligand.logFC = log2(1.25), receptor.logFC = NULL, receptor.pct.1 = 0.001)

gene.F_Fast_sx <- extractGeneSubsetFromPair(net.F_Fast_sx, ARH.Fst.F.v.M.CC)
gene.M_Fast_sx <- extractGeneSubsetFromPair(net.M_Fast_sx, ARH.Fst.F.v.M.CC)

# chord plots for female fasted and male fasted conditions for SOURCES
netVisual_chord_gene(object.list[[2]], slot.name = 'net',
                     net = net.F_Fast_sx, lab.cex = 1, small.gap = 3.5,
                     sources.use = c('Agrp','KNDy','DA','Pomc','Oligodendrocytes','Microglia'),
                     title.name = 'Greater in Fasted Females')


netVisual_chord_gene(object.list[[4]], slot.name = 'net',
                     net = net.M_Fast_sx, lab.cex = 1, small.gap = 3.5,
                     sources.use = c('Agrp','KNDy','DA','Pomc','Oligodendrocytes','Microglia'),
                     title.name = 'Greater in Fasted Males')




netR.F_Fast_sx <- subsetCommunication(ARH.Fst.F.v.M.CC, net = Fst_net, datasets = "F_Fast",ligand.logFC = NULL, receptor.logFC = -log2(1.25), receptor.pct.1 = 25, ligand.pct.1 = 5)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
netR.M_Fast_sx <- subsetCommunication(ARH.Fst.F.v.M.CC, net = Fst_net, datasets = "M_Fast",ligand.logFC = NULL, receptor.logFC = log2(1.25), receptor.pct.1 = 25, ligand.pct.1 = 5)


# chord plots for female fasted and male fasted for TARGETS
netVisual_chord_gene(object.list[[2]], slot.name = 'net',
                     net = netR.F_Fast_sx, lab.cex = 1, small.gap = 3.5, big.gap = 10,
                     targets.use = c('Agrp','KNDy','DA','Pomc','Oligodendrocytes','Microglia'),
                     title.name = 'Greater in Fasted Females',
                     link.visible = T)



netVisual_chord_gene(object.list[[4]], slot.name = 'net', 
                     net = netR.M_Fast_sx, lab.cex = 1, small.gap = 3.5,
                     targets.use = c('Agrp','KNDy','DA','Pomc','Oligodendrocytes','Microglia'),
                     title.name = 'Greater in Fasted Males')




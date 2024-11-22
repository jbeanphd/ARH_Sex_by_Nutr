libs <- c( 'gplots','stringi','reshape2','cowplot','RColorBrewer',
           'sctransform','stringr','org.Mm.eg.db','AnnotationDbi',
           'IRanges','S4Vectors','Biobase','BiocGenerics','clusterProfiler',
           'biomaRt','Matrix','DESeq2','RcppThread', 'extrafont', 'openxlsx',
           'Seurat','dplyr','tidyr','ggplot2','harmony','ggalluvial',
           'scDblFinder','SoupX','UpSetR','ComplexUpset','CellChat','NeuronChat','NMF','ggalluvial')

lapply(libs, require, character.only = TRUE)

#ARH_CellChat <- readRDS('data3/ARH_CellChat.rds')


#library(NMF)
#library(ggalluvial)


#create basic clustering for sample
ARH_Sex_by_Nutr <- readRDS('data/ARH_Sex_by_Nutr.rds')



####

F_Fed_ARH <- subset(ARH_Sex_by_Nutr, subset = sexXnutr == 'F_Fed')



F_Fed_ARH_CC <- createCellChat(F_Fed_ARH, group.by = 'cell_type3', assay = 'SCT')
#CellChatDB <- CellChatDB.mouse
#showDatabaseCategory(CellChatDB)

F_Fed_ARH_CC@DB <- CellChatDB
F_Fed_ARH_CC <- subsetData(F_Fed_ARH_CC) # This step is necessary even if using the whole database


F_Fed_ARH_CC <- identifyOverExpressedGenes(F_Fed_ARH_CC)
F_Fed_ARH_CC <- identifyOverExpressedInteractions(F_Fed_ARH_CC)


options(future.globals.maxSize = 8000 * 1024^2)
F_Fed_ARH_CC <- computeCommunProb(F_Fed_ARH_CC)


F_Fed_ARH_CC <- computeCommunProbPathway(F_Fed_ARH_CC)
F_Fed_ARH_CC <- aggregateNet(F_Fed_ARH_CC)


F_Fed_ARH_CC <- netAnalysis_computeCentrality(F_Fed_ARH_CC)

netAnalysis_signalingRole_scatter(F_Fed_ARH_CC)
netAnalysis_signalingRole_network(F_Fed_ARH_CC, width = 10, height = 3, font.size = 8)





####


####

F_Fast_ARH <- subset(ARH_Sex_by_Nutr, subset = sexXnutr == 'F_Fast')



F_Fast_ARH_CC <- createCellChat(F_Fast_ARH, group.by = 'cell_type3', assay = 'SCT')


F_Fast_ARH_CC@DB <- CellChatDB
F_Fast_ARH_CC <- subsetData(F_Fast_ARH_CC) # This step is necessary even if using the whole database


F_Fast_ARH_CC <- identifyOverExpressedGenes(F_Fast_ARH_CC)
F_Fast_ARH_CC <- identifyOverExpressedInteractions(F_Fast_ARH_CC)


F_Fast_ARH_CC <- computeCommunProb(F_Fast_ARH_CC)


F_Fast_ARH_CC <- computeCommunProbPathway(F_Fast_ARH_CC)
F_Fast_ARH_CC <- aggregateNet(F_Fast_ARH_CC)



F_Fast_ARH_CC <- netAnalysis_computeCentrality(F_Fast_ARH_CC)

netAnalysis_signalingRole_scatter(F_Fast_ARH_CC)
netAnalysis_signalingRole_network(F_Fast_ARH_CC, width = 10, height = 3, font.size = 8)



####
####

M_Fed_ARH <- subset(ARH_Sex_by_Nutr, subset = sexXnutr == 'M_Fed')



M_Fed_ARH_CC <- createCellChat(M_Fed_ARH, group.by = 'cell_type3', assay = 'SCT')


M_Fed_ARH_CC@DB <- CellChatDB
M_Fed_ARH_CC <- subsetData(M_Fed_ARH_CC) # This step is necessary even if using the whole database


M_Fed_ARH_CC <- identifyOverExpressedGenes(M_Fed_ARH_CC)
M_Fed_ARH_CC <- identifyOverExpressedInteractions(M_Fed_ARH_CC)


M_Fed_ARH_CC <- computeCommunProb(M_Fed_ARH_CC)


M_Fed_ARH_CC <- computeCommunProbPathway(M_Fed_ARH_CC)
M_Fed_ARH_CC <- aggregateNet(M_Fed_ARH_CC)


M_Fed_ARH_CC <- netAnalysis_computeCentrality(M_Fed_ARH_CC)

netAnalysis_signalingRole_scatter(M_Fed_ARH_CC)
netAnalysis_signalingRole_network(M_Fed_ARH_CC, width = 10, height = 3, font.size = 8)




####

M_Fast_ARH <- subset(ARH_Sex_by_Nutr, subset = sexXnutr == 'M_Fast')



M_Fast_ARH_CC <- createCellChat(M_Fast_ARH, group.by = 'cell_type3', assay = 'SCT')


M_Fast_ARH_CC@DB <- CellChatDB
M_Fast_ARH_CC <- subsetData(M_Fast_ARH_CC) # This step is necessary even if using the whole database


M_Fast_ARH_CC <- identifyOverExpressedGenes(M_Fast_ARH_CC)
M_Fast_ARH_CC <- identifyOverExpressedInteractions(M_Fast_ARH_CC)


M_Fast_ARH_CC <- computeCommunProb(M_Fast_ARH_CC)


M_Fast_ARH_CC <- computeCommunProbPathway(M_Fast_ARH_CC)
M_Fast_ARH_CC <- aggregateNet(M_Fast_ARH_CC)


M_Fast_ARH_CC <- netAnalysis_computeCentrality(M_Fast_ARH_CC)

netAnalysis_signalingRole_scatter(M_Fast_ARH_CC)
netAnalysis_signalingRole_network(M_Fast_ARH_CC, width = 10, height = 3, font.size = 8)



#####

ARH.F.M.Fd.Fst.CC <- mergeCellChat(list(F_Fed = F_Fed_ARH_CC, 
                                        F_Fast = F_Fast_ARH_CC, 
                                        M_Fed = M_Fed_ARH_CC, 
                                        M_Fast = M_Fast_ARH_CC), 
                                   add.names = names(list(F_Fed = F_Fed_ARH_CC, 
                                                          F_Fast = F_Fast_ARH_CC, 
                                                          M_Fed = M_Fed_ARH_CC, 
                                                          M_Fast = M_Fast_ARH_CC)))


compareInteractions(ARH.F.M.Fd.Fst.CC, show.legend = F, group = c(1,2,3,4), color.use = c('#6a816a','#f9994f','#19552b','#ff6e00'))
ggsave(filename = 'figures/num_interactionsARH_cells_CC.tiff', width = 3, height = 3, dpi = 600)
compareInteractions(ARH.F.M.Fd.Fst.CC, show.legend = F, group = c(1,2,3,4), measure = "weight", color.use = c('#6a816a','#f9994f','#19552b','#ff6e00'))
ggsave(filename = 'figures/weight_interactionsARH_cells_CC.tiff', width = 3, height = 3, dpi = 600)

netVisual_diffInteraction(ARH.F.M.Fd.Fst.CC, weight.scale = T)
netVisual_diffInteraction(ARH.F.M.Fd.Fst.CC, weight.scale = T, measure = "weight")


netVisual_heatmap(ARH.F.M.Fd.Fst.CC)
netVisual_heatmap(ARH.F.M.Fd.Fst.CC, measure = 'weight')





saveRDS(F_Fed_ARH_CC, file = 'data/F_Fed_ARH_CC.rds')
saveRDS(F_Fast_ARH_CC, file = 'data/F_Fast_ARH_CC.rds')
saveRDS(M_Fed_ARH_CC, file = 'data/M_Fed_ARH_CC.rds')
saveRDS(M_Fast_ARH_CC, file = 'data/M_Fast_ARH_CC.rds')
saveRDS(ARH.F.M.Fd.Fst.CC, file = 'data/ARH.F.M.Fd.Fst.CC.rds')



#ARH.F.M.Fd.Fst.CC <- readRDS('data3/ARH.F.M.Fd.Fst.CC.rds')
F_Fed_ARH_CC <- readRDS('data3/F_Fed_ARH_CC.rds')
F_Fast_ARH_CC <- readRDS('data3/F_Fast_ARH_CC.rds')
M_Fed_ARH_CC <- readRDS('data3/M_Fed_ARH_CC.rds')
M_Fast_ARH_CC <- readRDS('data3/M_Fast_ARH_CC.rds')

object.list = list(F_Fed = F_Fed_ARH_CC, 
                   F_Fast = F_Fast_ARH_CC, 
                   M_Fed = M_Fed_ARH_CC, 
                   M_Fast = M_Fast_ARH_CC)




pathway.union <- union(object.list[[1]]@netP$pathways, c(object.list[[2]]@netP$pathways,object.list[[3]]@netP$pathways, object.list[[4]]@netP$pathways ))
netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "outgoing", signaling = pathway.union, 
                                  title = names(object.list)[1], width = 7, height = 12, font.size = 6, font.size.title = 8)

netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "outgoing", signaling = pathway.union, 
                                  title = names(object.list)[3], width = 7, height = 12, font.size = 6, font.size.title = 8)


netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing", signaling = pathway.union, 
                                  title = names(object.list)[2], width = 7, height = 12, font.size = 6, font.size.title = 8)

netAnalysis_signalingRole_heatmap(object.list[[4]], pattern = "outgoing", signaling = pathway.union, 
                                  title = names(object.list)[4], width = 7, height = 12, font.size = 6, font.size.title = 8)


netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = pathway.union, 
                                  title = names(object.list)[1], width = 7, height = 12, font.size = 6, font.size.title = 8, color.heatmap = "GnBu")

netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "incoming", signaling = pathway.union, 
                                  title = names(object.list)[3], width = 7, height = 12, font.size = 6, font.size.title = 8, color.heatmap = "GnBu")


netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = pathway.union, 
                                  title = names(object.list)[2], width = 7, height = 12, font.size = 6, font.size.title = 8, color.heatmap = "GnBu")

netAnalysis_signalingRole_heatmap(object.list[[4]], pattern = "incoming", signaling = pathway.union, 
                                  title = names(object.list)[4], width = 7, height = 12, font.size = 6, font.size.title = 8, color.heatmap = "GnBu")






ARH.F.Fd.v.Fst.CC <- mergeCellChat(list(F_Fed = F_Fed_ARH_CC, 
                                        F_Fast = F_Fast_ARH_CC), 
                                   add.names = names(list(F_Fed = F_Fed_ARH_CC, 
                                                          F_Fast = F_Fast_ARH_CC)))

ARH.F.Fd.v.Fst.CC <- identifyOverExpressedGenes(ARH.F.Fd.v.Fst.CC, group.dataset = 'datasets', 
                                                pos.dataset = 'F_Fast', only.pos = FALSE, thresh.pc = 0.25, thresh.fc = log2(1.25), thresh.p = 0.05)



F_net <- netMappingDEG(ARH.F.Fd.v.Fst.CC,features.name = 'features')
# extract the ligand-receptor pairs with upregulated ligands in LS
net.F_Fast <- subsetCommunication(ARH.F.Fd.v.Fst.CC, net = F_net, datasets = "F_Fast",ligand.logFC = log2(1.25), receptor.logFC = NULL, receptor.pct.1 = 0.001)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.F_Fed <- subsetCommunication(ARH.F.Fd.v.Fst.CC, net = F_net, datasets = "F_Fed",ligand.logFC = -log2(1.25), receptor.logFC = NULL, receptor.pct.1 = 0.001)

gene.F_Fast <- extractGeneSubsetFromPair(net.F_Fast, ARH.F.Fd.v.Fst.CC)
gene.F_Fed <- extractGeneSubsetFromPair(net.F_Fed, ARH.F.Fd.v.Fst.CC)

# Chord diagram
#par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[1]], slot.name = 'net', 
                     net = net.F_Fed, lab.cex = 1, small.gap = 3.5,
                     sources.use = c('Agrp','KNDy','DA','Ghrh/Chat','Oligodendrocytes','Microglia'),
                     title.name = 'Fasting Suppressed in Females')


netVisual_chord_gene(object.list[[2]], slot.name = 'net', 
                     net = net.F_Fast, lab.cex = 1, small.gap = 3.5,
                     sources.use = c('Agrp','KNDy','DA','Ghrh/Chat','Oligodendrocytes','Microglia'),
                     title.name = 'Fasting Induced in Females')




netR.F_Fast <- subsetCommunication(ARH.F.Fd.v.Fst.CC, net = F_net, datasets = "F_Fast",ligand.logFC = NULL, receptor.logFC = log2(1.25), receptor.pct.1 = 25, ligand.pct.1 = 5)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
netR.F_Fed <- subsetCommunication(ARH.F.Fd.v.Fst.CC, net = F_net, datasets = "F_Fed",ligand.logFC = NULL, receptor.logFC = -log2(1.25), receptor.pct.2 = 25, ligand.pct.2 = 5)


netVisual_chord_gene(object.list[[1]], slot.name = 'net', 
                     net = netR.F_Fed, lab.cex = 1, small.gap = 3.5,
                     targets.use = c('Agrp','KNDy','DA','Ghrh/Chat','Oligodendrocytes','Microglia'),
                     title.name = 'Fasting Suppressed in Females')


netVisual_chord_gene(object.list[[2]], slot.name = 'net', 
                     net = netR.F_Fast, lab.cex = 1, small.gap = 3.5,
                     targets.use = c('Agrp','KNDy','DA','Ghrh/Chat','Oligodendrocytes','Microglia'),
                     title.name = 'Fasting Induced in Females')



######





ARH.M.Fd.v.Fst.CC <- mergeCellChat(list(M_Fed = M_Fed_ARH_CC, 
                                        M_Fast = M_Fast_ARH_CC), 
                                   add.names = names(list(M_Fed = M_Fed_ARH_CC, 
                                                          M_Fast = M_Fast_ARH_CC)))

ARH.M.Fd.v.Fst.CC <- identifyOverExpressedGenes(ARH.M.Fd.v.Fst.CC, group.dataset = 'datasets', 
                                                     pos.dataset = 'M_Fast', only.pos = FALSE, thresh.pc = 0.25, thresh.fc = log2(1.25), thresh.p = 0.05)



M_net <- netMappingDEG(ARH.M.Fd.v.Fst.CC,features.name = 'features')
# extract the ligand-receptor pairs with upregulated ligands in LS
net.M_Fast <- subsetCommunication(ARH.M.Fd.v.Fst.CC, net = M_net, datasets = "M_Fast",ligand.logFC = log2(1.25), receptor.logFC = NULL, receptor.pct.1 = 0.001)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.M_Fed <- subsetCommunication(ARH.M.Fd.v.Fst.CC, net = M_net, datasets = "M_Fed",ligand.logFC = -log2(1.25), receptor.logFC = NULL, receptor.pct.1 = 0.001)

gene.M_Fast <- extractGeneSubsetFromPair(net.M_Fast, ARH.M.Fd.v.Fst.CC)
gene.M_Fed <- extractGeneSubsetFromPair(net.M_Fed, ARH.M.Fd.v.Fst.CC)




netVisual_chord_gene(object.list[[3]], slot.name = 'net', 
                     net = net.M_Fed, lab.cex = 1, small.gap = 3.5,
                     sources.use = c('Agrp','KNDy','DA','Ghrh/Chat','Oligodendrocytes','Microglia'),
                     title.name = 'Fasting Suppressed in Males')


netVisual_chord_gene(object.list[[4]], slot.name = 'net', 
                     net = net.M_Fast, lab.cex = 1, small.gap = 3.5,
                     sources.use = c('Agrp','KNDy','DA','Ghrh/Chat','Oligodendrocytes','Microglia'),
                     title.name = 'Fasting Induced in Males')








netR.M_Fast <- subsetCommunication(ARH.M.Fd.v.Fst.CC, net = M_net, datasets = "M_Fast",ligand.logFC = NULL, receptor.logFC = log2(1.25), receptor.pct.1 = 25, ligand.pct.1 = 5)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
netR.M_Fed <- subsetCommunication(ARH.M.Fd.v.Fst.CC, net = M_net, datasets = "M_Fed",ligand.logFC = NULL, receptor.logFC = -log2(1.25), receptor.pct.2 = 25, ligand.pct.2 = 5)


netVisual_chord_gene(object.list[[1]], slot.name = 'net', 
                     net = netR.M_Fed, lab.cex = 1, small.gap = 3.5,
                     targets.use = c('Agrp','KNDy','DA','Ghrh/Chat','Oligodendrocytes','Microglia'),
                     title.name = 'Fasting Suppressed in Males')


netVisual_chord_gene(object.list[[2]], slot.name = 'net', 
                     net = netR.M_Fast, lab.cex = 1, small.gap = 3.5,
                     targets.use = c('Agrp','KNDy','DA','Ghrh/Chat','Oligodendrocytes','Microglia'),
                     title.name = 'Fasting Induced in Males')







saveRDS(F_Fed_ARH_CC, file = 'data3/F_Fed_ARH_CC.rds')
saveRDS(F_Fast_ARH_CC, file = 'data3/F_Fast_ARH_CC.rds')
saveRDS(ARH.F.Fd.v.Fst.CC, file = 'data/ARH.F.Fd.v.Fst.CC.rds')

saveRDS(M_Fed_ARH_CC, file = 'data3/M_Fed_ARH_CC.rds')
saveRDS(M_Fast_ARH_CC, file = 'data3/M_Fast_ARH_CC.rds')
saveRDS(ARH.M.Fd.v.Fst.CC, file = 'data/ARH.M.Fd.v.Fst.CC.rds')

saveRDS(ARH.Fd.F.v.M.CC, file = 'data/ARH.Fd.F.v.M.CC.rds')
saveRDS(ARH.Fst.F.v.M.CC, file = 'data/ARH.Fst.F.v.M.CC.rds')


#### 20231117 ####
F_Fed_ARH_CC <- readRDS(file = 'data/F_Fed_ARH_CC.rds')
F_Fast_ARH_CC <- readRDS(file = 'data/F_Fast_ARH_CC.rds')
ARH.F.Fd.v.Fst.CC <- readRDS(file = 'data/ARH.F.Fd.v.Fst.CC.rds')

M_Fed_ARH_CC <- readRDS(file = 'data/M_Fed_ARH_CC.rds')
M_Fast_ARH_CC <- readRDS(file = 'data/M_Fast_ARH_CC.rds')
ARH.M.Fd.v.Fst.CC <- readRDS(file = 'data/ARH.M.Fd.v.Fst.CC.rds')

ARH.Fd.F.v.M.CC <- readRDS('data/ARH.Fd.F.v.M.CC.rds')
ARH.Fst.F.v.M.CC <- readRDS('data/ARH.Fst.F.v.M.CC.rds')


object.list = list(F_Fed = F_Fed_ARH_CC, 
                   F_Fast = F_Fast_ARH_CC, 
                   M_Fed = M_Fed_ARH_CC, 
                   M_Fast = M_Fast_ARH_CC)



ARH.Fd.F.v.M.CC <- mergeCellChat(list(F_Fed = F_Fed_ARH_CC, 
                                        M_Fed = M_Fed_ARH_CC), 
                                   add.names = names(list(F_Fed = F_Fed_ARH_CC, 
                                                          M_Fed = M_Fed_ARH_CC)))

ARH.Fd.F.v.M.CC <- identifyOverExpressedGenes(ARH.Fd.F.v.M.CC, group.dataset = 'datasets', 
                                                pos.dataset = 'M_Fed', only.pos = FALSE, thresh.pc = 0.25, thresh.fc = log2(1.25), thresh.p = 0.05)



Fd_net <- netMappingDEG(ARH.Fd.F.v.M.CC,features.name = 'features')
# extract the ligand-receptor pairs with upregulated ligands in LS
net.F_Fed_sx <- subsetCommunication(ARH.Fd.F.v.M.CC, net = Fd_net, datasets = "F_Fed",ligand.logFC = -log2(1.25), receptor.logFC = NULL, receptor.pct.1 = 0.001)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.M_Fed_sx <- subsetCommunication(ARH.Fd.F.v.M.CC, net = Fd_net, datasets = "M_Fed",ligand.logFC = log2(1.25), receptor.logFC = NULL, receptor.pct.1 = 0.001)

gene.F_Fed_sx <- extractGeneSubsetFromPair(net.F_Fed_sx, ARH.Fd.F.v.M.CC)
gene.M_Fed_sx <- extractGeneSubsetFromPair(net.M_Fed_sx, ARH.Fd.F.v.M.CC)


netVisual_chord_gene(object.list[[1]], slot.name = 'net',
                     net = net.F_Fed_sx, lab.cex = 1, small.gap = 3.5,
                     sources.use = c('Agrp','KNDy','DA','Ghrh/Chat','Oligodendrocytes','Microglia'),
                     title.name = 'Greater in Fed Females',
                     link.visible = T)



netVisual_chord_gene(object.list[[3]], slot.name = 'net', 
                     net = net.M_Fed_sx, lab.cex = 1, small.gap = 3.5,
                     sources.use = c('Agrp','KNDy','DA','Ghrh/Chat','Oligodendrocytes','Microglia'),
                     title.name = 'Greater in Fed Males')




netR.F_Fed_sx <- subsetCommunication(ARH.Fd.F.v.M.CC, net = Fd_net, datasets = "F_Fed",ligand.logFC = NULL, receptor.logFC = -log2(1.25), receptor.pct.1 = 25, ligand.pct.1 = 5)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
netR.M_Fed_sx <- subsetCommunication(ARH.Fd.F.v.M.CC, net = Fd_net, datasets = "M_Fed",ligand.logFC = NULL, receptor.logFC = log2(1.25), receptor.pct.1 = 25, ligand.pct.1 = 5)

netVisual_chord_gene(object.list[[1]], slot.name = 'net',
                     net = netR.F_Fed_sx, lab.cex = 1, small.gap = 3.5, big.gap = 10,
                     targets.use = c('Agrp','KNDy','DA','Ghrh/Chat','Oligodendrocytes','Microglia'),
                     title.name = 'Greater in Fed Females',
                     link.visible = T)



netVisual_chord_gene(object.list[[3]], slot.name = 'net', 
                     net = netR.M_Fed_sx, lab.cex = 1, small.gap = 3.5,
                     targets.use = c('Agrp','KNDy','DA','Ghrh/Chat','Oligodendrocytes','Microglia'),
                     title.name = 'Greater in Fed Males')




ARH.Fst.F.v.M.CC <- mergeCellChat(list(F_Fast = F_Fast_ARH_CC, 
                                      M_Fast = M_Fast_ARH_CC), 
                                 add.names = names(list(F_Fast = F_Fast_ARH_CC, 
                                                        M_Fast = M_Fast_ARH_CC)))

ARH.Fst.F.v.M.CC <- identifyOverExpressedGenes(ARH.Fst.F.v.M.CC, group.dataset = 'datasets', 
                                              pos.dataset = 'M_Fast', only.pos = FALSE, thresh.pc = 0.25, thresh.fc = log2(1.25), thresh.p = 0.05)



Fst_net <- netMappingDEG(ARH.Fst.F.v.M.CC,features.name = 'features')
# extract the ligand-receptor pairs with upregulated ligands in LS
net.F_Fast_sx <- subsetCommunication(ARH.Fst.F.v.M.CC, net = Fst_net, datasets = "F_Fast",ligand.logFC = -log2(1.25), receptor.logFC = NULL, receptor.pct.1 = 0.001)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.M_Fast_sx <- subsetCommunication(ARH.Fst.F.v.M.CC, net = Fst_net, datasets = "M_Fast",ligand.logFC = log2(1.25), receptor.logFC = NULL, receptor.pct.1 = 0.001)

gene.F_Fast_sx <- extractGeneSubsetFromPair(net.F_Fast_sx, ARH.Fst.F.v.M.CC)
gene.M_Fast_sx <- extractGeneSubsetFromPair(net.M_Fast_sx, ARH.Fst.F.v.M.CC)


netVisual_chord_gene(object.list[[2]], slot.name = 'net',
                     net = net.F_Fast_sx, lab.cex = 1, small.gap = 3.5,
                     sources.use = c('Agrp','KNDy','DA','Ghrh/Chat','Oligodendrocytes','Microglia'),
                     title.name = 'Greater in Fasted Females')


netVisual_chord_gene(object.list[[4]], slot.name = 'net',
                     net = net.M_Fast_sx, lab.cex = 1, small.gap = 3.5,
                     sources.use = c('Agrp','KNDy','DA','Ghrh/Chat','Oligodendrocytes','Microglia'),
                     title.name = 'Greater in Fasted Males')




netR.F_Fast_sx <- subsetCommunication(ARH.Fst.F.v.M.CC, net = Fst_net, datasets = "F_Fast",ligand.logFC = NULL, receptor.logFC = -log2(1.25), receptor.pct.1 = 25, ligand.pct.1 = 5)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
netR.M_Fast_sx <- subsetCommunication(ARH.Fst.F.v.M.CC, net = Fst_net, datasets = "M_Fast",ligand.logFC = NULL, receptor.logFC = log2(1.25), receptor.pct.1 = 25, ligand.pct.1 = 5)

netVisual_chord_gene(object.list[[2]], slot.name = 'net',
                     net = netR.F_Fast_sx, lab.cex = 1, small.gap = 3.5, big.gap = 10,
                     targets.use = c('Agrp','KNDy','DA','Ghrh/Chat','Oligodendrocytes','Microglia'),
                     title.name = 'Greater in Fasted Females',
                     link.visible = T)



netVisual_chord_gene(object.list[[4]], slot.name = 'net', 
                     net = netR.M_Fast_sx, lab.cex = 1, small.gap = 3.5,
                     targets.use = c('Agrp','KNDy','DA','Ghrh/Chat','Oligodendrocytes','Microglia'),
                     title.name = 'Greater in Fasted Males')








num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

netAnalysis_signalingRole_scatter(object.list[[1]], weight.MinMax = weight.MinMax) & 
  coord_cartesian(xlim = c(0,175), ylim = c(0,175))

netAnalysis_signalingRole_scatter(object.list[[2]], weight.MinMax = weight.MinMax) & 
  coord_cartesian(xlim = c(0,175), ylim = c(0,175))

netAnalysis_signalingRole_scatter(object.list[[3]], weight.MinMax = weight.MinMax) & 
  coord_cartesian(xlim = c(0,175), ylim = c(0,175))

netAnalysis_signalingRole_scatter(object.list[[4]], weight.MinMax = weight.MinMax) & 
  coord_cartesian(xlim = c(0,175), ylim = c(0,175))




tm = Sys.time()

ARH.F.Fd.v.Fst.CC <- computeNetSimilarityPairwise(ARH.F.Fd.v.Fst.CC, type = "functional")
#> Compute signaling network similarity for datasets 1 2
ARH.F.Fd.v.Fst.CC <- netEmbedding(ARH.F.Fd.v.Fst.CC, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
ARH.F.Fd.v.Fst.CC <- netClustering(ARH.F.Fd.v.Fst.CC, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(ARH.F.Fd.v.Fst.CC, type = "functional", label.size = 3, 
                            title = "Female: Fed vs Fasted")
#> 2D visualization of signaling networks from datasets 1 2
ggsave(filename = 'figures/female_fd_v_fst_signaling_umap.tiff', device = 'tiff', units = 'in', width = 7, height = 4, dpi = 300)

rankSimilarity(ARH.F.Fd.v.Fst.CC, type = "functional", font.size = 12)
ggsave(filename = 'figures/female_fd_v_fst_signaling_distance.tiff', device = 'tiff', units = 'in', width = 7, height = 14, dpi = 300)



tm = Sys.time()

ARH.M.Fd.v.Fst.CC <- computeNetSimilarityPairwise(ARH.M.Fd.v.Fst.CC, type = "functional")
#> Compute signaling network similarity for datasets 1 2
ARH.M.Fd.v.Fst.CC <- netEmbedding(ARH.M.Fd.v.Fst.CC, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
ARH.M.Fd.v.Fst.CC <- netClustering(ARH.M.Fd.v.Fst.CC, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(ARH.M.Fd.v.Fst.CC, type = "functional", label.size = 3, 
                            title = "Male: Fed vs Fasted")
#> 2D visualization of signaling networks from datasets 1 2
ggsave(filename = 'figures/male_fd_v_fst_signaling_umap.tiff', device = 'tiff', units = 'in', width = 7, height = 4, dpi = 300)

rankSimilarity(ARH.M.Fd.v.Fst.CC, type = "functional", font.size = 12)
ggsave(filename = 'figures/male_fd_v_fst_signaling_distance.tiff', device = 'tiff', units = 'in', width = 7, height = 14, dpi = 300)




tm = Sys.time()

ARH.Fd.F.v.M.CC <- computeNetSimilarityPairwise(ARH.Fd.F.v.M.CC, type = "functional")
#> Compute signaling network similarity for datasets 1 2
ARH.Fd.F.v.M.CC <- netEmbedding(ARH.Fd.F.v.M.CC, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
ARH.Fd.F.v.M.CC <- netClustering(ARH.Fd.F.v.M.CC, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(ARH.Fd.F.v.M.CC, type = "functional", label.size = 3, 
                            title = "Fed: Female vs Male")
#> 2D visualization of signaling networks from datasets 1 2
ggsave(filename = 'figures/fd_f_v_m_signaling_umap.tiff', device = 'tiff', units = 'in', width = 7, height = 4, dpi = 300)

rankSimilarity(ARH.Fd.F.v.M.CC, type = "functional", font.size = 12)
ggsave(filename = 'figures/fd_f_v_m_signaling_distance.tiff', device = 'tiff', units = 'in', width = 7, height = 14, dpi = 300)



tm = Sys.time()

ARH.Fst.F.v.M.CC <- computeNetSimilarityPairwise(ARH.Fst.F.v.M.CC, type = "functional")
#> Compute signaling network similarity for datasets 1 2
ARH.Fst.F.v.M.CC <- netEmbedding(ARH.Fst.F.v.M.CC, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
ARH.Fst.F.v.M.CC <- netClustering(ARH.Fst.F.v.M.CC, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(ARH.Fst.F.v.M.CC, type = "functional", label.size = 3, 
                            title = "Fasted: Female vs Male")
#> 2D visualization of signaling networks from datasets 1 2
ggsave(filename = 'figures/fst_f_v_m_signaling_umap.tiff', device = 'tiff', units = 'in', width = 7, height = 4, dpi = 300)

rankSimilarity(ARH.Fst.F.v.M.CC, type = "functional", font.size = 12)
ggsave(filename = 'figures/fst_f_v_m_signaling_distance.tiff', device = 'tiff', units = 'in', width = 7, height = 14, dpi = 300)




rankNet(ARH.F.Fd.v.Fst.CC, mode = "comparison", measure = "weight", 
        sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE,
        do.flip = F, font.size = 10, x.rotation = 90,color.use = c('#6a816a','#f9994f'))
ggsave(filename = 'figures/female_fd_v_fst_signaling_flow.tiff', device = 'tiff', units = 'in', width = 14, height = 3, dpi = 300)


rankNet(ARH.M.Fd.v.Fst.CC, mode = "comparison", measure = "weight", 
        sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE,
        do.flip = F, font.size = 10, x.rotation = 90,color.use = c('#19552b','#ff6e00'))
ggsave(filename = 'figures/male_fd_v_fst_signaling_flow.tiff', device = 'tiff', units = 'in', width = 14, height = 3, dpi = 300)




rankNet(ARH.Fd.F.v.M.CC, mode = "comparison", measure = "weight", 
        sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE,
        do.flip = F, font.size = 10, x.rotation = 90,color.use = c('pink','lightblue'))
ggsave(filename = 'figures/fd_f_v_m_signaling_flow.tiff', device = 'tiff', units = 'in', width = 14, height = 3, dpi = 300)



rankNet(ARH.Fst.F.v.M.CC, mode = "comparison", measure = "weight", 
        sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE,
        do.flip = F, font.size = 10, x.rotation = 90,color.use = c('red','blue'))
ggsave(filename = 'figures/fst_f_v_m_signaling_flow.tiff', device = 'tiff', units = 'in', width = 14, height = 3, dpi = 300)



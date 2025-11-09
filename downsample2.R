# read in required libraries 
libs <- c( 'gplots','stringi','reshape2','cowplot','RColorBrewer',
           'sctransform','stringr','org.Mm.eg.db','AnnotationDbi',
           'IRanges','S4Vectors','Biobase','BiocGenerics','clusterProfiler',
           'biomaRt','Matrix','DESeq2','RcppThread', 'extrafont', 'openxlsx',
           'Seurat','dplyr','tidyr','ggplot2','harmony','ggalluvial',
           'scDblFinder','SoupX','UpSetR','ComplexUpset','CellChat','NeuronChat','NMF','ggalluvial')

lapply(libs, require, character.only = TRUE)

# create dataframe counting DE genes for each conditions and cell-type combinations 
# while downsampling to 200 cells 
# 4 conditions X 42 cell-types = 168 rows
4*42

ds_DE_celltype3 <- tibble(cell_type3 = 'a', comparison = 'a', condition = 'a', Fed = 0, Fasted = 0, Female = 0, Male = 0, .rows = 168)
rownum = 0


# iterate one row
      
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "VMH.01"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
# MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "VMH.01", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      #iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'VMH.01'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200     
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'VMH.01', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'VMH.01'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'VMH.01', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      #iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'VMH.01'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'VMH.01', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
#track progress      
      print(rownum/168*100)

  # iterate one row 
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Pomc"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Pomc", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
     # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Pomc'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Pomc', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
     # iterate one row 
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Pomc'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Pomc', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Pomc'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Pomc', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)

      
# iterate one row    
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Astrocytes"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Astrocytes", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Astrocytes'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Astrocytes', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Astrocytes'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Astrocytes', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Astrocytes'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Astrocytes', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
# iterate one row    
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Htr3b"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Htr3b", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row 
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Htr3b'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Htr3b', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      #iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Htr3b'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Htr3b', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Htr3b'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Htr3b', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      #track progress
      print(rownum/168*100)
      
      
      
# iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Gad2/Htr2c"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Gad2/Htr2c", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Gad2/Htr2c'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Gad2/Htr2c', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Gad2/Htr2c'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Gad2/Htr2c', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Gad2/Htr2c'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Gad2/Htr2c', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
# iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "PVp.03"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "PVp.03", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'PVp.03'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'PVp.03', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'PVp.03'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'PVp.03', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'PVp.03'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'PVp.03', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
# iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "α-Tanycytes"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "α-Tanycytes", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      #iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'α-Tanycytes'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'α-Tanycytes', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'α-Tanycytes'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'α-Tanycytes', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'α-Tanycytes'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'α-Tanycytes', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
   
      
      
      
 # iterate one row     
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "MM.01"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "MM.01", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'MM.01'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'MM.01', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'MM.01'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'MM.01', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'MM.01'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'MM.01', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
    
 # iterate one row     
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Agrp"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Agrp", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Agrp'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Agrp', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Agrp'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Agrp', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Agrp'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Agrp', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
    
    # iterate one row  
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "β-Tanycytes"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "β-Tanycytes", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'β-Tanycytes'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'β-Tanycytes', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      #iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'β-Tanycytes'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'β-Tanycytes', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'β-Tanycytes'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'β-Tanycytes', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
    
    # iterate one row  
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Lamp5/Npy5r"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Lamp5/Npy5r", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Lamp5/Npy5r'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Lamp5/Npy5r', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Lamp5/Npy5r'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Lamp5/Npy5r', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Lamp5/Npy5r'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Lamp5/Npy5r', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      #track progress
      print(rownum/168*100)
      
      
      
      
     
# iterate one row      
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "MM.03"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "MM.03", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'MM.03'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'MM.03', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'MM.03'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'MM.03', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'MM.03'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'MM.03', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
 # iterate one row     
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "MM.02"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "MM.02", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'MM.02'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'MM.02', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'MM.02'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'MM.02', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'MM.02'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'MM.02', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
     
# iterate one row      
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Oligodendrocytes"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Oligodendrocytes", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row 
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Oligodendrocytes'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Oligodendrocytes', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Oligodendrocytes'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Oligodendrocytes', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Oligodendrocytes'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Oligodendrocytes', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "PVp.01"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "PVp.01", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'PVp.01'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'PVp.01', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'PVp.01'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'PVp.01', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'PVp.01'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'PVp.01', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "DA"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "DA", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'DA'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'DA', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'DA'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'DA', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'DA'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'DA', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
  
# iterate one row      
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Tbx19"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Tbx19", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Tbx19'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Tbx19', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Tbx19'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Tbx19', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Tbx19'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Tbx19', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
    
    # iterate one row  
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Slc17a6/Alk"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Slc17a6/Alk", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Slc17a6/Alk'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Slc17a6/Alk', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Slc17a6/Alk'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Slc17a6/Alk', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Slc17a6/Alk'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Slc17a6/Alk', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Sst/Unc13c"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Sst/Unc13c", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Sst/Unc13c'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Sst/Unc13c', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Sst/Unc13c'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Sst/Unc13c', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Sst/Unc13c'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Sst/Unc13c', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
    
    # iterate one row  
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "KNDy"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "KNDy", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'KNDy'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'KNDy', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'KNDy'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'KNDy', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'KNDy'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'KNDy', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
       # track progress
      print(rownum/168*100)
     
      
      
      
    
    # iterate one row  
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Tu"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Tu", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Tu'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Tu', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Tu'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Tu', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Tu'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Tu', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
      
     # iterate one row 
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Klhl1/Ebf3"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Klhl1/Ebf3", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Klhl1/Ebf3'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Klhl1/Ebf3', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Klhl1/Ebf3'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Klhl1/Ebf3', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Klhl1/Ebf3'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Klhl1/Ebf3', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      #track progress
      print(rownum/168*100)
   
      
      
      
     
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Ghrh/Chat"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Ghrh/Chat", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Ghrh/Chat'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Ghrh/Chat', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Ghrh/Chat'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Ghrh/Chat', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Ghrh/Chat'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Ghrh/Chat', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Ependymal"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Ependymal", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Ependymal'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Ependymal', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Ependymal'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Ependymal', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Ependymal'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Ependymal', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
      
    
    # iterate one row  
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Microglia"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Microglia", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Microglia'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Microglia', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Microglia'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Microglia', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Microglia'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Microglia', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
     
     # iterate one row 
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Lef1"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Lef1", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Lef1'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Lef1', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Lef1'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Lef1', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Lef1'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Lef1', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
     
      
      
      
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "SCN"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "SCN", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'SCN'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'SCN', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'SCN'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'SCN', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'SCN'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'SCN', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
     
     # iterate one row 
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "VLMC"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "VLMC", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'VLMC'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'VLMC', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'VLMC'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'VLMC', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'VLMC'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'VLMC', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
      
      
      # iterate one row 
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "PVp.02"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "PVp.02", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'PVp.02'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'PVp.02', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'PVp.02'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'PVp.02', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'PVp.02'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'PVp.02', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
   
   # iterate one row   
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "TOP"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "TOP", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'TOP'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'TOP', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'TOP'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'TOP', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      #iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'TOP'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'TOP', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "PVp.04"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "PVp.04", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'PVp.04'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'PVp.04', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'PVp.04'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'PVp.04', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'PVp.04'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'PVp.04', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
    
    # iterate one row  
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Endothelial"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Endothelial", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Endothelial'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Endothelial', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Endothelial'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Endothelial', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Endothelial'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Endothelial', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Ebf3/Htr2c"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Ebf3/Htr2c", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Ebf3/Htr2c'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Ebf3/Htr2c', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Ebf3/Htr2c'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Ebf3/Htr2c', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Ebf3/Htr2c'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Ebf3/Htr2c', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
    
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Satb2/Slc18a2"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Satb2/Slc18a2", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Satb2/Slc18a2'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Satb2/Slc18a2', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Satb2/Slc18a2'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Satb2/Slc18a2', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Satb2/Slc18a2'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Satb2/Slc18a2', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Erg/Lepr"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Erg/Lepr", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Erg/Lepr'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Erg/Lepr', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Erg/Lepr'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Erg/Lepr', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Erg/Lepr'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Erg/Lepr', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Pars Tuberalis"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Pars Tuberalis", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Pars Tuberalis'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Pars Tuberalis', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Pars Tuberalis'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Pars Tuberalis', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Pars Tuberalis'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Pars Tuberalis', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
    
 # iterate one row     
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Ros1/Alk"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Ros1/Alk", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Ros1/Alk'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Ros1/Alk', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Ros1/Alk'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Ros1/Alk', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Ros1/Alk'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Ros1/Alk', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
    
   # iterate one row   
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "VMH.02"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "VMH.02", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'VMH.02'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'VMH.02', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'VMH.02'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'VMH.02', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'VMH.02'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'VMH.02', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
      
     
     # iterate one row 
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "OPC"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "OPC", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'OPC'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'OPC', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'OPC'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'OPC', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'OPC'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'OPC', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
      
    
    # iterate one row  
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Tbx15"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Tbx15", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Tbx15'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Tbx15', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Tbx15'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Tbx15', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Tbx15'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Tbx15', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
      
     
     # iterate one row 
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Tac1/Reln"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Tac1/Reln", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Tac1/Reln'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Tac1/Reln', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Tac1/Reln'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Tac1/Reln', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Tac1/Reln'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Tac1/Reln', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
      
      
      
      
      
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = "Coch/Slc18a2"
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Female'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = 200,
                       subset.ident = "Coch/Slc18a2", min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Coch/Slc18a2'
      ds_DE_celltype3$comparison[rownum] = 'Nutr'
      ds_DE_celltype3$condition[rownum] = 'Male'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Coch/Slc18a2', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Fed[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Fasted[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Coch/Slc18a2'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fed'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', max.cells.per.ident = 200,
                       subset.ident = 'Coch/Slc18a2', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      # iterate one row
      rownum = rownum + 1
      
      
      ds_DE_celltype3$cell_type3[rownum] = 'Coch/Slc18a2'
      ds_DE_celltype3$comparison[rownum] = 'Sex'
      ds_DE_celltype3$condition[rownum] = 'Fasted'
      
      # MAST test with latent variable Sample_ID, min pct 0.25, logfc threshold log2(1.25), max cells 200
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', max.cells.per.ident = 200,
                       subset.ident = 'Coch/Slc18a2', min.pct = 0.25, logfc.threshold = log2(1.25), test.use = 'MAST',
                       latent.vars = 'Sample_ID') |> 
        filter(p_val_adj < 0.05)
      
      ds_DE_celltype3$Female[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      ds_DE_celltype3$Male[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      # track progress
      print(rownum/168*100)
 
 ds_DE_celltype3$Fed <- ds_DE_celltype3$Fed %>% as.numeric()     
 ds_DE_celltype3$Fasted <- ds_DE_celltype3$Fasted %>% as.numeric() 
 ds_DE_celltype3$Female <- ds_DE_celltype3$Female %>% as.numeric() 
 ds_DE_celltype3$Male <- ds_DE_celltype3$Male %>% as.numeric() 
      
 # custom graphic plotting DE genes when downsampled to 200 
 # nutritional regulation for females
      ds_DE_celltype3 %>% filter(condition == 'Female') %>% 
        ggplot(aes(Fed, Fasted)) + geom_point() + 
        geom_label(aes(label = cell_type3), size = 8/.pt) +
        lims(x = c(-5,35), y = c(0,30)) +
        labs(title = 'Females: Nutritionally Regulated') +
        theme_classic() +
        theme(text = element_text(family = "Arial", size = 6, color = 'black'),
              axis.text = element_text(family = 'Arial', size = 6, color = 'black'),
              plot.title = element_text(hjust = 0.5))
      
      ggsave(filename = 'figures/ds_DE_Females.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.8, dpi = 600)
      
    
      # custom graphic plotting DE genes when downsampled to 200 
      # nutritional regulation for males  
      ds_DE_celltype3 %>% filter(condition == 'Male') %>% 
        ggplot(aes(Fed, Fasted)) + geom_point() + 
        geom_label(aes(label = cell_type3), size = 8/.pt) +
        lims(x = c(-2,8.5), y = c(0,35)) +
        labs(title = 'Males: Nutritionally Regulated') +
        theme_classic() +
        theme(text = element_text(family = "Arial", size = 6, color = 'black'),
              axis.text = element_text(family = 'Arial', size = 6, color = 'black'),
              plot.title = element_text(hjust = 0.5))
      
      ggsave(filename = 'figures/ds_DE_Males.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.8, dpi = 600)
      
      
      # custom graphic plotting DE genes when downsampled to 200 
      # sexual regulation for fed
      ds_DE_celltype3 %>% filter(condition == 'Fed') %>% 
        ggplot(aes(Female, Male)) + geom_point() + 
        geom_label(aes(label = cell_type3), size = 8/.pt) +
        lims(x = c(-5,35), y = c(0,20)) +
        labs(title = 'Fed: Sexually Dimorphic') +
        theme_classic() +
        theme(text = element_text(family = "Arial", size = 6, color = 'black'),
              axis.text = element_text(family = 'Arial', size = 6, color = 'black'),
              plot.title = element_text(hjust = 0.5))
      
      ggsave(filename = 'figures/ds_DE_Fed.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.8, dpi = 600)
      
      
      # custom graphic plotting DE genes when downsampled to 200 
      # sexual regulation for fasted
      ds_DE_celltype3 %>% filter(condition == 'Fasted') %>% 
        ggplot(aes(Female, Male)) + geom_point() + 
        geom_label(aes(label = cell_type3), size = 8/.pt) +
        lims(x = c(-5,45), y = c(0, 25)) +
        labs(title = 'Fasted: Sexually Dimorphic') +
        theme_classic() +
        theme(text = element_text(family = "Arial", size = 6, color = 'black'),
              axis.text = element_text(family = 'Arial', size = 6, color = 'black'),
              plot.title = element_text(hjust = 0.5))
      
      ggsave(filename = 'figures/ds_DE_Fasted.tiff', device = 'tiff', units = 'in', width = 3.5, height = 1.8, dpi = 600)
      
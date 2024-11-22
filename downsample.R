########

12*42*3*5

downsampled_DE_celltype3 <- tibble(cell_type3 = 0, comparison = 'a', condition = 'a', downsampled = 0, randomseed = 0, DE = 0, .rows = 7560)
rownum = 0


for(cl in as.list(unique(ARH_Sex_by_Nutr@meta.data$cell_type3))){
  for(i in c(50, 100, 200)){
    for(j in c(43445,746774,411735,275672,957057)){
      
      
      rownum = rownum + 1
      
      
      downsampled_DE_celltype3$cell_type3[rownum] = cl
      downsampled_DE_celltype3$comparison[rownum] = 'Nutr'
      downsampled_DE_celltype3$condition1[rownum] = 'Female'
      downsampled_DE_celltype3$condition2[rownum] = 'Fed'
      downsampled_DE_celltype3$downsampled[rownum] = i
      downsampled_DE_celltype3$randomseed[rownum] = j
      
      t1 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'F_Fast', max.cells.per.ident = i, random.seed = j,
                       subset.ident = cl, pseudocount.use = 0, min.pct = 0.25, logfc.threshold = log2(1.25)) |> 
        filter(p_val_adj < 0.05)
      t1 = t1 |> mutate(gene = row.names(t1))
      
      downsampled_DE_celltype3$DE[rownum] = t1 |> filter(avg_log2FC > 0) |> count() 
      
      
      rownum = rownum + 1
      
      
      downsampled_DE_celltype3$cell_type3[rownum] = cl
      downsampled_DE_celltype3$comparison[rownum] = 'Nutr'
      downsampled_DE_celltype3$condition1[rownum] = 'Female'
      downsampled_DE_celltype3$condition2[rownum] = 'Fasted'
      downsampled_DE_celltype3$downsampled[rownum] = i
      downsampled_DE_celltype3$randomseed[rownum] = j
      
      downsampled_DE_celltype3$DE[rownum] = t1 |> filter(avg_log2FC < 0) |> count() 
      
      
      rownum = rownum + 1
      
      downsampled_DE_celltype3$cell_type3[rownum] = cl
      downsampled_DE_celltype3$comparison[rownum] = 'Nutr'
      downsampled_DE_celltype3$condition1[rownum] = 'Male'
      downsampled_DE_celltype3$condition2[rownum] = 'Fed'
      downsampled_DE_celltype3$downsampled[rownum] = i
      downsampled_DE_celltype3$randomseed[rownum] = j
      
      t2 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'M_Fed', ident.2 = 'M_Fast', max.cells.per.ident = i, random.seed = j,
                       subset.ident = cl, pseudocount.use = 0, min.pct = 0.25, logfc.threshold = log2(1.25)) |> 
        filter(p_val_adj < 0.05)
      t2 = t2 |> mutate(gene = row.names(t2))
      
      downsampled_DE_celltype3$DE[rownum] = t2 |> filter(avg_log2FC > 0) |> count()
      
      
      rownum = rownum + 1
      
      downsampled_DE_celltype3$cell_type3[rownum] = cl
      downsampled_DE_celltype3$comparison[rownum] = 'Nutr'
      downsampled_DE_celltype3$condition1[rownum] = 'Male'
      downsampled_DE_celltype3$condition2[rownum] = 'Fasted'
      downsampled_DE_celltype3$downsampled[rownum] = i
      downsampled_DE_celltype3$randomseed[rownum] = j
      
      downsampled_DE_celltype3$DE[rownum] = t2 |> filter(avg_log2FC < 0) |> count()
      
      rownum = rownum + 1
      
      
      downsampled_DE_celltype3$cell_type3[rownum] = cl
      downsampled_DE_celltype3$comparison[rownum] = 'Nutr'
      downsampled_DE_celltype3$condition1[rownum] = 'F.M.Ovlp'
      downsampled_DE_celltype3$condition2[rownum] = 'Fed'
      downsampled_DE_celltype3$downsampled[rownum] = i
      downsampled_DE_celltype3$randomseed[rownum] = j
      
      downsampled_DE_celltype3$DE[rownum] = inner_join(filter(t1, avg_log2FC > 0),filter(t2, avg_log2FC > 0), by = 'gene') |> count()
      
      
      rownum = rownum + 1
      
      
      downsampled_DE_celltype3$cell_type3[rownum] = cl
      downsampled_DE_celltype3$comparison[rownum] = 'Nutr'
      downsampled_DE_celltype3$condition1[rownum] = 'F.M.Ovlp'
      downsampled_DE_celltype3$condition2[rownum] = 'Fasted'
      downsampled_DE_celltype3$downsampled[rownum] = i
      downsampled_DE_celltype3$randomseed[rownum] = j
      
      downsampled_DE_celltype3$DE[rownum] = inner_join(filter(t1, avg_log2FC < 0),filter(t2, avg_log2FC < 0), by = 'gene') |> count()
      
      rownum = rownum + 1
      
      downsampled_DE_celltype3$cell_type3[rownum] = cl
      downsampled_DE_celltype3$comparison[rownum] = 'Sex'
      downsampled_DE_celltype3$condition1[rownum] = 'Fed'
      downsampled_DE_celltype3$condition2[rownum] = 'Female'
      downsampled_DE_celltype3$downsampled[rownum] = i
      downsampled_DE_celltype3$randomseed[rownum] = j
      
      t3 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fed', ident.2 = 'M_Fed', 
                       max.cells.per.ident = i, random.seed = j,
                       subset.ident = cl, pseudocount.use = 0, min.pct = 0.25, logfc.threshold = log2(1.25)) |> 
        filter(p_val_adj < 0.05)
      
      t3 = t3 |> mutate(gene = row.names(t3))
      
      downsampled_DE_celltype3$DE[rownum] = t3 |> filter(avg_log2FC > 0) |> count()
      
      
      rownum = rownum + 1
      
      downsampled_DE_celltype3$cell_type3[rownum] = cl
      downsampled_DE_celltype3$comparison[rownum] = 'Sex'
      downsampled_DE_celltype3$condition1[rownum] = 'Fed'
      downsampled_DE_celltype3$condition2[rownum] = 'Male'
      downsampled_DE_celltype3$downsampled[rownum] = i
      downsampled_DE_celltype3$randomseed[rownum] = j
      
      downsampled_DE_celltype3$DE[rownum] = t3 |> filter(avg_log2FC < 0) |> count()
      
      rownum = rownum + 1
      
      downsampled_DE_celltype3$cell_type3[rownum] = cl
      downsampled_DE_celltype3$comparison[rownum] = 'Sex'
      downsampled_DE_celltype3$condition1[rownum] = 'Fasted'
      downsampled_DE_celltype3$condition2[rownum] = 'Female'
      downsampled_DE_celltype3$downsampled[rownum] = i
      downsampled_DE_celltype3$randomseed[rownum] = j
      
      t4 = FindMarkers(ARH_Sex_by_Nutr, 
                       group.by = 'sexXnutr', ident.1 = 'F_Fast', ident.2 = 'M_Fast', 
                       max.cells.per.ident = i, random.seed = j,
                       subset.ident = cl, pseudocount.use = 0, min.pct = 0.25, logfc.threshold = log2(1.25)) |> 
        filter(p_val_adj < 0.05)
      
      t4 = t4 |> mutate(gene = row.names(t4))
      
      downsampled_DE_celltype3$DE[rownum] = t4 |> filter(avg_log2FC > 0) |> count()
      
      
      rownum = rownum + 1
      
      downsampled_DE_celltype3$cell_type3[rownum] = cl
      downsampled_DE_celltype3$comparison[rownum] = 'Sex'
      downsampled_DE_celltype3$condition1[rownum] = 'Fasted'
      downsampled_DE_celltype3$condition2[rownum] = 'Male'
      downsampled_DE_celltype3$downsampled[rownum] = i
      downsampled_DE_celltype3$randomseed[rownum] = j
      
      downsampled_DE_celltype3$DE[rownum] = t4 |> filter(avg_log2FC < 0) |> count()
      
      
      rownum = rownum + 1
      
      downsampled_DE_celltype3$cell_type3[rownum] = cl
      downsampled_DE_celltype3$comparison[rownum] = 'Sex'
      downsampled_DE_celltype3$condition1[rownum] = 'Fd.Fst.Ovlp'
      downsampled_DE_celltype3$condition2[rownum] = 'Female'
      downsampled_DE_celltype3$downsampled[rownum] = i
      downsampled_DE_celltype3$randomseed[rownum] = j
      
      
      downsampled_DE_celltype3$DE[rownum] = inner_join(filter(t3, avg_log2FC > 0),filter(t4, avg_log2FC > 0), by = 'gene') |> count()
      
      
      rownum = rownum + 1
      
      downsampled_DE_celltype3$cell_type3[rownum] = cl
      downsampled_DE_celltype3$comparison[rownum] = 'Sex'
      downsampled_DE_celltype3$condition1[rownum] = 'Fd.Fst.Ovlp'
      downsampled_DE_celltype3$condition2[rownum] = 'Male'
      downsampled_DE_celltype3$downsampled[rownum] = i
      downsampled_DE_celltype3$randomseed[rownum] = j
      
      
      downsampled_DE_celltype3$DE[rownum] = inner_join(filter(t3, avg_log2FC < 0),filter(t4, avg_log2FC < 0), by = 'gene') |> count()
      
      
      
      
      print(rownum/7560*100)
      
    }
  }
  
}



downsampled_DE_celltype3$DE <- downsampled_DE_celltype3$DE |> as.numeric()

for(i in seq(from = 1, to = 7560, by = 6)){
  
  downsampled_DE_celltype3$DE2[i] = downsampled_DE_celltype3$DE[i] - downsampled_DE_celltype3$DE[i+4]
  downsampled_DE_celltype3$DE2[i+1] = downsampled_DE_celltype3$DE[i+1] - downsampled_DE_celltype3$DE[i+5]
  downsampled_DE_celltype3$DE2[i+2] = downsampled_DE_celltype3$DE[i+2] - downsampled_DE_celltype3$DE[i+4]
  downsampled_DE_celltype3$DE2[i+3] = downsampled_DE_celltype3$DE[i+3] - downsampled_DE_celltype3$DE[i+5]
  downsampled_DE_celltype3$DE2[i+4] = downsampled_DE_celltype3$DE[i+4]
  downsampled_DE_celltype3$DE2[i+5] = downsampled_DE_celltype3$DE[i+5]
}





#####


for(i in seq(from = 1, to = 7560, by = 6)){
  
  downsampled_DE_celltype3$DE3[i] = downsampled_DE_celltype3$DE[i] + downsampled_DE_celltype3$DE[i+2] - downsampled_DE_celltype3$DE[i+4]
  downsampled_DE_celltype3$DE3[i+1] = downsampled_DE_celltype3$DE[i+1] + downsampled_DE_celltype3$DE[i+3] - downsampled_DE_celltype3$DE[i+5]
  downsampled_DE_celltype3$DE3[i+2] = downsampled_DE_celltype3$DE[i+2] - downsampled_DE_celltype3$DE[i+4]
  downsampled_DE_celltype3$DE3[i+3] = downsampled_DE_celltype3$DE[i+3] - downsampled_DE_celltype3$DE[i+5]
  downsampled_DE_celltype3$DE3[i+4] = downsampled_DE_celltype3$DE[i+2]
  downsampled_DE_celltype3$DE3[i+5] = downsampled_DE_celltype3$DE[i+3]
}


########



# Create an empty tibble
sum_ds200_ct3 <- tibble(cell_type3 = character(),nutrDE = double(),sexDE = double())

# Loop through unique values in downsampled_DE_celltype3$cell_type3
for (i in unique(downsampled_DE_celltype3$cell_type3)) {
  # Filter the data and calculate means
  nutr_mean <- downsampled_DE_celltype3 %>%
    filter(cell_type3 == i, downsampled == 200, comparison == 'Nutr', condition1 == 'Female', condition2 == 'Fed') %>%
    summarize(mean_nutr_DE2 = mean(DE2)) %>%
    pull(mean_nutr_DE2) + downsampled_DE_celltype3 %>%
    filter(cell_type3 == i, downsampled == 200, comparison == 'Nutr', condition1 == 'Female', condition2 == 'Fasted') %>%
    summarize(mean_nutr_DE2 = mean(DE2)) %>%
    pull(mean_nutr_DE2) + downsampled_DE_celltype3 %>%
    filter(cell_type3 == i, downsampled == 200, comparison == 'Nutr', condition1 == 'Male', condition2 == 'Fed') %>%
    summarize(mean_nutr_DE2 = mean(DE2)) %>%
    pull(mean_nutr_DE2) + downsampled_DE_celltype3 %>%
    filter(cell_type3 == i, downsampled == 200, comparison == 'Nutr', condition1 == 'Male', condition2 == 'Fasted') %>%
    summarize(mean_nutr_DE2 = mean(DE2)) %>%
    pull(mean_nutr_DE2) + downsampled_DE_celltype3 %>%
    filter(cell_type3 == i, downsampled == 200, comparison == 'Nutr', condition1 == 'F.M.Ovlp', condition2 == 'Fed') %>%
    summarize(mean_nutr_DE2 = mean(DE2)) %>%
    pull(mean_nutr_DE2) + downsampled_DE_celltype3 %>%
    filter(cell_type3 == i, downsampled == 200, comparison == 'Nutr', condition1 == 'F.M.Ovlp', condition2 == 'Fasted') %>%
    summarize(mean_nutr_DE2 = mean(DE2)) %>%
    pull(mean_nutr_DE2)
  
  sex_mean <- downsampled_DE_celltype3 %>%
    filter(cell_type3 == i, downsampled == 200, comparison == 'Sex', condition1 == 'Fed', condition2 == 'Female') %>%
    summarize(mean_sex_DE2 = mean(DE2)) %>%
    pull(mean_sex_DE2) + downsampled_DE_celltype3 %>%
    filter(cell_type3 == i, downsampled == 200, comparison == 'Sex', condition1 == 'Fed', condition2 == 'Male') %>%
    summarize(mean_sex_DE2 = mean(DE2)) %>%
    pull(mean_sex_DE2) + downsampled_DE_celltype3 %>%
    filter(cell_type3 == i, downsampled == 200, comparison == 'Sex', condition1 == 'Fasted', condition2 == 'Female') %>%
    summarize(mean_sex_DE2 = mean(DE2)) %>%
    pull(mean_sex_DE2) + downsampled_DE_celltype3 %>%
    filter(cell_type3 == i, downsampled == 200, comparison == 'Sex', condition1 == 'Fasted', condition2 == 'Male') %>%
    summarize(mean_sex_DE2 = mean(DE2)) %>%
    pull(mean_sex_DE2) + downsampled_DE_celltype3 %>%
    filter(cell_type3 == i, downsampled == 200, comparison == 'Sex', condition1 == 'Fd.Fst.Ovlp', condition2 == 'Female') %>%
    summarize(mean_sex_DE2 = mean(DE2)) %>%
    pull(mean_sex_DE2) + downsampled_DE_celltype3 %>%
    filter(cell_type3 == i, downsampled == 200, comparison == 'Sex', condition1 == 'Fd.Fst.Ovlp', condition2 == 'Male') %>%
    summarize(mean_sex_DE2 = mean(DE2)) %>%
    pull(mean_sex_DE2)
  
  # Add a row to the summary tibble
  sum_ds200_ct3 <- sum_ds200_ct3 %>%
    add_row(cell_type3 = i, nutrDE = nutr_mean, sexDE = sex_mean)
}



#sum_ds200_ct3$nutrDE[6] <- 200
#sum_ds200_ct3$sexDE[12] <- 200

sum_ds200_ct3 |>
  ggplot(aes(x = nutrDE, y = sexDE)) + 
  geom_point(shape = 21, size = 7, fill = c(
    '#cdcdcd',
    '#b5d2c2',
    '#b5d2c2',
    '#b5d2c2',
    '#c6bcc9',
    '#cdcdcd',
    '#d2c2b5',
    '#c6bcc9',
    
    '#f6ac47',
    '#bac7aa',
    '#84dbac',
    '#c6bcc9',
    '#c6bcc9',
    '#d0a965',
    '#c6bcc9',
    '#5989ca',
    '#b5d2c2',
    '#cdcdcd',
    '#bac7aa',
    '#09f374',
    '#cdcdcd',
    '#cdcdcd',
    '#8fe19a',
    '#b5d2c2',
    '#bac7aa',
    '#b5d2c2',
    '#b5d2c2',
    '#c6bcc9',
    '#cdcdcd',
    '#cdcdcd',
    '#c6bcc9',
    '#c6bcc9',
    '#cdcdcd',
    '#b5d2c2',
    '#bebc91',
    '#a6c186',
    '#cdcdcd',
    '#c6bcc9',
    '#b5d2c2',
    '#cdcdcd',
    '#cdcdcd',
    '#c0cc7d'
    
  )) +
  geom_label(label = sum_ds200_ct3$cell_type3) +
  coord_cartesian(xlim = c(0,200), ylim = c(0,200)) +
  labs(title = 'Number of DE Genes When Downsampled to 200 Nuclei per Condition', x = 'Nutritionally Regulated', y = 'Sexually Regulated') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,family = "Arial", size = 8, color = 'black'),
        plot.subtitle =  element_text(hjust = 0.5,family = "Arial", size = 6, color = 'black'),
        text = element_text(family = "Arial", size = 8, color = 'black'),
        axis.text.x= element_text(family = 'Arial', color = 'black', size = 6),
        axis.text.y= element_text(family = 'Arial', color = 'black', size = 6),
        legend.position = 'none')
ggsave(filename = 'figures/celltype3_number_of_DE_DS200_color_code.svg', device = 'svg', units = 'in', width = 4, height = 3, dpi = 600)







DimPlot(ARH_Sex_by_Nutr, label = F, group.by = 'cell_type3', 
        repel = T, shuffle = F,
        label.size = 12/.pt,
        cols = c(
          
          '#f6ac47',
          '#b5d2c2',
          '#c0cc7d',
          '#5989ca',
          '#cdcdcd',
          '#c6bcc9',
          '#b5d2c2',
          '#bebc91',
          '#c6bcc9',
          '#8fe19a',
          '#b5d2c2',
          '#cdcdcd',
          '#09f374',
          '#84dbac',
          '#b5d2c2',
          '#bac7aa',
          '#c6bcc9',
          '#c6bcc9',
          '#c6bcc9',
          '#d0a965',
          '#b5d2c2',
          '#a6c186',
          '#b5d2c2',
          '#c6bcc9',
          '#cdcdcd',
          '#cdcdcd',
          '#c6bcc9',
          '#cdcdcd',
          '#b5d2c2',
          '#b5d2c2',
          '#cdcdcd',
          '#bac7aa',
          '#cdcdcd',
          '#cdcdcd',
          '#b5d2c2',
          '#cdcdcd',
          '#cdcdcd',
          '#c6bcc9',
          '#c6bcc9',
          '#c6bcc9',
          '#d2c2b5',
          '#bac7aa'  
        )
) &
 
  labs(title = "", x = '', y = '') &
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(family = 'Arial', size = 10)) & NoLegend()


ggsave('figures/cell_type_dimplot_color_coded.tiff', device = 'tiff', units = 'in', width = 7,height = 7,dpi = 600)
ggsave('figures/cell_type_dimplot_color_coded.svg', device = 'svg', units = 'in', width = 7,height = 7,dpi = 600)

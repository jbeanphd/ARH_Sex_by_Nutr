
Sex_by_Nutr_DE_numsXY <- tibble(cell_type = 'a', comparison = 'a', condition = 'a', 
                                DE = 0, XY = 0, Auto = 0, XY_Perc = 0, Auto_Perc = 0,
                                .rows = 168)
row.number <- 0
{
  {
    #Agrp
   
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Agrp' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Agrp.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Agrp.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Agrp.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Agrp.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Agrp' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Agrp.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Agrp.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Agrp.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Agrp.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Agrp' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Agrp.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Agrp.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Agrp.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05,
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Agrp.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Agrp' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    #Pomc
    
  
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Pomc' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Pomc.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                   (Pomc.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Pomc.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Pomc.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Pomc' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Pomc.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                   (Pomc.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Pomc.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Pomc.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Pomc' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Pomc.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                    (Pomc.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Pomc.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Pomc.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Pomc' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    #KNDy
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'KNDy' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((KNDy.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (KNDy.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((KNDy.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (KNDy.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'KNDy' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((KNDy.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (KNDy.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((KNDy.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (KNDy.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'KNDy' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((KNDy.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (KNDy.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((KNDy.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (KNDy.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'KNDy' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    #DA
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'DA' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((DA.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (DA.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((DA.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (DA.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'DA' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((DA.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (DA.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((DA.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (DA.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'DA' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((DA.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (DA.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((DA.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (DA.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'DA' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #Ghrh/Chat
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Ghrh/Chat' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Ghrh.Chat.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Ghrh.Chat.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Ghrh.Chat.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Ghrh.Chat.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Ghrh/Chat' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Ghrh.Chat.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Ghrh.Chat.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Ghrh.Chat.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Ghrh.Chat.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Ghrh/Chat' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Ghrh.Chat.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Ghrh.Chat.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Ghrh.Chat.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Ghrh.Chat.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Ghrh/Chat' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    
    #Sst/Unc13c
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Sst/Unc13c' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Sst.Unc13c.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Sst.Unc13c.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Sst.Unc13c.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Sst.Unc13c.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Sst/Unc13c' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Sst.Unc13c.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Sst.Unc13c.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Sst.Unc13c.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Sst.Unc13c.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Sst/Unc13c' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Sst.Unc13c.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Sst.Unc13c.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Sst.Unc13c.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Sst.Unc13c.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Sst/Unc13c' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #Lamp5/Npy5r
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Lamp5/Npy5r' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Lamp5.Npy5r.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Lamp5.Npy5r.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Lamp5.Npy5r.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Lamp5.Npy5r.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Lamp5/Npy5r' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Lamp5.Npy5r.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Lamp5.Npy5r.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Lamp5.Npy5r.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Lamp5.Npy5r.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Lamp5/Npy5r' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Lamp5.Npy5r.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Lamp5.Npy5r.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Lamp5.Npy5r.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Lamp5.Npy5r.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Lamp5/Npy5r' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    #Lef1
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Lef1' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Lef1.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Lef1.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Lef1.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Lef1.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Lef1' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Lef1.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Lef1.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Lef1.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Lef1.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Lef1' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Lef1.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Lef1.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Lef1.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Lef1.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Lef1' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    #Htr3b
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Htr3b' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Htr3b.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Htr3b.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Htr3b.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Htr3b.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Htr3b' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Htr3b.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Htr3b.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Htr3b.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Htr3b.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Htr3b' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Htr3b.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Htr3b.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Htr3b.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Htr3b.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Htr3b' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    #Gad2/Htr2c
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Gad2/Htr2c' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Gad2.Htr2c.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Gad2.Htr2c.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Gad2.Htr2c.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Gad2.Htr2c.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Gad2/Htr2c' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Gad2.Htr2c.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Gad2.Htr2c.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Gad2.Htr2c.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Gad2.Htr2c.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Gad2/Htr2c' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Gad2.Htr2c.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Gad2.Htr2c.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Gad2.Htr2c.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Gad2.Htr2c.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Gad2/Htr2c' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #Ros1/Alk
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Ros1/Alk' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Ros1.Alk.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Ros1.Alk.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Ros1.Alk.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Ros1.Alk.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Ros1/Alk' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Ros1.Alk.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Ros1.Alk.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Ros1.Alk.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Ros1.Alk.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Ros1/Alk' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Ros1.Alk.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Ros1.Alk.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Ros1.Alk.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Ros1.Alk.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Ros1/Alk' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    #Satb2/Slc18a2
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Satb2/Slc18a2' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Satb2.Slc18a2.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Satb2.Slc18a2.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Satb2.Slc18a2.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Satb2.Slc18a2.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Satb2/Slc18a2' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Satb2.Slc18a2.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Satb2.Slc18a2.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Satb2.Slc18a2.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Satb2.Slc18a2.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Satb2/Slc18a2' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Satb2.Slc18a2.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Satb2.Slc18a2.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Satb2.Slc18a2.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Satb2.Slc18a2.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Satb2/Slc18a2' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #Coch/Slc18a2
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Coch/Slc18a2' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Coch.Slc18a2.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Coch.Slc18a2.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Coch.Slc18a2.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Coch.Slc18a2.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Coch/Slc18a2' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Coch.Slc18a2.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Coch.Slc18a2.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Coch.Slc18a2.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Coch.Slc18a2.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Coch/Slc18a2' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Coch.Slc18a2.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Coch.Slc18a2.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Coch.Slc18a2.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Coch.Slc18a2.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Coch/Slc18a2' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #Tbx19
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Tbx19' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Tbx19.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Tbx19.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Tbx19.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Tbx19.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Tbx19' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Tbx19.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Tbx19.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Tbx19.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Tbx19.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Tbx19' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Tbx19.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Tbx19.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Tbx19.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Tbx19.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Tbx19' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    #Tbx15
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Tbx15' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Tbx15.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Tbx15.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Tbx15.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Tbx15.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Tbx15' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Tbx15.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Tbx15.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Tbx15.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Tbx15.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Tbx15' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Tbx15.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Tbx15.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Tbx15.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Tbx15.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Tbx15' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #Ebf3/Htr2c
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Ebf3/Htr2c' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Ebf3.Htr2c.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Ebf3.Htr2c.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Ebf3.Htr2c.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Ebf3.Htr2c.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Ebf3/Htr2c' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Ebf3.Htr2c.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Ebf3.Htr2c.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Ebf3.Htr2c.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Ebf3.Htr2c.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Ebf3/Htr2c' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Ebf3.Htr2c.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Ebf3.Htr2c.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Ebf3.Htr2c.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Ebf3.Htr2c.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Ebf3/Htr2c' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #Slc17a6/Alk
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Slc17a6/Alk' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Slc17a6.Alk.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Slc17a6.Alk.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Slc17a6.Alk.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Slc17a6.Alk.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Slc17a6/Alk' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Slc17a6.Alk.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Slc17a6.Alk.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Slc17a6.Alk.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Slc17a6.Alk.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Slc17a6/Alk' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Slc17a6.Alk.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Slc17a6.Alk.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Slc17a6.Alk.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Slc17a6.Alk.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Slc17a6/Alk' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    #Erg/Lepr
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Erg/Lepr' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Erg.Lepr.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Erg.Lepr.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Erg.Lepr.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Erg.Lepr.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Erg/Lepr' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Erg.Lepr.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Erg.Lepr.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Erg.Lepr.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Erg.Lepr.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Erg/Lepr' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Erg.Lepr.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Erg.Lepr.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Erg.Lepr.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Erg.Lepr.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Erg/Lepr' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #Tac1/Reln
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Tac1/Reln' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Tac1.Reln.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Tac1.Reln.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Tac1.Reln.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Tac1.Reln.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Tac1/Reln' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Tac1.Reln.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Tac1.Reln.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Tac1.Reln.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Tac1.Reln.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Tac1/Reln' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Tac1.Reln.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Tac1.Reln.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Tac1.Reln.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Tac1.Reln.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Tac1/Reln' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #Klhl1/Ebf3
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Klhl1/Ebf3' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Klhl1.Ebf3.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Klhl1.Ebf3.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Klhl1.Ebf3.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Klhl1.Ebf3.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Klhl1/Ebf3' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Klhl1.Ebf3.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Klhl1.Ebf3.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Klhl1.Ebf3.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Klhl1.Ebf3.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Klhl1/Ebf3' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Klhl1.Ebf3.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Klhl1.Ebf3.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Klhl1.Ebf3.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Klhl1.Ebf3.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Klhl1/Ebf3' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #SCN
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'SCN' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((SCN.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (SCN.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((SCN.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (SCN.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'SCN' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((SCN.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (SCN.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((SCN.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (SCN.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'SCN' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((SCN.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (SCN.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((SCN.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (SCN.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'SCN' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    #VMH.01
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'VMH.01' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((VMH.01.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (VMH.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((VMH.01.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (VMH.01.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'VMH.01' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((VMH.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (VMH.01.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((VMH.01.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (VMH.01.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'VMH.01' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((VMH.01.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (VMH.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((VMH.01.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (VMH.01.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'VMH.01' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    #VMH.02
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'VMH.02' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((VMH.02.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (VMH.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((VMH.02.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (VMH.02.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'VMH.02' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((VMH.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (VMH.02.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((VMH.02.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (VMH.02.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'VMH.02' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((VMH.02.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (VMH.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((VMH.02.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (VMH.02.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'VMH.02' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #PVp.01
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'PVp.01' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((PVp.01.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (PVp.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((PVp.01.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (PVp.01.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'PVp.01' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((PVp.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (PVp.01.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((PVp.01.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (PVp.01.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'PVp.01' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((PVp.01.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (PVp.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((PVp.01.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (PVp.01.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'PVp.01' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #PVp.02
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'PVp.02' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((PVp.02.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (PVp.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((PVp.02.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (PVp.02.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'PVp.02' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((PVp.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (PVp.02.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((PVp.02.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (PVp.02.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'PVp.02' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((PVp.02.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (PVp.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((PVp.02.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (PVp.02.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'PVp.02' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #PVp.03
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'PVp.03' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((PVp.03.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (PVp.03.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((PVp.03.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (PVp.03.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'PVp.03' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((PVp.03.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (PVp.03.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((PVp.03.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (PVp.03.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'PVp.03' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((PVp.03.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (PVp.03.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((PVp.03.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (PVp.03.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'PVp.03' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #PVp.04
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'PVp.04' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((PVp.04.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (PVp.04.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((PVp.04.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (PVp.04.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'PVp.04' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((PVp.04.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (PVp.04.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((PVp.04.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (PVp.04.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'PVp.04' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((PVp.04.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (PVp.04.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((PVp.04.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (PVp.04.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'PVp.04' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #MM.01
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'MM.01' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((MM.01.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (MM.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((MM.01.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (MM.01.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'MM.01' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((MM.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (MM.01.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((MM.01.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (MM.01.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'MM.01' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((MM.01.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (MM.01.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((MM.01.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (MM.01.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'MM.01' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    #MM.02
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'MM.02' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((MM.02.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (MM.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((MM.02.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (MM.02.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'MM.02' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((MM.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (MM.02.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((MM.02.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (MM.02.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'MM.02' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((MM.02.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (MM.02.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((MM.02.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (MM.02.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'MM.02' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #MM.03
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'MM.03' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((MM.03.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (MM.03.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((MM.03.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (MM.03.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'MM.03' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((MM.03.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (MM.03.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((MM.03.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (MM.03.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'MM.03' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((MM.03.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (MM.03.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((MM.03.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (MM.03.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'MM.03' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    #Tu
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Tu' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Tu.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Tu.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Tu.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Tu.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Tu' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Tu.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Tu.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Tu.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Tu.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Tu' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Tu.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Tu.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Tu.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Tu.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Tu' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    #Pars Tuberalis
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Pars Tuberalis' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((ParsTub.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (ParsTub.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((ParsTub.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (ParsTub.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Pars Tuberalis' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((ParsTub.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (ParsTub.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((ParsTub.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (ParsTub.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Pars Tuberalis' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((ParsTub.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (ParsTub.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((ParsTub.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (ParsTub.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Pars Tuberalis' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    #Astrocytes
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Astrocytes' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Astrocytes.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Astrocytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Astrocytes.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Astrocytes.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Astrocytes' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Astrocytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Astrocytes.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Astrocytes.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Astrocytes.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Astrocytes' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Astrocytes.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Astrocytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Astrocytes.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Astrocytes.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Astrocytes' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]  
    
    
    
    #alpha-Tanycytes α.Tanycytes
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'α-Tanycytes' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((α.Tanycytes.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (α.Tanycytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((α.Tanycytes.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (α.Tanycytes.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'α-Tanycytes' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((α.Tanycytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (α.Tanycytes.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((α.Tanycytes.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (α.Tanycytes.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'α-Tanycytes' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((α.Tanycytes.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (α.Tanycytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((α.Tanycytes.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (α.Tanycytes.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'α-Tanycytes' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]
    #Beta-Tanycytes β.Tanycytes
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'β-Tanycytes' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((β.Tanycytes.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (β.Tanycytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((β.Tanycytes.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (β.Tanycytes.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'β-Tanycytes' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((β.Tanycytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (β.Tanycytes.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((β.Tanycytes.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (β.Tanycytes.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'β-Tanycytes' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((β.Tanycytes.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (β.Tanycytes.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((β.Tanycytes.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (β.Tanycytes.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'β-Tanycytes' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]
    
    #Ependymal
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Ependymal' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Ependymal.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Ependymal.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Ependymal.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Ependymal.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Ependymal' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Ependymal.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Ependymal.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Ependymal.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Ependymal.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Ependymal' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Ependymal.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Ependymal.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Ependymal.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Ependymal.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Ependymal' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]
    
    #oligodendrocytes
    
    
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Oligodendrocytes' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Oligo.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Oligo.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Oligo.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Oligo.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Oligodendrocytes' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Oligo.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Oligo.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Oligo.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Oligo.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Oligodendrocytes' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Oligo.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Oligo.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Oligo.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Oligo.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Oligodendrocytes' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]
    
    
    #TOP
    
    
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'TOP' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((TOP.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (TOP.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((TOP.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (TOP.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'TOP' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((TOP.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (TOP.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((TOP.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (TOP.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'TOP' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((TOP.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (TOP.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((TOP.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (TOP.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'TOP' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]
    
    
    
    #OPC
    
    
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'OPC' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((OPC.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (OPC.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((OPC.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (OPC.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'OPC' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((OPC.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (OPC.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((OPC.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (OPC.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'OPC' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((OPC.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (OPC.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((OPC.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (OPC.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'OPC' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]
    
    #Microglia
    
    
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Microglia' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Microglia.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Microglia.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Microglia.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Microglia.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Microglia' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Microglia.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Microglia.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Microglia.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Microglia.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Microglia' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Microglia.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Microglia.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Microglia.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Microglia.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Microglia' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]
    
    #Endothelial
    
    
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Endothelial' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Endothelial.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (Endothelial.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Endothelial.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (Endothelial.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Endothelial' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((Endothelial.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (Endothelial.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((Endothelial.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (Endothelial.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Endothelial' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((Endothelial.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (Endothelial.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((Endothelial.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (Endothelial.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'Endothelial' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]
    
    #VLMC
    
    
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'VLMC' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fed'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((VLMC.Fd.F.v.M |> filter(p_val_adj < 0.05)),
                                                     (VLMC.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((VLMC.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))),
                                                     (VLMC.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'VLMC' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fasted'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = anti_join((VLMC.Fst.F.v.M |> filter(p_val_adj < 0.05)), 
                                                     (VLMC.Fd.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = anti_join((VLMC.Fst.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                     (VLMC.Fd.F.v.M |> 
                                                        filter(p_val_adj < 0.05, 
                                                               (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'VLMC' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Fd.Fst.Ovlp'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = inner_join((VLMC.Fd.F.v.M |> filter(p_val_adj < 0.05)), 
                                                      (VLMC.Fst.F.v.M |> filter(p_val_adj < 0.05)), by = 'gene') |> count() |> as.numeric()
    Sex_by_Nutr_DE_numsXY$XY[row.number] = inner_join((VLMC.Fd.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), 
                                                      (VLMC.Fst.F.v.M |> 
                                                         filter(p_val_adj < 0.05, 
                                                                (chromosome_name == 'X' | chromosome_name == 'Y'))), by = 'gene') |> count() |> as.numeric()
    
    row.number = row.number + 1
    
    Sex_by_Nutr_DE_numsXY$cell_type[row.number] = 'VLMC' 
    Sex_by_Nutr_DE_numsXY$comparison[row.number] = 'Sex'
    Sex_by_Nutr_DE_numsXY$condition[row.number] = 'Total'
    Sex_by_Nutr_DE_numsXY$DE[row.number] = Sex_by_Nutr_DE_numsXY$DE[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$DE[row.number -2] + Sex_by_Nutr_DE_numsXY$DE[row.number - 1]  
    Sex_by_Nutr_DE_numsXY$XY[row.number] = Sex_by_Nutr_DE_numsXY$XY[row.number - 3] + 
      Sex_by_Nutr_DE_numsXY$XY[row.number -2] + Sex_by_Nutr_DE_numsXY$XY[row.number - 1]
    
    
    
  }
  
  
}

Sex_by_Nutr_DE_numsXY$DE <- Sex_by_Nutr_DE_numsXY$DE |> as.numeric()
Sex_by_Nutr_DE_numsXY$XY <- Sex_by_Nutr_DE_numsXY$XY |> as.numeric()

Sex_by_Nutr_DE_numsXY$Auto <- Sex_by_Nutr_DE_numsXY$DE - Sex_by_Nutr_DE_numsXY$XY

Sex_by_Nutr_DE_numsXY$XY_Perc = Sex_by_Nutr_DE_numsXY$XY / Sex_by_Nutr_DE_numsXY$DE * 100
Sex_by_Nutr_DE_numsXY$Auto_Perc = Sex_by_Nutr_DE_numsXY$Auto / Sex_by_Nutr_DE_numsXY$DE * 100






Sex_by_Nutr_DE_numsXY |> filter(comparison == 'Sex', condition == 'Total') |> 
  ggplot(aes(x = cell_type, y = XY)) + 
  geom_col(aes(y = XY), color = 'black', size = 0.25) + 
  scale_fill_manual(values = c('#C39BD3','black'), name = '') +
  scale_x_discrete(limits = c('KNDy',
                              'Agrp',
                              'Lamp5/Npy5r',
                              'Pomc',
                              'Astrocytes',
                              'β-Tanycytes',
                              'Ghrh/Chat',
                              'Oligodendrocytes',
                              'Sst/Unc13c',
                              'MM.01',
                              'α-Tanycytes',
                              'VMH.01',
                              'Gad2/Htr2c',
                              'DA',
                              'Tbx19',
                              'Htr3b',
                              'Coch/Slc18a2',
                              'Slc17a6/Alk',
                              'MM.02',
                              'PVp.01',
                              'PVp.02',
                              'Erg/Lepr',
                              'Pars Tuberalis',
                              'Lef1',
                              'OPC',
                              'Ependymal',
                              'Microglia',
                              'PVp.03',
                              'SCN',
                              'PVp.04',
                              'Satb2/Slc18a2',
                              'MM.03',
                              'Endothelial',
                              'VMH.02',
                              'Klhl1/Ebf3',
                              'Ros1/Alk',
                              'Ebf3/Htr2c',
                              'VLMC',
                              'Tu',
                              'Tac1/Reln',
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
  labs(x = '', y= 'DE Genes on XY Chromosomes') +
  scale_y_reverse() +
  #coord_cartesian(ylim = c(1250, 50)) +
  theme_classic() +
  theme(text = element_text(family = "Arial", size = 6, color = 'black'),
        axis.text.x= element_text(family = 'Arial', color = 'black', size = 6,angle = 90, hjust = 1, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y= element_text(family = 'Arial', color = 'black', size = 6),
        legend.key.height = unit(.1, 'in'),
        legend.direction = 'horizontal',
        legend.position = c(0.5,0.25))
ggsave(filename = 'figures/celltype_number_of_Sex_DE_on_XY2.tiff', device = 'tiff', units = 'in', width = 5, height = 1.5, dpi = 600)


#!/bin/bash

ulimit -n 16000 && \
STAR \
--runThreadN 20 \
--soloType CB_UMI_Simple \
--soloCBwhitelist None \
--soloBarcodeReadLength 0 \
--soloFeatures GeneFull_Ex50pAS \
--soloCellFilter EmptyDrops_CR 10000 0.99 10 NULL NULL 500 0.01 20000 0.01 10000 \
--soloCellReadStats Standard \
--genomeDir /home/jonathan/STAR_genome \
--readFilesCommand zcat \
--readFilesIn /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B1/01.RawData/F_Fast_B1/F_Fast_B1_CKDL220027065-1A_HGK3YDSX5_S1_L001_R2_001.fastq.gz \
/home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B1/01.RawData/F_Fast_B1/F_Fast_B1_CKDL220027065-1A_HGK3YDSX5_S1_L001_R1_001.fastq.gz \
--outFileNamePrefix /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/STARsolo/F_Fast_B1; \
STAR \
--runThreadN 20 \
--soloType CB_UMI_Simple \
--soloCBwhitelist None \
--soloBarcodeReadLength 0 \
--soloFeatures GeneFull_Ex50pAS \
--soloCellFilter EmptyDrops_CR 10000 0.99 10 NULL NULL 500 0.01 20000 0.01 10000 \
--soloCellReadStats Standard \
--genomeDir /home/jonathan/STAR_genome \
--readFilesCommand zcat \
--readFilesIn /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B1/01.RawData/F_Fed_B1/F_Fed_B1_CKDL220027064-1A_HGK3YDSX5_S4_L001_R2_001.fastq.gz \
//home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B1/01.RawData/F_Fed_B1/F_Fed_B1_CKDL220027064-1A_HGK3YDSX5_S4_L001_R1_001.fastq.gz \
--outFileNamePrefix /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/STARsolo/F_Fed_B1; \
STAR \
--runThreadN 20 \
--soloType CB_UMI_Simple \
--soloCBwhitelist None \
--soloBarcodeReadLength 0 \
--soloFeatures GeneFull_Ex50pAS \
--soloCellFilter EmptyDrops_CR 10000 0.99 10 NULL NULL 500 0.01 20000 0.01 10000 \
--soloCellReadStats Standard \
--genomeDir /home/jonathan/STAR_genome \
--readFilesCommand zcat \
--readFilesIn /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B1/01.RawData/M_Fast_B1/M_Fast_B1_CKDL220027067-1A_HGK3YDSX5_S2_L001_R2_001.fastq.gz \
/home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B1/01.RawData/M_Fast_B1/M_Fast_B1_CKDL220027067-1A_HGK3YDSX5_S2_L001_R1_001.fastq.gz \
--outFileNamePrefix /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/STARsolo/M_Fast_B1; \
STAR \
--runThreadN 20 \
--soloType CB_UMI_Simple \
--soloCBwhitelist None \
--soloBarcodeReadLength 0 \
--soloFeatures GeneFull_Ex50pAS \
--soloCellFilter EmptyDrops_CR 10000 0.99 10 NULL NULL 500 0.01 20000 0.01 10000 \
--soloCellReadStats Standard \
--genomeDir /home/jonathan/STAR_genome \
--readFilesCommand zcat \
--readFilesIn /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B1/01.RawData/M_Fed_B1/M_Fed_B1_CKDL220027066-1A_HGK3YDSX5_S3_L001_R2_001.fastq.gz \
/home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B1/01.RawData/M_Fed_B1/M_Fed_B1_CKDL220027066-1A_HGK3YDSX5_S3_L001_R1_001.fastq.gz \
--outFileNamePrefix /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/STARsolo/M_Fed_B1; \
STAR \
--runThreadN 20 \
--soloType CB_UMI_Simple \
--soloCBwhitelist None \
--soloBarcodeReadLength 0 \
--soloFeatures GeneFull_Ex50pAS \
--soloCellFilter EmptyDrops_CR 10000 0.99 10 NULL NULL 500 0.01 20000 0.01 10000 \
--soloCellReadStats Standard \
--genomeDir /home/jonathan/STAR_genome \
--readFilesCommand zcat \
--readFilesIn /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B2-3/01.RawData/F_Fast_B2/F_Fast_B2_CKDL220034057-1A_HT3MHDSX5_S4_L004_R2_001.fastq.gz \
/home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B2-3/01.RawData/F_Fast_B2/F_Fast_B2_CKDL220034057-1A_HT3MHDSX5_S4_L004_R1_001.fastq.gz \
--outFileNamePrefix /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/STARsolo/F_Fast_B2; \
STAR \
--runThreadN 20 \
--soloType CB_UMI_Simple \
--soloCBwhitelist None \
--soloBarcodeReadLength 0 \
--soloFeatures GeneFull_Ex50pAS \
--soloCellFilter EmptyDrops_CR 10000 0.99 10 NULL NULL 500 0.01 20000 0.01 10000 \
--soloCellReadStats Standard \
--genomeDir /home/jonathan/STAR_genome \
--readFilesCommand zcat \
--readFilesIn /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B2-3/01.RawData/F_Fast_B3/F_Fast_B3_CKDL220034056-1A_HT3MHDSX5_S4_L003_R2_001.fastq.gz \
/home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B2-3/01.RawData/F_Fast_B3/F_Fast_B3_CKDL220034056-1A_HT3MHDSX5_S4_L003_R1_001.fastq.gz \
--outFileNamePrefix /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/STARsolo/F_Fast_B3; \
STAR \
--runThreadN 20 \
--soloType CB_UMI_Simple \
--soloCBwhitelist None \
--soloBarcodeReadLength 0 \
--soloFeatures GeneFull_Ex50pAS \
--soloCellFilter EmptyDrops_CR 10000 0.99 10 NULL NULL 500 0.01 20000 0.01 10000 \
--soloCellReadStats Standard \
--genomeDir /home/jonathan/STAR_genome \
--readFilesCommand zcat \
--readFilesIn /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B2-3/01.RawData/F_Fed_B2/F_Fed_B2_CKDL220034055-1A_HT3MHDSX5_S2_L004_R2_001.fastq.gz \
/home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B2-3/01.RawData/F_Fed_B2/F_Fed_B2_CKDL220034055-1A_HT3MHDSX5_S2_L004_R1_001.fastq.gz \
--outFileNamePrefix /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/STARsolo/F_Fed_B2; \
STAR \
--runThreadN 20 \
--soloType CB_UMI_Simple \
--soloCBwhitelist None \
--soloBarcodeReadLength 0 \
--soloFeatures GeneFull_Ex50pAS \
--soloCellFilter EmptyDrops_CR 10000 0.99 10 NULL NULL 500 0.01 20000 0.01 10000 \
--soloCellReadStats Standard \
--genomeDir /home/jonathan/STAR_genome \
--readFilesCommand zcat \
--readFilesIn /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B2-3/01.RawData/F_Fed_B3/F_Fed_B3_CKDL220034054-1A_HT3MHDSX5_S2_L003_R2_001.fastq.gz \
/home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B2-3/01.RawData/F_Fed_B3/F_Fed_B3_CKDL220034054-1A_HT3MHDSX5_S2_L003_R1_001.fastq.gz \
--outFileNamePrefix /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/STARsolo/F_Fed_B3; \
STAR \
--runThreadN 20 \
--soloType CB_UMI_Simple \
--soloCBwhitelist None \
--soloBarcodeReadLength 0 \
--soloFeatures GeneFull_Ex50pAS \
--soloCellFilter EmptyDrops_CR 10000 0.99 10 NULL NULL 500 0.01 20000 0.01 10000 \
--soloCellReadStats Standard \
--genomeDir /home/jonathan/STAR_genome \
--readFilesCommand zcat \
--readFilesIn /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B2-3/01.RawData/M_Fast_B2/M_Fast_B2_CKDL220034053-1A_HT3MHDSX5_S1_L004_R2_001.fastq.gz \
/home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B2-3/01.RawData/M_Fast_B2/M_Fast_B2_CKDL220034053-1A_HT3MHDSX5_S1_L004_R1_001.fastq.gz \
--outFileNamePrefix /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/STARsolo/M_Fast_B2; \
STAR \
--runThreadN 20 \
--soloType CB_UMI_Simple \
--soloCBwhitelist None \
--soloBarcodeReadLength 0 \
--soloFeatures GeneFull_Ex50pAS \
--soloCellFilter EmptyDrops_CR 10000 0.99 10 NULL NULL 500 0.01 20000 0.01 10000 \
--soloCellReadStats Standard \
--genomeDir /home/jonathan/STAR_genome \
--readFilesCommand zcat \
--readFilesIn /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B2-3/01.RawData/M_Fast_B3/M_Fast_B3_CKDL220034052-1A_HT3MHDSX5_S1_L003_R2_001.fastq.gz \
/home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B2-3/01.RawData/M_Fast_B3/M_Fast_B3_CKDL220034052-1A_HT3MHDSX5_S1_L003_R1_001.fastq.gz \
--outFileNamePrefix /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/STARsolo/M_Fast_B3; \
STAR \
--runThreadN 20 \
--soloType CB_UMI_Simple \
--soloCBwhitelist None \
--soloBarcodeReadLength 0 \
--soloFeatures GeneFull_Ex50pAS \
--soloCellFilter EmptyDrops_CR 10000 0.99 10 NULL NULL 500 0.01 20000 0.01 10000 \
--soloCellReadStats Standard \
--genomeDir /home/jonathan/STAR_genome \
--readFilesCommand zcat \
--readFilesIn /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B2-3/01.RawData/M_Fed_B2/M_Fed_B2_CKDL220034051-1A_HT3MHDSX5_S3_L004_R2_001.fastq.gz \
/home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B2-3/01.RawData/M_Fed_B2/M_Fed_B2_CKDL220034051-1A_HT3MHDSX5_S3_L004_R1_001.fastq.gz \
--outFileNamePrefix /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/STARsolo/M_Fed_B2; \
STAR \
--runThreadN 20 \
--soloType CB_UMI_Simple \
--soloCBwhitelist None \
--soloBarcodeReadLength 0 \
--soloFeatures GeneFull_Ex50pAS \
--soloCellFilter EmptyDrops_CR 10000 0.99 10 NULL NULL 500 0.01 20000 0.01 10000 \
--soloCellReadStats Standard \
--genomeDir /home/jonathan/STAR_genome \
--readFilesCommand zcat \
--readFilesIn /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B2-3/01.RawData/M_Fed_B3/M_Fed_B3_CKDL220034050-1A_HT3MHDSX5_S3_L003_R2_001.fastq.gz \
/home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/B2-3/01.RawData/M_Fed_B3/M_Fed_B3_CKDL220034050-1A_HT3MHDSX5_S3_L003_R1_001.fastq.gz \
--outFileNamePrefix /home/jonathan/storage_SSD/RNA_Seq/ARH_Sex_by_Nutr/STARsolo/M_Fed_B3

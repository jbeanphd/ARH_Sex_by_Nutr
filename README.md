# snRNA-seq analysis pipeline

#This repository contains the scripts used for the snRNA-seq analyses described in this study.
The code can be used to reproduce the results either starting from raw fastq files or from the processed Seurat object (.rds).

#The data generated from this study can be downloaded from:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282955

#Data can be downloaded as individual fastq files or as processed rds file. 


#From fastq files run:
#Input
- Paired-end fastq files generated from the 10x Genomics platform.

#Output
- Processed seurat object saved as an '.rds' file.

starsolo.sh - Aligns fastq files using STARsolo and generates 'barcodes.tsv', 'features.tsv', and 'matrix.mtx' files for each sample.


SoupX_scDblFnd.R - The purpose of this script is to filter out ambient RNA using SoupX and exclude possible doublets using scDblFnd. Input: barcodes.tsv, features.tsv, matrix.mtx Output：pre-proccessed Seurat object (.rds).

Annotation.R - Performs clustering and cell-type annotation and saves the proccessed Seurat object (.rds).

#By running these three scripts you will generate the same rds file that can be downloaded directly.



#From processed rds file run:

#Input
- Processed seurat object (.rds).

#Output
- Figures and result tables corresponding to main and supplementary figures in the manuscript.

Annotation.R - This will generate tSNE plots found in figures 1c, d, supplementary fig 1e.

Differential_expression.R - The purpose of this script is to find differentially expressed genes between conditions within cell-types. This script will generate bar chart found in figure 1e.

XY_Chromosome.R - The purpose of this script is to count the number of differentially expressed genes that are on autosomes or on sex chromosomes. This script will generate the inset bar chart found in figure 1e.

downsample2.R - The purpose of this script is to count the number of differentially expressed genes when downsampling to 200 nuclei per condition. This script will generate plots in figure 1f.

isolate_cell_types.R - The purpose of this script is to subcluster individual cell-types and investigate them in finer detail. This script will generate plots seen in figures 2-5, supplementary figures 15-17.

CellChat.R - The purpose of this script is to determine cell-cell communication. This script will generate plots seen in figures 6 and 7. 

NicheNet.R - The purpose of this script is to validate potential ligand–receptor interactions predicted by CellChat using a complementary analysis with NicheNet. This script will generate plots seen in supplementary figure 18.

hdWGCNA.R - The purpose of this script is to run hdWGCNA analysis. This script will generate plots seen in supplementary figure 19.

RPCA.R - This script mitigates potential batch effects in some cell types across conditions and integrates cell clusters from the current study with those from Campbell et al., 2017. This script will generate plots seen in supplementary figures 4, 8 and 10.  
  
#Each script contains detailed, line-by-line explanations of the code.

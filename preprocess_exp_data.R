if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')

library(monocle3)
library(readxl)
log2TPM_matrix <- read.delim("~/Spring20/Genomics/project/GSM4307111_GEO_processed_BC159-T_3_log2TPM_matrix_final.txt", row.names=1)
cell_metadata <- read_excel("Spring20/Genomics/project/GSE145137_BC159-T_3_All_SC_final_QC_Celltype_Information.xlsx")
gene_metadata <- data.frame(rownames(log2TPM_matrix))

# making things work for cds
log2TPM_matrix <- data.matrix(log2TPM_matrix)
rownames(cell_metadata) <- cell_metadata$Samples
rownames(gene_metadata) <- rownames(log2TPM_matrix)

# load data into cds
cds <- new_cell_data_set(as(log2TPM_matrix, "sparseMatrix"),
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)

cds_2d <- reduce_dimension(cds)
plot_cells(cds_2d, color_cells_by="cell_type")

saveRDS(cds, file = "cds.rds")
saveRDS(cell_metadata, file = "cell_met.rds")
saveRDS(gene_metadata, file = "gene_met.rds")

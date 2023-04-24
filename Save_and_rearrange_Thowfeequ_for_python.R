library(Seurat)
library(SeuratDisk)
library(data.table)
library(ggplot2)
library(pheatmap)
options(future.globals.maxSize = 16000 * 1024^2)

# Read in Normalized Thowfeequ data seperated by timepoints
# Define new labels comprising information about both Timepoint and Cell type
# E5.5
Convert("./Data/Thowfeequ_dataset/Thowfeequ_norm_counts/norm_data_55.gzip.h5ad", dest = "h5seurat", overwrite = TRUE)
Thowfeequ_E5.5_2 <- LoadH5Seurat("./Data/Thowfeequ_dataset/Thowfeequ_norm_counts/norm_data_55.gzip.h5seurat")
meta.data <- fread("./Data/Thowfeequ_dataset/Thowfeequ_raw_counts/metadata_AVE_modified_2.csv", data.table=FALSE)
rownames(meta.data) <- meta.data[,1]
meta.data[,c(1,3,4)] <- NULL
Thowfeequ_E5.5_2 <- AddMetaData(Thowfeequ_E5.5_2, metadata = meta.data)
Thowfeequ_E5.5_2 <- FindVariableFeatures(Thowfeequ_E5.5_2)
Thowfeequ_E5.5_2 <- RunPCA(Thowfeequ_E5.5_2, npcs = 30, verbose = FALSE)
Thowfeequ_E5.5_2 <- RunUMAP(Thowfeequ_E5.5_2, reduction = "pca", dims = 1:12, verbose = FALSE)
Idents(Thowfeequ_E5.5_2) <- Thowfeequ_E5.5_2$annotation
Thowfeequ_E5.5_2 <- RenameIdents(Thowfeequ_E5.5_2, 'Epi-VE'='E5.5 Epi-VE', 'anterior visceral endoderm'='E5.5 anterior visceral endoderm','ExE-VE'='E5.5 ExE-VE','unassigned'='E5.5 unassigned','epiblast'='E5.5 epiblast','extraembryonic ectoderm'='E5.5 extraembryonic ectoderm')
Thowfeequ_E5.5_2$Embryo_CellType <- Thowfeequ_E5.5_2@active.ident
DimPlot(Thowfeequ_E5.5_2, group.by = "annotation")
# E6.25
Convert("./Data/Thowfeequ_dataset/Thowfeequ_norm_counts/norm_data_625.gzip.h5ad", dest = "h5seurat", overwrite = TRUE)
Thowfeequ_E6.25_2 <- LoadH5Seurat("./Data/Thowfeequ_dataset/Thowfeequ_norm_counts/norm_data_625.gzip.h5seurat")
Thowfeequ_E6.25_2 <- AddMetaData(Thowfeequ_E6.25_2, metadata = meta.data)

Thowfeequ_E6.25_2 <- FindVariableFeatures(Thowfeequ_E6.25_2)
Thowfeequ_E6.25_2 <- RunPCA(Thowfeequ_E6.25_2, npcs = 30, verbose = FALSE)
Thowfeequ_E6.25_2 <- RunUMAP(Thowfeequ_E6.25_2, reduction = "pca", dims = 1:12, verbose = FALSE)
Idents(Thowfeequ_E6.25_2) <- Thowfeequ_E6.25_2$annotation
Thowfeequ_E6.25_2 <- RenameIdents(Thowfeequ_E6.25_2, 'Epi-VE'='E6.25 Epi-VE', 'anterior visceral endoderm'='E6.25 anterior visceral endoderm','ExE-VE'='E6.25 ExE-VE','unassigned'='E6.25 unassigned','epiblast'='E6.25 epiblast','extraembryonic ectoderm'='E6.25 extraembryonic ectoderm')
Thowfeequ_E6.25_2$Embryo_CellType <- Thowfeequ_E6.25_2@active.ident
DimPlot(Thowfeequ_E6.25_2, group.by = "annotation")

Thowfeequ <- merge(Thowfeequ_E5.5_2, Thowfeequ_E6.25_2)

# Add metadata about cell origin
Thowfeequ$orig_ident <- "Embryo"

# Save Data for opening in python 
SaveH5Seurat(Thowfeequ, filename = "./Data/Thowfeequ_dataset/Thowfeequ.h5Seurat")
Convert("./Data/Thowfeequ_dataset/Thowfeequ.h5Seurat", dest = "h5ad")

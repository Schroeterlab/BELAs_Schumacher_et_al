Analysis of the in vitro scRNAseq dataset
================

### Including integration for batch correction, clustering and DEG identification and visualization

Import required packages

``` r
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(SeuratDisk)
  library(pheatmap)
})
options(future.globals.maxSize = 8000 * 1024^2)
```

Load in vitro data (CellRanger Output)

``` r
read10x_filter_seurat <- function(matrix_path, sample_id){
  raw_counts <- Read10X(matrix_path)
  raw_counts <- CreateSeuratObject(counts=raw_counts, project=sample_id) 
  raw_counts$percent.MT <- PercentageFeatureSet(raw_counts, pattern="^mt-" )
  raw_counts <- subset(raw_counts, subset = nFeature_RNA > 4000 & percent.MT < 10)
  return(raw_counts)
}
path_BELAs = "./Data/BELAs/"
path_EpiCysts ="./Data/EpiCysts/"
path_PrECysts = "./Data/VECysts/"

BELAs = read10x_filter_seurat(path_BELAs, "BELAs")
```

    ## as(<dgTMatrix>, "dgCMatrix") is deprecated since Matrix 1.5-0; do as(., "CsparseMatrix") instead

``` r
EpiCysts = read10x_filter_seurat(path_EpiCysts, "Epi Cysts")
VECysts = read10x_filter_seurat(path_PrECysts, "VE Cysts")
```

To remove batch effects between sample we integrated the sample with
RPCA based integration in Seurat.

``` r
data.list <- list(BELAs, EpiCysts, VECysts)
data.list <- lapply(X = data.list, FUN = SCTransform, verbose = FALSE)
data.list <- lapply(X = data.list, FUN = RunPCA, verbose = FALSE)

data.list.features <- SelectIntegrationFeatures(object.list = data.list, 
                                                nfeatures = 3000, verbose = FALSE)
data.list <- PrepSCTIntegration(object.list = data.list, 
                                anchor.features = data.list.features, verbose = FALSE)
anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", 
                                  dims = 1:30, anchor.features = data.list.features, 
                                  reduction = "rpca", verbose = FALSE)
in_vitro <- IntegrateData(anchorset = anchors, normalization.method = "SCT", 
                          dims = 1:30, verbose = FALSE)
```

Run dimensionality reduction and plot UMAP annotated by sample of
origin.

``` r
DefaultAssay(in_vitro) <- "integrated"
in_vitro <- RunPCA(in_vitro, npcs = 30, verbose = FALSE)
in_vitro <- RunUMAP(in_vitro, reduction = "pca", dims = 1:12, verbose = FALSE)
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

``` r
DimPlot(in_vitro, reduction = "umap", group.by = "orig.ident", pt.size = 1,
        cols = c("BELAs"="#F8766D",  "Epi Cysts"="#D3D30B", "VE Cysts"="#619CFF")) + 
  theme(aspect.ratio = 1, axis.text= element_blank(), axis.ticks = element_blank())
```

![](in_vitro_Integration+FeaturePlots+Clustering+AVE_DEGs_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

The expression of known VE and Epi marker genes was visualized as
FeaturePlots.

``` r
DefaultAssay(in_vitro) <- "SCT"
markers <- c("Gata6", "Sox17", "Dab2", "Cubn", "Pou5f1", "Sox2", "Nanog", "Fgf4")
cutoffs<- c(NA, NA, NA, NA, NA, "2", "2", "1.5")
for (i in 1:length(markers)){
  print(FeaturePlot(in_vitro, features = markers[i], slot = "data", 
                    max.cutoff=cutoffs[i], pt.size = 1) + 
          scale_colour_gradientn(colours = viridis::cividis(100)) + theme(aspect.ratio = 1) +
          theme(axis.text= element_blank(), axis.ticks = element_blank()))
}
```

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](in_vitro_Integration+FeaturePlots+Clustering+AVE_DEGs_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](in_vitro_Integration+FeaturePlots+Clustering+AVE_DEGs_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](in_vitro_Integration+FeaturePlots+Clustering+AVE_DEGs_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](in_vitro_Integration+FeaturePlots+Clustering+AVE_DEGs_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](in_vitro_Integration+FeaturePlots+Clustering+AVE_DEGs_files/figure-gfm/unnamed-chunk-5-5.png)<!-- -->

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](in_vitro_Integration+FeaturePlots+Clustering+AVE_DEGs_files/figure-gfm/unnamed-chunk-5-6.png)<!-- -->

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](in_vitro_Integration+FeaturePlots+Clustering+AVE_DEGs_files/figure-gfm/unnamed-chunk-5-7.png)<!-- -->![](in_vitro_Integration+FeaturePlots+Clustering+AVE_DEGs_files/figure-gfm/unnamed-chunk-5-8.png)<!-- -->

Clustering and visualization in UMAP space.

``` r
DefaultAssay(in_vitro) <- "integrated"
in_vitro <- FindNeighbors(in_vitro, reduction = "pca", verbose = FALSE)
in_vitro <- FindClusters(in_vitro, resolution = 0.2, verbose = FALSE)
in_vitro <- RenameIdents(in_vitro, '0' = '3', '1' = '1', '2' = '2', '3' = '4')
levels(x = in_vitro) <- c('1', '2', '3', '4')
in_vitro$Cluster <- in_vitro@active.ident

# Umap by Cluster
DimPlot(in_vitro, reduction = "umap", group.by = "Cluster", pt.size = 1,
        cols = c("#8ABB93", "#0F5C35", "#C2A3C1", "#F6A348")) + 
  theme(aspect.ratio = 1, axis.text= element_blank(), axis.ticks = element_blank())
```

![](in_vitro_Integration+FeaturePlots+Clustering+AVE_DEGs_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

The proportions of cells from each sample in each cluster are shown as a
heatmap, using the custom pl_cell_frac_pheatmap_v2 function.

``` r
in_vitro$Sample <- factor(in_vitro$orig.ident, levels = c("BELAs","VE Cysts","Epi Cysts"))
source("./func_cell_fraction_heatmap.R")
pl_cell_frac_pheatmap_v2(in_vitro,
                         column_data = "Cluster",
                         row_data = "Sample",
                         include_absolute_values = TRUE,
                         ratio = "column")
```

![](in_vitro_Integration+FeaturePlots+Clustering+AVE_DEGs_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Find and visualize the differentially expressed genes between the AVE
cluster and all other VE-like cells, including cells from BELAs and VE
cysts.

``` r
####
### DOES NOT REPRODUCE EXACT HEATMAP ###
#####
# Display differentially expressed genes of AVE as heatmap
# Problem seems to be: SCTransform before integration results in multiple models so one has to run it again. This results in small changes in list, not reproducing list of poster.
# Prepare Data for subsetting
in_vitro$orig_ident_CellType <- paste(in_vitro$Sample, in_vitro$Cluster, sep = ": ")
Idents(in_vitro) <- "orig_ident_CellType"

# DE Genes for AVE cluster
in_vitro <- PrepSCTFindMarkers(in_vitro)
```

    ## Found 3 SCT models. Recorrecting SCT counts using minimum median counts: 58554

``` r
AVE.dif_AVE <- FindMarkers(in_vitro, assay = "SCT", ident.1 = "BELAs: 4", 
                           ident.2 = c("BELAs: 3", "VE Cysts: 3"), only.pos = TRUE, 
                           verbose = FALSE)
AVE.dif_PrE <- FindMarkers(in_vitro, ident.1 = c("BELAs: 3", "VE Cysts: 3"), 
                           ident.2 = "BELAs: 4", only.pos = TRUE, verbose = FALSE)


heat_data <- subset(in_vitro, idents = c("BELAs: 4", "BELAs: 3", "VE Cysts: 3"))
heat_data$orig_ident_CellType <- factor(heat_data@active.ident, 
                                        levels = c("BELAs: 4", "BELAs: 3", "VE Cysts: 3"),
                                        ordered = TRUE)

heat_data <- SCTransform(heat_data, return.only.var.genes = FALSE, verbose = FALSE) # To get SCT value for all genes

genes <- c(rownames(AVE.dif_AVE[order(-AVE.dif_AVE$avg_log2FC), ])[1:30], 
           rownames(AVE.dif_PrE[order(-AVE.dif_PrE$avg_log2FC), ])[1:30])
heat_frame <- heat_data@assays[["SCT"]]@scale.data[genes,]
heatmap_col_anno <- data.frame("Cluster" = heat_data$orig_ident_CellType)
heat_frame <- heat_frame[,order(heat_data$orig_ident_CellType)]
pheatmap(heat_frame, scale = "none", cluster_rows = FALSE, 
         col=viridis::cividis(100), show_colnames = FALSE,
         cluster_cols = FALSE, annotation_col = heatmap_col_anno, gaps_col = c(62,532), 
         cellwidth = 0.2, cellheight = 10, gaps_row = 30)
```

![](in_vitro_Integration+FeaturePlots+Clustering+AVE_DEGs_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
# For horizontal heatmap
# pheatmap(t(heat_frame), scale = "none", cluster_rows = FALSE, 
#          col=viridis::cividis(100), show_rownames = FALSE,
#          cluster_cols = FALSE, annotation_row = heatmap_col_anno, gaps_row = c(62,532), 
#          cellwidth = 10, cellheight = 0.2, gaps_col = 30, angle_col = 90)
```

For the visualization of single genes of interest, or example Nodal,
FeaturePlots were used.

``` r
DefaultAssay(in_vitro) <- "SCT"
FeaturePlot(in_vitro, feature = "Nodal") + theme(aspect.ratio = 1) +
  theme(axis.text= element_blank(), axis.ticks = element_blank())
```

![](in_vitro_Integration+FeaturePlots+Clustering+AVE_DEGs_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Export dataset (as h5ad file) for ingest integration with in vivo
datasets performed in python. Only raw counts are used for this.

``` r
in_vitro$orig_ident <- in_vitro$orig.ident
DefaultAssay(in_vitro) <- "RNA"
in_vitro[['integrated']] <- NULL
in_vitro[['SCT']] <- NULL

# SaveH5Seurat(in_vitro, filename = "./Data/in_vitro.h5Seurat")
# Convert("./Data/in_vitro.h5Seurat", dest = "h5ad")
```

single cell RNA-sequencing analysis of BELAs v20220223
================
Max Fernkorn
2/23/2022

-   [1 BELAs, PrE Cysts and Epi Cysts dataset
    preprocessing](#1-belas-pre-cysts-and-epi-cysts-dataset-preprocessing)
    -   [1.1 Data loading and filtering](#11-data-loading-and-filtering)
    -   [1.2 Basic data analysis](#12-basic-data-analysis)
-   [2. Visualisation and Annotation of the in vitro
    dataset](#2-visualisation-and-annotation-of-the-in-vitro-dataset)
-   [3. Loading and normalization of the embryo dataset from Nowotschin
    et
    al](#3-loading-and-normalization-of-the-embryo-dataset-from-nowotschin-et-al)
-   [4. Integration analysis](#4-integration-analysis)
    -   [4.1 General Integration
        approach](#41-general-integration-approach)
    -   [4.2 Integration with all embryonic
        lineages](#42-integration-with-all-embryonic-lineages)
        -   [4.2.1 Integration and
            Quantification](#421-integration-and-quantification)
        -   [4.2.2 Ectodermal marker gene
            expression](#422-ectodermal-marker-gene-expression)
    -   [4.3 Integration of the Epiblast
        lineage](#43-integration-of-the-epiblast-lineage)
    -   [4.4 Integration of the primitive Endoderm
        lineage](#44-integration-of-the-primitive-endoderm-lineage)
        -   [4.4.1 Integration and
            Quantification](#441-integration-and-quantification)
        -   [4.4.2 Analysis of embryonic time point specific genes in
            BELAs](#442-analysis-of-embryonic-time-point-specific-genes-in-belas)
        -   [4.4.3 AVE marker gene
            determination](#443-ave-marker-gene-determination)

# 1 BELAs, PrE Cysts and Epi Cysts dataset preprocessing

## 1.1 Data loading and filtering

Demultiplexing of the sequencing data, alignment to the mouse genome
(mm10) and read quantification was performed with CellRanger (10x
Genomics, v4.0.0) for each sample of origin (BELAs, Epi cysts and PrE
Cysts) seperatly. After loading into R and creation of a Seurat Object,
cells with more than 10% of the reads aligned to mitochondrial genes as
well as less than 4000 different detected features were filtered out.

``` r
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(pheatmap)
  library(data.table)
  library(SeuratDisk)
  library(dendsort)
})
options(future.globals.maxSize = 8000 * 1024^2)

# Load in vitro data and filter for mt% and nFeature
read10x_filter_seurat <- function(matrix_path, sample_id){
  raw_counts <- Read10X(matrix_path)
  raw_counts <- CreateSeuratObject(counts=raw_counts, project=sample_id) 
  raw_counts$percent.MT <- PercentageFeatureSet(raw_counts, pattern="^mt-" )
  raw_counts <- subset(raw_counts, subset = nFeature_RNA > 4000 & percent.MT < 10)
  return(raw_counts)
}
path_BELAs = "./Data/BELAs/filtered_feature_bc_matrix/"
path_EpiCysts ="./Data/EpiCysts/filtered_feature_bc_matrix/"
path_PrECysts = "./Data/PrECysts/filtered_feature_bc_matrix/"

BELAs = read10x_filter_seurat(path_BELAs, "BELAs")
EpiCysts = read10x_filter_seurat(path_EpiCysts, "Epi Cysts")
PrECysts = read10x_filter_seurat(path_PrECysts, "PrE Cysts")
```

## 1.2 Basic data analysis

Datasets from all three sample of origins were merged and normalization
was perfomred using `SCTransform()` with default parameters. Principal
component analysis was followed by Uniform Manifold Approximation and
Projection (UMAP) for visualization using `RunPCA()` and `RunUMAP()`.

``` r
# Clusters based on total dataset
in_vitro <- merge(BELAs, y = c(EpiCysts, PrECysts))
in_vitro <- SCTransform(in_vitro, verbose = FALSE)
in_vitro <- RunPCA(in_vitro, verbose = FALSE)
in_vitro <- RunUMAP(in_vitro, dims = 1:12, verbose = FALSE)
```

# 2. Visualisation and Annotation of the in vitro dataset

To provide an overview about the dataset a UMAP Plot with the cells
colored according to the sample of origin was created.

``` r
# Umap plot for Figure 4B
DimPlot(in_vitro, reduction = "umap", group.by = "orig.ident", 
        cols = c( "BELAs"="#F8766D",  "Epi Cysts"="#D3D30B", "PrE Cysts"="#619CFF"), pt.size = 2) + 
  theme(aspect.ratio = 1, axis.text= element_blank(), axis.ticks = element_blank())
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

With Louvain clustering analysis using `FindNeighbors()`and
`FindClusters()` two groups of cells were identified within the BELA
cells. The resolution for clustering was set in order to obtain four
clusters to determine, which sample of origin divides into two seperate
clusters first. The identified clusters were renamed and shown in a UMAP
plot colored by the four clusters.

``` r
# Rename seurat clusters and save as Named CellTypes for later reference
in_vitro <- FindNeighbors(in_vitro, dims = 1:12, verbose = FALSE)
in_vitro <- FindClusters(in_vitro, resolution = 0.2, verbose = FALSE) # 4 clusters for resolutions between 0.07 and 0.27
in_vitro <- RenameIdents(in_vitro, '0' = '1', '1' = '2', '2' = '4', '3' = '3')
in_vitro$CellType <- in_vitro@active.ident

# Umap plot for Figure 4C
DimPlot(in_vitro, reduction = "umap", cols = c("1"="#01663A","2"= "#05CED8", "3"="#4C37B5", "4"="#BC0ABC"),
        pt.size = 2, group.by = "CellType") +
  theme(aspect.ratio = 1, axis.text= element_blank(), axis.ticks = element_blank())  
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

To annotate the identified clusters based on sample of origin a heatmap
with the relative proportion of each sample of origin per cluster was
used.

``` r
heatmap_prop = data.frame("Sample" = in_vitro@meta.data[["orig.ident"]],
                          "CellType" = in_vitro@meta.data[["CellType"]])
heatmap_prop <- prop.table(t(table(heatmap_prop)), margin = 2)[c(1,2,4,3),c(3,2,1)]
pheatmap(heatmap_prop, scale = "none", cluster_rows = FALSE,
         cluster_cols = FALSE, col=viridis::cividis(100), display_numbers =  TRUE, number_color = "black",
         gaps_col = 3, cellwidth = 50, cellheight = 50, fontsize_number = 14, breaks = seq(0,1,length = 100))
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

To annotate the identified clusters based on the cell lineage they
resemble expression of known marker genes for Epiblast and PrE lineages
were shown as a `FeaturePlot()` with the total dataset. We emphasized
changes in genes expression for selected marker genes by thresholding
the maximum value on the logarithmic scale.

``` r
# Feature Plots for Figure 4E
markers <- c("Gata6", "Sox17", "Dab2", "Cubn", "Pou5f1", "Sox2", "Nanog", "Fgf4")
cutoffs<- c("2", NA, NA, NA, NA, "2", NA, "1.5")

for (i in 1:length(markers)){
  print(FeaturePlot(in_vitro, features = markers[i], slot = "data", max.cutoff=cutoffs[i], pt.size = 1) +
    scale_colour_gradientn(colours = viridis::cividis(100)) + theme(aspect.ratio = 1) +
    theme(axis.text= element_blank(), axis.ticks = element_blank()))
}
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-6-5.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-6-6.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-6-7.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-6-8.png)<!-- -->

Based on this new labels were assigned to the identified clusters
indicating the main sample of origin characterizing the cluster (BELAs,
Epi cysts or PrE cysts) together with the information about the lineage
(Epi or PrE). For later comparison with an embryo dataset addtional
labels for the in vitro cells were assigned.

``` r
# Annotate in_vitro dataset
in_vitro <- RenameIdents(in_vitro, '1' = 'BELA-Epi', '2' = 'cyst-Epi', '4' = 'BELA-PrE', '3' = 'cyst-PrE')
in_vitro$CellType_Named <- in_vitro@active.ident
in_vitro$Origin <- "mESC"
in_vitro$Lineage <- in_vitro$orig.ident
in_vitro$Lineage_CellType <- in_vitro$CellType_Named
```

# 3. Loading and normalization of the embryo dataset from Nowotschin et al

Count and metadata tables of the embryonic dataset by Nowotschin et al
can be downloaded as text files [here](https://endoderm-explorer.com)
(Downloads tab). After loading single time points between E4.5 and E6.5
were selected to add Timepoint specific CellType annotations. After
merging the time points, annotation of the three major embryonic
lineages was performed based on the cell type annotations from
Nowotschin et al. All following plots are generated based on only one
replicate from the embryo dataset, since batch effects complicated the
analysis. Which replicate is selected does not change any conclusions
drawn. For normalization `SCTransform()` was used to ensure
comparability with the in vitro dataset.

``` r
counts<-fread("./Data/sc_endoderm_all_cells_counts.csv", data.table=FALSE)
meta.data <- fread("./Data/sc_endoderm_all_cells_metadata.csv", data.table=FALSE)

rownames(counts) <- counts[,1]
Cellnames <- counts[,1]
counts[,1] <- NULL
counts <- data.frame(t(counts))
colnames(counts) <- Cellnames
rownames(meta.data) <- meta.data[,1]
meta.data[,1] <- NULL

sc_endo <- CreateSeuratObject(counts, meta.data = meta.data)
Idents(object = sc_endo) <- "Timepoint"
sc_endo$Origin <- "Embryo"

E4.5 <- subset(x = sc_endo, idents = c("E4.5"))
Idents(E4.5) <- E4.5$CellType
new.Sample <- c("E4.5 PrE",  "E4.5 EPI",  "E4.5 TE")
names(new.Sample) <- levels(E4.5)
E4.5 <- RenameIdents(E4.5, new.Sample)
E4.5$Timepoint_CellType <- E4.5@active.ident

E5.5 <- subset(x = sc_endo, idents = c("E5.5"))
Idents(E5.5) <- E5.5$CellType
new.Sample <- c("E5.5 ExE",  "E5.5 EPI",  "E5.5 emVE", "E5.5 exVE", "E5.5 VE")
names(new.Sample) <- levels(E5.5)
E5.5 <- RenameIdents(E5.5, new.Sample)
E5.5$Timepoint_CellType <- E5.5@active.ident

E6.5 <- subset(x = sc_endo, idents = c("E6.5"))
Idents(E6.5) <- E6.5$CellType
new.Sample <- c("E6.5 exVE","E6.5 emVE","E6.5 EPI", "E6.5 Mes", "E6.5 TE" )
names(new.Sample) <- levels(E6.5)
E6.5 <- RenameIdents(E6.5, new.Sample)
E6.5$Timepoint_CellType <- E6.5@active.ident

sc_endo <- merge(E4.5, y = c(E5.5, E6.5), project = "Embryo")

# Create Lineage Classification
sc_endo$Lineage <- sc_endo$CellType
Idents(sc_endo) <- sc_endo$Lineage
new.Sample <- c("Endodermal Lineages","Epiblast Lineages","Ectodermal Lineages","Ectodermal Lineages",
                "Endodermal Lineages","Endodermal Lineages","Endodermal Lineages","Epiblast Lineages")
names(new.Sample) <- levels(sc_endo)
sc_endo <- RenameIdents(sc_endo, new.Sample)
sc_endo$Lineage <- sc_endo@active.ident
sc_endo$Lineage_CellType <- sc_endo$Lineage

Idents(object = sc_endo) <- "orig.ident"
#Subset and normalize Replicate
sc_endo_rep1 <- subset(sc_endo, idents = c("Lib1-1","Lib2-3","Lib3-1"))
sc_endo_rep1 <- SCTransform(sc_endo_rep1, verbose = FALSE)
```

# 4. Integration analysis

## 4.1 General Integration approach

For further dataset integration-based analysis the following function
was used. It was based on the recommended approach for SCTransformed
from the Seurat Package.

``` r
Seurat_integration_SCT <- function(Embryo_dataset, invitro_dataset){
  data.list <- list(Embryo_dataset, invitro_dataset)
  data.list.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
  data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = data.list.features, 
                                  verbose = FALSE)
  anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", dims = 1:30, 
                                    anchor.features = data.list.features)
  integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
  
  DefaultAssay(integrated) <- "integrated"
  integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
  integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:12, verbose = FALSE)
  return(integrated)
}
```

## 4.2 Integration with all embryonic lineages

### 4.2.1 Integration and Quantification

We integrated the total in vitro dataset with E4.5 to E6.5 embryo with
the integration algorithm of the Seurat package to determine which
embryonic celltypes and time points are resembled by the in vitro
differentiated cells. Clustering analsysis was performed in order to
seperate the three main embryonic lineages.

``` r
integrated_rep1 <- Seurat_integration_SCT(sc_endo_rep1, in_vitro)
integrated_rep1 <- FindNeighbors(integrated_rep1, dims = 1:12)
integrated_rep1 <- FindClusters(integrated_rep1, resolution = 0.02, verbose = FALSE) # 3 clusters for res = 0.016 to 0.039
```

A UMAP plot was used to visualize the annotated cell lineage of the
embryo as well as the in vitro annotations based on clustering and
marker gene expression from Figure 4D.

``` r
DimPlot(integrated_rep1, reduction = "umap", group.by = "Lineage_CellType", pt.size = 0.9, shape.by = "Origin",
        cols=c('Endodermal Lineages'='#d1b3d1','Epiblast Lineages'='#94c39b','Ectodermal Lineages' = '#D8D8D8',
               'BELA-Epi'='#01663A','BELA-PrE'='#BC0ABC','cyst-PrE'='#4C37B5', 'cyst-Epi'='#05CED8')) + 
  scale_shape_manual(values = c("Embryo" = 16, "mESC" = 17)) + theme(aspect.ratio = 1) +
  theme(axis.text= element_blank(), axis.ticks = element_blank())
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Louvain clustering grouped the cells into three clusters, which was
visualized in a UMAP plot.

``` r
DimPlot(integrated_rep1, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.9, shape.by = "Origin",
        cols=c('0'='#d1b3d1','1'='#94c39b','2' = '#D8D8D8')) + 
  scale_shape_manual(values = c("Embryo" = 16, "mESC" = 17)) + theme(aspect.ratio = 1) +
  theme(axis.text= element_blank(), axis.ticks = element_blank())
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Again, relative proportions of each embryonic lineage in every cluster
were shown in a heatmap. Thereby each cluster could be linked to a
embryonic lineage. Additionally, proportions of the cell types
differentiated in vitro (see Figure 4D) for each cluster were shown in
the same heatmap.

``` r
heatmap_prop = data.frame("Cluster" = integrated_rep1@active.ident, 
                          "Lineage_CellType" = integrated_rep1@meta.data[["Lineage_CellType"]])
heatmap_prop <- prop.table(table(heatmap_prop), margin = 2)[,c(5,7,6,1,2,3,4)]
pheatmap(heatmap_prop, scale = "none", cluster_rows = FALSE,
         cluster_cols = FALSE, col=viridis::cividis(100), display_numbers =  TRUE,
         number_color = "black", fontsize_number = 14, cellwidth = 50, cellheight = 50,
         gaps_col = 3, width = 5, height = 3.5, breaks = seq(0,1,length = 100))
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

### 4.2.2 Ectodermal marker gene expression

The expression of ectordermal marker genes was shown in a
`FeaturePlot()` in the integrated datasets, as well as in the cluster
defined by the embryonic ectodermal lineage divided by the the origin of
the cells (embryo or in vitro).

``` r
DefaultAssay(integrated_rep1) <- "SCT"
Idents(object = integrated_rep1) <- "seurat_clusters"
Ecto <- subset(integrated_rep1, idents = "2")

Idents(object = Ecto) <- "Origin"
Ecto_Embryo <- subset(Ecto, idents = "Embryo")
Idents(object = Ecto) <- "orig.ident"
Ecto_mESC <- subset(Ecto, idents = "BELAs")

marker_genes <- c("Cdx2", "Gata3", "Elf5", "Krt7", "Tead4", "Tfap2c", "Ascl2", "Bmp4", "Eomes", "Esrrb")
for (i in 1:length(marker_genes)){
  max = max(integrated_rep1@assays[["SCT"]]@data[marker_genes[i],],
            Ecto_mESC@assays[["SCT"]]@data[marker_genes[i],], Ecto_Embryo@assays[["SCT"]]@data[marker_genes[i],])
  print(FeaturePlot(integrated_rep1, features = marker_genes[i], slot = "data") +
    scale_colour_gradientn(colours = viridis::cividis(100), limits = c(0,max)) +
    theme(aspect.ratio = 1, axis.text= element_blank(), axis.ticks = element_blank()))
  print(FeaturePlot(Ecto_mESC, features = marker_genes[i], slot = "data", pt.size = 2) +ylim(6.5,12) +xlim(0,6)+
    scale_colour_gradientn(colours = viridis::cividis(100), limits = c(0,max)) + 
    theme(aspect.ratio = 1, axis.text= element_blank(), axis.ticks = element_blank()))
  print(FeaturePlot(Ecto_Embryo, features = marker_genes[i], slot = "data", pt.size = 2) +ylim(6.5,12)+xlim(0,6)+
    scale_colour_gradientn(colours = viridis::cividis(100), limits = c(0,max)) + 
    theme(aspect.ratio = 1, axis.text= element_blank(), axis.ticks = element_blank()))
}
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-3.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-4.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-5.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-6.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-7.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-8.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-9.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-10.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-11.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-12.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-13.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-14.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-15.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-16.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-17.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-18.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-19.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-20.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-21.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-22.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-23.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-24.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-25.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-26.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-27.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-28.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-29.png)<!-- -->![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-14-30.png)<!-- -->

## 4.3 Integration of the Epiblast lineage

For a more detailed analysis of the epiblast lineage of the embryo
(Epiblast and Mesoderm cell types) together with the in vitro cells
originating from BELAs or Epi cysts and clustering within the two Epi
clusters integration was performed after specifically selecting these
cells from both datasets.

``` r
# Integration of Epi Lineage (Fig 5, Suppl. 1)
# Subset BELA dataset
in_vitro$Timepoint_CellType <- in_vitro$CellType_Named
Idents(in_vitro) <- in_vitro$orig.ident
in_vitro_Epi <- subset(in_vitro, idents = c("BELAs", "Epi Cysts"))
Idents(in_vitro_Epi) <- in_vitro_Epi$CellType_Named
in_vitro_Epi <- subset(in_vitro_Epi, idents = c("BELA-Epi", "cyst-Epi"))
in_vitro_Epi <- SCTransform(in_vitro_Epi, verbose = FALSE)
# Subset Embryo dataset
Idents(object = sc_endo_rep1) <- "CellType"
sc_endo_rep1_Epi <- subset(sc_endo_rep1, idents = c("EPI", "Mes"))
sc_endo_rep1_Epi <- SCTransform(sc_endo_rep1_Epi, verbose = FALSE)
# Integration
integrated_rep1_Epi <- Seurat_integration_SCT(sc_endo_rep1_Epi, in_vitro_Epi)
```

To represent the selected cells after integration a UMAP plot was used
coloring the cells either by their time point and cell type (for the
embryonic cells) or by the annotated cluster from Fig 4D for the in
vitro differentiated cells.

``` r
# Fig 5 Supp 1 A
DimPlot(integrated_rep1_Epi, reduction = "umap", shape.by = "Origin", group.by = "Timepoint_CellType", pt.size=2,
        cols=c('E4.5 EPI'='#66CBE2', 'E5.5 EPI'='#95C39C', 'E6.5 EPI'='#5E7673', 
               'E6.5 Mes'='#979962','BELA-Epi'='#03673B','cyst-Epi'='#1D5E9A')) + 
  scale_shape_manual(values = c("Embryo" = 16, "mESC" = 17)) +
  theme(aspect.ratio = 1, axis.text= element_blank(), axis.ticks = element_blank())
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

The same cells were reclustered with Louvain clustering with a
resolution resulting in three clusters. These two UMAP plots revealed
differences between Epi cysts and the Epi-like cells of BELAs.

``` r
# Fig 5 Supp 1 B
integrated_rep1_Epi <- FindNeighbors(integrated_rep1_Epi, reduction = "pca", dims = 1:30)
integrated_rep1_Epi <- FindClusters(integrated_rep1_Epi, resolution = 0.03, verbose = FALSE)
integrated_rep1_Epi <- RenameIdents(integrated_rep1_Epi, '0' = '2', '1' = '1', '2' = '3')
DimPlot(integrated_rep1_Epi, reduction = "umap", shape.by = "Origin", pt.size = 2, 
        cols=c('1'='#7CAA00', '2'='#F4766D','3'= '#C37CFF')) +
  scale_shape_manual(values = c("Embryo" = 16, "mESC" = 17)) +
  theme(aspect.ratio = 1, axis.text= element_blank(), axis.ticks = element_blank())
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

In order to quantify these differences, the relative proportion of cells
per cluster were shown in a heatmap based on either their time point and
cell type (for the embryonic cells) or by the annotated cluster from Fig
4D for the in vitro differentiated cells.

``` r
# Fig 5 Supp 1 C
heatmap_prop = data.frame("Cluster" = integrated_rep1_Epi@active.ident, 
                          "Timepoint_CellType" = integrated_rep1_Epi@meta.data[["Timepoint_CellType"]])
heatmap_prop <- prop.table(table(heatmap_prop), margin = 2)[c(2,1,3),c(3,4,5,6,1,2)]
pheatmap(heatmap_prop, scale = "none", cluster_rows = FALSE, cellwidth = 50, cellheight = 50,
         cluster_cols = FALSE, col=viridis::cividis(100), display_numbers =  TRUE, number_color = "black", 
         fontsize_number = 14, gaps_col = 4, width = 5, height = 3.5, breaks = seq(0,1,length = 100)) 
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

## 4.4 Integration of the primitive Endoderm lineage

### 4.4.1 Integration and Quantification

For the endodermal lineage of both the embryonic as well as the in vitro
cells a similar approach was used to find differences between PrE cysts
and PrE-like cells of BELAs. Integration was performed after
specifically selecting all PrE-derived or PrE-like cells from both
datasets. Only cells annotated as visceral endoderm were excluded from
the analysis as described in Materials and Methods.

``` r
# Integration of PrE Lineage (Fig 5 A-C)
Idents(in_vitro) <- in_vitro$orig.ident
in_vitro_PrE <- subset(in_vitro, idents = c("BELAs", "PrE Cysts"))
Idents(in_vitro_PrE) <- in_vitro_PrE$CellType_Named
in_vitro_PrE <- subset(in_vitro_PrE, idents = c("BELA-PrE", "cyst-PrE"))
in_vitro_PrE <- SCTransform(in_vitro_PrE, verbose = FALSE)
# Subset Embryo dataset
Idents(object = sc_endo_rep1) <- "CellType"
sc_endo_rep1_PrE <- subset(sc_endo_rep1, idents = c("PrE", "emVE", "exVE")) # deleted "VE"
sc_endo_rep1_PrE <- SCTransform(sc_endo_rep1_PrE, verbose = FALSE)
# Integration
integrated_rep1_PrE <- Seurat_integration_SCT(sc_endo_rep1_PrE, in_vitro_PrE)
```

The integrateded dataset was shown in a UMAP plot with cells colored by
either by their time point and cell type (for the embryonic cells) or by
the annotated cluster from Fig 4D for the in vitro differentiated cells.

``` r
# Fig 5A
DimPlot(integrated_rep1_PrE, reduction = "umap", group.by = "Timepoint_CellType", pt.size=2, shape.by = "Origin",
        cols=c('E4.5 PrE'='#F7D0F7', 'E5.5 emVE'='#D3ABF4', 'E5.5 exVE'='#C1B7F7', 'E6.5 emVE'='#7BA2FF', 
               'E6.5 exVE'='#B5CEFF', 'BELA-PrE'='#BC0ABC','cyst-PrE'='#4C37B5')) + 
  scale_shape_manual(values = c("Embryo" = 16, "mESC" = 17)) + theme(aspect.ratio = 1) +
  theme(axis.text= element_blank(), axis.ticks = element_blank())
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

Louvain clustering on the integrated dataset with a resolution to obtain
four clusters was represented in another UMAP plot.

``` r
# Fig 5B
# Clustering
integrated_rep1_PrE <- FindNeighbors(integrated_rep1_PrE, reduction = "pca", dims = 1:30)
integrated_rep1_PrE <- FindClusters(integrated_rep1_PrE, resolution = 0.15, verbose = FALSE)
integrated_rep1_PrE <- RenameIdents(integrated_rep1_PrE, '0' = '1', '1' = '2', '2' = '3', '3' = '4')
integrated_rep1_PrE$seurat_clusters <- integrated_rep1_PrE@active.ident

DimPlot(integrated_rep1_PrE, reduction = "umap", pt.size=2, shape.by = "Origin", group.by = "seurat_clusters",
        cols=c('1'='#F4766D', '2'='#7CAA00','3' = '#C37CFF', '4' = '#00BBC0')) + 
  scale_shape_manual(values = c("Embryo" = 16, "mESC" = 17)) + 
  theme(aspect.ratio = 1, axis.text= element_blank(), axis.ticks = element_blank())
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

Since the UMAP plots above, together with the clustering analysis,
suggest a difference between BELA cells and PrE cyst cells in terms of
their similarity to either E6.5 emVE or E6.5 exVE cells, the proportion
of cells clustering with one of these clusters were compared in a
heatmap construced as described above.

``` r
# Fig 5C
# plot Heatmap for 3 clusters
heatmap_prop = data.frame("Cluster" = integrated_rep1_PrE@active.ident, 
                          "Timepoint_CellType" = integrated_rep1_PrE@meta.data[["Timepoint_CellType"]])
heatmap_prop <- prop.table(table(heatmap_prop), margin = 2)[,c(3,4,5,6,7,1,2)]
pheatmap(heatmap_prop, scale = "none", cluster_rows = FALSE,cellwidth = 50, cellheight = 50,
         cluster_cols = FALSE, col=viridis::cividis(100), display_numbers =  TRUE, number_color = "black", 
         fontsize_number = 14, gaps_col = 5, width = 5, height = 3.5, breaks = seq(0,1,length = 100)) 
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

### 4.4.2 Analysis of embryonic time point specific genes in BELAs

The PrE-like BELA cells were selected from the intgerateion analysis of
the endodermal lineage in order to investigate if the differences
between cells mapping onto E4.5 cells and cells mapping onto E6.5 cells
do resemble differences in developmental progression. Therefore, the
differentially expressed genes between these two groups of cells
(cluster 1 and 2) were identified by the `FindMarker()`function. The 50
most upregulated genes sorted by their log1p-transformed foldchange were
shown in a heatmap. For that cells were clustered by their scaled mean
expression in embryonic cells divided by time point (see annotation
row). Moreover, the scaled single cell expression values are shown from
embryonic cells, annotaed by time point. Epop and Dipk1a were excluded
from the heatmap, since there were no expression values for these genes
in the embyonic dataset.

``` r
# Fig 5 Supp 2 A
# Find DE genes between Agg E4.5 and Agg E6.5 exVE cells
BELA_afin <- subset(integrated_rep1_PrE, orig.ident == "BELAs")
Idents(BELA_afin) <- BELA_afin$seurat_clusters
BELA_afin_markers <- FindMarkers(BELA_afin, ident.1 = "1", ident.2 = "2", assay = "SCT")
BELA_afin_markers$Gene <- rownames(BELA_afin_markers)
BELA_afin_markers <- subset(BELA_afin_markers, Gene!="Epop")
BELA_afin_markers <- subset(BELA_afin_markers, Gene!="Dipk1a") 
up <- rownames(BELA_afin_markers[order(BELA_afin_markers$avg_log2FC), ])[1:50]
down <- rownames(BELA_afin_markers[order(-BELA_afin_markers$avg_log2FC), ])[1:50]

# Make Heatmap of DE genes in Agg in Embryo, comparing timepoints
sc_endo_rep1_PrE <- SCTransform(sc_endo_rep1_PrE, return.only.var.genes = FALSE, verbose = FALSE) # To get SCT for all genes
Idents(sc_endo_rep1_PrE) <- sc_endo_rep1_PrE$Timepoint
heat_frame <- sc_endo_rep1_PrE@assays[["SCT"]]@scale.data[up,]
order_for_heat <- data.frame(sc_endo_rep1_PrE$Timepoint)
colnames(order_for_heat) <- "Timepoint"

# Make a properly scaled heatmap with mean Expression and clustering based on this
heat_frame_order <- log1p(data.frame(AverageExpression(sc_endo_rep1_PrE, assay = "SCT", features = up)))
heat_frame_order <- data.frame(t(scale(t(heat_frame_order))))[,c(3,2,1)]
colnames(heat_frame_order) <- c("Mean(E6.5)", "Mean(E5.5)","Mean(E4.5)")
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...), isReverse = TRUE))
mat_cluster_rows <- sort_hclust(hclust(dist(heat_frame_order)))

ann_colors = list("Mean(E4.5)" = viridis::plasma(100), 
                  "Mean(E5.5)" = viridis::plasma(100), 
                  "Mean(E6.5)" = viridis::plasma(100))

pheatmap(heat_frame, scale = "none", cluster_rows = mat_cluster_rows, border_color = NA,
         annotation_col = order_for_heat, gaps_col = c(60,536), breaks = seq(-2,2,length = 100),
         annotation_row = heat_frame_order, annotation_colors = ann_colors, cluster_cols = FALSE, 
         col=viridis::cividis(100), number_color = "black", cellwidth = 0.05, cellheight = 8, 
         fontsize_number = 14,show_colnames = FALSE)
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

These differentially expressed genes, both up- and downregulated, were
used to investgate the correlation between embryonic timepoints. Pearson
correlation coefficients were calculated on the average expression in
the endoderm lineage of the embryo divided by time point. Results were
shown as a heatmap.

``` r
# Fig 5 Supp 2 B
Idents(object = sc_endo) <- "CellType"
sc_endo_PrE<- subset(sc_endo, idents = c("PrE", "emVE", "exVE"))
sc_endo_PrE <- SCTransform(sc_endo_PrE, verbose = FALSE)
Idents(sc_endo_PrE) <- sc_endo_rep1_PrE$Timepoint
av.exp <- AverageExpression(sc_endo_PrE, features = c(down,up), assays = "SCT")$SCT
cor.exp <- as.data.frame(cor(av.exp))
pheatmap(cor.exp, scale = "none", cluster_rows = FALSE, border_color = NA, display_numbers = TRUE, 
         cluster_cols = FALSE, col=viridis::cividis(200), number_color = "black", cellwidth = 50, 
         cellheight = 50, fontsize_number = 16, show_colnames = TRUE)
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

### 4.4.3 AVE marker gene determination

To determine which genes drive the BELA-PrE cells more towards the E6.5
emVE cells during integration in comparison to the cyst-PrE cells,
marker genes distinguishing the BELA cells which cluster together with
the E6.5 emVE from the BELA cells which cluster together with the E6.5
exVE cells were identified using `FindMarkers()`. The 50 most
upregulated genes were selected for plotting. Changing the
`pseudocount.use` parameter slightly changed the list of genes as
expected for lowly expressed genes. We therefore included a second
version of these heatmaps in Fig 5 Suppl. 3.

``` r
# Fig 5 D (Heatmap for AVE marker genes)
Idents(integrated_rep1_PrE) <- integrated_rep1_PrE$seurat_clusters
integrated_rep1_PrE <- RenameIdents(integrated_rep1_PrE, '1' = 'exVE', '2' = 'emVE', '3' = 'E5.5', '4' = 'E4.5')
integrated_rep1_PrE$seurat_clusters_named <- integrated_rep1_PrE@active.ident

BELAs_PrE <- subset(integrated_rep1_PrE, orig.ident =="BELAs")
BELAs_PrE <- RenameIdents(BELAs_PrE, 'exVE' = 'exVE (BELAs)', 'emVE' = 'emVE (BELAs)')
PrECysts_PrE <- subset(integrated_rep1_PrE, orig.ident =="PrE Cysts")
PrECysts_PrE <- RenameIdents(PrECysts_PrE, 'exVE' = 'exVE (PrE Cysts)', 'emVE' = 'emVE (PrE Cysts)')

DefaultAssay(object = BELAs_PrE) <- "SCT"
AVE.markers <- FindMarkers(BELAs_PrE, ident.1 = "emVE (BELAs)", ident.2 = "exVE (BELAs)", verbose = FALSE) # pseudocount.use = 0, min.diff.pct = 0.15 for Fig Sup 3
up <- rownames(AVE.markers[order(-AVE.markers$avg_log2FC), ])[1:50]
```

Heatmaps were construced showing the scaled and thresholded single cell
expression values for each cluster of cells from BELAs and PrE cysts.
Cells were sorted by hierachical clustering per single cluster. Again,
Dipk1a was excluded, since it was not detected in the original dataset.

``` r
Breaks <- seq(-3,6,length = 100)
# Do pheatmap on BELA emVE with clustering
BELAs_PrE_emVE <- subset(BELAs_PrE, idents = "emVE (BELAs)")
heat_frame <- BELAs_PrE_emVE@assays[["SCT"]]@scale.data[subset(up, up!="Dipk1a"),]
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...), isReverse = TRUE))
mat_cluster_cols <- sort_hclust(hclust(dist(t(heat_frame)))) # For sorting dendrogram
annotation <- data.frame("Cluster" = rep("emVE (BELAs)", times = 139), row.names = colnames(heat_frame)) # to get annotation row
pheatmap(heat_frame, scale = "none", cluster_rows = FALSE, border_color = NA,
         annotation_col = annotation, annotation_colors = list(Cluster = c('emVE (BELAs)' = "#80B25F")), 
         breaks = Breaks, cluster_cols = mat_cluster_cols, col=viridis::cividis(100), number_color = "black", 
         cellwidth = 0.5, cellheight = 10, fontsize_number = 14,show_colnames = FALSE)
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
# Do pheatmap on BELA exVE with clustering
BELAs_PrE_exVE <- subset(BELAs_PrE, idents = "exVE (BELAs)")
heat_frame <- BELAs_PrE_exVE@assays[["SCT"]]@scale.data[subset(up, up!="Dipk1a"),]
mat_cluster_cols <- sort_hclust(hclust(dist(t(heat_frame)))) # For sorting dendrogram
annotation <- data.frame("Cluster" = rep("exVE (BELAs)", times = 283), row.names = colnames(heat_frame)) # to get annotation row
pheatmap(heat_frame, scale = "none", cluster_rows = FALSE, border_color = NA,
         annotation_col = annotation, annotation_colors = list(Cluster = c('exVE (BELAs)' = "#F49089")), 
         breaks = Breaks, cluster_cols = mat_cluster_cols, col=viridis::cividis(100), number_color = "black", 
         cellwidth = 0.5, cellheight = 10, fontsize_number = 14,show_colnames = FALSE)
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-26-2.png)<!-- -->

``` r
# Do pheatmap on PrE Cyst emVE with clustering
PrECysts_PrE_emVE <- subset(PrECysts_PrE, idents = "emVE (PrE Cysts)")
heat_frame <- PrECysts_PrE_emVE@assays[["SCT"]]@scale.data[subset(up, up!="Dipk1a"),]
mat_cluster_cols <- sort_hclust(hclust(dist(t(heat_frame)))) # For sorting dendrogram
annotation <- data.frame("Cluster" = rep("emVE (PrE Cysts)", times = 6), row.names = colnames(heat_frame)) # to get annotation row
pheatmap(heat_frame, scale = "none", cluster_rows = FALSE, border_color = NA,
         annotation_col = annotation, annotation_colors = list(Cluster = c('emVE (PrE Cysts)' = "#086B30")), 
         breaks = Breaks, cluster_cols = mat_cluster_cols, col=viridis::cividis(100), number_color = "black", 
         cellwidth = 0.5, cellheight = 10, fontsize_number = 14,show_colnames = FALSE)
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-26-3.png)<!-- -->

``` r
# Do pheatmap on PrE Cyst exVE with clustering
PrECysts_PrE_exVE <- subset(PrECysts_PrE, idents = "exVE (PrE Cysts)")
heat_frame <- PrECysts_PrE_exVE@assays[["SCT"]]@scale.data[subset(up, up!="Dipk1a"),]
mat_cluster_cols <- sort_hclust(hclust(dist(t(heat_frame)))) # For sorting dendrogram
annotation <- data.frame("Cluster" = rep("exVE (PrE Cysts)", times = 505), row.names = colnames(heat_frame)) # to get annotation row
pheatmap(heat_frame, scale = "none", cluster_rows = FALSE, border_color = NA,
         annotation_col = annotation, annotation_colors = list(Cluster = c('exVE (PrE Cysts)' = "#F05A28")),
         breaks = Breaks, cluster_cols = mat_cluster_cols, col=viridis::cividis(100), number_color = "black", 
         cellwidth = 0.5, cellheight = 10, fontsize_number = 14,show_colnames = FALSE)
```

![](20220223_scRNAseq_analysis_files/figure-gfm/unnamed-chunk-26-4.png)<!-- -->

``` r
sessionInfo()
```

    ## R version 4.0.4 (2021-02-15)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] dendsort_0.3.4        SeuratDisk_0.0.0.9019 data.table_1.14.2    
    ## [4] pheatmap_1.0.12       ggplot2_3.3.5         SeuratObject_4.0.2   
    ## [7] Seurat_4.0.5         
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.15            colorspace_2.0-2      deldir_1.0-6         
    ##   [4] ellipsis_0.3.2        ggridges_0.5.3        rstudioapi_0.13      
    ##   [7] spatstat.data_2.1-0   farver_2.1.0          leiden_0.3.9         
    ##  [10] listenv_0.8.0         ggrepel_0.9.1         bit64_4.0.5          
    ##  [13] RSpectra_0.16-0       fansi_0.5.0           codetools_0.2-18     
    ##  [16] splines_4.0.4         knitr_1.36            polyclip_1.10-0      
    ##  [19] jsonlite_1.7.2        ica_1.0-2             cluster_2.1.0        
    ##  [22] png_0.1-7             uwot_0.1.10           shiny_1.7.1          
    ##  [25] sctransform_0.3.2     spatstat.sparse_2.0-0 compiler_4.0.4       
    ##  [28] httr_1.4.2            assertthat_0.2.1      Matrix_1.3-4         
    ##  [31] fastmap_1.1.0         lazyeval_0.2.2        limma_3.46.0         
    ##  [34] cli_3.0.1             later_1.3.0           htmltools_0.5.2      
    ##  [37] tools_4.0.4           igraph_1.2.7          gtable_0.3.0         
    ##  [40] glue_1.4.2            RANN_2.6.1            reshape2_1.4.4       
    ##  [43] dplyr_1.0.7           Rcpp_1.0.7            scattermore_0.7      
    ##  [46] vctrs_0.3.8           nlme_3.1-152          lmtest_0.9-38        
    ##  [49] xfun_0.27             stringr_1.4.0         globals_0.14.0       
    ##  [52] mime_0.12             miniUI_0.1.1.1        lifecycle_1.0.1      
    ##  [55] irlba_2.3.3           goftest_1.2-3         future_1.22.1        
    ##  [58] MASS_7.3-53           zoo_1.8-9             scales_1.1.1         
    ##  [61] spatstat.core_2.3-0   promises_1.2.0.1      spatstat.utils_2.2-0 
    ##  [64] parallel_4.0.4        RColorBrewer_1.1-2    yaml_2.3.5           
    ##  [67] reticulate_1.22       pbapply_1.5-0         gridExtra_2.3        
    ##  [70] rpart_4.1-15          stringi_1.7.5         highr_0.9            
    ##  [73] rlang_0.4.12          pkgconfig_2.0.3       matrixStats_0.61.0   
    ##  [76] evaluate_0.14         lattice_0.20-41       ROCR_1.0-11          
    ##  [79] purrr_0.3.4           tensor_1.5            labeling_0.4.2       
    ##  [82] patchwork_1.1.1       htmlwidgets_1.5.4     bit_4.0.4            
    ##  [85] cowplot_1.1.1         tidyselect_1.1.1      parallelly_1.28.1    
    ##  [88] RcppAnnoy_0.0.19      plyr_1.8.6            magrittr_2.0.1       
    ##  [91] R6_2.5.1              generics_0.1.1        DBI_1.1.1            
    ##  [94] pillar_1.6.4          withr_2.4.2           mgcv_1.8-33          
    ##  [97] fitdistrplus_1.1-6    survival_3.2-7        abind_1.4-5          
    ## [100] tibble_3.1.5          future.apply_1.8.1    crayon_1.4.1         
    ## [103] hdf5r_1.3.4           KernSmooth_2.23-18    utf8_1.2.2           
    ## [106] spatstat.geom_2.3-0   plotly_4.10.0         rmarkdown_2.11       
    ## [109] viridis_0.6.2         grid_4.0.4            digest_0.6.28        
    ## [112] xtable_1.8-4          tidyr_1.1.4           httpuv_1.6.3         
    ## [115] munsell_0.5.0         viridisLite_0.4.0
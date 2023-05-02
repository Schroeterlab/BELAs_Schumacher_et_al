# BELAs_Schumacher_et_al

This repository contains all code used for the analysis of single cell RNA sequencing data in:

>**Bilayered embryo-like stem cell aggregates reveal an autonomous potential of the primitive endoderm for morphogenesis and differentiation**  
>Sina Schumacher, Max Fernkorn, Michelle Marten, Christian Schröter  
>*[Insert Link/Publication details here]*

This includes the normalization, clustering, the identification of differentially expressed genes and cell-cell-communication analysis based on the dataset of BELAs, Epi cysts and VE cysts as well as the integration with the public embryo datasets.
The single cell sequencing data from the in vitro system described in this study can be found under [GSE198780](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198780).

The analysis of only the in vitro dataset as well as all visualizations were performed in R, the integration and label transfer based on public embryo datasets in python. The structure of the analysis and which files contain which anlysis steps is described here:

![Overview_Analysis_Workflow](https://user-images.githubusercontent.com/88881773/233935939-d5b703ea-1cdb-4540-9dd4-b67338935a35.png)

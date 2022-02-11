### Overview

**scRegulocity** is publicly available from https://github.com/akdess/scRegulocity

This vignette shows the basic steps for running scram.


### Installation

Install latest version from GitHub (requires [devtools](https://github.com/hadley/devtools) package):

```r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("akdess/scRegulocity", dependencies = TRUE, build_vignettes = FALSE)
```

### Input data

1. Generation of the intron exon read count statistics from RNA-Seq BAM files:

The input to scRegulocity consists of bam file. We first use our c++ code, IntrExtract, for extracting intron exon read count statistics from bam file. Intron and exon read count statistics is extracted directly from the mapped RNA-Seq reads using our algorithm called IntrExtract. You can download  IntrExtract C++ source or binary code from [here](https://github.com/harmancilab/IntrExtract).

For extracting Intron and exon read count statistics from RNA-Seq bam files download BAFExtract C++ source or binary code from here.

After downloading source code type the following:

```
cd IntrExtract
make clean
make
```
The executable is located under directory /bin.
IntrExtract output files are read using the following functions:

```r
intex_stats <- readIntrexStats (file.dir="./intrex_stats_ST.txt/", file.name="intrex_stats_ST.txt", emb=emb) 
stats <- normalizeAndSmoothIntrexStats (emat=intex_stats$emat, nmat=intex_stats$nmat, smat=intex_stats$smat)
```
emb: contains the reduced dimension projection (PCA/tSNE/UMAP) coordinates of the cells. 


2. Creating scRegulocity object

The scRegulocity object is required for performing velocity analysis on single-cell data. scRegulocity object is created using following functions


```r
motif_ref <- "hg19-tss-centered-10kb-7species.mc9nr.feather"

scRegulo_obj <- CreateSCReguloCityObject(emat = stats$emat, nmat = stats$nmat, smat =NULL,
 emb=intex_stats$emb,  species="mouse" , motif_ref=motif_ref, cell.types=cell.types)
scRegulo_obj@vel<- scRegulo_obj@vel[which(apply(scRegulo_obj@vel, 1, sd)>1), ]
```

motif_ref: file can be downloaded from [here](https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/)

3. Running scRegulocity object

After creating the scRegulocity object we perform our analysis using the below function: 

```r
kNNdistplot(scRegulo_obj@emb[, 1:2], k = 20)
neigh <- 6 
scRegulo<- applyReguloCity  (object=scRegulo_obj, neigh=neigh, minsize=10)
```

We extract regulotory network using velocity changes using the below function: 
```r
scReguloTF <- extractRegulatoryNetwork (object=scRegulo,  minNumGenesInPattern=5, motifAnnotations=motifAnnotations)
```

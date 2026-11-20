# TB-Gene-Expression-Shiny-App
Interactive Shiny app for TB gene-expression clustering
# TB Gene Expression Clustering App

This Shiny application explores gene-expression patterns in patients with active tuberculosis (TB) and latent tuberculosis infection (LTBI).

## Features

- Choose between PCA, UMAP, and NMF for dimensionality reduction.
- Adjust the number of K-means clusters.
- Visualize patient samples colored by cluster and shaped by TB status.
- Display K-means cluster centers.
- Compare clusters with actual TB status using a contingency table and Adjusted Rand Index (ARI).

## How to Run

1. Download the repository files into the same folder.
2. Open `app.R` in RStudio.
3. Install the required R packages if necessary:

```r
install.packages(c("shiny", "tidyverse", "uwot", "NMF", "mclust"))
```

4. Make sure `TBnanostring.rds` is in the same folder as `app.R`.
5. Click **Run App** in RStudio.

## Dataset

The application uses TB Nanostring gene-expression data containing 179 patient samples and 107 genes.

The dataset must be available locally for the application to run.

## Purpose

This application was developed for an Introduction to Machine Learning assignment to explore dimensionality reduction and unsupervised clustering of gene-expression data.

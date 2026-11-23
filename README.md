## Intro

TODO


## Setup & Installation

This repository uses **Git LFS** (Large File Storage) to store large Seurat `.rds` datasets. 

### 1. Install Git LFS (Required before cloning)
If you do not already have Git LFS installed, install and initialize it once on your system:

  ```bash
  brew install git-lfs
  git lfs install
  ```

### 2. Clone the Repository
Once Git LFS is enabled, clone the repository normally

> **Note:** If you cloned the repository *before* installing Git LFS, or if you encounter `Error in readRDS(): unknown input format`, run `git lfs pull` inside the repository folder to fetch the actual `.rds` files.


### Installing Packages: BiocManager
Some packages require BiocManager to install, use this command to install every required package before running the file. 

    install.packages("BiocManager") 

    BiocManager::install(c(
    "DESeq2",
    "limma",
    "Seurat",
    "ggplot2",
    "ggrepel",
    "ggforce",
    "dplyr",
    "readr",
    "matrixStats"
    ), update = FALSE, ask = FALSE)
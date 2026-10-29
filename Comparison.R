## =============================================================================
## Cross-Dataset Comparison
## =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggplot2)
  library(ggrepel)
  library(limma)          # for removeBatchEffect
  library(ggforce)        # for ellipse in PCA 
})

proj_lab   <- "/Users/avlyn/Desktop/Shariati/RNAseq_data/Cas12a_RNASeq_hESC"
proj_jang  <- "/Users/avlyn/Desktop/Shariati/Jang_Kim_Induction_May_2022"
proj_soncin    <- "/Users/avlyn/Desktop/Shariati/Soncin_Derivation_May_2022"

out_dir_lab  <- file.path(proj_lab,  "02.Results/matrices")
out_dir_jang <- file.path(proj_jang, "02.Results/matrices")
out_dir_soncin   <- file.path(proj_soncin,   "02.Results/matrices")

# Where to save outputs from this comparison script
compare_out <- file.path("/Users/avlyn/Desktop/Shariati/Comparison")
dir.create(compare_out, recursive = TRUE, showWarnings = FALSE)

# Marker gene panels (reused from your existing scripts)
Hypoblast   <- c("GATA6", "SOX17", "PDGFRA", "FGFR2", "LAMA1")
Epiblast    <- c("POU5F1", "SOX2", "NANOG", "KLF17", "PRDM14", "DPPA4", "UTF1", "ZSCAN10", "TDGF1")
Amnion      <- c("ISL1", "GABRP",  "HEY1", "HLA-A", "DIRAS2", "TMEM54", "DLX5")
Trophoblast <- c( "LRP2", "VGLL1", "TACSTD2", "TEAD3", "TEAD4", "CDX2", "TP63", "DLX3")
TB_Amnion <- c("TFAP2A", "GATA3", "TFAP2C", "BMP4", "WNT6", "KRT7", "VTCN1", "DLX5", "HAND1", "GRHL1", "GRHL3", "ENPEP")
Cytotrophoblast <- c("GCM1", "TP63", "HSD3B1", "HSD11B2", "PEG10", "ITGA5", "CDH1")
Syncytiotrophoblast <- c("GATA2",  "XAGE2", "EGFR", "PAGE4", "SPINT1", "CGB1", "NCF4", "TCAF2C", "CGB2")
EVT <- c("HLA-G", "PDE6H", "TNNI2", "IFIH1", "PADI1", "PDCD1LG2")

## =============================================================================
## 2) Joint PCA across all three datasets with batch correction
## =============================================================================
## Merges raw counts from all three projects (on shared genes), builds one
## DESeqDataSet with a `dataset` covariate, VST-transforms, then removes the
## dataset-level batch effect before PCA — so you can see whether biologically
## similar groups from different datasets land near each other once technical
## (sequencing/batch) variation is regressed out.

build_joint_pca <- function(counts_paths, meta_list, dataset_names,
                            out_png_prefix, group_colors = NULL,
                            ntop = 1000) {
  
  stopifnot(length(counts_paths) == length(meta_list),
            length(counts_paths) == length(dataset_names))
  
  # --- Load and align raw count matrices on shared genes ---
  count_mats <- purrr::map(counts_paths, function(p) {
    df <- readr::read_csv(p, show_col_types = FALSE)
    mat <- as.matrix(df[,-1])
    rownames(mat) <- sub("\\..*$", "", df[[1]])   # strip Ensembl version
    mat
  })
  
  common_genes <- Reduce(intersect, lapply(count_mats, rownames))
  message("Shared genes across all datasets: ", length(common_genes))
  if (length(common_genes) < 5000) {
    warning("Fewer than 5000 shared genes — check that all datasets used the same GTF/annotation build.")
  }
  
  count_mats <- lapply(count_mats, function(m) m[common_genes, , drop = FALSE])
  
  # Keep only groups that have a color assigned
  if (!is.null(group_colors)) {
    meta_list <- lapply(meta_list, function(m) dplyr::filter(m, group %in% names(group_colors)))
  }
  count_mats <- purrr::map2(count_mats, meta_list, function(mat, m) {
    missing <- setdiff(m$sample, colnames(mat))
    if (length(missing)) stop("Samples not in counts: ", paste(missing, collapse = ", "))
    mat[, m$sample, drop = FALSE]
  })
  
  merged_counts <- do.call(cbind, count_mats)
  merged_counts <- round(merged_counts, 0)
  storage.mode(merged_counts) <- "integer"
  
  # --- Build combined coldata ---
  coldata_list <- purrr::map2(meta_list, dataset_names, function(m, dname) {
    m %>% mutate(dataset = dname)
  })
  coldata <- bind_rows(coldata_list) %>%
    column_to_rownames("sample")
  
  stopifnot(identical(colnames(merged_counts), rownames(coldata)))
  
  coldata$dataset <- factor(coldata$dataset)
  coldata$group <- factor(coldata$group, levels = names(group_colors))
  
  # --- DESeq2 + VST (design ~ dataset here; not testing DE across datasets,
  #     just need a valid design for size-factor/dispersion estimation) ---
  dds_joint <- DESeqDataSetFromMatrix(countData = merged_counts, colData = coldata, design = ~ dataset)
  
  keep <- rowSums(counts(dds_joint) >= 10) >= 3
  dds_joint <- dds_joint[keep, ]
  
  vsd_joint <- tryCatch(
    vst(dds_joint, blind = TRUE),
    error = function(e) varianceStabilizingTransformation(dds_joint, blind = TRUE)
  )
  
  # --- Remove dataset-level batch effect ---
  mat_corrected <- limma::removeBatchEffect(assay(vsd_joint), batch = colData(vsd_joint)$dataset)
  
  # --- PCA on top-variable genes, before and after correction, for comparison ---
  run_pca <- function(mat, coldata, tag) {
    rv <- matrixStats::rowVars(mat)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(mat[select, ]))
    pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2))
    
    df <- as.data.frame(pca$x[, 1:2])
    df$sample  <- rownames(df)
    df$group   <- coldata$group
    df$dataset <- coldata$dataset
    
    df <- as.data.frame(pca$x[, 1:2])
    df$sample  <- rownames(df)
    df$group   <- coldata$group
    df$dataset <- coldata$dataset
    
    # Fixed shape per dataset, and the matching shape for each group's legend key
    shape_vals <- setNames(c(16, 17, 15, 18)[seq_along(levels(df$dataset))],
                           levels(df$dataset))
    grp_levels <- levels(droplevels(df$group))
    grp_shapes <- unname(shape_vals[as.character(df$dataset[match(grp_levels, df$group)])])
    
    p <- ggplot(df, aes(PC1, PC2)) +
      ggforce::geom_mark_ellipse(
        aes(group = group),     # still one ellipse per group, but not colored by it
        colour = "red",
        fill = NA,              # no shading
        linewidth = 0.25,       # thinner line (was 0.4)
        expand = unit(2, "mm"),
        show.legend = FALSE
      ) +
      geom_point(aes(colour = group, shape = dataset), size = 3, alpha = 0.9) +
      { if (!is.null(group_colors)) scale_colour_manual(values = group_colors) } +
      { if (!is.null(group_colors)) scale_fill_manual(values = group_colors) } +
      xlab(paste0("PC1: ", pct[1], "% variance")) +
      ylab(paste0("PC2: ", pct[2], "% variance")) +
      labs(title = tag, colour = NULL, shape = "Dataset") +
      scale_shape_manual(values = shape_vals) +
      guides(colour = guide_legend(override.aes = list(shape = grp_shapes, alpha = 1))) +
      theme_classic()
    
    list(plot = p, data = df)
  }
  
  pca_before <- run_pca(assay(vsd_joint), coldata, "Before batch correction (dataset = shape)")
  pca_after  <- run_pca(mat_corrected,    coldata, "After batch correction (dataset = shape)")
  
  ggsave(paste0(out_png_prefix, "_before_correction.png"), pca_before$plot, width = 9, height = 6, dpi = 300)
  ggsave(paste0(out_png_prefix, "_after_correction.png"),  pca_after$plot,  width = 9, height = 6, dpi = 300)
  
  readr::write_csv(pca_after$data, paste0(out_png_prefix, "_after_correction_data.csv"))
  
  message("Saved: ", paste0(out_png_prefix, "_before_correction.png"))
  message("Saved: ", paste0(out_png_prefix, "_after_correction.png"))
  
  invisible(list(dds = dds_joint, vsd = vsd_joint, corrected_matrix = mat_corrected,
                 pca_before = pca_before, pca_after = pca_after))
}

meta_lab  <- readr::read_csv(file.path(proj_lab,  "metadata.csv"), show_col_types = FALSE)
meta_jang <- readr::read_csv(file.path(proj_jang, "metadata.csv"), show_col_types = FALSE)
meta_soncin   <- readr::read_csv(file.path(proj_soncin,   "metadata.csv"), show_col_types = FALSE)

joint_group_colors <- c("grhl1_3" = "#9043B7", "gata2_3 tfap2a_c" = "#C180D2", 
                        "Control_gRNA" = "#D3D3D3", "No_gRNA" = "#9E9E9E", "H9_hESC_d0" = "#696969", # controls
                        "ESTSLC"="#A0D2E7", "H9_TSC_p6"  = "#81B1D5",
                        "TSC"="#3D60A7", "TSC_1048"   = "#26408B", "TSC_1049"   = "#0F084B",
                        "CTB_1stE"   = "#C5DE83",
                        "ST"="#F7B0DB", "ESST"="#FACFEC", 
                        "EVT"="#23BB98", "ESEVT"="#82DFC6",
                        "H9_BMP4_d4" = "#B52471"
)
build_joint_pca(
  counts_paths   = c(file.path(out_dir_lab, "gene_counts_raw.csv"),
                    file.path(out_dir_jang, "gene_counts_raw.csv"),
                    file.path(out_dir_soncin,   "gene_counts_raw.csv")),
  meta_list      = list(meta_lab, meta_jang, meta_soncin),
  dataset_names  = c("Lab", "Jang", "Soncin"),
  out_png_prefix = file.path(compare_out, "joint_PCA"),
  group_colors   = joint_group_colors
)


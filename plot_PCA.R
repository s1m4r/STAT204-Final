# ==============================================================================
# plot_PCA.R Output Files
#   1) PCA_Shariati.png           - Shariati lab samples only
#   2) PCA_Shariati_Jang_Soncin.png - joint PCA, dataset batch removed (limma)
#   3) PCA_Shariati_Zhao.png      - Shariati samples projected onto the Zhao
#                                   embryo reference (pseudobulk) PCA
#
# Tree
# ├── gencode.v49.primary_assembly.annotation.gtf
# ├── Jang
# │   ├── gene_counts_raw.csv
# │   ├── gene_map_human_gencode.csv
# │   └── metadata.csv
# ├── plot_PCA.R
# ├── Shariati
# │   ├── gene_counts_raw.csv
# │   ├── gene_map_human_gencode.csv
# │   └── metadata.csv
# ├── Soncin
# │   ├── gene_counts_raw.csv
# │   ├── gene_map_human_gencode.csv
# │   └── metadata.csv
# └── Zhao
# └── psd.R3.6.em.seurat.ob.rds

suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(ggrepel); library(dplyr)
})

root_dir <- "."
out_dir  <- file.path(root_dir, "PCA_outputs")
dir.create(out_dir, showWarnings = FALSE)

## ---- Colors -----------------------------------------------------------------
lab_colors <- c("grhl1_3" = "#9043B7", "gata2_3 tfap2a_c" = "#C180D2",
                "Control_gRNA" = "#D3D3D3", "No_gRNA" = "#9E9E9E")

joint_colors <- c(lab_colors, "H9_hESC_d0" = "#696969",
                  "ESTSLC" = "#A0D2E7", "H9_TSC_p6" = "#81B1D5",
                  "TSC" = "#3D60A7", "TSC_1048" = "#26408B", "TSC_1049" = "#0F084B",
                  "CTB_1stE" = "#C5DE83", "ST" = "#F7B0DB", "ESST" = "#FACFEC",
                  "EVT" = "#23BB98", "ESEVT" = "#82DFC6", "H9_BMP4_d4" = "#B52471")

lineage_colors <- c("Zygote"= "#999999", "2-4 cell"= "#f39b7f", "8 cell"= "#f8766d",
                    "Morula"= "#ff7f0e", "Prelineage"="#c49c94", "ICM"= "#bc81ff",
                    "Epiblast"= "#00b813", "PriS"= "#f863df", "Mesoderm"= "#00c1a7",
                    "AdvMes"= "#e9842c", "Axial Mes"= "#094334", "Amnion"= "#e71f19",
                    "ExE_Mes"= "#0085ee", "Hypoblast"= "#7c96ff", "DE"= "#a4c9dd",
                    "YSE"= "#ea628e", "TE"= "#00b5ee", "CTB"= "#9ca700",
                    "STB"= "#f7b6d2", "EVT"= "#54b986", "Erythroblasts"= "#393b79",
                    "HEP"= "#00ebed", "EPI.PrE.INT"= "#8ccc6e")

## ---- Helpers ----------------------------------------------------------------
# Raw counts (first column = Ensembl ID) + metadata (sample, group), keeping colored groups only
load_dataset <- function(name, colors) {
  df   <- readr::read_csv(file.path(root_dir, name, "gene_counts_raw.csv"), show_col_types = FALSE)
  meta <- readr::read_csv(file.path(root_dir, name, "metadata.csv"), show_col_types = FALSE)
  meta <- meta[meta$group %in% names(colors), ]
  counts <- round(as.matrix(df[, meta$sample]))
  rownames(counts) <- sub("\\..*$", "", df[[1]])  # strip Ensembl version
  list(counts = counts, meta = data.frame(sample = meta$sample, group = meta$group, dataset = name))
}

run_vst <- function(counts, coldata) {
  storage.mode(counts) <- "integer"
  dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ 1)
  assay(vst(dds[rowSums(counts(dds) >= 10) >= 3, ], blind = TRUE))
}

pca_top <- function(mat, ntop) {
  sel <- order(matrixStats::rowVars(mat), decreasing = TRUE)[1:min(ntop, nrow(mat))]
  pca <- prcomp(t(mat[sel, ]))
  pca$pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2))
  pca$genes <- rownames(mat)[sel]
  pca
}

axis_labs <- function(pct) labs(x = paste0("PC1: ", pct[1], "% variance"),
                                y = paste0("PC2: ", pct[2], "% variance"))

shariati <- load_dataset("Shariati", lab_colors)


# ==============================================================================
# 1) Shariati only
# ==============================================================================
pca1 <- pca_top(run_vst(shariati$counts, shariati$meta), 500)
df1  <- data.frame(pca1$x[, 1:2], shariati$meta)

p1 <- ggplot(df1, aes(PC1, PC2, colour = group, label = sample)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, colour = "grey20", box.padding = 0.6, segment.color = "grey50",
                  max.overlaps = Inf) +
  scale_colour_manual(values = lab_colors, breaks = names(lab_colors)) +
  axis_labs(pca1$pct) + labs(title = "Shariati lab", colour = NULL) +
  theme_classic()
ggsave(file.path(out_dir, "PCA_Shariati.png"), p1, width = 7, height = 5, dpi = 300)


# ==============================================================================
# 2) Shariati + Jang + Soncin, dataset effect removed
# ==============================================================================
sets   <- list(shariati, load_dataset("Jang", joint_colors), load_dataset("Soncin", joint_colors))
common <- Reduce(intersect, lapply(sets, function(d) rownames(d$counts)))
counts2 <- do.call(cbind, lapply(sets, function(d) d$counts[common, ]))
meta2   <- do.call(rbind, lapply(sets, `[[`, "meta"))

mat2 <- limma::removeBatchEffect(run_vst(counts2, meta2), batch = meta2$dataset)
pca2 <- pca_top(mat2, 1000)
df2  <- data.frame(pca2$x[, 1:2], meta2)
df2$group <- factor(df2$group, levels = intersect(names(joint_colors), df2$group))

shapes <- c(Jang = 16, Shariati = 17, Soncin = 15)
p2 <- ggplot(df2, aes(PC1, PC2)) +
  ggforce::geom_mark_ellipse(aes(group = group), colour = "red", linewidth = 0.25,
                             expand = unit(2, "mm")) +
  geom_point(aes(colour = group, shape = dataset), size = 3, alpha = 0.9) +
  scale_colour_manual(values = joint_colors) +
  scale_shape_manual(values = shapes) +
  guides(colour = guide_legend(override.aes = list(
    shape = unname(shapes[df2$dataset[match(levels(df2$group), df2$group)]])))) +
  axis_labs(pca2$pct) +
  labs(title = "Shariati + Jang + Soncin (dataset effect removed)", colour = NULL, shape = "Dataset") +
  theme_classic()
ggsave(file.path(out_dir, "PCA_Shariati_Jang_Soncin.png"), p2, width = 9, height = 6, dpi = 300)


# ==============================================================================
# 3) Shariati projected onto the Zhao embryo reference
# ==============================================================================
# Pseudobulk = summed counts per lineage x source study. Slow, so cached after run 1.
pb_cache <- file.path(root_dir, "Zhao", "zhao_pseudobulk.rds")
if (!file.exists(pb_cache)) {
  library(Seurat)
  ob <- UpdateSeuratObject(readRDS(file.path(root_dir, "Zhao", "psd.R3.6.em.seurat.ob.rds")))
  study <- intersect(c("dataset", "orig.ident", "sample", "batch"), colnames(ob@meta.data))[1]
  ob <- subset(ob, cells = colnames(ob)[ob$rename_EML %in% names(lineage_colors)])
  # AggregateExpression mangles "_" in names (ExE_Mes), so group on a safe ID
  key <- paste(ob$rename_EML, ob@meta.data[[study]], sep = "||")
  ob$pb_id <- paste0("pb", as.integer(factor(key)))
  pb <- AggregateExpression(ob, assays = "RNA", group.by = "pb_id")$RNA
  colnames(pb) <- sub("^g", "", colnames(pb))
  pb_meta <- data.frame(pb_id = colnames(pb),
                        lineage = ob$rename_EML[match(colnames(pb), ob$pb_id)],
                        study   = ob@meta.data[[study]][match(colnames(pb), ob$pb_id)],
                        n_cells = as.vector(table(ob$pb_id)[colnames(pb)]))
  saveRDS(list(pb = pb, pb_meta = pb_meta), pb_cache)
}
zb <- readRDS(pb_cache)
ref_meta   <- subset(zb$pb_meta, n_cells >= 10) # drops Zygote (too few cells)
ref_counts <- as.matrix(zb$pb[, ref_meta$pb_id])

# Shariati: Ensembl ID -> gene symbol, to match the reference
gmap <- readr::read_csv(file.path(root_dir, "Shariati", "gene_map_human_gencode.csv"), show_col_types = FALSE)
sym  <- gmap$gene_name[match(rownames(shariati$counts), sub("\\..*$", "", gmap$gene_id))]
lab_counts <- rowsum(shariati$counts[!is.na(sym), ], sym[!is.na(sym)])

common3 <- intersect(rownames(ref_counts), rownames(lab_counts))
is_ref  <- c(rep(TRUE, ncol(ref_counts)), rep(FALSE, ncol(lab_counts)))
vmat <- run_vst(round(cbind(ref_counts[common3, ], lab_counts[common3, ])),
                data.frame(row.names = c(ref_meta$pb_id, colnames(lab_counts)), is_ref))

# PCA fit on reference only, Shariati samples projected onto it 
# No batch correction since reference and Shariati don't share groups
pca3 <- pca_top(vmat[, is_ref], 500)
rdf  <- data.frame(pca3$x[, 1:2], group = factor(ref_meta$lineage, levels = names(lineage_colors)))
ldf  <- data.frame(predict(pca3, t(vmat[pca3$genes, !is_ref]))[, 1:2], shariati$meta)
ldf$group <- factor(ldf$group, levels = names(lab_colors))
centers <- aggregate(cbind(PC1, PC2) ~ group, rdf, median)

p3 <- ggplot(mapping = aes(PC1, PC2)) +
  geom_point(data = rdf, aes(colour = group), size = 2, alpha = 0.6) +
  geom_text_repel(data = centers, aes(label = group), size = 3, colour = "grey35",
                  fontface = "italic", max.overlaps = Inf, seed = 1) +
  geom_point(data = ldf, aes(fill = group), shape = 22, colour = "transparent", size = 4) +
  geom_text_repel(data = ldf, aes(label = sample), size = 3, box.padding = 0.6,
                  segment.color = "grey50", max.overlaps = Inf, seed = 1) +
  scale_colour_manual(values = lineage_colors, name = "Reference lineage", drop = TRUE) +
  scale_fill_manual(values = lab_colors, name = "Lab group") +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1))) +
  axis_labs(pca3$pct) +
  labs(title = "Shariati samples projected onto Zhao embryo reference",
       subtitle = "PCs from reference only; % variance refers to the reference") +
  theme_classic()
ggsave(file.path(out_dir, "PCA_Shariati_Zhao.png"), p3, width = 10, height = 7, dpi = 300)

message("Saved 3 plots to ", out_dir)
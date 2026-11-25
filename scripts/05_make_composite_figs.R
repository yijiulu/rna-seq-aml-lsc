#!/usr/bin/env Rscript
# Assemble panel-tagged Fig 1 (A-G: QC, platform stratification, lineage markers & CQN diagnostic)
# and Fig S1 (A-B: library QC & Euclidean sample distance) from saved Module 1 objects.

suppressPackageStartupMessages({
  library(DESeq2)
  library(cqn)
  library(splines)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(matrixStats)
  library(RColorBrewer)
})

cmd_args   <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", grep("^--file=", cmd_args, value = TRUE))
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(script_arg[1])) else file.path(getwd(), "scripts")
source(file.path(dirname(script_dir), "R", "utils.R"))

ROOT      <- find_root(file.path("data", GSE))
DATA_DIR  <- file.path(ROOT, "data", GSE)
QC_DIR    <- file.path(ROOT, "results", "qc")
CELL_COLS <- get_cell_palette()

vsd            <- readRDS(file.path(QC_DIR, "vsd_blind.rds"))
dds            <- readRDS(file.path(QC_DIR, "dds_filtered.rds"))
qc             <- fread(file.path(QC_DIR, "qc_metrics.tsv"))
gc_dt          <- fread(file.path(DATA_DIR, "gene_gc_length.tsv"))
samp_gc_diag   <- fread(file.path(QC_DIR, "sample_gc_bias_diagnostic.tsv"))
tech_axis_diag <- fread(file.path(QC_DIR, "technical_axis_gc_r2.tsv"))
setkey(gc_dt, gene)

vmat <- assay(vsd)
cd   <- as.data.frame(colData(vsd))
cd$cell_type <- factor(as.character(cd$cell_type), levels = CELL_LEVELS)
sample_order <- rownames(cd)[order(cd$cell_type, cd$donor)]
top_var      <- head(order(rowVars(vmat), decreasing = TRUE), N_TOP_VAR)
theme_set(theme_bw(base_size = 10))

pca_plot <- function(samples, title, colour = "cell_type", show_legend = TRUE) {
  m   <- vmat[, samples]
  sel <- head(order(rowVars(m), decreasing = TRUE), N_TOP_VAR)
  pc  <- prcomp(t(m[sel, ]))
  ve  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  df  <- cbind(cd[samples, ], pc$x[, 1:2])
  g <- ggplot(df, aes(PC1, PC2, colour = .data[[colour]], shape = platform)) +
    geom_point(size = 2.3, alpha = 0.9) +
    labs(title = title, x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
         colour = if (colour == "cell_type") "Cell type" else "Platform", shape = "Platform") +
    theme(plot.title = element_text(face = "bold", size = 10.5))
  if (colour == "cell_type") {
    g <- g + scale_colour_manual(values = CELL_COLS, drop = FALSE) +
      guides(colour = guide_legend(ncol = 4, override.aes = list(size = 2.2)), shape = guide_legend(ncol = 2))
  }
  if (!show_legend) g <- g + theme(legend.position = "none")
  g
}

s_all  <- rownames(cd)
s_hspc <- rownames(cd)[cd$compartment != "Mature"]
s_plat <- rownames(cd)[cd$cell_type %in% c("pHSC", "Blast")]

p_pca_all  <- pca_plot(s_all, "All samples (n = 81)", show_legend = TRUE) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8), legend.key.size = unit(3.5, "mm"))
p_pca_hspc <- pca_plot(s_hspc, "HSPC + AML (n = 58)", show_legend = FALSE)
p_pca_plat <- pca_plot(s_plat, "pHSC + Blast: platform check (n = 24)", colour = "platform") +
  aes(shape = cell_type) + labs(colour = "Platform", shape = "Cell type") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8), legend.key.size = unit(3.5, "mm"))

mk <- intersect(names(MARKERS), rownames(vmat))
z  <- t(scale(t(vmat[mk, sample_order])))
zl <- melt(as.data.table(z, keep.rownames = "gene"), id.vars = "gene", variable.name = "sample_id", value.name = "z")
zl[, `:=`(gene = factor(gene, levels = rev(mk)), sample_id = factor(sample_id, levels = sample_order))]
zl[, cell_type := cd[as.character(sample_id), "cell_type"]]

p_mk <- ggplot(zl, aes(sample_id, gene, fill = pmax(pmin(z, 2.5), -2.5))) +
  geom_tile() + facet_grid(~ cell_type, scales = "free_x", space = "free_x") +
  scale_fill_distiller(palette = "RdBu", limits = c(-2.5, 2.5)) +
  scale_y_discrete(labels = function(g) sprintf("%s (%s)", g, MARKERS[g])) +
  labs(title = "Canonical haematopoietic lineage markers (all 81 samples, blind VST)", x = NULL, y = NULL, fill = "z (VST)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 7.5),
        strip.text = element_text(size = 7, angle = 90), panel.spacing = unit(0.6, "mm"),
        plot.title = element_text(face = "bold", size = 10.5))

ns_dds <- dds[, dds$platform == "NextSeq500"]
ns_dds$platform <- droplevels(ns_dds$platform)
common_gc    <- intersect(rownames(ns_dds), gc_dt$gene)
ns_gc_diag   <- estimateSizeFactors(ns_dds[common_gc, ])
raw_norm_ns  <- log2(counts(ns_gc_diag, normalized = TRUE) + 1)
cq_ns        <- cqn(counts(ns_gc_diag), x = gc_dt[common_gc, gc_median], lengths = gc_dt[common_gc, length_median],
                    sizeFactors = colSums(counts(ns_gc_diag)), verbose = FALSE)
cqn_norm_ns  <- cq_ns$y + cq_ns$offset
diag_genes   <- common_gc[rowSums(counts(ns_gc_diag) >= MIN_COUNT) >= 6]
gc_d         <- gc_dt[diag_genes, gc_median]
hsc_others   <- c("5852-HSC", "6792-HSC", "7653-HSC")

df_ax <- data.table(
  gc        = gc_d,
  hsc_pre   = raw_norm_ns[diag_genes, "7256-HSC"] - rowMeans(raw_norm_ns[diag_genes, hsc_others]),
  lmpp_pre  = rowMeans(raw_norm_ns[diag_genes, c("5483-LMPP", "7256-LMPP")]) - raw_norm_ns[diag_genes, "7653-LMPP"],
  hsc_post  = cqn_norm_ns[diag_genes, "7256-HSC"] - rowMeans(cqn_norm_ns[diag_genes, hsc_others]),
  lmpp_post = rowMeans(cqn_norm_ns[diag_genes, c("5483-LMPP", "7256-LMPP")]) - cqn_norm_ns[diag_genes, "7653-LMPP"]
)

focus_cts <- c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "CLP", "Mono", "pHSC", "LSC", "Blast")
p_samp_gc <- ggplot(
  melt(samp_gc_diag[cell_type %in% focus_cts], id.vars = c("sample_id", "cell_type"), measure.vars = c("r_gc_pre", "r_gc_post")),
  aes(factor(cell_type, levels = focus_cts), value, colour = variable)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(position = position_dodge(width = 0.5), size = 2.3, alpha = 0.85) +
  scale_colour_manual(values = c(r_gc_pre = "#D73027", r_gc_post = "#4575B4"),
                      labels = c(r_gc_pre = "Before CQN (DESeq2 SF)", r_gc_post = "After CQN")) +
  labs(title = "Per-sample deviation from cohort median vs transcript GC%",
       subtitle = sprintf("NextSeq 500 HSPC + Mono + AML samples (n = 52; %d expressed genes)", length(diag_genes)),
       x = NULL, y = "Pearson r (sample deviation vs GC%)", colour = NULL) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title = element_text(face = "bold", size = 10.5))

p_ax_pre <- ggplot(df_ax, aes(hsc_pre, lmpp_pre, colour = gc)) +
  geom_point(size = 0.6, alpha = 0.5) + scale_colour_viridis_c(option = "plasma", name = "GC%") +
  geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") +
  labs(title = sprintf("Before CQN: within-HSC vs within-LMPP (r = %+.3f)", cor(df_ax$hsc_pre, df_ax$lmpp_pre)),
       subtitle = sprintf("GC spline R2 = %.2f (HSC), %.2f (LMPP)",
                          tech_axis_diag[axis == "within_HSC_pre", R2_gc_spline],
                          tech_axis_diag[axis == "within_LMPP_pre", R2_gc_spline]),
       x = "7256-HSC vs other 3 HSCs (log2)", y = "5483/7256-LMPP vs 7653-LMPP (log2)") +
  theme(plot.title = element_text(face = "bold", size = 10))

p_ax_post <- ggplot(df_ax, aes(hsc_post, lmpp_post, colour = gc)) +
  geom_point(size = 0.6, alpha = 0.5) + scale_colour_viridis_c(option = "plasma", name = "GC%") +
  geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") +
  labs(title = sprintf("After CQN: within-HSC vs within-LMPP (r = %+.3f)", cor(df_ax$hsc_post, df_ax$lmpp_post)),
       subtitle = sprintf("GC spline R2 = %.3f (HSC), %.3f (LMPP)",
                          tech_axis_diag[axis == "within_HSC_post", R2_gc_spline],
                          tech_axis_diag[axis == "within_LMPP_post", R2_gc_spline]),
       x = "7256-HSC vs other 3 HSCs (CQN log2)", y = "5483/7256-LMPP vs 7653-LMPP (CQN log2)") +
  theme(plot.title = element_text(face = "bold", size = 10))

fig1 <- ((p_pca_all | p_pca_hspc | p_pca_plat) + plot_layout(widths = c(1.05, 1, 1))) /
  p_mk / p_samp_gc / ((p_ax_pre | p_ax_post) + plot_layout(guides = "collect")) +
  plot_layout(heights = c(1.05, 1.05, 0.85, 0.95)) +
  plot_annotation(
    title = "Fig 1. Quality control, platform stratification, lineage markers, and CQN GC-bias correction (GSE74246)",
    tag_levels = "A", theme = theme(plot.title = element_text(face = "bold", size = 12.5))
  )
ggsave(file.path(QC_DIR, "Fig1_qc_platform_cqn.png"), fig1, width = 16.5, height = 19.5, dpi = 180)

qc[, cell_type := factor(cell_type, levels = CELL_LEVELS)]
qc_long <- melt(qc, id.vars = c("sample_id", "cell_type"), measure.vars = c("lib_size", "detected", "pct_ribo", "pct_top100"))
qc_long[, sample_id := factor(sample_id, levels = sample_order)]
p_lib <- ggplot(qc_long, aes(sample_id, value, fill = cell_type)) +
  geom_col() +
  facet_wrap(~ variable, ncol = 1, scales = "free_y",
             labeller = as_labeller(c(lib_size = "Library size (reads)", detected = "Genes detected (>0)",
                                      pct_ribo = "% ribosomal protein reads", pct_top100 = "% reads in top 100 genes"))) +
  scale_fill_manual(values = CELL_COLS) +
  labs(title = "Per-sample sequencing QC metrics across all 81 libraries", x = NULL, y = NULL, fill = "Cell type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5), plot.title = element_text(face = "bold", size = 11))

d  <- dist(t(vmat[top_var, ]))
hc <- hclust(d, method = "average")
dm <- melt(as.data.table(as.matrix(d), keep.rownames = "s1"), id.vars = "s1", variable.name = "s2", value.name = "distance")
dm[, `:=`(s1 = factor(s1, levels = hc$labels[hc$order]), s2 = factor(s2, levels = hc$labels[hc$order]))]
p_dist <- ggplot(dm, aes(s1, s2, fill = distance)) +
  geom_tile() + scale_fill_distiller(palette = "Blues", direction = -1) +
  labs(title = "Hierarchically clustered sample-to-sample Euclidean distance (top 500 variable genes)",
       x = NULL, y = NULL, fill = "Euclidean\ndistance") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), axis.text.y = element_text(size = 5),
        plot.title = element_text(face = "bold", size = 11))

figS1 <- p_lib / p_dist +
  plot_layout(heights = c(0.92, 1.08)) +
  plot_annotation(
    title = "Fig S1. Sequencing quality control and global transcriptome distance structure (GSE74246, n = 81)",
    tag_levels = "A", theme = theme(plot.title = element_text(face = "bold", size = 12.5))
  )
ggsave(file.path(QC_DIR, "FigS1_qc_metrics_distance.png"), figS1, width = 13.5, height = 18.5, dpi = 180)
cat("Generated Fig1_qc_platform_cqn.png and FigS1_qc_metrics_distance.png successfully.\n")

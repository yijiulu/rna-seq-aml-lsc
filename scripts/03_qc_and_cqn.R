#!/usr/bin/env Rscript
# Module 1: GSE74246 data import, sequencing QC, platform stratification, lineage markers,
# and NextSeq 500 GC/length bias diagnostic (pre- vs post-CQN).
#
# Input : data/GSE74246/{counts.tsv.gz, samples.tsv, gene_gc_length.tsv}
# Output: results/qc/{qc_metrics.tsv, sample_gc_bias_diagnostic.tsv, technical_axis_gc_r2.tsv,
#                     dds_filtered.rds, vsd_blind.rds, 01..05 QC figures, sessionInfo.txt}

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

ROOT     <- find_root(file.path("data", GSE))
DATA_DIR <- file.path(ROOT, "data", GSE)
QC_DIR   <- file.path(ROOT, "results", "qc")
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
CELL_COLS <- get_cell_palette()

# ---- 1. Import ----------------------------------------------------------------

counts_path  <- file.path(DATA_DIR, "counts.tsv.gz")
samples_path <- file.path(DATA_DIR, "samples.tsv")
gc_path      <- file.path(DATA_DIR, "gene_gc_length.tsv")
if (!file.exists(counts_path) || !file.exists(samples_path)) {
  stop("Missing processed files in ", DATA_DIR, ". Run `python scripts/01_fetch_geo.py ", GSE, "` first.", call. = FALSE)
}
if (!file.exists(gc_path)) {
  stop("Missing ", gc_path, ". Run `python scripts/02_build_gene_gc.py ", GSE, "` first.", call. = FALSE)
}

counts_raw <- read_counts(counts_path)
meta       <- fread(samples_path, data.table = FALSE)
gc_dt      <- fread(gc_path)
gc_dt      <- gc_dt[is.finite(gc_median) & is.finite(length_median) & length_median > 0]
setkey(gc_dt, gene)

# ---- 2. Validate sample alignment and build colData ---------------------------

if (anyDuplicated(meta$sample_id)) stop("Duplicate sample_id in samples.tsv", call. = FALSE)
if (!identical(colnames(counts_raw), meta$sample_id)) {
  missing_meta   <- setdiff(colnames(counts_raw), meta$sample_id)
  missing_counts <- setdiff(meta$sample_id, colnames(counts_raw))
  if (length(missing_meta) || length(missing_counts)) {
    stop(
      "Unmatched samples.\n  In counts only: ", toString(missing_meta),
      "\n  In samples.tsv only: ", toString(missing_counts),
      call. = FALSE
    )
  }
  meta <- meta[match(colnames(counts_raw), meta$sample_id), ]
}

coldata <- data.frame(
  sample_id   = meta$sample_id,
  gsm         = meta$geo_accession,
  cell_type   = factor(meta$source_name_ch1, levels = CELL_LEVELS),
  donor       = factor(meta$donorid),
  platform    = factor(unname(PLATFORMS[meta$platform_id]), levels = PLATFORMS),
  status      = factor(ifelse(grepl("leukemia", meta$sample_type), "AML", "Normal"), levels = c("Normal", "AML")),
  orig_column = if ("count_column_original" %in% names(meta)) meta$count_column_original else meta$sample_id,
  row.names   = meta$sample_id
)
coldata$compartment <- factor(
  ifelse(coldata$status == "AML", "AML", ifelse(coldata$cell_type %in% MATURE, "Mature", "HSPC")),
  levels = c("HSPC", "Mature", "AML")
)

if (anyNA(coldata$cell_type)) stop("Unexpected cell types: ", toString(unique(meta$source_name_ch1[is.na(coldata$cell_type)])), call. = FALSE)
if (anyNA(coldata$platform))  stop("Unexpected platform IDs", call. = FALSE)
stopifnot(identical(colnames(counts_raw), rownames(coldata)))

message(sprintf("Imported %d genes x %d samples; GC/length table: %d genes", nrow(counts_raw), ncol(counts_raw), nrow(gc_dt)))
print(table(coldata$cell_type, coldata$platform))

# ---- 3. Sample-level QC metrics -----------------------------------------------

ribo <- grepl("^RP[LS]\\d", rownames(counts_raw))
lib  <- colSums(counts_raw)
top100_share <- apply(counts_raw, 2, function(x) sum(sort(x, decreasing = TRUE)[1:100])) / lib

qc <- data.frame(
  coldata[, c("sample_id", "cell_type", "status", "platform", "donor")],
  lib_size   = lib,
  detected   = colSums(counts_raw > 0),
  pct_ribo   = 100 * colSums(counts_raw[ribo, ]) / lib,
  pct_top100 = 100 * top100_share
)
qc$flag_low_libsize  <- low_mad(log10(qc$lib_size))
qc$flag_low_detected <- low_mad(qc$detected)

# ---- 4. Gene filtering and DESeq2 object --------------------------------------

keep_samples <- setdiff(colnames(counts_raw), EXCLUDE)
counts <- counts_raw[, keep_samples]
cd     <- droplevels(coldata[keep_samples, ])

keep_genes <- rowSums(counts >= MIN_COUNT) >= MIN_SAMPLES
message(sprintf(
  "Gene filter (>= %d reads in >= %d samples): kept %d of %d genes (%d all-zero)",
  MIN_COUNT, MIN_SAMPLES, sum(keep_genes), nrow(counts), sum(rowSums(counts) == 0)
))

dds <- DESeqDataSetFromMatrix(counts[keep_genes, ], cd, design = ~ cell_type)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
vmat <- assay(vsd)
qc$size_factor <- NA_real_
qc[keep_samples, "size_factor"] <- sizeFactors(dds)

# ---- 5. Leave-one-out centroid identity check ---------------------------------

top_var <- head(order(rowVars(vmat), decreasing = TRUE), N_TOP_VAR)
nearest_group <- vapply(seq_len(ncol(vmat)), function(i) {
  r <- vapply(levels(cd$cell_type), function(g) {
    idx <- setdiff(which(cd$cell_type == g), i)
    if (!length(idx)) return(NA_real_)
    cor(vmat[top_var, i], rowMeans(vmat[top_var, idx, drop = FALSE]))
  }, numeric(1))
  names(which.max(r))
}, character(1))

qc$nearest_group <- NA_character_
qc[keep_samples, "nearest_group"] <- nearest_group
qc$flag_identity <- qc$status == "Normal" & !is.na(qc$nearest_group) &
  qc$nearest_group != as.character(qc$cell_type)

# ---- 6. QC Figures & NextSeq 500 GC-bias diagnostic ---------------------------

theme_set(theme_bw(base_size = 10))
sample_order <- rownames(cd)[order(cd$cell_type, cd$donor)]

qc_long <- melt(
  as.data.table(qc[keep_samples, ]),
  id.vars = c("sample_id", "cell_type"),
  measure.vars = c("lib_size", "detected", "pct_ribo", "pct_top100")
)
qc_long[, sample_id := factor(sample_id, levels = sample_order)]
p_lib <- ggplot(qc_long, aes(sample_id, value, fill = cell_type)) +
  geom_col() +
  facet_wrap(
    ~ variable, ncol = 1, scales = "free_y",
    labeller = as_labeller(c(
      lib_size = "Library size (reads)", detected = "Genes detected (>0)",
      pct_ribo = "% ribosomal protein reads", pct_top100 = "% reads in top 100 genes"
    ))
  ) +
  scale_fill_manual(values = CELL_COLS) +
  labs(x = NULL, y = NULL, fill = "Cell type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
ggsave(file.path(QC_DIR, "01_library_qc.png"), p_lib, width = 12, height = 9, dpi = 150)

pca_plot <- function(samples, title, colour = "cell_type") {
  m   <- vmat[, samples]
  sel <- head(order(rowVars(m), decreasing = TRUE), N_TOP_VAR)
  pc  <- prcomp(t(m[sel, ]))
  ve  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  df  <- cbind(cd[samples, ], pc$x[, 1:2])
  g <- ggplot(df, aes(PC1, PC2, colour = .data[[colour]], shape = platform)) +
    geom_point(size = 2.5, alpha = 0.9) +
    labs(title = title, x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]))
  if (colour == "cell_type") g + scale_colour_manual(values = CELL_COLS) else g
}
s_all  <- rownames(cd)
s_hspc <- rownames(cd)[cd$compartment != "Mature"]
s_plat <- rownames(cd)[cd$cell_type %in% c("pHSC", "Blast")]
p_pca <- (pca_plot(s_all, "All samples") | pca_plot(s_hspc, "HSPC + AML")) /
  (pca_plot(s_plat, "pHSC + Blast: platform check", colour = "platform") + aes(shape = cell_type) | plot_spacer()) +
  plot_annotation(caption = sprintf("Blind VST, top %d variable genes per panel", N_TOP_VAR))
ggsave(file.path(QC_DIR, "02_pca.png"), p_pca, width = 13, height = 10, dpi = 150)

d  <- dist(t(vmat[top_var, ]))
hc <- hclust(d, method = "average")
dm <- melt(as.data.table(as.matrix(d), keep.rownames = "s1"), id.vars = "s1", variable.name = "s2", value.name = "distance")
dm[, `:=`(s1 = factor(s1, levels = hc$labels[hc$order]), s2 = factor(s2, levels = hc$labels[hc$order]))]
p_dist <- ggplot(dm, aes(s1, s2, fill = distance)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", direction = -1) +
  labs(x = NULL, y = NULL, fill = "Euclidean\ndistance") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), axis.text.y = element_text(size = 5))
ggsave(file.path(QC_DIR, "03_sample_distance.png"), p_dist, width = 11, height = 10, dpi = 150)

mk <- intersect(names(MARKERS), rownames(vmat))
z  <- t(scale(t(vmat[mk, sample_order])))
zl <- melt(as.data.table(z, keep.rownames = "gene"), id.vars = "gene", variable.name = "sample_id", value.name = "z")
zl[, `:=`(gene = factor(gene, levels = rev(mk)), sample_id = factor(sample_id, levels = sample_order))]
zl[, cell_type := cd[as.character(sample_id), "cell_type"]]
p_mk <- ggplot(zl, aes(sample_id, gene, fill = pmax(pmin(z, 2.5), -2.5))) +
  geom_tile() +
  facet_grid(~ cell_type, scales = "free_x", space = "free_x") +
  scale_fill_distiller(palette = "RdBu", limits = c(-2.5, 2.5)) +
  scale_y_discrete(labels = function(g) sprintf("%s (%s)", g, MARKERS[g])) +
  labs(x = NULL, y = NULL, fill = "z (VST)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text = element_text(size = 6, angle = 90), panel.spacing = unit(0.5, "mm"))
ggsave(file.path(QC_DIR, "04_marker_heatmap.png"), p_mk, width = 13, height = 6, dpi = 150)

# NextSeq 500 GC/length bias diagnostic (justifies CQN offsets in Module 2)
ns_dds <- dds[, dds$platform == "NextSeq500"]
ns_dds$platform <- droplevels(ns_dds$platform)

common_gc    <- intersect(rownames(ns_dds), gc_dt$gene)
ns_gc_diag   <- estimateSizeFactors(ns_dds[common_gc, ])
gc_vec_diag  <- gc_dt[common_gc, gc_median]
len_vec_diag <- gc_dt[common_gc, length_median]

raw_norm_ns <- log2(counts(ns_gc_diag, normalized = TRUE) + 1)
cq_ns <- cqn(counts(ns_gc_diag), x = gc_vec_diag, lengths = len_vec_diag,
             sizeFactors = colSums(counts(ns_gc_diag)), verbose = FALSE)
cqn_norm_ns <- cq_ns$y + cq_ns$offset

diag_genes <- common_gc[rowSums(counts(ns_gc_diag) >= MIN_COUNT) >= 6]
gc_d <- gc_dt[diag_genes, gc_median]
ll_d <- log2(gc_dt[diag_genes, length_median])

hsc_others       <- c("5852-HSC", "6792-HSC", "7653-HSC")
within_hsc_pre   <- raw_norm_ns[diag_genes, "7256-HSC"] - rowMeans(raw_norm_ns[diag_genes, hsc_others])
within_lmpp_pre  <- rowMeans(raw_norm_ns[diag_genes, c("5483-LMPP", "7256-LMPP")]) - raw_norm_ns[diag_genes, "7653-LMPP"]
within_hsc_post  <- cqn_norm_ns[diag_genes, "7256-HSC"] - rowMeans(cqn_norm_ns[diag_genes, hsc_others])
within_lmpp_post <- rowMeans(cqn_norm_ns[diag_genes, c("5483-LMPP", "7256-LMPP")]) - cqn_norm_ns[diag_genes, "7653-LMPP"]

tech_axis_diag <- rbindlist(lapply(
  list(
    within_HSC_pre = within_hsc_pre, within_LMPP_pre = within_lmpp_pre,
    within_HSC_post = within_hsc_post, within_LMPP_post = within_lmpp_post
  ),
  function(v) as.data.table(r2_spline(v, gc_d, ll_d))
), idcol = "axis")
fwrite(tech_axis_diag, file.path(QC_DIR, "technical_axis_gc_r2.tsv"), sep = "\t")

med_pre  <- rowMedians(raw_norm_ns[diag_genes, ])
med_post <- rowMedians(cqn_norm_ns[diag_genes, ])
samp_gc_diag <- data.table(
  sample_id  = colnames(ns_gc_diag),
  cell_type  = as.character(ns_gc_diag$cell_type),
  donor      = as.character(ns_gc_diag$donor),
  status     = as.character(ns_gc_diag$status),
  lib_size   = colSums(counts(ns_gc_diag)),
  r_gc_pre   = apply(raw_norm_ns[diag_genes, ] - med_pre,  2, cor, y = gc_d),
  r_gc_post  = apply(cqn_norm_ns[diag_genes, ] - med_post, 2, cor, y = gc_d),
  r_len_pre  = apply(raw_norm_ns[diag_genes, ] - med_pre,  2, cor, y = ll_d),
  r_len_post = apply(cqn_norm_ns[diag_genes, ] - med_post, 2, cor, y = ll_d)
)
fwrite(samp_gc_diag, file.path(QC_DIR, "sample_gc_bias_diagnostic.tsv"), sep = "\t")

focus_cts <- c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "CLP", "Mono", "pHSC", "LSC", "Blast")
p_samp_gc <- ggplot(
  melt(samp_gc_diag[cell_type %in% focus_cts], id.vars = c("sample_id", "cell_type"), measure.vars = c("r_gc_pre", "r_gc_post")),
  aes(factor(cell_type, levels = focus_cts), value, colour = variable)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(position = position_dodge(width = 0.5), size = 2.2, alpha = 0.85) +
  scale_colour_manual(values = c(r_gc_pre = "#D73027", r_gc_post = "#4575B4"),
                      labels = c(r_gc_pre = "Before CQN (DESeq2 SF)", r_gc_post = "After CQN")) +
  labs(title = "Per-sample deviation from cohort median vs transcript GC%",
       subtitle = sprintf("NextSeq 500 HSPC + Mono + AML samples (n = %d expressed genes)", length(diag_genes)),
       x = NULL, y = "Pearson r (sample deviation vs GC%)", colour = NULL) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))

df_ax <- data.table(gc = gc_d, hsc_pre = within_hsc_pre, lmpp_pre = within_lmpp_pre,
                    hsc_post = within_hsc_post, lmpp_post = within_lmpp_post)
p_ax_pre <- ggplot(df_ax, aes(hsc_pre, lmpp_pre, colour = gc)) +
  geom_point(size = 0.6, alpha = 0.5) + scale_colour_viridis_c(option = "plasma", name = "GC%") +
  geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") +
  labs(title = sprintf("Before CQN: within-HSC vs within-LMPP (r = %+.3f)", cor(df_ax$hsc_pre, df_ax$lmpp_pre)),
       subtitle = sprintf("GC spline R2 = %.2f (HSC), %.2f (LMPP)",
                          tech_axis_diag[axis == "within_HSC_pre", R2_gc_spline],
                          tech_axis_diag[axis == "within_LMPP_pre", R2_gc_spline]),
       x = "7256-HSC vs other 3 HSCs (log2)", y = "5483/7256-LMPP vs 7653-LMPP (log2)")
p_ax_post <- ggplot(df_ax, aes(hsc_post, lmpp_post, colour = gc)) +
  geom_point(size = 0.6, alpha = 0.5) + scale_colour_viridis_c(option = "plasma", name = "GC%") +
  geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") +
  labs(title = sprintf("After CQN: within-HSC vs within-LMPP (r = %+.3f)", cor(df_ax$hsc_post, df_ax$lmpp_post)),
       subtitle = sprintf("GC spline R2 = %.3f (HSC), %.3f (LMPP)",
                          tech_axis_diag[axis == "within_HSC_post", R2_gc_spline],
                          tech_axis_diag[axis == "within_LMPP_post", R2_gc_spline]),
       x = "7256-HSC vs other 3 HSCs (CQN log2)", y = "5483/7256-LMPP vs 7653-LMPP (CQN log2)")
ggsave(file.path(QC_DIR, "05_gc_bias_diagnostic.png"), p_samp_gc / (p_ax_pre | p_ax_post), width = 14, height = 10, dpi = 150)

# ---- 7. Save QC outputs -------------------------------------------------------

fwrite(qc, file.path(QC_DIR, "qc_metrics.tsv"), sep = "\t")
saveRDS(dds, file.path(QC_DIR, "dds_filtered.rds"))
saveRDS(vsd, file.path(QC_DIR, "vsd_blind.rds"))
writeLines(capture.output(sessionInfo()), file.path(QC_DIR, "sessionInfo.txt"))

message("GC/length technical axis R2 (before vs after CQN):")
print(tech_axis_diag)
message("Module 1 QC outputs written to ", QC_DIR)

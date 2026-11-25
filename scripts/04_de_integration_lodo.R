#!/usr/bin/env Rscript
# Module 2: NextSeq-only CQN + DESeq2 + apeglm differential expression across Models A, B, C, and C2,
# GO:BP fGSEA + ORA enrichment, cross-model integration, matched-N H1/H3 specificity,
# integrated Figure 3, and 6-patient Leave-One-Donor-Out (LODO) sensitivity analysis.
#
# Input : results/qc/dds_filtered.rds, data/GSE74246/gene_gc_length.tsv
# Output: results/de/{de_*.tsv, fgsea_*.tsv, ora_*.tsv, summary TSVs, 01..07 & Fig3 PNGs, de_results.rds}

suppressPackageStartupMessages({
  library(DESeq2)
  library(cqn)
  library(apeglm)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(ggforce)
  library(fgsea)
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
})

cmd_args   <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", grep("^--file=", cmd_args, value = TRUE))
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(script_arg[1])) else file.path(getwd(), "scripts")
source(file.path(dirname(script_dir), "R", "utils.R"))

ROOT     <- find_root(file.path("data", GSE))
DATA_DIR <- file.path(ROOT, "data", GSE)
QC_DIR   <- file.path(ROOT, "results", "qc")
DE_DIR   <- file.path(ROOT, "results", "de")
dir.create(DE_DIR, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 10))

dds_path <- file.path(QC_DIR, "dds_filtered.rds")
gc_path  <- file.path(DATA_DIR, "gene_gc_length.tsv")
if (!file.exists(dds_path)) stop("Missing ", dds_path, ". Run `Rscript scripts/03_qc_and_cqn.R` first.", call. = FALSE)
if (!file.exists(gc_path))  stop("Missing ", gc_path, ". Run `python scripts/02_build_gene_gc.py` first.", call. = FALSE)

dds   <- readRDS(dds_path)
gc_dt <- fread(gc_path)
gc_dt <- gc_dt[is.finite(gc_median) & is.finite(length_median) & length_median > 0]
setkey(gc_dt, gene)

ns_dds <- dds[, dds$platform == "NextSeq500"]
ns_dds$platform <- droplevels(ns_dds$platform)
message(sprintf("Module 2: NextSeq 500 subset (%d samples) with CQN + apeglm", ncol(ns_dds)))

fit_deseq_model <- function(dds_pool, sample_mask, design_formula, label, gc_table = gc_dt) {
  sub <- dds_pool[, sample_mask]
  colData(sub) <- droplevels(colData(sub))
  min_grp   <- min(table(sub$cell_type))
  keep_expr <- rowSums(counts(sub) >= MIN_COUNT) >= min_grp
  keep_gc   <- rownames(sub) %in% gc_table$gene
  n_expr    <- sum(keep_expr)
  sub       <- sub[keep_expr & keep_gc, ]

  # Centre per-gene CQN normalisation factors to geometric mean 1 so baseMean retains count scale.
  gc_v  <- gc_table[rownames(sub), gc_median]
  len_v <- gc_table[rownames(sub), length_median]
  cq <- cqn(counts(sub), x = gc_v, lengths = len_v, sizeFactors = colSums(counts(sub)), verbose = FALSE)
  nf <- exp(cq$glm.offset)
  nf <- nf / exp(rowMeans(log(nf)))
  normalizationFactors(sub) <- nf
  design(sub) <- design_formula
  message(sprintf("Fitting CQN %s: %d samples, %d/%d expressed genes with GC/len (%.1f%%)",
                  label, ncol(sub), nrow(sub), n_expr, 100 * nrow(sub) / n_expr))
  DESeq(sub, quiet = TRUE)
}

extract_contrast <- function(dds_sub, contrast_name, num, den) {
  dds_re <- dds_sub
  dds_re$cell_type <- relevel(dds_re$cell_type, ref = den)
  dds_re <- nbinomWaldTest(dds_re, quiet = TRUE)
  coef_name <- paste0("cell_type_", num, "_vs_", den)
  res <- results(dds_re, name = coef_name, alpha = PADJ_CUTOFF)
  lfc <- lfcShrink(dds_re, coef = coef_name, res = res, type = "apeglm", quiet = TRUE)
  dt <- data.table(
    gene          = rownames(res),
    contrast      = contrast_name,
    baseMean      = res$baseMean,
    log2FC_raw    = res$log2FoldChange,
    log2FC_shrunk = lfc$log2FoldChange,
    lfcSE         = res$lfcSE,
    stat          = res$stat,
    pvalue        = res$pvalue,
    padj          = res$padj
  )
  dt[, direction := fifelse(!is.na(padj) & padj < PADJ_CUTOFF & log2FC_shrunk > LFC_CUTOFF, "Up",
                    fifelse(!is.na(padj) & padj < PADJ_CUTOFF & log2FC_shrunk < -LFC_CUTOFF, "Down", "NS"))]
  dt[order(padj, pvalue)]
}

build_gobp_pathways <- function(genes, min_size = 15, max_size = 500) {
  go_tab <- unique(AnnotationDbi::select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", columns = c("GOALL", "ONTOLOGYALL")))
  go_bp  <- go_tab[go_tab$ONTOLOGYALL == "BP" & !is.na(go_tab$GOALL), c("GOALL", "SYMBOL")]
  pw     <- split(go_bp$SYMBOL, go_bp$GOALL)
  pw     <- pw[lengths(pw) >= min_size & lengths(pw) <= max_size]
  terms  <- AnnotationDbi::select(GO.db, keys = names(pw), keytype = "GOID", columns = "TERM")
  list(pathways = pw, term_map = setNames(terms$TERM, terms$GOID))
}

run_enrichment <- function(de_dt, gobp) {
  universe  <- de_dt[!is.na(padj), gene]
  ranked_dt <- de_dt[!is.na(stat)]
  ranks     <- sort(setNames(ranked_dt$stat, ranked_dt$gene), decreasing = TRUE)

  set.seed(42)
  fg <- fgseaMultilevel(pathways = gobp$pathways, stats = ranks, minSize = 15, maxSize = 500, eps = 1e-10)
  fg[, term := unname(gobp$term_map[pathway])]
  sig_fg <- fg[!is.na(padj) & padj < PADJ_CUTOFF][order(pval)]
  main_ids <- if (nrow(sig_fg) >= 2) collapsePathways(head(sig_fg, 250), gobp$pathways, ranks)$mainPathways else sig_fg$pathway
  fg[, is_main := pathway %in% main_ids]
  fg[, leadingEdge_genes := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
  fg_out <- fg[, .(contrast = de_dt$contrast[1], pathway, term, pval, padj, log2err, ES, NES,
                   size, is_main, leadingEdge = leadingEdge_genes)][order(padj, pval)]

  ora_parts <- lapply(c("Up", "Down"), function(dir) {
    sig_genes <- de_dt[direction == dir, gene]
    if (length(sig_genes) < 5) return(NULL)
    o <- fora(pathways = gobp$pathways, genes = sig_genes, universe = universe, minSize = 15, maxSize = 500)
    o[, `:=`(
      contrast      = de_dt$contrast[1],
      direction     = dir,
      term          = unname(gobp$term_map[pathway]),
      n_sig_genes   = length(sig_genes),
      n_universe    = length(universe),
      gene_ratio    = overlap / length(sig_genes),
      overlap_genes = vapply(overlapGenes, paste, collapse = ";", FUN.VALUE = character(1))
    )]
    o[, .(contrast, direction, pathway, term, pval, padj, overlap, size,
          n_sig_genes, n_universe, gene_ratio, overlap_genes)][order(padj, pval)]
  })
  ora_out <- rbindlist(ora_parts, use.names = TRUE, fill = TRUE)
  if (!nrow(ora_out)) {
    ora_out <- data.table(
      contrast = character(), direction = character(), pathway = character(), term = character(),
      pval = numeric(), padj = numeric(), overlap = integer(), size = integer(),
      n_sig_genes = integer(), n_universe = integer(), gene_ratio = numeric(), overlap_genes = character()
    )
  }
  list(fgsea = fg_out, ora = ora_out)
}

# ---- 1. Fit Models A, B, C, C2 and run GO:BP fGSEA + ORA ---------------------

trios6   <- c("SU353", "SU444", "SU496", "SU575", "SU583", "SU654")
donors_c <- c("5852", "6792", "7256", "7653")

dds_A  <- fit_deseq_model(ns_dds, ns_dds$status == "AML" & ns_dds$donor %in% trios6,
                          ~ donor + cell_type, "Model A (6 NextSeq AML trios)")
dds_B  <- fit_deseq_model(ns_dds, ns_dds$compartment != "Mature",
                          ~ cell_type, "Model B (NextSeq HSPC + AML)")
dds_C  <- fit_deseq_model(ns_dds, ns_dds$donor %in% donors_c & ns_dds$cell_type %in% c("HSC", "MPP", "GMP", "Mono"),
                          ~ donor + cell_type, "Model C (4 paired normal donors)")
dds_C2 <- fit_deseq_model(ns_dds, ns_dds$compartment == "HSPC" | ns_dds$cell_type == "Mono",
                          ~ cell_type, "Model C2 (NextSeq normal HSPC + Mono)")

contrast_specs <- list(
  M1_LSC_vs_pHSC  = list(model = "A",  dds = dds_A,  num = "LSC",   den = "pHSC", label = "M1: LSC vs pHSC"),
  M2_Blast_vs_LSC = list(model = "A",  dds = dds_A,  num = "Blast", den = "LSC",  label = "M2: Blast vs LSC"),
  LSC_vs_GMP      = list(model = "B",  dds = dds_B,  num = "LSC",   den = "GMP",  label = "LSC vs GMP"),
  LSC_vs_LMPP     = list(model = "B",  dds = dds_B,  num = "LSC",   den = "LMPP", label = "LSC vs LMPP"),
  pHSC_vs_HSC     = list(model = "B",  dds = dds_B,  num = "pHSC",  den = "HSC",  label = "pHSC vs HSC"),
  GMP_vs_HSC      = list(model = "C",  dds = dds_C,  num = "GMP",   den = "HSC",  label = "GMP vs HSC"),
  Mono_vs_GMP     = list(model = "C",  dds = dds_C,  num = "Mono",  den = "GMP",  label = "Mono vs GMP"),
  LMPP_vs_HSC     = list(model = "C2", dds = dds_C2, num = "LMPP",  den = "HSC",  label = "LMPP vs HSC"),
  GMP_vs_LMPP     = list(model = "C2", dds = dds_C2, num = "GMP",   den = "LMPP", label = "GMP vs LMPP")
)

de_res <- lapply(names(contrast_specs), function(nm) {
  sp <- contrast_specs[[nm]]
  dt <- extract_contrast(sp$dds, nm, sp$num, sp$den)
  fwrite(dt, file.path(DE_DIR, sprintf("de_%s.tsv", nm)), sep = "\t")
  dt
})
names(de_res) <- names(contrast_specs)

gobp <- build_gobp_pathways(rownames(ns_dds))
enrich_res <- lapply(names(de_res), function(nm) {
  er <- run_enrichment(de_res[[nm]], gobp)
  fwrite(er$fgsea, file.path(DE_DIR, sprintf("fgsea_%s.tsv", nm)), sep = "\t")
  fwrite(er$ora,   file.path(DE_DIR, sprintf("ora_%s.tsv", nm)),   sep = "\t")
  er
})
names(enrich_res) <- names(de_res)

# ---- 2. Plotting helpers & per-model composite panels -------------------------

plot_volcano <- function(de_dt, title, n_label = 12, xlim = c(-11.2, 11.2), y_cap = NULL, ylim = NULL) {
  df <- rank_volcano_degs(de_dt)
  df[, capped := if (!is.null(y_cap)) neg_log10_padj > y_cap else FALSE]
  df[, y_plot := if (!is.null(y_cap)) pmin(neg_log10_padj, y_cap) else neg_log10_padj]
  labs_dt  <- head(df[direction != "NS"][order(-score)], n_label)
  x_rng    <- range(df$log2FC_shrunk, na.rm = TRUE)
  xlim_use <- c(min(xlim[1], x_rng[1] - 0.5), max(xlim[2], x_rng[2] + 0.5))
  sub_txt  <- sprintf("Up: %d | Down: %d (padj < %.2f, |apeglm LFC| > %.1f)",
                      sum(df$direction == "Up"), sum(df$direction == "Down"), PADJ_CUTOFF, LFC_CUTOFF)
  if (sum(df$capped) > 0) sub_txt <- paste0(sub_txt, sprintf(" | %d capped at y=%.0f", sum(df$capped), y_cap))

  g <- ggplot(df, aes(log2FC_shrunk, y_plot, colour = direction, shape = capped)) +
    geom_point(size = 0.85, alpha = 0.65) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", linewidth = 0.35, colour = "grey45") +
    geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", linewidth = 0.35, colour = "grey45") +
    scale_colour_manual(values = DIR_COLS, drop = FALSE) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), guide = "none") +
    coord_cartesian(xlim = xlim_use, ylim = ylim) +
    labs(title = title, subtitle = sub_txt,
         x = "Shrunk log2 fold change (apeglm)", y = "-log10(adjusted p)", colour = NULL) +
    theme(legend.position = "bottom")
  if (nrow(labs_dt)) {
    g <- g + geom_text_repel(data = labs_dt, aes(x = log2FC_shrunk, y = y_plot, label = gene),
                             size = 2.4, max.overlaps = 15, box.padding = 0.25,
                             show.legend = FALSE, inherit.aes = FALSE)
  }
  g
}

plot_fgsea_panel <- function(contrast_keys, top_n = 5, wrap_width = 52) {
  picked <- rbindlist(lapply(contrast_keys, function(k) {
    fg <- enrich_res[[k]]$fgsea[!is.na(padj) & padj < PADJ_CUTOFF]
    if (!nrow(fg)) return(NULL)
    fg_main <- if (nrow(fg[is_main == TRUE])) fg[is_main == TRUE] else fg
    rbind(head(fg_main[NES > 0][order(-NES)], top_n), head(fg_main[NES < 0][order(NES)], top_n))
  }))
  if (!nrow(picked)) {
    return(ggplot() + annotate("text", x = 0, y = 0, label = "No fGSEA GO:BP terms at padj < 0.05") + theme_void() + labs(title = "fGSEA GO:BP"))
  }
  picked_ids <- unique(picked$pathway)
  id_labels  <- setNames(make.unique(wrap_term(gobp$term_map[picked_ids], wrap_width)), picked_ids)
  tab <- rbindlist(lapply(contrast_keys, function(k) {
    enrich_res[[k]]$fgsea[pathway %in% picked_ids,
                          .(clabel = contrast_specs[[k]]$label, term_label = unname(id_labels[pathway]),
                            NES, padj = pmax(padj, 1e-50))]
  }))
  tab[, clabel := factor(clabel, levels = vapply(contrast_keys, function(k) contrast_specs[[k]]$label, character(1)))]
  ord <- tab[, .(m = max(abs(NES))), by = term_label][order(m), term_label]
  tab[, term_label := factor(term_label, levels = ord)]
  ggplot(tab, aes(clabel, term_label, colour = NES, size = pmin(-log10(padj), 15))) +
    geom_point() +
    scale_colour_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B", midpoint = 0) +
    scale_size_continuous(range = c(1.5, 5), name = "-log10(padj)") +
    labs(title = "fGSEA (main non-redundant GO:BP)", x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 8),
          axis.text.y = element_text(size = 6.6, lineheight = 0.85))
}

plot_ora_panel <- function(contrast_keys, top_n = 4, wrap_width = 52) {
  picked <- rbindlist(lapply(contrast_keys, function(k) {
    o <- enrich_res[[k]]$ora[!is.na(padj) & padj < PADJ_CUTOFF]
    if (!nrow(o)) return(NULL)
    o[order(padj), head(.SD, top_n), by = direction]
  }))
  if (!nrow(picked)) {
    return(ggplot() + annotate("text", x = 0, y = 0, label = "No ORA GO:BP terms at padj < 0.05") + theme_void() + labs(title = "ORA GO:BP"))
  }
  grp_lvls <- as.vector(outer(vapply(contrast_keys, function(k) contrast_specs[[k]]$label, character(1)),
                              c("Up", "Down"), function(a, b) paste0(a, "\n(", b, ")")))
  picked[, group := factor(paste0(vapply(contrast, function(k) contrast_specs[[k]]$label, character(1)), "\n(", direction, ")"), levels = grp_lvls)]
  picked_ids <- unique(picked$pathway)
  id_labels  <- setNames(make.unique(wrap_term(gobp$term_map[picked_ids], wrap_width)), picked_ids)
  picked[, term_label := unname(id_labels[pathway])]
  ord <- picked[, .(m = max(-log10(pmax(padj, 1e-50)))), by = term_label][order(m), term_label]
  picked[, term_label := factor(term_label, levels = ord)]
  ggplot(picked, aes(group, term_label, colour = direction, size = gene_ratio)) +
    geom_point(aes(alpha = pmin(-log10(pmax(padj, 1e-50)), 15))) +
    scale_colour_manual(values = DIR_COLS[c("Up", "Down")]) +
    scale_alpha_continuous(range = c(0.45, 1), name = "-log10(padj)") +
    scale_size_continuous(range = c(1.8, 5), name = "Gene ratio") +
    labs(title = "ORA (GO:BP by direction)", x = NULL, y = NULL, colour = NULL) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 7.5),
          axis.text.y = element_text(size = 6.6, lineheight = 0.85))
}

plot_lfc_scatter <- function(k_x, k_y, title = NULL, n_label = 10, shared_note = NULL) {
  dx <- de_res[[k_x]][!is.na(padj), .(gene, x = log2FC_shrunk, dir_x = direction)]
  dy <- de_res[[k_y]][!is.na(padj), .(gene, y = log2FC_shrunk, dir_y = direction)]
  m  <- merge(dx, dy, by = "gene")
  m[, quad := fifelse(dir_x == "Up" & dir_y == "Up", "Concordant Up",
              fifelse(dir_x == "Down" & dir_y == "Down", "Concordant Down",
              fifelse((dir_x == "Up" & dir_y == "Down") | (dir_x == "Down" & dir_y == "Up"), "Discordant", "Single / NS")))]
  labs_dt <- head(m[quad != "Single / NS"][order(-(abs(x) + abs(y)))], n_label)
  sub_txt <- sprintf("Pearson r = %.3f | Spearman rho = %.3f (n = %d genes)", cor(m$x, m$y), cor(m$x, m$y, method = "spearman"), nrow(m))
  if (!is.null(shared_note)) sub_txt <- paste0(sub_txt, "\n", shared_note)

  g <- ggplot(m, aes(x, y, colour = quad)) +
    geom_point(data = m[quad == "Single / NS"], size = 0.7, alpha = 0.4) +
    geom_point(data = m[quad != "Single / NS"], size = 1.1, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey45") +
    scale_colour_manual(values = QUAD_COLS) +
    labs(title = title %||% sprintf("%s vs %s", contrast_specs[[k_y]]$label, contrast_specs[[k_x]]$label),
         subtitle = sub_txt,
         x = paste0(contrast_specs[[k_x]]$label, " (apeglm log2FC)"),
         y = paste0(contrast_specs[[k_y]]$label, " (apeglm log2FC)"), colour = NULL) +
    theme(legend.position = "bottom", plot.subtitle = element_text(size = 8))
  if (nrow(labs_dt)) {
    g <- g + geom_text_repel(data = labs_dt, aes(label = gene), size = 2.3, max.overlaps = 14, box.padding = 0.25, show.legend = FALSE)
  }
  g
}

pick_heatmap_degs <- function(contrast_keys, n_each = 8) {
  unique(unlist(lapply(contrast_keys, function(k) {
    dt <- de_res[[k]]
    c(head(dt[direction == "Up", gene], n_each), head(dt[direction == "Down", gene], n_each))
  })))
}

plot_deg_heatmap <- function(dds_sub, genes, title, facet_by = c("donor", "cell_type"), cell_levels = NULL, z_cap = 2.5) {
  facet_by <- match.arg(facet_by)
  vsd_sub  <- vst(dds_sub, blind = FALSE)
  genes    <- intersect(genes, rownames(vsd_sub))
  z        <- t(scale(t(assay(vsd_sub)[genes, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  gene_ord <- if (nrow(z) >= 3) rownames(z)[hclust(dist(z))$order] else rev(genes)
  cd_sub   <- as.data.table(as.data.frame(colData(dds_sub)), keep.rownames = "sample_id")
  if (!is.null(cell_levels)) {
    cd_sub <- cd_sub[as.character(cell_type) %in% cell_levels]
    cd_sub[, cell_type := factor(as.character(cell_type), levels = cell_levels)]
  } else {
    cd_sub[, cell_type := factor(as.character(cell_type))]
  }
  zl <- melt(as.data.table(z[, cd_sub$sample_id, drop = FALSE], keep.rownames = "gene"),
             id.vars = "gene", variable.name = "sample_id", value.name = "z")
  zl[, sample_id := as.character(sample_id)]
  zl <- merge(zl, cd_sub[, .(sample_id, donor = as.character(donor), cell_type)], by = "sample_id")
  zl[, gene := factor(gene, levels = gene_ord)]

  if (facet_by == "donor") {
    ggplot(zl, aes(cell_type, gene, fill = pmax(pmin(z, z_cap), -z_cap))) +
      geom_tile() + facet_grid(~ donor) +
      scale_fill_distiller(palette = "RdBu", limits = c(-z_cap, z_cap)) +
      labs(title = title, x = NULL, y = NULL, fill = "z (CQN VST)") +
      theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7.5), axis.text.y = element_text(size = 6.5))
  } else {
    samp_ord <- cd_sub[order(cell_type, donor, sample_id), sample_id]
    zl[, sample_id := factor(sample_id, levels = samp_ord)]
    ggplot(zl, aes(sample_id, gene, fill = pmax(pmin(z, z_cap), -z_cap))) +
      geom_tile() + facet_grid(~ cell_type, scales = "free_x", space = "free_x") +
      scale_fill_distiller(palette = "RdBu", limits = c(-z_cap, z_cap)) +
      labs(title = title, x = NULL, y = NULL, fill = "z (CQN VST)") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 6.5))
  }
}

fig_A <- (plot_volcano(de_res$M1_LSC_vs_pHSC, "M1: LSC vs pHSC", y_cap = 55, ylim = c(0, 60)) |
          plot_volcano(de_res$M2_Blast_vs_LSC, "M2: Blast vs LSC", y_cap = 55, ylim = c(0, 60))) /
         (plot_deg_heatmap(dds_A, pick_heatmap_degs(c("M1_LSC_vs_pHSC", "M2_Blast_vs_LSC"), 8),
                           "Top M1 / M2 DEGs per NextSeq trio (CQN)", facet_by = "donor", cell_levels = c("pHSC", "LSC", "Blast")) |
          plot_lfc_scatter("M1_LSC_vs_pHSC", "M2_Blast_vs_LSC", "Partitioning pHSC -> LSC -> Blast",
                           shared_note = "Note: contrasts share LSC with opposite signs (mechanical negative coupling)")) /
         (plot_fgsea_panel(c("M1_LSC_vs_pHSC", "M2_Blast_vs_LSC")) | plot_ora_panel(c("M1_LSC_vs_pHSC", "M2_Blast_vs_LSC"))) +
  plot_layout(heights = c(1, 1.05, 1.25)) +
  plot_annotation(title = "Model A: Paired AML progression (6 NextSeq trios, CQN + ~ donor + cell_type, apeglm)")
ggsave(file.path(DE_DIR, "01_modelA_AML_progression.png"), fig_A, width = 18.5, height = 16, dpi = 150)

fig_B <- (plot_volcano(de_res$LSC_vs_GMP, "LSC vs GMP", y_cap = 35, ylim = c(0, 38)) |
          plot_volcano(de_res$LSC_vs_LMPP, "LSC vs LMPP", y_cap = 35, ylim = c(0, 38)) |
          plot_volcano(de_res$pHSC_vs_HSC, "pHSC vs HSC", y_cap = 35, ylim = c(0, 38))) /
         ((plot_deg_heatmap(dds_B, pick_heatmap_degs(c("LSC_vs_GMP", "LSC_vs_LMPP", "pHSC_vs_HSC"), 6),
                            "Top Model B DEGs across NextSeq HSPC + AML (CQN)", facet_by = "cell_type",
                            cell_levels = c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "CLP", "pHSC", "LSC", "Blast")) |
           plot_lfc_scatter("LSC_vs_GMP", "LSC_vs_LMPP", "LSC vs LMPP vs LSC vs GMP",
                            shared_note = "Note: contrasts share LSC with same sign (mechanical positive coupling)")) + plot_layout(widths = c(1.35, 1))) /
         (plot_fgsea_panel(c("LSC_vs_GMP", "LSC_vs_LMPP", "pHSC_vs_HSC"), top_n = 4) |
          plot_ora_panel(c("LSC_vs_GMP", "LSC_vs_LMPP", "pHSC_vs_HSC"), top_n = 3)) +
  plot_layout(heights = c(1, 1.05, 1.3)) +
  plot_annotation(title = "Model B: AML vs healthy counterparts (48 NextSeq HSPC + AML samples, CQN + ~ cell_type, apeglm)")
ggsave(file.path(DE_DIR, "02_modelB_AML_vs_normal.png"), fig_B, width = 19.5, height = 16.5, dpi = 150)

fig_C <- (plot_volcano(de_res$GMP_vs_HSC, "GMP vs HSC", y_cap = 61, ylim = c(0, 66)) |
          plot_volcano(de_res$Mono_vs_GMP, "Mono vs GMP", y_cap = 61, ylim = c(0, 66))) /
         (plot_deg_heatmap(dds_C, pick_heatmap_degs(c("GMP_vs_HSC", "Mono_vs_GMP"), 8),
                           "Top Model C DEGs per healthy donor (CQN)", facet_by = "donor", cell_levels = c("HSC", "MPP", "GMP", "Mono")) |
          plot_lfc_scatter("GMP_vs_HSC", "Mono_vs_GMP", "Normal myeloid transitions",
                           shared_note = "Note: contrasts share GMP with opposite signs (mechanical negative coupling)")) /
         (plot_fgsea_panel(c("GMP_vs_HSC", "Mono_vs_GMP")) | plot_ora_panel(c("GMP_vs_HSC", "Mono_vs_GMP"))) +
  plot_layout(heights = c(1, 1.05, 1.25)) +
  plot_annotation(title = "Model C: Paired normal differentiation (4 healthy donors x 4 populations, CQN + ~ donor + cell_type, apeglm)")
ggsave(file.path(DE_DIR, "03_modelC_normal_paired.png"), fig_C, width = 18.5, height = 16, dpi = 150)

fig_C2 <- (plot_volcano(de_res$LMPP_vs_HSC, "LMPP vs HSC", y_cap = 35, ylim = c(0, 38)) |
           plot_volcano(de_res$GMP_vs_LMPP, "GMP vs LMPP", y_cap = 35, ylim = c(0, 38))) /
          (plot_deg_heatmap(dds_C2, pick_heatmap_degs(c("LMPP_vs_HSC", "GMP_vs_LMPP"), 8),
                            "Top Model C2 DEGs across NextSeq normal HSPC + Mono (CQN)", facet_by = "cell_type",
                            cell_levels = c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "CLP", "Mono")) |
           plot_lfc_scatter("LMPP_vs_HSC", "GMP_vs_LMPP", "Partitioning HSC -> LMPP -> GMP",
                            shared_note = "Note: contrasts share LMPP with opposite signs (mechanical negative coupling)")) /
          (plot_fgsea_panel(c("LMPP_vs_HSC", "GMP_vs_LMPP")) | plot_ora_panel(c("LMPP_vs_HSC", "GMP_vs_LMPP"))) +
  plot_layout(heights = c(1, 1.05, 1.25)) +
  plot_annotation(title = "Model C2: Normal LMPP priming steps (30 NextSeq normal HSPC + Mono samples, CQN + ~ cell_type, apeglm)")
ggsave(file.path(DE_DIR, "04_modelC2_normal_LMPP.png"), fig_C2, width = 18.5, height = 16, dpi = 150)

# ---- 3. Cross-model overlap, size-matched H1/H3 specificity & Fig 3 ----------

common_univ <- Reduce(intersect, lapply(de_res, function(dt) dt[!is.na(padj), gene]))
cc_genes <- unique(AnnotationDbi::select(org.Hs.eg.db, keys = "GO:0007049", keytype = "GOALL", columns = "SYMBOL")$SYMBOL)
stopifnot(length(cc_genes) > 1000)
ribo_genes      <- grep("^RP[LS]\\d", common_univ, value = TRUE)
non_prolif_univ <- setdiff(common_univ, union(cc_genes, ribo_genes))

sets_hyp <- list(
  HSC_gt_GMP   = intersect(de_res$GMP_vs_HSC[direction == "Down", gene], common_univ),
  HSC_gt_LMPP  = intersect(de_res$LMPP_vs_HSC[direction == "Down", gene], common_univ),
  LSC_gt_Blast = intersect(de_res$M2_Blast_vs_LSC[direction == "Down", gene], common_univ),
  LSC_gt_LMPP  = intersect(de_res$LSC_vs_LMPP[direction == "Up", gene], common_univ),
  LSC_gt_GMP   = intersect(de_res$LSC_vs_GMP[direction == "Up", gene], common_univ),
  M1_up        = intersect(de_res$M1_LSC_vs_pHSC[direction == "Up", gene], common_univ),
  LMPP_gt_HSC  = intersect(de_res$LMPP_vs_HSC[direction == "Up", gene], common_univ),
  GMP_gt_HSC   = intersect(de_res$GMP_vs_HSC[direction == "Up", gene], common_univ)
)

overlap_pairs <- list(
  list("HSC>GMP vs LSC>Blast",  "HSC_gt_GMP",  "LSC_gt_Blast"),
  list("HSC>GMP vs LSC>LMPP",   "HSC_gt_GMP",  "LSC_gt_LMPP"),
  list("HSC>LMPP vs LSC>Blast", "HSC_gt_LMPP", "LSC_gt_Blast"),
  list("HSC>LMPP vs LSC>GMP",   "HSC_gt_LMPP", "LSC_gt_GMP"),
  list("M1 Up vs LMPP>HSC",     "M1_up",       "LMPP_gt_HSC"),
  list("M1 Up vs GMP>HSC",      "M1_up",       "GMP_gt_HSC")
)
overlap_stats <- rbindlist(lapply(overlap_pairs, function(p) {
  rbind(
    fisher_overlap(p[[1]], sets_hyp[[p[[2]]]], sets_hyp[[p[[3]]]], common_univ, p[[2]], p[[3]], "All tested genes"),
    fisher_overlap(p[[1]], sets_hyp[[p[[2]]]], sets_hyp[[p[[3]]]], non_prolif_univ, p[[2]], p[[3]], "No cell-cycle/ribo")
  )
}))
overlap_stats[, padj := p.adjust(p_value, method = "BH")]
fwrite(overlap_stats, file.path(DE_DIR, "stemness_overlap_stats.tsv"), sep = "\t")

core_stem <- data.table(gene = sort(unique(c(intersect(sets_hyp$HSC_gt_GMP, sets_hyp$LSC_gt_Blast),
                                             intersect(sets_hyp$HSC_gt_LMPP, sets_hyp$LSC_gt_Blast)))))
core_stem[, `:=`(
  in_HSC_gt_GMP         = gene %in% sets_hyp$HSC_gt_GMP,
  in_HSC_gt_LMPP        = gene %in% sets_hyp$HSC_gt_LMPP,
  in_LSC_gt_Blast       = gene %in% sets_hyp$LSC_gt_Blast,
  in_LSC_gt_LMPP        = gene %in% sets_hyp$LSC_gt_LMPP,
  in_LSC_gt_GMP         = gene %in% sets_hyp$LSC_gt_GMP,
  is_cell_cycle_or_ribo = gene %in% union(cc_genes, ribo_genes),
  LFC_GMP_vs_HSC        = de_res$GMP_vs_HSC[match(core_stem$gene, de_res$GMP_vs_HSC$gene), log2FC_shrunk],
  LFC_Blast_vs_LSC      = de_res$M2_Blast_vs_LSC[match(core_stem$gene, de_res$M2_Blast_vs_LSC$gene), log2FC_shrunk],
  LFC_LSC_vs_LMPP       = de_res$LSC_vs_LMPP[match(core_stem$gene, de_res$LSC_vs_LMPP$gene), log2FC_shrunk],
  LFC_LSC_vs_GMP        = de_res$LSC_vs_GMP[match(core_stem$gene, de_res$LSC_vs_GMP$gene), log2FC_shrunk]
)]
fwrite(core_stem, file.path(DE_DIR, "stemness_shared_genes.tsv"), sep = "\t")

stat_of <- function(contrast, res_list = de_res, univ = common_univ) {
  d <- res_list[[contrast]][gene %in% univ & !is.na(stat)]
  setNames(d$stat, d$gene)
}
S_M1 <- stat_of("M1_LSC_vs_pHSC"); S_LMPP <- stat_of("LMPP_vs_HSC"); S_GMP <- stat_of("GMP_vs_HSC")
N_GRID <- c(150, 300, 500)

spec_H1 <- rbindlist(lapply(N_GRID, function(N) rbindlist(lapply(c("Up", "Down"), function(dir) {
  L <- top_genes(S_LMPP, N, dir); G <- top_genes(S_GMP, N, dir)
  x <- S_M1[setdiff(L, G)]; y <- S_M1[setdiff(G, L)]
  a <- if (dir == "Up") x else y; b <- if (dir == "Up") y else x
  set.seed(42)
  bs <- replicate(1000, auc_rank(sample(a, replace = TRUE), sample(b, replace = TRUE)))
  data.table(
    N = N, direction = dir, n_shared = length(intersect(L, G)),
    n_LMPP_specific = length(x), n_GMP_specific = length(y),
    median_M1_stat_LMPP_specific = median(x), median_M1_stat_GMP_specific = median(y),
    AUC_LMPP_like = auc_rank(a, b),
    AUC_lo = unname(quantile(bs, 0.025)), AUC_hi = unname(quantile(bs, 0.975)),
    p_wilcox_descriptive = wilcox.test(a, b, alternative = "greater", exact = FALSE)$p.value
  )
}))))
fwrite(spec_H1, file.path(DE_DIR, "specificity_H1_matchedN.tsv"), sep = "\t")

reg_fit <- lm(rk(S_M1, common_univ) ~ rk(S_LMPP, common_univ) + rk(S_GMP, common_univ))
reg_tab <- data.table(
  term     = c("LMPP_vs_HSC", "GMP_vs_HSC"),
  std_beta = unname(coef(reg_fit)[-1]),
  se       = summary(reg_fit)$coefficients[-1, "Std. Error"],
  R2       = summary(reg_fit)$r.squared,
  n_genes  = length(common_univ)
)
fwrite(reg_tab, file.path(DE_DIR, "specificity_H1_regression.tsv"), sep = "\t")

h3_universes <- list("All tested genes" = common_univ, "No cell-cycle/ribo" = non_prolif_univ)
h3_matched <- rbindlist(lapply(N_GRID, function(N) rbindlist(lapply(c("GMP", "LMPP"), function(route) {
  hsc_set <- top_genes(if (route == "GMP") S_GMP else S_LMPP, N, "Down")
  rbindlist(lapply(names(h3_universes), function(u) {
    as.data.table(c(list(N = N, route = paste0("HSC>", route), filter = u),
                    fisher_ci(hsc_set, sets_hyp$LSC_gt_Blast, h3_universes[[u]])))
  }))
}))))
fwrite(h3_matched, file.path(DE_DIR, "stemness_H3_matchedN.tsv"), sep = "\t")

nes_pair <- function(k_x, k_y) {
  fx <- enrich_res[[k_x]]$fgsea[, .(pathway, term, NES_x = NES, padj_x = padj, main_x = is_main)]
  fy <- enrich_res[[k_y]]$fgsea[, .(pathway, NES_y = NES, padj_y = padj, main_y = is_main)]
  m  <- merge(fx, fy, by = "pathway")[is.finite(NES_x) & is.finite(NES_y)]
  sig_both <- !is.na(m$padj_x) & !is.na(m$padj_y) & m$padj_x < PADJ_CUTOFF & m$padj_y < PADJ_CUTOFF
  m[, quad := fifelse(sig_both & NES_x > 0 & NES_y > 0, "Concordant Up",
              fifelse(sig_both & NES_x < 0 & NES_y < 0, "Concordant Down",
              fifelse(sig_both, "Discordant", "Single / NS")))]
  m[, `:=`(x_contrast = k_x, y_contrast = k_y)]
  m
}

plot_nes_scatter <- function(k_x, k_y, title, n_label = 5, min_sep = 1.0) {
  m    <- nes_pair(k_x, k_y)
  rho  <- cor(m$NES_x, m$NES_y, method = "spearman")
  cand <- m[quad %in% c("Concordant Up", "Concordant Down") & main_x & main_y][order(-(abs(NES_x) + abs(NES_y)))]
  keep <- integer(0)
  for (i in seq_len(nrow(cand))) {
    if (length(keep) >= n_label) break
    d <- sqrt((cand$NES_x[i] - cand$NES_x[keep])^2 + (cand$NES_y[i] - cand$NES_y[keep])^2)
    if (!length(keep) || all(d > min_sep)) keep <- c(keep, i)
  }
  labs_dt <- cand[keep]
  labs_dt[, `:=`(lab = wrap_term(term, 24), nx = fifelse(NES_x > 0, 0.5, -0.5), ny = fifelse(NES_y > 0, 0.9, -0.9))]
  g <- ggplot(m, aes(NES_x, NES_y, colour = quad)) +
    geom_point(data = m[quad == "Single / NS"], size = 0.7, alpha = 0.4) +
    geom_point(data = m[quad != "Single / NS"], size = 1.1, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey45") +
    scale_colour_manual(values = QUAD_COLS) +
    labs(title = title,
         subtitle = sprintf("GO:BP NES | Spearman rho = %.3f (n = %d sets) | %d concordant, %d discordant (padj < %.2f in both)",
                            rho, nrow(m), sum(grepl("^Concordant", m$quad)), sum(m$quad == "Discordant"), PADJ_CUTOFF),
         x = paste0(contrast_specs[[k_x]]$label, " (NES)"), y = paste0(contrast_specs[[k_y]]$label, " (NES)"), colour = NULL) +
    theme(legend.position = "bottom", plot.subtitle = element_text(size = 7))
  if (nrow(labs_dt)) {
    g <- g + geom_label_repel(data = labs_dt, aes(label = lab), size = 2.1, lineheight = 0.9,
                              nudge_x = labs_dt$nx, nudge_y = labs_dt$ny, fill = alpha("white", 0.85),
                              label.size = 0, label.padding = 0.15, max.overlaps = Inf, box.padding = 0.6,
                              point.padding = 0.3, force = 4, min.segment.length = 0, segment.colour = "grey50", show.legend = FALSE) +
      scale_x_continuous(expand = expansion(mult = 0.18)) + scale_y_continuous(expand = expansion(mult = 0.18))
  }
  g
}

FIG3_PAIRS <- list(
  list("LMPP_vs_HSC", "M1_LSC_vs_pHSC",  "M1 vs LMPP priming"),
  list("GMP_vs_HSC",  "M1_LSC_vs_pHSC",  "M1 vs GMP commitment"),
  list("GMP_vs_LMPP", "M2_Blast_vs_LSC", "M2 vs LMPP->GMP"),
  list("Mono_vs_GMP", "M2_Blast_vs_LSC", "M2 vs GMP->Mono")
)
nes_all <- rbindlist(lapply(FIG3_PAIRS, function(p) nes_pair(p[[1]], p[[2]])[, pair := p[[3]]]))
fwrite(nes_all[, .(pair, x_contrast, y_contrast, pathway, term, NES_x, padj_x, NES_y, padj_y, main_x, main_y, quad)],
       file.path(DE_DIR, "pathway_crossmodel_nes.tsv"), sep = "\t")
nes_summary <- nes_all[, .(n_sets = .N, rho_NES = cor(NES_x, NES_y, method = "spearman"),
                           n_concordant = sum(grepl("^Concordant", quad)), n_discordant = sum(quad == "Discordant")), by = pair]
fwrite(nes_summary, file.path(DE_DIR, "pathway_crossmodel_summary.tsv"), sep = "\t")

p_h1 <- ggplot(spec_H1, aes(factor(N), AUC_LMPP_like, colour = direction)) +
  geom_pointrange(aes(ymin = AUC_lo, ymax = AUC_hi), position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey45") +
  scale_colour_manual(values = c(Up = "#B2182B", Down = "#2166AC")) +
  labs(title = "H1: M1 vs LMPP-specific or GMP-specific genes (equal size)",
       subtitle = "AUC > 0.5 = M1 more LMPP-like; gene bootstrap 95% CI (descriptive)",
       x = "Top-N genes per normal contrast", y = "AUC (LMPP-like)", colour = "Axis direction") +
  theme(legend.position = "bottom")

h3_plot_dt <- copy(h3_matched[filter == "All tested genes"])
h3_plot_dt[, `:=`(OR_plot = pmax(OR, 0.05), lo_plot = pmax(OR_lo, 0.05), hi_plot = pmin(OR_hi, 100), zero_or = OR == 0)]
p_h3 <- ggplot(h3_plot_dt, aes(factor(N), OR_plot, colour = route, shape = zero_or)) +
  geom_pointrange(aes(ymin = lo_plot, ymax = hi_plot), position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
  scale_y_log10() +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 25), guide = "none") +
  scale_colour_manual(values = c("HSC>GMP" = "#D95F02", "HSC>LMPP" = "#1B9E77")) +
  labs(title = "H3: overlap of HSC>X (top-N) with LSC>Blast (CQN + apeglm)",
       subtitle = "OR with 95% CI (log scale, descriptive; triangle = OR 0 floored at 0.05)",
       x = "Top-N genes in HSC>X set", y = "Odds ratio", colour = NULL) +
  theme(legend.position = "bottom")

p_reg <- ggplot(reg_tab, aes(factor(term, levels = c("LMPP_vs_HSC", "GMP_vs_HSC")), std_beta, fill = term)) +
  geom_col(width = 0.55) + geom_hline(yintercept = 0, colour = "grey45") +
  scale_x_discrete(labels = c(LMPP_vs_HSC = "LMPP vs HSC", GMP_vs_HSC = "GMP vs HSC")) +
  scale_fill_manual(values = c(LMPP_vs_HSC = "#1B9E77", GMP_vs_HSC = "#D95F02")) +
  labs(title = "H1: joint rank regression of M1 Wald stat",
       subtitle = sprintf("R2 = %.2f; descriptive only (genes not independent)", reg_tab$R2[1]),
       x = NULL, y = "Standardised beta") + theme(legend.position = "none")

p_fisher <- ggplot(overlap_stats, aes(reorder(comparison, odds_ratio), odds_ratio, fill = filter)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_text(aes(label = sprintf("padj=%.1e", padj)), position = position_dodge(width = 0.7), hjust = -0.08, size = 2.5) +
  coord_flip() +
  scale_fill_manual(values = c("All tested genes" = "#4575B4", "No cell-cycle/ribo" = "#D73027")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(title = "Directional overlap Fisher's test (descriptive)",
       subtitle = sprintf("Universe: %d common (%d non-prolif); genes correlated", length(common_univ), length(non_prolif_univ)),
       x = NULL, y = "Odds ratio (one-sided Fisher)", fill = NULL) +
  theme(legend.position = "bottom")

fig3 <- (plot_lfc_scatter("LMPP_vs_HSC", "M1_LSC_vs_pHSC", "M1 vs LMPP priming") |
         plot_lfc_scatter("GMP_vs_HSC",  "M1_LSC_vs_pHSC", "M1 vs GMP commitment") |
         plot_lfc_scatter("GMP_vs_LMPP", "M2_Blast_vs_LSC", "M2 vs LMPP->GMP") |
         plot_lfc_scatter("Mono_vs_GMP", "M2_Blast_vs_LSC", "M2 vs GMP->Mono")) /
        (plot_nes_scatter(FIG3_PAIRS[[1]][[1]], FIG3_PAIRS[[1]][[2]], "M1 vs LMPP priming (pathways)") |
         plot_nes_scatter(FIG3_PAIRS[[2]][[1]], FIG3_PAIRS[[2]][[2]], "M1 vs GMP commitment (pathways)") |
         plot_nes_scatter(FIG3_PAIRS[[3]][[1]], FIG3_PAIRS[[3]][[2]], "M2 vs LMPP->GMP (pathways)", n_label = 4, min_sep = 0.3) |
         plot_nes_scatter(FIG3_PAIRS[[4]][[1]], FIG3_PAIRS[[4]][[2]], "M2 vs GMP->Mono (pathways)", n_label = 4, min_sep = 0.3)) /
        (p_h1 | p_h3 | p_reg | p_fisher) +
  plot_layout(heights = c(1, 1, 0.85)) +
  plot_annotation(
    title = "Fig 3. AML progression (Model A) against normal differentiation (Models C, C2), NextSeq 500 + CQN + apeglm",
    caption = paste("A-D gene level; E-H pathway level; I-K size-matched specificity (threshold-free); L threshold-based overlap.",
                    "AML-derived and normal-derived contrasts share no samples. Gene-level CIs and p-values are descriptive (genes correlated)."),
    tag_levels = "A"
  )
ggsave(file.path(DE_DIR, "Fig3_cross_model_integrated.png"), fig3, width = 24, height = 18, dpi = 150)

# ---- 4. Model A leave-one-donor-out (LODO) sensitivity check -----------------

message("Running Model A leave-one-donor-out (LODO) sensitivity check...")
lodo_tab <- rbindlist(lapply(c("ALL_6", trios6), function(drop_d) {
  keep_d <- if (drop_d == "ALL_6") trios6 else setdiff(trios6, drop_d)
  dds_l  <- if (drop_d == "ALL_6") dds_A else fit_deseq_model(
    ns_dds, ns_dds$status == "AML" & ns_dds$donor %in% keep_d, ~ donor + cell_type, paste0("LODO drop ", drop_d)
  )
  m1_l <- if (drop_d == "ALL_6") de_res$M1_LSC_vs_pHSC else extract_contrast(dds_l, "M1", "LSC", "pHSC")
  m2_l <- if (drop_d == "ALL_6") de_res$M2_Blast_vs_LSC else extract_contrast(dds_l, "M2", "Blast", "LSC")

  u    <- intersect(common_univ, intersect(m1_l[!is.na(padj), gene], m2_l[!is.na(padj), gene]))
  s_m1 <- setNames(m1_l[match(u, gene), stat], u)
  s_lh <- S_LMPP[u]; s_gh <- S_GMP[u]

  L_up <- top_genes(s_lh, 300, "Up"); G_up <- top_genes(s_gh, 300, "Up")
  L_dn <- top_genes(s_lh, 300, "Down"); G_dn <- top_genes(s_gh, 300, "Down")
  betas <- coef(lm(rk(s_m1, u) ~ rk(s_lh, u) + rk(s_gh, u)))[-1]
  lsc_gt_blast_l <- intersect(m2_l[direction == "Down", gene], u)

  data.table(
    dropped       = drop_d,
    M1_up         = sum(m1_l$direction == "Up"),
    M1_dn         = sum(m1_l$direction == "Down"),
    M2_up         = sum(m2_l$direction == "Up"),
    M2_dn         = sum(m2_l$direction == "Down"),
    r_M1_LMPP     = cor(m1_l[match(u, gene), log2FC_shrunk], de_res$LMPP_vs_HSC[match(u, gene), log2FC_shrunk]),
    r_M1_GMP      = cor(m1_l[match(u, gene), log2FC_shrunk], de_res$GMP_vs_HSC[match(u, gene), log2FC_shrunk]),
    r_M2_GMPvLMPP = cor(m2_l[match(u, gene), log2FC_shrunk], de_res$GMP_vs_LMPP[match(u, gene), log2FC_shrunk]),
    r_M2_MonovGMP = cor(m2_l[match(u, gene), log2FC_shrunk], de_res$Mono_vs_GMP[match(u, gene), log2FC_shrunk]),
    H1_AUC_Up_300 = auc_rank(s_m1[setdiff(L_up, G_up)], s_m1[setdiff(G_up, L_up)]),
    H1_AUC_Dn_300 = auc_rank(s_m1[setdiff(G_dn, L_dn)], s_m1[setdiff(L_dn, G_dn)]),
    H1_beta_LMPP  = unname(betas[1]),
    H1_beta_GMP   = unname(betas[2]),
    H3_OR_GMP_300 = fisher_ci(top_genes(s_gh, 300, "Down"), lsc_gt_blast_l, u)$OR,
    H3_OR_LMPP_300= fisher_ci(top_genes(s_lh, 300, "Down"), lsc_gt_blast_l, u)$OR
  )
}))
fwrite(lodo_tab, file.path(DE_DIR, "modelA_lodo_sensitivity.tsv"), sep = "\t")

lodo_lvls <- c("ALL_6", trios6)
plot_lodo <- function(cols, title, ylab, hlines) {
  m <- melt(lodo_tab, id.vars = "dropped", measure.vars = cols)
  m[, dropped := factor(dropped, levels = lodo_lvls)]
  g <- ggplot(m, aes(dropped, value, colour = variable, group = variable)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.4) +
    geom_hline(yintercept = hlines, linetype = "dashed", colour = "grey55") +
    labs(title = title, x = "Dropped AML patient", y = ylab, colour = NULL) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))
  g
}
p_lodo_cor <- plot_lodo(c("r_M1_LMPP", "r_M1_GMP", "r_M2_GMPvLMPP", "r_M2_MonovGMP"),
                        "Model A LODO: cross-model gene LFC correlation", "Pearson r", 0)
p_lodo_h1  <- plot_lodo(c("H1_AUC_Up_300", "H1_AUC_Dn_300", "H1_beta_LMPP", "H1_beta_GMP"),
                        "Model A LODO: H1 specificity (N = 300) & joint rank beta", "Metric value", c(0, 0.5))
lodo_h3_m  <- melt(lodo_tab, id.vars = "dropped", measure.vars = c("H3_OR_GMP_300", "H3_OR_LMPP_300"))
lodo_h3_m[, dropped := factor(dropped, levels = lodo_lvls)]
p_lodo_h3  <- ggplot(lodo_h3_m, aes(dropped, pmax(value, 0.1), colour = variable, group = variable)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2.4) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  scale_y_log10() +
  labs(title = "Model A LODO: H3 matched-N=300 odds ratio", x = "Dropped AML patient", y = "Odds ratio (log scale)", colour = NULL) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(DE_DIR, "Fig4_modelA_lodo_sensitivity.png"), p_lodo_cor | p_lodo_h1 | p_lodo_h3, width = 16, height = 5, dpi = 150)

saveRDS(
  list(models = list(A = dds_A, B = dds_B, C = dds_C, C2 = dds_C2),
       de = de_res, enrichment = enrich_res, overlap_stats = overlap_stats,
       shared_genes = core_stem, matched = list(H1 = spec_H1, H1_regression = reg_tab, H3 = h3_matched),
       pathway_crossmodel = nes_summary, lodo = lodo_tab),
  file.path(DE_DIR, "de_results.rds")
)
writeLines(capture.output(sessionInfo()), file.path(DE_DIR, "sessionInfo.txt"))
message("Module 2 DE & integration outputs written to ", DE_DIR)

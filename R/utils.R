# Shared constants and pure statistical / I/O helpers for the GSE74246 AML CQN-DESeq2 pipeline.

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

GSE         <- "GSE74246"
MIN_COUNT   <- 10L
MIN_SAMPLES <- 3L
N_TOP_VAR   <- 500L
MAD_CUTOFF  <- 3
EXCLUDE     <- character(0)

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1

CELL_LEVELS <- c(
  "HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "CLP",
  "Mono", "Ery", "Bcell", "CD4Tcell", "CD8Tcell", "NKcell",
  "pHSC", "LSC", "Blast"
)
MATURE    <- c("Mono", "Ery", "Bcell", "CD4Tcell", "CD8Tcell", "NKcell")
PLATFORMS <- c(GPL11154 = "HiSeq2000", GPL18573 = "NextSeq500")

MARKERS <- c(
  HLF = "HSC", CRHBP = "HSC", AVP = "HSC", MECOM = "HSC",
  CD34 = "HSPC", PROM1 = "HSPC", FLT3 = "MPP/LMPP", CD38 = "Progenitor",
  DNTT = "CLP", MPO = "GMP", ELANE = "GMP",
  GATA1 = "MEP/Ery", GYPA = "Ery", HBB = "Ery",
  CD14 = "Mono", LYZ = "Mono", MS4A1 = "B", CD79A = "B",
  CD3E = "T", CD4 = "CD4 T", CD8A = "CD8 T", NCAM1 = "NK", GNLY = "NK"
)

DIR_COLS <- c(Up = "#B2182B", Down = "#2166AC", NS = "grey80")
QUAD_COLS <- c(
  "Concordant Up"   = "#B2182B",
  "Concordant Down" = "#2166AC",
  "Discordant"      = "#7570B3",
  "Single / NS"     = "grey82"
)

get_cell_palette <- function() {
  stats::setNames(c(RColorBrewer::brewer.pal(12, "Paired"), RColorBrewer::brewer.pal(4, "Dark2")), CELL_LEVELS)
}

find_root <- function(marker = file.path("data", GSE)) {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  candidates <- normalizePath(c(getwd(), dirname(script)), mustWork = FALSE)
  for (d in candidates) {
    repeat {
      if (dir.exists(file.path(d, marker))) return(d)
      if (dirname(d) == d) break
      d <- dirname(d)
    }
  }
  stop("Cannot find '", marker, "'. Run from the repository root or generate data first.", call. = FALSE)
}

read_counts <- function(path) {
  dt <- data.table::fread(path)
  m  <- as.matrix(dt[, -1])
  rownames(m) <- dt[[1]]
  if (anyNA(m))                   stop("Counts contain NA", call. = FALSE)
  if (any(m < 0))                 stop("Counts contain negative values", call. = FALSE)
  if (any(m != round(m)))         stop("Counts are not integers; DESeq2 needs raw counts", call. = FALSE)
  if (anyDuplicated(rownames(m))) stop("Duplicate gene IDs in count matrix", call. = FALSE)
  storage.mode(m) <- "integer"
  m
}

low_mad <- function(x, cutoff = MAD_CUTOFF) {
  x < stats::median(x) - cutoff * stats::mad(x)
}

r2_spline <- function(y, g, l) {
  list(
    r_gc           = stats::cor(y, g),
    r_loglen       = stats::cor(y, l),
    R2_gc_spline   = summary(stats::lm(y ~ splines::ns(g, df = 5)))$r.squared,
    R2_len_spline  = summary(stats::lm(y ~ splines::ns(l, df = 5)))$r.squared,
    R2_gc_plus_len = summary(stats::lm(y ~ splines::ns(g, df = 5) + splines::ns(l, df = 5)))$r.squared
  )
}

wrap_term <- function(x, width = 52) {
  vapply(x, function(s) {
    if (is.na(s) || !nzchar(s)) return("")
    paste(strwrap(s, width = width), collapse = "\n")
  }, character(1))
}

rank_volcano_degs <- function(de_dt) {
  df <- data.table::copy(de_dt[!is.na(padj)])
  df[, neg_log10_padj := -log10(pmax(padj, 1e-300))]
  df[, score := abs(log2FC_shrunk) * neg_log10_padj]
  df
}

top_genes <- function(stat_vec, n, dir = c("Up", "Down")) {
  dir <- match.arg(dir)
  names(utils::head(sort(stat_vec, decreasing = (dir == "Up")), n))
}

auc_rank <- function(x, y) {
  r <- rank(c(x, y))
  (sum(r[seq_along(x)]) - length(x) * (length(x) + 1) / 2) / (length(x) * length(y))
}

rk <- function(stat_vec, univ) {
  as.numeric(scale(rank(stat_vec[univ])))
}

fisher_overlap <- function(name, set1, set2, universe, set1_label, set2_label, filter_label) {
  s1  <- intersect(set1, universe)
  s2  <- intersect(set2, universe)
  n11 <- length(intersect(s1, s2))
  n10 <- length(setdiff(s1, s2))
  n01 <- length(setdiff(s2, s1))
  n00 <- length(universe) - n11 - n10 - n01
  ft  <- stats::fisher.test(matrix(c(n11, n10, n01, n00), nrow = 2), alternative = "greater")
  data.table::data.table(
    comparison = name,
    filter     = filter_label,
    set1       = set1_label,
    set2       = set2_label,
    n_universe = length(universe),
    n_set1     = length(s1),
    n_set2     = length(s2),
    n_overlap  = n11,
    odds_ratio = unname(ft$estimate),
    p_value    = ft$p.value
  )
}

fisher_ci <- function(s1, s2, U) {
  s1 <- intersect(s1, U)
  s2 <- intersect(s2, U)
  a  <- length(intersect(s1, s2))
  b  <- length(s1) - a
  c_ <- length(s2) - a
  d  <- length(U) - a - b - c_
  m  <- matrix(c(a, b, c_, d), 2)
  ft2 <- stats::fisher.test(m)
  list(
    n1                    = length(s1),
    n2                    = length(s2),
    n_overlap             = a,
    OR                    = unname(ft2$estimate),
    OR_lo                 = ft2$conf.int[1],
    OR_hi                 = ft2$conf.int[2],
    p_greater_descriptive = stats::fisher.test(m, alternative = "greater")$p.value
  )
}

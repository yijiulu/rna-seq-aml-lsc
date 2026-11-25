#!/usr/bin/env Rscript
# Unit tests for pure R statistical/I/O utilities and syntax validation across all R scripts.

suppressPackageStartupMessages(library(data.table))

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
test_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
repo_root <- normalizePath(file.path(test_dir, ".."))

# 1. Verify every .R file across R/, scripts/, and tests/ parses without syntax errors
r_files <- sort(c(
  list.files(file.path(repo_root, "R"), pattern = "\\.R$", full.names = TRUE),
  list.files(file.path(repo_root, "scripts"), pattern = "\\.R$", full.names = TRUE)
))
stopifnot(length(r_files) >= 4L)
for (f in r_files) {
  parse(file = f)
}

source(file.path(repo_root, "R", "utils.R"))

# 2. Rank-sum AUC exact values
stopifnot(
  abs(auc_rank(c(10, 20, 30), c(1, 2, 3)) - 1.0) < 1e-12,
  abs(auc_rank(c(1, 2, 3), c(10, 20, 30)) - 0.0) < 1e-12,
  abs(auc_rank(c(1, 4), c(2, 3)) - 0.5) < 1e-12
)

# 3. top_genes selection and rank standardisation
s <- c(G1 = 5.0, G2 = -4.0, G3 = 2.0, G4 = -1.0, G5 = 0.5)
stopifnot(
  identical(top_genes(s, 2, "Up"), c("G1", "G3")),
  identical(top_genes(s, 2, "Down"), c("G2", "G4"))
)
r_std <- rk(s, names(s))
stopifnot(abs(mean(r_std)) < 1e-12, abs(sd(r_std) - 1.0) < 1e-12)

# 4. Spline GC bias diagnostic recovers synthetic GC curve
set.seed(42)
g <- seq(0.3, 0.7, length.out = 200)
l <- runif(200, 8, 13)
y <- -6 * (g - 0.5) + rnorm(200, sd = 0.1)
diag_res <- r2_spline(y, g, l)
stopifnot(diag_res$r_gc < -0.95, diag_res$R2_gc_spline > 0.90)

# 5. Fisher overlap and CI helpers
univ <- paste0("G", seq_len(100))
s1 <- paste0("G", 1:20)
s2 <- paste0("G", 1:15)
fo <- fisher_overlap("test", s1, s2, univ, "S1", "S2", "All")
fc <- fisher_ci(s1, s2, univ)
stopifnot(
  fo$n_overlap == 15L,
  fo$odds_ratio > 10,
  fc$n_overlap == 15L,
  fc$OR_lo > 1
)

# 6. low_mad outlier detection & wrap_term
vals <- c(rep(100, 10), 10)
stopifnot(low_mad(vals)[11], !any(low_mad(vals)[1:10]))
stopifnot(grepl("\n", wrap_term("very long gene ontology biological process term description", width = 20)))

# 7. read_counts input validation (valid integer table vs malformed tables)
tmp_ok <- tempfile(fileext = ".tsv")
tmp_bad <- tempfile(fileext = ".tsv")
on.exit(unlink(c(tmp_ok, tmp_bad)), add = TRUE)

writeLines(c("gene\tS1\tS2", "HLF\t10\t20", "MPO\t0\t50"), tmp_ok)
m_ok <- read_counts(tmp_ok)
stopifnot(is.integer(m_ok), identical(rownames(m_ok), c("HLF", "MPO")), m_ok["MPO", "S2"] == 50L)

writeLines(c("gene\tS1\tS2", "HLF\t10.5\t20"), tmp_bad)
caught_nonint <- tryCatch({ read_counts(tmp_bad); FALSE }, error = function(e) grepl("integers", conditionMessage(e)))
stopifnot(isTRUE(caught_nonint))

cat("[OK] All 7 R utility and syntax unit tests passed.\n")

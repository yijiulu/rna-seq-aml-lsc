#!/usr/bin/env bash
# End-to-end reproducible runner for the GSE74246 AML CQN-DESeq2 transcriptomics pipeline.
# Usage:
#   bash scripts/run_pipeline.sh [--skip-fetch]
set -euo pipefail

SKIP_FETCH=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-fetch) SKIP_FETCH=1; shift ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

PYTHON=${PYTHON:-python3}
RSCRIPT=${RSCRIPT:-Rscript}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"

if [[ "$SKIP_FETCH" -eq 0 ]]; then
  echo "=== [1/5] Fetching GSE74246 counts and metadata from NCBI GEO ==="
  "$PYTHON" scripts/01_fetch_geo.py GSE74246
  echo "=== [2/5] Building UCSC hg19 refGene GC-content and transcript length table ==="
  "$PYTHON" scripts/02_build_gene_gc.py GSE74246
fi

echo "=== [3/5] Running Module 1: QC, platform check & CQN GC-bias diagnostics ==="
"$RSCRIPT" scripts/03_qc_and_cqn.R

echo "=== [4/5] Running Module 2: CQN + DESeq2 + apeglm, fGSEA/ORA, integration & LODO ==="
"$RSCRIPT" scripts/04_de_integration_lodo.R

echo "=== [5/5] Assembling composite QC figures (Fig 1 & Fig S1) ==="
"$RSCRIPT" scripts/05_make_composite_figs.R

echo "[OK] Full GSE74246 pipeline completed cleanly."

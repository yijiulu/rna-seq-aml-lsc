#!/usr/bin/env python3
"""Build hg19 UCSC refGene per-gene transcript GC content and length table for GSE74246 CQN.

Corces et al. (2016, PMID: 27526324) quantified GSE74246 RNA-seq counts with HTSeq against the
UCSC hg19 refGene transcriptome and applied Conditional Quantile Normalisation (CQN) in DESeq2.
This script reproduces per-gene spliced mRNA GC fraction and transcript length from UCSC hg19
`refGene.txt.gz` and `refMrna.fa.gz`, resolving legacy 2015 HGNC symbols via local `org.Hs.eg.db`.

Usage:
  python scripts/02_build_gene_gc.py [GSE74246] [--cache-dir /path/to/ucsc_hg19]

Output:
  data/<GSE>/gene_gc_length.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import statistics
import subprocess
import tempfile
import urllib.request
from collections import defaultdict
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
REFGENE_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"
REFMRNA_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz"


def ensure_file(path: Path, url: str) -> Path:
    if path.exists() and path.stat().st_size > 0:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(path)
    return path


def load_counts_genes(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        next(fh)
        return [line.split("\t", 1)[0] for line in fh]


def load_refgene_transcripts(path: Path) -> dict[str, set[str]]:
    sym_to_accs: dict[str, set[str]] = defaultdict(set)
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 13:
                acc = parts[1].split(".")[0]
                sym = parts[12]
                sym_to_accs[sym].add(acc)
    return sym_to_accs


def resolve_aliases_via_r(counts_path: Path, refgene_path: Path) -> dict[str, str]:
    """Use local Bioconductor org.Hs.eg.db to resolve legacy hg19 symbols unambiguously."""
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as tf:
        out_tsv = Path(tf.name)
    r_code = f"""
suppressPackageStartupMessages({{library(org.Hs.eg.db); library(data.table)}})
cnt_genes <- fread("{counts_path}", select = 1)[[1]]
rg <- fread("{refgene_path}", header = FALSE, select = c(2, 13))
rg_syms <- unique(rg[[2]])
missing <- cnt_genes[!(cnt_genes %in% rg_syms)]
al <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = missing, keytype = "ALIAS", columns = c("ENTREZID", "SYMBOL")))
al_dt <- unique(as.data.table(al)[!is.na(SYMBOL) & SYMBOL %in% rg_syms, .(old = ALIAS, new = SYMBOL)])
al_dt[, in_cnt := new %in% cnt_genes]
n1 <- al_dt[, .N, by = old]
m1 <- al_dt[old %in% n1[N == 1, old], .(old, new)]
rem <- al_dt[!old %in% m1$old & !in_cnt]
n2 <- rem[, .N, by = old]
m2 <- rem[old %in% n2[N == 1, old], .(old, new)]
fwrite(rbind(m1, m2), "{out_tsv}", sep = "\\t")
"""
    try:
        subprocess.run(["Rscript", "-e", r_code], check=True)
        with open(out_tsv, encoding="utf-8") as fh:
            return {r["old"]: r["new"] for r in csv.DictReader(fh, delimiter="\t")}
    finally:
        if out_tsv.exists():
            out_tsv.unlink()


def parse_refmrna_gc_len(path: Path, needed_accs: set[str]) -> dict[str, tuple[float, int]]:
    stats: dict[str, tuple[float, int]] = {}
    cur_acc = None
    gc = 0
    acgt = 0

    def flush() -> None:
        nonlocal cur_acc, gc, acgt
        if cur_acc is not None and acgt > 0:
            stats[cur_acc] = (gc / acgt, acgt)

    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                header = line[1:].split()[0].split(".")[0]
                cur_acc = header if header in needed_accs else None
                gc = 0
                acgt = 0
            elif cur_acc is not None:
                seq = line.strip().upper()
                g = seq.count("G")
                c = seq.count("C")
                a = seq.count("A")
                t = seq.count("T")
                gc += g + c
                acgt += a + c + g + t
    flush()
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gse", nargs="?", default="GSE74246")
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Optional directory containing cached refGene.txt.gz and refMrna.fa.gz.",
    )
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR, help="Base data directory (default: data/).")
    args = parser.parse_args()

    gse_dir = args.data_dir / args.gse.upper()
    counts_path = gse_dir / "counts.tsv.gz"
    if not counts_path.exists():
        raise FileNotFoundError(f"Missing {counts_path}; run scripts/01_fetch_geo.py first.")

    with tempfile.TemporaryDirectory() as tmpdir:
        cache = args.cache_dir if args.cache_dir is not None else Path(tmpdir)
        refgene_path = ensure_file(cache / "refGene.txt.gz", REFGENE_URL)
        refmrna_path = ensure_file(cache / "refMrna.fa.gz", REFMRNA_URL)

        genes = load_counts_genes(counts_path)
        sym_to_accs = load_refgene_transcripts(refgene_path)
        alias_map = resolve_aliases_via_r(counts_path, refgene_path)

        gene_to_target: dict[str, tuple[str, str]] = {}
        needed_accs: set[str] = set()
        for g in genes:
            if g in sym_to_accs:
                gene_to_target[g] = (g, "direct")
                needed_accs.update(sym_to_accs[g])
            elif g in alias_map and alias_map[g] in sym_to_accs:
                tgt = alias_map[g]
                gene_to_target[g] = (tgt, "alias")
                needed_accs.update(sym_to_accs[tgt])

        acc_stats = parse_refmrna_gc_len(refmrna_path, needed_accs)

    out_path = gse_dir / "gene_gc_length.tsv"
    tmp_path = out_path.with_suffix(".tsv.part")
    rows = []
    for g in genes:
        if g not in gene_to_target:
            continue
        tgt, mtype = gene_to_target[g]
        accs = [a for a in sym_to_accs[tgt] if a in acc_stats]
        if not accs:
            continue
        nm_accs = [a for a in accs if a.startswith("NM_")]
        use_accs = nm_accs if nm_accs else accs
        gcs = [acc_stats[a][0] for a in use_accs]
        lens = [acc_stats[a][1] for a in use_accs]
        longest_idx = max(range(len(use_accs)), key=lambda i: lens[i])
        rows.append((
            g,
            tgt,
            mtype,
            len(use_accs),
            f"{statistics.median(gcs):.6f}",
            f"{statistics.median(lens):.1f}",
            f"{gcs[longest_idx]:.6f}",
            str(lens[longest_idx]),
        ))

    with open(tmp_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "gene", "refgene_symbol", "mapping_type", "n_isoforms",
            "gc_median", "length_median", "gc_longest", "length_longest",
        ])
        writer.writerows(rows)
    tmp_path.replace(out_path)
    print(f"Wrote {len(rows)}/{len(genes)} genes ({100 * len(rows) / len(genes):.2f}%) to {out_path}")


if __name__ == "__main__":
    main()

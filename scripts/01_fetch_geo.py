#!/usr/bin/env python3
"""Fetch a GEO RNA-seq series (processed counts + multi-platform metadata) for downstream analysis.

Usage:
  python scripts/01_fetch_geo.py [GSE74246] [--suppl FILE] [--outdir data]

Output (in <outdir>/<GSE>/):
  counts.tsv.gz  genes x samples integer count matrix
  samples.tsv    sample metadata aligned 1:1 with counts columns
  manifest.json  source URLs, MD5 checksums, retrieval timestamp, and dimensions
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import requests

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

# Some GEO submitters label count-matrix columns differently from !Sample_title.
# Each regex rule harmonises count-column headers to match GEO sample titles exactly.
COLUMN_ALIASES: Dict[str, List[Tuple[str, str]]] = {
    "GSE74246": [(r"-rHSC$", "-pHSC"), (r"-Blasts$", "-Blast")],
}

DEFAULT_SUPPL: Dict[str, str] = {
    "GSE74246": "GSE74246_RNAseq_All_Counts.txt.gz",
}


def geo_base(gse: str) -> str:
    """Return the NCBI FTP base URL for a GEO series accession (e.g. GSE74246 -> GSE74nnn)."""
    bucket = re.sub(r"\d{1,3}$", "nnn", gse)
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{gse}"


def list_files(url: str, pattern: str) -> List[str]:
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    return sorted(set(re.findall(rf'href="({pattern})"', resp.text)))


def download_and_hash(url: str, dest: Path) -> Tuple[str, int]:
    """Stream-download url to dest and return (md5_hex, n_bytes)."""
    digest = hashlib.md5()
    nbytes = 0
    with requests.get(url, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        with open(dest, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
                digest.update(chunk)
                nbytes += len(chunk)
    return digest.hexdigest(), nbytes


def parse_series_matrix(path: Path) -> pd.DataFrame:
    """Extract one metadata row per sample from the !Sample_ header lines of a GEO series matrix."""
    key_counts: Counter[str] = Counter()
    meta_dict: Dict[str, List[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_"):
                parts = line.rstrip("\r\n").split("\t")
                raw_key = parts[0][len("!Sample_"):]
                key_counts[raw_key] += 1
                key = raw_key if key_counts[raw_key] == 1 else f"{raw_key}_{key_counts[raw_key]}"
                meta_dict[key] = [x.strip().strip('"') for x in parts[1:]]
    meta = pd.DataFrame(meta_dict)

    # Expand "characteristics_ch1*" key: value columns (e.g. "cell type: HSC") into dedicated fields
    for col in [c for c in meta.columns if c.startswith("characteristics_ch1")]:
        keys = meta[col].str.split(":", n=1).str[0].str.strip().unique()
        if len(keys) == 1 and keys[0]:
            name = re.sub(r"\W+", "_", keys[0]).strip("_").lower()
            meta[name] = meta[col].str.split(":", n=1).str[1].str.strip()
    meta["series_matrix_file"] = path.name
    return meta


def fetch_geo(gse: str, suppl_file: Optional[str] = None) -> Tuple[pd.DataFrame, pd.DataFrame, List[dict]]:
    base = geo_base(gse)

    if suppl_file is None:
        tables = list_files(f"{base}/suppl/", r'[^"/]+\.(?:txt|tsv|csv)\.gz')
        if len(tables) != 1:
            raise ValueError(f"Choose --suppl explicitly; found {tables or 'no tables'} in {base}/suppl/")
        suppl_file = tables[0]

    # Multi-platform series provide one matrix file per GPL platform (e.g. HiSeq2000 + NextSeq500)
    matrix_files = list_files(f"{base}/matrix/", r'[^"/]+_series_matrix\.txt\.gz')
    if not matrix_files:
        raise FileNotFoundError(f"No series matrix found in {base}/matrix/")

    urls = {suppl_file: f"{base}/suppl/{suppl_file}"}
    urls.update({m: f"{base}/matrix/{m}" for m in matrix_files})

    sources: List[dict] = []
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        for name, url in urls.items():
            digest, nbytes = download_and_hash(url, tmp / name)
            sources.append({"file": name, "url": url, "md5": digest, "bytes": nbytes})

        sep = "," if suppl_file.endswith(".csv.gz") else "\t"
        counts = pd.read_csv(tmp / suppl_file, sep=sep, index_col=0)
        meta = pd.concat([parse_series_matrix(tmp / m) for m in matrix_files], ignore_index=True)

    return counts, meta, sources


def harmonise(gse: str, counts: pd.DataFrame, meta: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Match every count column to exactly one GEO sample; fail loudly rather than silently drop samples."""
    renamed: Dict[str, str] = {}
    for col in counts.columns:
        new = col
        for pattern, repl in COLUMN_ALIASES.get(gse, []):
            new = re.sub(pattern, repl, new)
        renamed[col] = new

    titles = meta["title"]
    if titles.duplicated().any():
        raise ValueError(f"Duplicate GEO titles: {titles[titles.duplicated()].tolist()}")
    unmatched_counts = sorted(set(renamed.values()) - set(titles))
    unmatched_geo = sorted(set(titles) - set(renamed.values()))
    if unmatched_counts or unmatched_geo:
        raise ValueError(
            f"Unmatched count columns: {unmatched_counts}\n"
            f"Unmatched GEO samples: {unmatched_geo}\n"
            f"Add a rule to COLUMN_ALIASES['{gse}']."
        )

    counts = counts.rename(columns=renamed)
    counts.index.name = "gene"
    meta = meta.set_index("title").loc[counts.columns].reset_index()
    meta.insert(1, "count_column_original", list(renamed.keys()))
    meta = meta.rename(columns={"title": "sample_id"})
    return counts, meta


def save_geo(gse: str, counts: pd.DataFrame, meta: pd.DataFrame, sources: List[dict], out_dir: Path = DATA_DIR) -> Path:
    if counts.isna().any().any():
        raise ValueError("Counts contain missing values")
    if counts.index.duplicated().any():
        raise ValueError(f"Duplicate gene IDs, e.g. {counts.index[counts.index.duplicated()][:5].tolist()}")
    if not all(pd.api.types.is_integer_dtype(t) for t in counts.dtypes):
        print("Warning: counts are not all integers; DESeq2 expects raw integer counts.")

    gse_dir = out_dir / gse
    gse_dir.mkdir(parents=True, exist_ok=True)

    def atomic_tsv(df: pd.DataFrame, dest: Path, **kw) -> None:
        tmp = dest.with_name(dest.name + ".part")
        df.to_csv(tmp, sep="\t", compression="gzip" if dest.suffix == ".gz" else None, **kw)
        tmp.replace(dest)

    atomic_tsv(counts, gse_dir / "counts.tsv.gz")
    front = [
        c for c in [
            "sample_id", "geo_accession", "source_name_ch1", "cell_type",
            "donorid", "sample_type", "platform_id", "count_column_original",
        ]
        if c in meta.columns
    ]
    atomic_tsv(meta[front + [c for c in meta.columns if c not in front]], gse_dir / "samples.tsv", index=False)

    manifest = {
        "accession": gse,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sources": sources,
        "n_genes": int(counts.shape[0]),
        "n_samples": int(counts.shape[1]),
        "outputs": ["counts.tsv.gz", "samples.tsv"],
        "pandas_version": pd.__version__,
    }
    (gse_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return gse_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gse", nargs="?", default="GSE74246", help="GEO series accession (default: GSE74246).")
    parser.add_argument("--suppl", default=None, help="Explicit supplementary count-table filename.")
    parser.add_argument("--outdir", type=Path, default=DATA_DIR, help="Output base directory (default: data/).")
    args = parser.parse_args()

    gse = args.gse.strip().upper()
    suppl = args.suppl or DEFAULT_SUPPL.get(gse)

    counts, meta, sources = fetch_geo(gse, suppl_file=suppl)
    counts, meta = harmonise(gse, counts, meta)
    out = save_geo(gse, counts, meta, sources, out_dir=args.outdir)

    print(f"Saved {counts.shape[0]} genes x {counts.shape[1]} samples to {out}")
    group_cols = [c for c in ["source_name_ch1", "platform_id"] if c in meta.columns]
    if group_cols:
        print(meta.groupby(group_cols).size().to_string())


if __name__ == "__main__":
    main()

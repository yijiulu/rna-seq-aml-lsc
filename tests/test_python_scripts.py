"""Unit tests for GEO fetching/harmonisation and UCSC hg19 refGene GC/length parsing."""

from __future__ import annotations

import ast
import gzip
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestAmlPythonScripts(unittest.TestCase):
    def test_all_python_scripts_parse_cleanly(self) -> None:
        py_files = sorted(SCRIPTS_DIR.glob("*.py"))
        self.assertGreaterEqual(len(py_files), 2)
        for p in py_files:
            ast.parse(p.read_text(encoding="utf-8"), filename=str(p))

    def test_refgene_and_refmrna_gc_parsing(self) -> None:
        mod = _load_module("build_gene_gc", SCRIPTS_DIR / "02_build_gene_gc.py")
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            rg = td_path / "refGene.txt.gz"
            fa = td_path / "refMrna.fa.gz"
            cnt = td_path / "counts.tsv.gz"

            with gzip.open(rg, "wt", encoding="utf-8") as fh:
                fh.write("1\tNM_000546.6\tchr17\t-\t1\t10\t1\t10\t1\t1,\t10,\t0\tTP53\tcmpl\tcmpl\t0,\n")
                fh.write("2\tNR_001234.1\tchr1\t+\t1\t10\t1\t10\t1\t1,\t10,\t0\tLINC1\tunk\tunk\t-1,\n")

            with gzip.open(fa, "wt", encoding="utf-8") as fh:
                # 6 GC out of 10 ACGT (N ignored) -> GC = 0.6, length = 10
                fh.write(">NM_000546.6 some description\nGCGCGCATATNN\n")
                fh.write(">NM_999999.1 unused\nAAAAAA\n")

            with gzip.open(cnt, "wt", encoding="utf-8") as fh:
                fh.write("gene\tS1\nTP53\t100\nLINC1\t20\n")

            self.assertEqual(mod.load_counts_genes(cnt), ["TP53", "LINC1"])
            sym_map = mod.load_refgene_transcripts(rg)
            self.assertEqual(sym_map["TP53"], {"NM_000546"})
            self.assertEqual(sym_map["LINC1"], {"NR_001234"})

            stats = mod.parse_refmrna_gc_len(fa, {"NM_000546"})
            self.assertEqual(set(stats.keys()), {"NM_000546"})
            gc_frac, length = stats["NM_000546"]
            self.assertAlmostEqual(gc_frac, 0.6, places=6)
            self.assertEqual(length, 10)

    def test_fetch_geo_helpers_if_pandas_available(self) -> None:
        if importlib.util.find_spec("pandas") is None or importlib.util.find_spec("requests") is None:
            self.skipTest("pandas/requests not installed in base interpreter")
        import pandas as pd

        mod = _load_module("fetch_geo", SCRIPTS_DIR / "01_fetch_geo.py")
        self.assertEqual(mod.geo_base("GSE74246"), "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74246")

        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            mat_path = td_path / "GSE74246-GPL18573_series_matrix.txt.gz"
            with gzip.open(mat_path, "wt", encoding="utf-8") as fh:
                fh.write('!Sample_title\t"SU353-pHSC"\t"SU353-Blast"\n')
                fh.write('!Sample_geo_accession\t"GSM1"\t"GSM2"\n')
                fh.write('!Sample_characteristics_ch1\t"cell type: pHSC"\t"cell type: Blast"\n')
                fh.write('!Sample_characteristics_ch1\t"donorid: SU353"\t"donorid: SU353"\n')
                fh.write("!series_matrix_table_begin\n")

            meta = mod.parse_series_matrix(mat_path)
            self.assertEqual(list(meta["title"]), ["SU353-pHSC", "SU353-Blast"])
            self.assertEqual(list(meta["cell_type"]), ["pHSC", "Blast"])
            self.assertEqual(list(meta["donorid"]), ["SU353", "SU353"])

            # Verify COLUMN_ALIASES (-rHSC -> -pHSC, -Blasts -> -Blast)
            counts = pd.DataFrame({"SU353-rHSC": [10, 20], "SU353-Blasts": [30, 40]}, index=["HLF", "MPO"])
            h_counts, h_meta = mod.harmonise("GSE74246", counts, meta)
            self.assertEqual(list(h_counts.columns), ["SU353-pHSC", "SU353-Blast"])
            self.assertEqual(list(h_meta["sample_id"]), ["SU353-pHSC", "SU353-Blast"])

            # Unmatched column must fail loudly
            bad_counts = pd.DataFrame({"WRONG": [1, 2]}, index=["HLF", "MPO"])
            with self.assertRaises(ValueError):
                mod.harmonise("GSE74246", bad_counts, meta)

            out_dir = mod.save_geo("GSE74246", h_counts, h_meta, [{"file": "test.gz", "md5": "abc", "bytes": 10}], out_dir=td_path)
            manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["n_genes"], 2)
            self.assertEqual(manifest["n_samples"], 2)


if __name__ == "__main__":
    unittest.main()

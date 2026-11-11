#!/usr/bin/env python3
"""Self-contained tests using SYNTHETIC data only. Run: python3 test_respondent_tools.py"""
import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import respondent_core as core  # noqa: E402


def make_dir(n=5, seed=3):
    d = Path(tempfile.mkdtemp(prefix="synthetic_respondents_"))
    subprocess.run([sys.executable, str(HERE / "make_synthetic_data.py"), "--out-dir", str(d), "--n", str(n),
                    "--seed", str(seed)], check=True, capture_output=True)
    return d


class Tests(unittest.TestCase):
    def test_refuses_paths_inside_git(self):
        with self.assertRaises(SystemExit):
            core.assert_outside_git(HERE / "x.csv", "test file")
        with self.assertRaises(SystemExit):
            core.assert_outside_git(core.STAGE_DIR / "local", "test dir")
        core.assert_outside_git(tempfile.gettempdir(), "tmp")      # outside: no error

    def test_construction_and_median_of_products(self):
        d = make_dir(); cfg = core.load_config(d / "config.json")
        r = core.effective_r(core.load_respondents(cfg))
        self.assertTrue(np.allclose(r.r_prep, r.q_prep * (1 - r.gamma_prep)))
        self.assertTrue(np.allclose(r.r_vs, r.q_art * (1 - r.gamma_art)))
        self.assertTrue(np.allclose(r.r_prep_alt, r.q_prep))
        mp, mv = core.median_point(r)
        self.assertAlmostEqual(mp, np.median(r.r_prep)); self.assertAlmostEqual(mv, np.median(r.r_vs))
        pts = core.p_points(r)
        self.assertEqual(pts["P1"], (0.0, mv)); self.assertEqual(pts["P2"], (mp, 0.0)); self.assertEqual(pts["P3"], (mp, mv))

    def test_median_of_products_differs_from_product_of_medians(self):
        df = pd.DataFrame({"q": [0.9, 0.1, 0.5, 0.2, 0.8], "g": [0.9, 0.1, 0.5, 0.2, 0.8]})
        self.assertNotAlmostEqual(np.median(df.q * (1 - df.g)), np.median(df.q) * (1 - np.median(df.g)))

    def test_rejects_bad_data(self):
        d = make_dir(); cfg_p = d / "config.json"
        csv = d / "respondents_synthetic.csv"
        df = pd.read_csv(csv); df.iloc[0, 1] = 150; df.to_csv(csv, index=False)
        with self.assertRaises(SystemExit):
            core.load_respondents(core.load_config(cfg_p))

    def test_one_anonymous_pdf_per_respondent_and_no_other_values(self):
        import fitz
        d = make_dir(n=5, seed=7); out = Path(tempfile.mkdtemp(prefix="synthetic_pdfs_"))
        subprocess.run([sys.executable, str(HERE / "respondent_validation.py"), "--config", str(d / "config.json"),
                        "--out-dir", str(out)], check=True, capture_output=True)
        pdfs = sorted(out.glob("respondent_*.pdf")); self.assertEqual(len(pdfs), 5)
        df = pd.read_csv(d / "respondents_synthetic.csv")
        for pdf, (_, row) in zip(pdfs, df.iterrows()):
            text = fitz.open(pdf)[0].get_text()
            self.assertIn("SYNTHETIC TEST DATA", text)
            self.assertIn("Question for you", text)
            for key in ("P1", "P2", "P3"):
                self.assertIn(key, text)
            own = {int(v) for v in row.iloc[1:5]}
            others = df.drop(row.name)
            for _, o in others.iterrows():          # another respondent's id must not appear
                self.assertNotIn(f"Respondent {o.iloc[0]}:", text)

    def test_refuses_output_inside_git(self):
        d = make_dir()
        res = subprocess.run([sys.executable, str(HERE / "respondent_validation.py"), "--config", str(d / "config.json"),
                              "--out-dir", str(core.STAGE_DIR / "local_pdfs")], capture_output=True, text=True)
        self.assertNotEqual(res.returncode, 0); self.assertIn("REFUSED", res.stderr + res.stdout)

    def test_no_network_imports(self):
        banned = re.compile(r"^\s*(import|from)\s+(requests|urllib|http|httpx|aiohttp|ftplib|smtplib)\b", re.M)
        for f in HERE.glob("*.py"):
            if f.name != Path(__file__).name:
                self.assertIsNone(banned.search(f.read_text()), f.name)


if __name__ == "__main__":
    unittest.main(verbosity=2)

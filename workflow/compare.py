"""File comparison rules used by `make verify`.

A regenerated file is compared with the committed version of the same file:

  identical      byte-for-byte equal
  numeric-match  same structure, every number within rtol/atol (csv, tsv, txt, md, json)
  image-match    same pixel size (within 2%) and mean absolute difference below a threshold (png, jpg); the numeric tables carry the
                 strict checks, images only catch gross changes
  skipped        formats that are not compared (vector graphics, html, binary models/arrays)
  DIFFERENT      anything else; the detail says why
"""
from pathlib import Path
import filecmp
import re

import numpy as np

NUM = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")
SKIP_SUFFIX = {".svg", ".pdf", ".eps", ".html", ".pth", ".pt", ".h5", ".npz", ".pkl", ".gz", ".tar", ".dat", ".log", ".pyc"}
TEXT_SUFFIX = {".txt", ".md", ".json", ".fasta", ".fa", ".tsv", ".csv"}


def _numbers_close(a, b, rtol, atol):
    """Compare two strings token by token: identical non-numeric text, numbers within tolerance."""
    ta, tb = NUM.split(a), NUM.split(b)
    if ta != tb:
        return False, "text differs outside the numbers"
    na = [float(x) for x in NUM.findall(a)]
    nb = [float(x) for x in NUM.findall(b)]
    if len(na) != len(nb):
        return False, "different count of numbers"
    if not na:
        return True, 0.0
    x, y = np.array(na), np.array(nb)
    ok = np.isclose(x, y, rtol=rtol, atol=atol, equal_nan=True)
    if ok.all():
        denom = np.maximum(np.abs(y), atol)
        return True, float(np.max(np.abs(x - y) / denom))
    i = int(np.argmin(ok))
    return False, f"number {i}: {x[i]} vs {y[i]}"


def compare(new: Path, ref: Path, rtol=1e-6, atol=1e-9, png_tol=12.0):
    """Return (status, detail)."""
    if not ref.exists():
        return "new-file", "no committed counterpart"
    if filecmp.cmp(new, ref, shallow=False):
        return "identical", ""
    suf = new.suffix.lower()
    if new.name.lower().endswith(".csv.gz"):          # compressed tables are compared as text (the compressed bytes differ with the gzip level and time stamp)
        import gzip
        ta = gzip.open(new, "rt", errors="replace").read().replace("\r\n", "\n").rstrip()
        tb = gzip.open(ref, "rt", errors="replace").read().replace("\r\n", "\n").rstrip()
        if ta == tb:
            return "identical", "same text, different compression"
        ok, d = _numbers_close(ta, tb, rtol, atol)
        return ("numeric-match", f"max relative diff {d:.1e}") if ok else ("DIFFERENT", d)
    if suf in SKIP_SUFFIX:
        return "skipped", f"{suf} files are not compared"
    if suf in {".png", ".jpg", ".jpeg"}:
        from PIL import Image
        a, b = Image.open(new).convert("L"), Image.open(ref).convert("L")
        if abs(a.width - b.width) > 0.02 * b.width or abs(a.height - b.height) > 0.02 * b.height:
            return "DIFFERENT", f"image size {a.size} vs {b.size}"
        if a.size != b.size:                      # a few pixels of canvas difference come from font / layout builds
            a = a.resize(b.size, Image.LANCZOS)
        d = float(np.abs(np.asarray(a, float) - np.asarray(b, float)).mean())
        if d <= png_tol:
            return "image-match", f"mean abs pixel diff {d:.2f}/255 (rendering noise tolerated up to {png_tol:g})"
        return "DIFFERENT", f"mean abs pixel diff {d:.2f}/255 > {png_tol:g}"
    if suf in TEXT_SUFFIX:
        ta = new.read_text(errors="replace").replace("\r\n", "\n").rstrip()
        tb = ref.read_text(errors="replace").replace("\r\n", "\n").rstrip()
        if ta == tb:
            return "identical", "whitespace only"
        ok, d = _numbers_close(ta, tb, rtol, atol)
        if ok:
            return "numeric-match", f"max relative diff {d:.1e}"
        return "DIFFERENT", d
    return "skipped", f"{suf or 'no suffix'} not compared"

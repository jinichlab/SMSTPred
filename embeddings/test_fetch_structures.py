#!/usr/bin/env python3
"""Offline tests of the sequence check of fetch_structures.py. Run: python test_fetch_structures.py"""
import tempfile
from pathlib import Path

import fetch_structures as f

S = "MKTAYIAKQRQISFVK"


def pdb_text(seq, chain="A"):
    three = {v: k for k, v in f.THREE.items() if k not in ("SEC", "PYL", "MSE")}
    lines = [f"ATOM  {i:5d}  CA  {three[a]} {chain}{i:4d}    {i:8.3f}{0:8.3f}{0:8.3f}  1.00 50.00           C" for i, a in enumerate(seq, 1)]
    return "\n".join(lines) + "\nEND\n"


def test_same_sequence():
    assert f.same_sequence(S, S)
    assert not f.same_sequence(S, S[:-1] + "R")                          # one real substitution (a changed UniProt sequence, same length)
    assert not f.same_sequence(S, S + "A")                               # a longer one
    assert f.same_sequence(S, S[:5] + "X" + S[6:])                       # an unknown residue in UniProt matches any residue of a full-length model
    assert f.same_sequence(S[:5] + S[6:], S[:5] + "X" + S[6:])           # a model that left the unknown residue out (ColabFold)
    assert f.same_sequence(S[:5] + S[7:], S[:5] + "XX" + S[7:])          # two of them
    assert not f.same_sequence(S[:5] + S[6:], S[:5] + "A" + S[6:])       # a model that is simply one residue short is NOT accepted
    assert not f.same_sequence(S[:5] + S[7:], S[:5] + "X" + S[6:])       # shorter by more than the unknown residues


def test_pdb_sequence():
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "x.pdb"
        p.write_text(pdb_text(S))
        assert f.pdb_sequence(p) == S
        p.write_text(pdb_text(S) + pdb_text("GGGG", "B"))              # only the first chain counts
        assert f.pdb_sequence(p) == S


if __name__ == "__main__":
    for t in (test_same_sequence, test_pdb_sequence):
        t()
        print("ok", t.__name__)

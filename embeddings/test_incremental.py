#!/usr/bin/env python3
"""Offline tests of incremental.py (no model, no GPU): which proteins get embedded again, which are kept, and that one file never mixes two pipelines.
Run:  python test_incremental.py   (or pytest)"""
import os
import tempfile
from pathlib import Path

import h5py
import numpy as np

import incremental as inc

S1, S2, S3 = "MKTAYIAKQR", "MKTAYIAKQK", "GGGGSGGGGS"  # S2 differs from S1 at the last residue


def make_h5(path, ids, stamps=None):
    with h5py.File(path, "w") as f:
        for sid in ids:
            ds = f.create_dataset(sid, data=np.zeros(4, dtype=np.float16))
            if stamps and sid in stamps:
                ds.attrs[inc.ATTR] = stamps[sid]


def test_fingerprint_ignores_case_gaps_whitespace():
    assert inc.md5_text("mktay-iak qr") == inc.md5_text(S1)
    assert inc.md5_text(S1) != inc.md5_text(S2)


def test_plan_every_case():
    with tempfile.TemporaryDirectory() as d:
        h5 = Path(d) / "a.h5"
        make_h5(h5, ["same", "edited", "bare_ok", "bare_old", "bare_nolog", "gone"],
                {"same": inc.md5_text(S1), "edited": inc.md5_text(S1)})
        have = inc.stored_fingerprints(h5)
        wanted = {"same": inc.md5_text(S1), "edited": inc.md5_text(S2), "bare_ok": inc.md5_text(S1), "bare_old": inc.md5_text(S2),
                  "bare_nolog": inc.md5_text(S3), "fresh": inc.md5_text(S3)}
        ledger = {"bare_ok": inc.md5_text(S1), "bare_old": inc.md5_text(S1)}  # what the unstamped vectors were made from
        pl = inc.plan(wanted, have, ledger)
        assert pl["todo"] == {"edited": "changed", "bare_old": "changed", "bare_nolog": "unverified", "fresh": "new"}
        assert sorted(pl["keep"]) == ["bare_ok", "same"] and pl["stamp"] == {"bare_ok": inc.md5_text(S1)}
        assert pl["extra"] == ["gone"]
        # without a ledger, an unstamped vector is embedded again, unless the user vouches for the file
        assert inc.plan(wanted, have)["todo"]["bare_ok"] == "unverified"
        trust = inc.plan(wanted, have, trust_existing=True)
        assert "bare_ok" in trust["keep"] and "bare_nolog" in trust["keep"] and trust["stamp"] == {}  # kept, but not stamped as verified
        # a second run after stamping keeps the verified one with no ledger at all
        inc.apply_stamps(h5, pl["stamp"])
        assert "bare_ok" in inc.plan(wanted, inc.stored_fingerprints(h5))["keep"]


def test_nothing_to_do_and_everything_new():
    assert inc.plan({"a": "x"}, {"a": "x"})["todo"] == {}
    assert inc.plan({"a": "x", "b": "y"}, {})["todo"] == {"a": "new", "b": "new"}


def test_ledgers():
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "l.fasta").write_text(f">P1 some name\n{S1[:5]}\n{S1[5:]}\n>Q.2/x\n{S3}\n")
        (d / "l.csv").write_text(f"Entry,Sequence,Other\nP1,{S1},z\nP9,{S2},z\n")
        (d / "pdbs").mkdir()
        (d / "pdbs" / "P1.pdb").write_text("ATOM\n")
        (d / "pdbs" / "Q5_unrelaxed_rank_001_model_3.pdb").write_text("ATOM2\n")
        fa = inc.read_ledger([d / "l.fasta"], norm=lambda s: s.replace("/", "_").replace(".", "_"))
        assert fa == {"P1": inc.md5_text(S1), "Q_2_x": inc.md5_text(S3)}  # wrapped lines joined; the id normalised like the embedder does
        assert inc.read_ledger([d / "l.csv"]) == {"P1": inc.md5_text(S1), "P9": inc.md5_text(S2)}
        pdb = inc.read_ledger([d / "pdbs"])
        assert set(pdb) == {"P1", "Q5"} and pdb["P1"] == inc.md5_file(d / "pdbs" / "P1.pdb")  # ColabFold names reduce to the Entry
        assert set(inc.read_ledger([d / "l.fasta", d / "l.csv"])) == {"P1", "Q.2/x", "P9"}  # several ledgers merge


def test_groups_and_replacement_and_config_guard():
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "g.h5"
        with h5py.File(p, "w") as f:  # ProstT5 layout: a group per protein
            g = f.create_group("P1")
            g.create_dataset("mean", data=np.zeros(3))
            g.attrs[inc.ATTR] = inc.md5_text(S1)
        assert inc.stored_fingerprints(p) == {"P1": inc.md5_text(S1)}
        with h5py.File(p, "a") as f:
            inc.replace_slot(f, "P1")
            f.create_group("P1").create_dataset("mean", data=np.ones(3))
            assert f["P1/mean"][0] == 1 and "P2" not in f
            inc.check_config(f, {"model": "A", "precision": "half"})      # a file without a record gets one ...
            inc.check_config(f, {"model": "A", "precision": "half"})      # ... and the same settings pass
            try:
                inc.check_config(f, {"model": "A", "precision": "full"})  # ... but other settings are refused
            except SystemExit as e:
                assert "mix two pipelines" in str(e)
            else:
                raise AssertionError("different settings must be refused")
            inc.check_config(f, {"model": "A", "precision": "full"}, allow=True)


if __name__ == "__main__":
    for t in (test_fingerprint_ignores_case_gaps_whitespace, test_plan_every_case, test_nothing_to_do_and_everything_new, test_ledgers, test_groups_and_replacement_and_config_guard):
        t()
        print("ok", t.__name__)

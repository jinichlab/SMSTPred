"""
Incremental embedding: embed only what the output file does not hold yet (or holds for a different input), and add it to the file.

Every vector carries a fingerprint of the input it was made from, the HDF5 attribute `input_md5`: the md5 of the (normalised) sequence for ProtT5 / ProstT5, of the
structure file for SaProt / ESM-IF. A run compares the fingerprint of every wanted protein with the one stored in the output file:

  not in the file                          -> new: embed
  stored fingerprint equal                 -> keep (nothing is computed)
  stored fingerprint different             -> changed (UniProt edited the sequence, the structure was re-predicted): embed again and replace the vector
  no stored fingerprint (the archived files, made before this existed):
      --ledger says what it was made from  -> equal: keep and stamp it; different: embed again and replace
      --trust-existing                     -> keep, by accession alone (use only when you know the file belongs to the same inputs)
      otherwise                            -> unverified: embed again (slow but correct)
Vectors in the file for proteins that are not wanted are left alone. The output file also records the model and settings it was made with (file attributes), and a run with
different settings refuses to append: one file never mixes two pipelines.

--ledger is what the vectors already in the file were made from: a FASTA, a CSV / TSV (Entry, Sequence), or a folder of <Entry>.pdb files (for the structure embedders).
Pure Python (h5py, hashlib; pandas only to read a CSV / TSV ledger).
"""
import gzip
import hashlib
import re
import sys
from pathlib import Path

import h5py

ATTR = "input_md5"


def normalize_sequence(seq):
    """The sequence as it is embedded: upper case, no whitespace, no gaps."""
    return re.sub(r"[\s-]+", "", str(seq)).upper()


def md5_text(seq):
    return hashlib.md5(normalize_sequence(seq).encode()).hexdigest()


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_ledger(paths, norm=lambda x: x):
    """{id: fingerprint} of the inputs the vectors already in a file were made from. A FASTA, a CSV / TSV with Entry and Sequence columns, or a folder of <Entry>.pdb files."""
    ledger = {}
    for path in paths or []:
        p = Path(path)
        if p.is_dir():
            for f in p.glob("*.pdb"):
                ledger[norm(f.stem.split("_")[0])] = md5_file(f)
        elif re.search(r"\.(csv|tsv)(\.gz)?$", p.name):
            import pandas as pd
            df = pd.read_csv(p, sep="\t" if ".tsv" in p.name else ",", usecols=["Entry", "Sequence"], low_memory=False).drop_duplicates("Entry")
            ledger.update({norm(e): md5_text(s) for e, s in zip(df["Entry"], df["Sequence"]) if isinstance(s, str)})
        else:
            opener = gzip.open if p.suffix == ".gz" else open
            sid, buf = None, []
            with opener(p, "rt") as fh:
                for line in fh:
                    line = line.strip()
                    if line.startswith(">"):
                        if sid is not None:
                            ledger[norm(sid)] = md5_text("".join(buf))
                        sid, buf = line[1:].split()[0], []
                    elif line:
                        buf.append(line)
            if sid is not None:
                ledger[norm(sid)] = md5_text("".join(buf))
    return ledger


def stored_fingerprints(h5_path):
    """{id: stored fingerprint or None} for every vector (dataset or group) in an HDF5 file; {} when the file does not exist."""
    if not Path(h5_path).exists():
        return {}
    with h5py.File(h5_path, "r") as f:
        return {k: f[k].attrs.get(ATTR) for k in f.keys()}


def plan(wanted, have, ledger=None, trust_existing=False):
    """
    wanted: {id: fingerprint of the current input}; have: stored_fingerprints(); ledger: read_ledger().
    -> {"todo": {id: reason}, "keep": [...], "stamp": {id: fingerprint} (verified through the ledger: write the fingerprint down), "extra": [ids in the file nobody asked for]}
    """
    ledger = ledger or {}
    todo, keep, stamp = {}, [], {}
    for sid, fp in wanted.items():
        if sid not in have:
            todo[sid] = "new"
            continue
        stored = have[sid]
        if isinstance(stored, bytes):
            stored = stored.decode()
        if stored is not None:
            if stored == fp:
                keep.append(sid)
            else:
                todo[sid] = "changed"
        elif sid in ledger:
            if ledger[sid] == fp:
                keep.append(sid)
                stamp[sid] = fp
            else:
                todo[sid] = "changed"
        elif trust_existing:
            keep.append(sid)
        else:
            todo[sid] = "unverified"
    return {"todo": todo, "keep": keep, "stamp": stamp, "extra": sorted(set(have) - set(wanted))}


def report(name, pl, n_wanted):
    reasons = {r: sum(1 for v in pl["todo"].values() if v == r) for r in ("new", "changed", "unverified")}
    print(f"[plan] {name}: {n_wanted} wanted | keep {len(pl['keep'])} (verified through the ledger and stamped: {len(pl['stamp'])}) | "
          f"embed {len(pl['todo'])} (new {reasons['new']}, changed {reasons['changed']}, unverified {reasons['unverified']}) | in the file but not wanted: {len(pl['extra'])}", flush=True)
    if reasons["unverified"]:
        print(f"[plan] WARNING: {reasons['unverified']} vectors have no fingerprint and no ledger entry, so they are embedded again; give --ledger (what the file was made from) "
              f"or --trust-existing to keep them", flush=True)
    changed = [s for s, r in pl["todo"].items() if r == "changed"]
    if changed:
        print(f"[plan] changed input, vector replaced: {changed[:12]}{' ...' if len(changed) > 12 else ''}", flush=True)


def check_config(h5, config, allow=False):
    """The file records how its vectors were made; refuse to add vectors made another way. A file without a record (the archived ones) gets the current one."""
    for k, v in config.items():
        old = h5.attrs.get(f"config_{k}")
        if old is None:
            h5.attrs[f"config_{k}"] = str(v)
        elif str(old) != str(v) and not allow:
            sys.exit(f"{h5.filename} was made with {k}={old}, this run uses {k}={v}: appending would mix two pipelines in one file "
                     f"(use another --out, or --overwrite to start the file again)")


def apply_stamps(h5_path, stamps):
    if stamps:
        with h5py.File(h5_path, "a") as f:
            for sid, fp in stamps.items():
                f[sid].attrs[ATTR] = fp


def replace_slot(h5, sid):
    """A changed input replaces its old vector: remove it so that the new one can be written under the same id."""
    if sid in h5:
        del h5[sid]

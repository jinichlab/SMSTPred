#!/usr/bin/env python3
"""
Stage 2 -- the classifier bake-off (Fig 2A / Table 1 of a refreshed release): ONE cluster-disjoint split per run, 16 runs = the 16 splits.

  python bakeoff.py --iteration 3            # <release>/benchmark/iterations/iter03_*
  python bakeoff.py --iteration 1 --smoke    # tiny end-to-end test (20 epochs, 600 / 150 / 150 proteins) -> <release>/benchmark/_smoke/
  BAKEOFF_METHODS=MLP_ProtT5,MMseqs2_1NN python bakeoff.py --iteration 1             # only some methods
  sbatch run_bakeoff.sh (16-task array) ; then sbatch run_bakeoff_report.sh

Inputs (all in the release folder, nothing else is needed): labels_SAM.csv (the in_benchmark proteins and their four labels), sam_mtases.tsv.gz (sequences),
embeddings/embeddings_SAM_{ProtT5,ProstT5,SaProt,ESM-IF}.h5, structures/<Entry>.pdb (Foldseek), clusters/splits_tvt_*.json (1_clustering/tvt_splits.py:
16 random splits of the 30 %-identity clusters into 65 % train / 15 % validation / 20 % test). Settings: ../refresh_config.json.

THE PROTOCOL (every method gets the same three roles):
  train       the only data a method is built from: the MLP weights, the MMseqs2 / Foldseek search database, the profile HMMs
  validation  early stopping + the best weights (MLPs) and the per-label threshold of the profile HMM (it has no natural cut-off)
  test        scored once, never used for any choice
MLPs (the six embedding conditions): in -> 512 -> 256 -> 4 outputs, plain BCE on all four labels, Adam (lr 0.001), one step per epoch on a random 90 % of the train proteins,
trained for ALL 2,000 epochs so that the train / validation loss can be plotted to the end. Early stopping only chooses the weights: those of the epoch with the lowest total
validation loss (sum over the four labels), taking only epochs before the loss has failed to improve for 80 epochs in a row (patience); the test proteins are predicted with
those weights. THE DECISION THRESHOLD IS A FIXED 0.5 for the MLPs and for the MMseqs2 / Foldseek vote (a majority vote of the top-10 hits' labels, weighted by bit score).
The profile HMM scores are bit scores, which have no 0.5: its per-label threshold is the MCC-optimal one on the validation scores.
Every method also gets mcc_valtuned: the MCC at per-label thresholds tuned on the validation predictions (a sensitivity check; on the 2026_03 release it was lower than 0.5
for every method, because a validation set holds only 9-33 DNA and 60+ Protein positives).

Writes per iteration: iterNN_detailed.csv, iterNN_per_label.csv (one row per method / label; columns auprc, auroc, mcc, f1, accuracy at the primary thresholds + mcc_valtuned),
iterNN_curves.csv.gz (train / validation loss of every epoch, embedding and label), iterNN_training.csv (the epoch kept, the epoch early stopping fired),
iterNN_scores.npz (validation and test scores of every method), and once config_used.json. Resumable: a finished iteration is skipped (--force redoes it).
Needs on PATH: mmseqs, mafft, hmmbuild, hmmsearch (the conda environment); Foldseek via $FOLDSEEK_BIN, <repo>/tools/foldseek/bin (make tools) or PATH.
"""
import argparse
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from sklearn.metrics import accuracy_score, auc, f1_score, matthews_corrcoef, precision_recall_curve, roc_auc_score

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
CFG = json.loads((REPO / "refresh_config.json").read_text())
LABELS = ["Metabolite", "Protein", "RNA", "DNA"]
SEED = CFG["seed"]
EPOCHS = int(os.environ.get("BAKEOFF_EPOCHS", CFG["mlp"]["epochs"]))
PATIENCE = int(os.environ.get("BAKEOFF_PATIENCE", CFG["mlp"]["patience"]))
LEARNING_RATE = CFG["mlp"]["learning_rate"]
BATCH_FRACTION = CFG["mlp"]["batch_fraction"]
THRESHOLD = CFG["threshold"]
HMM_MAX_SEQS = CFG["hmm"]["max_sequences_per_label"]
MAFFT_THREADS = CFG["hmm"]["mafft_threads"]          # 1: multi-threaded MAFFT may align slightly differently from run to run
THREADS = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# the six MLP conditions (name -> embedding files, concatenated) and the H5 layouts of the embedding files
SOURCES = {"prott5": ("embeddings_SAM_ProtT5.h5", "flat"), "prostt5": ("embeddings_SAM_ProstT5.h5", "mean_group"),
           "esmif": ("embeddings_SAM_ESM-IF.h5", "flat"), "saprot": ("embeddings_SAM_SaProt.h5", "flat")}
MLP_SOURCES = {"MLP_ProstT5": ["prostt5"], "MLP_SaProt": ["saprot"], "MLP_ProtT5": ["prott5"], "MLP_ESM-IF": ["esmif"],
               "MLP_ProtT5+ESM-IF": ["prott5", "esmif"], "MLP_ProstT5+ESM-IF": ["prostt5", "esmif"]}
SEARCH_METHODS = ["MMseqs2_1NN", "Profile_HMM", "Foldseek_1NN"]
ALL_METHODS = list(MLP_SOURCES) + SEARCH_METHODS


# ============================================================================
# release folder -> data
# ============================================================================

def release_dir():
    """$RELEASE_DIR, else the newest release_* of $UNIPROT_PULL_DIR (default <repo>/data/uniprot_pull)."""
    if os.environ.get("RELEASE_DIR"):
        return Path(os.environ["RELEASE_DIR"])
    pull = Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull"))
    rel = sorted(p for p in pull.glob("release_*") if p.is_dir() and not p.name.startswith("_"))
    if not rel:
        sys.exit(f"no release_* folder in {pull} (set RELEASE_DIR or UNIPROT_PULL_DIR)")
    return rel[-1]


def load_source(path, layout):
    out = {}
    with h5py.File(path, "r") as f:
        for k in f.keys():
            out[k] = np.asarray(f[k][:] if layout == "flat" else f[k]["mean"][:], dtype=np.float32)
    return out


def load_data(rel):
    """(table with Entry, Sequence and the four labels, {embedding source: matrix aligned to the table}) for the in_benchmark proteins that have a vector in every embedding file."""
    lab = pd.read_csv(rel / "labels_SAM.csv", usecols=["Entry", "in_benchmark", *LABELS], low_memory=False)
    lab = lab[lab.in_benchmark == True].drop_duplicates("Entry")                              # noqa: E712
    seq = pd.read_csv(rel / "sam_mtases.tsv.gz", sep="\t", usecols=["Entry", "Sequence"], low_memory=False).drop_duplicates("Entry").set_index("Entry")["Sequence"]
    df = lab.assign(Sequence=lab.Entry.map(seq)).reset_index(drop=True)
    maps = {s: load_source(rel / "embeddings" / f, lay) for s, (f, lay) in SOURCES.items()}
    n0 = len(df)
    df = df[df.Entry.apply(lambda e: all(e in m for m in maps.values())) & df.Sequence.notna()].reset_index(drop=True)
    if len(df) < n0:
        print(f"WARNING: {n0 - len(df)} of the {n0} in_benchmark proteins have no sequence or no vector in some embedding file and are left out (run embeddings/run_embeddings.sh)", flush=True)
    return df, {s: np.vstack([m[e] for e in df.Entry]).astype(np.float32) for s, m in maps.items()}


def foldseek_bin():
    for c in (os.environ.get("FOLDSEEK_BIN"), str(REPO / "tools" / "foldseek" / "bin" / "foldseek"), shutil.which("foldseek")):
        if c and (os.path.exists(c) or shutil.which(c)):
            return c
    return None


def sh(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(f"command failed ({r.returncode}): {' '.join(map(str, cmd))[:300]}\n{r.stderr[-800:]}")
    return r


def check_tools(methods, rel):
    missing = []
    if any(m in methods for m in ("MMseqs2_1NN", "Profile_HMM")):
        missing += [t for t in ("mmseqs", "mafft", "hmmbuild", "hmmsearch") if shutil.which(t) is None]
    if "Foldseek_1NN" in methods:
        if foldseek_bin() is None:
            missing.append("foldseek (make tools, or set FOLDSEEK_BIN)")
        if not (rel / "structures").is_dir():
            missing.append(f"{rel / 'structures'} (embeddings/fetch_structures.py)")
    if missing:
        sys.exit(f"Missing required tools/inputs: {missing}")


# ============================================================================
# metrics (macro over the labels that have both classes in the set)
# ============================================================================

def _macro(fn, y_true, y_pred):
    vals = [fn(y_true[:, i], y_pred[:, i]) for i in range(y_true.shape[1]) if len(np.unique(y_true[:, i])) > 1]
    return float(np.mean(vals)) if vals else 0.0


def all_metrics(y, probs, thr):
    pred = (probs >= thr).astype(int)
    return {"auprc": _macro(lambda t, s: auc(*precision_recall_curve(t, s)[1::-1]), y, probs), "auroc": _macro(roc_auc_score, y, probs),
            "mcc": _macro(matthews_corrcoef, y, pred), "f1": _macro(f1_score, y, pred),
            "accuracy": float(np.mean([accuracy_score(y[:, i], pred[:, i]) for i in range(y.shape[1])]))}


def per_label_metrics(y, probs, thr):
    thr = np.broadcast_to(np.asarray(thr, dtype=float), (len(LABELS),))
    pred = (probs >= thr).astype(int)
    out = {}
    for i, lab in enumerate(LABELS):
        if len(np.unique(y[:, i])) > 1:
            pr, rc, _ = precision_recall_curve(y[:, i], probs[:, i])
            out[lab] = {"auprc": auc(rc, pr), "auroc": roc_auc_score(y[:, i], probs[:, i]), "mcc": matthews_corrcoef(y[:, i], pred[:, i]),
                        "f1": f1_score(y[:, i], pred[:, i]), "accuracy": accuracy_score(y[:, i], pred[:, i])}
        else:
            out[lab] = {k: np.nan for k in ("auprc", "auroc", "mcc", "f1", "accuracy")}
    return out


def _mcc_over_thresholds(y, p, cand):
    """MCC of (p >= t) against y for every t of cand at once (a zero denominator gives 0, as in scikit-learn)."""
    pred = p[:, None] >= cand[None, :]
    pos = (y > 0.5)[:, None]
    tp, fp = (pred & pos).sum(0).astype(np.float64), (pred & ~pos).sum(0).astype(np.float64)
    fn, tn = (~pred & pos).sum(0).astype(np.float64), (~pred & ~pos).sum(0).astype(np.float64)
    den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return np.where(den > 0, (tp * tn - fp * fn) / np.where(den > 0, den, 1.0), 0.0)


def tune_thresholds(y, scores, valid=None):
    """Per-label MCC-optimal thresholds on (y, scores): 200 percentile candidates, the first maximum; 0.5 for a label with one class. `valid` (optional) masks usable scores per label."""
    thr = np.full(y.shape[1], 0.5)
    for j in range(y.shape[1]):
        ok = np.ones(len(y), bool) if valid is None else valid[:, j]
        if ok.sum() == 0 or len(np.unique(y[ok, j])) < 2:
            continue
        cand = np.unique(np.percentile(scores[ok, j], np.linspace(0, 100, 200)))
        thr[j] = cand[int(np.argmax(_mcc_over_thresholds(y[ok, j], scores[ok, j], cand)))]
    return thr


# ============================================================================
# MLP: trained through all epochs, the weights of the best validation epoch are kept
# ============================================================================

class MultilabelNN(nn.Module):
    def __init__(self, input_size, output_size):
        super().__init__()
        self.fc1, self.relu, self.fc2, self.fc3 = nn.Linear(input_size, 512), nn.ReLU(), nn.Linear(512, 256), nn.Linear(256, output_size)

    def forward(self, x):
        return self.fc3(self.relu(self.fc2(self.relu(self.fc1(x)))))      # raw logits


class EarlyStopper:
    """Which epoch's weights are kept. update(epoch, loss) -> True when this epoch is a new best that is kept. The best is only updated until the loss has failed to improve for
    `patience` epochs in a row (stop_epoch = that epoch); after that the run goes on (the curves continue) but nothing more can be kept."""
    def __init__(self, patience):
        self.patience, self.best, self.best_epoch, self.stop_epoch, self.since = patience, float("inf"), 0, None, 0

    def update(self, epoch, loss):
        if self.stop_epoch is not None:
            return False
        if loss < self.best:
            self.best, self.best_epoch, self.since = loss, epoch, 0
            return True
        self.since += 1
        if self.since >= self.patience:
            self.stop_epoch = epoch
        return False


def label_losses(logits, y):
    """[n_labels] mean BCE of every label."""
    return torch.stack([F.binary_cross_entropy_with_logits(logits[:, j], y[:, j]) for j in range(y.shape[1])])


def train_early_stop(Xtr, ytr, Xva, yva, epochs, patience, lr, seed):
    """-> (model with the best weights, histories {train, val}: [epochs, labels], EarlyStopper)."""
    torch.manual_seed(seed)
    model = MultilabelNN(Xtr.shape[1], ytr.shape[1]).to(DEVICE)
    opt = optim.Adam(model.parameters(), lr=lr)
    tt = lambda a: torch.tensor(a, dtype=torch.float32, device=DEVICE)
    Xtr_t, ytr_t, Xva_t, yva_t = tt(Xtr), tt(ytr), tt(Xva), tt(yva)
    n = Xtr_t.shape[0]
    k = max(1, int(BATCH_FRACTION * n))
    hist = {key: torch.zeros(epochs, ytr.shape[1], device=DEVICE) for key in ("train", "val")}
    stopper, best_state = EarlyStopper(patience), None
    for ep in range(1, epochs + 1):
        idx = torch.randperm(n, device=DEVICE)[:k]
        model.train()
        opt.zero_grad()
        tr = label_losses(model(Xtr_t[idx]), ytr_t[idx])
        tr.sum().backward()
        opt.step()
        model.eval()
        with torch.no_grad():
            va = label_losses(model(Xva_t), yva_t)
        hist["train"][ep - 1], hist["val"][ep - 1] = tr.detach(), va
        if stopper.update(ep, float(va.sum())):
            best_state = {key: v.detach().clone() for key, v in model.state_dict().items()}
    model.load_state_dict(best_state)
    return model, {key: v.cpu().numpy() for key, v in hist.items()}, stopper


def predict(model, X):
    model.eval()
    with torch.no_grad():
        return torch.sigmoid(model(torch.tensor(X, dtype=torch.float32, device=DEVICE))).cpu().numpy()


# ============================================================================
# search baselines: scores of the query proteins (validation + test) against the TRAIN proteins only
# ============================================================================

def _fasta(path, entries, seqs):
    with open(path, "w") as f:
        for e, s in zip(entries, seqs):
            f.write(f">{e}\n{s}\n")


def _vote(res, query_to_idx, n_query, train_label_map):
    """Top-10 hits by bit score; the score of a label = the bit-score-weighted fraction of those hits that carry it."""
    y = np.zeros((n_query, len(LABELS)), dtype=np.float32)
    for query, group in res.groupby("query"):
        if str(query) not in query_to_idx:
            continue
        top = group.nlargest(10, "bits")
        total = top["bits"].values.clip(min=0).sum()
        idx = query_to_idx[str(query)]
        if total > 0:
            for tgt, bits in zip(top["target"].astype(str), top["bits"]):
                if tgt in train_label_map and bits > 0:
                    y[idx] += (bits / total) * train_label_map[tgt]
        elif str(top.iloc[0]["target"]) in train_label_map:
            y[idx] = train_label_map[str(top.iloc[0]["target"])]
    return y


def mmseqs_scores(train_entries, train_seqs, train_y, q_entries, q_seqs):
    with tempfile.TemporaryDirectory() as tmp:
        tr_fa, q_fa = os.path.join(tmp, "train.fasta"), os.path.join(tmp, "query.fasta")
        _fasta(tr_fa, train_entries, train_seqs)
        _fasta(q_fa, q_entries, q_seqs)
        trainDB, qDB, resDB = (os.path.join(tmp, x) for x in ("trainDB", "qDB", "resDB"))
        sh(["mmseqs", "createdb", tr_fa, trainDB])
        sh(["mmseqs", "createdb", q_fa, qDB])
        sh(["mmseqs", "search", qDB, trainDB, resDB, os.path.join(tmp, "tmp"), "--max-seqs", "10", "-e", "100", "--threads", str(THREADS)])
        res_tsv = os.path.join(tmp, "res.tsv")
        sh(["mmseqs", "convertalis", qDB, trainDB, resDB, res_tsv, "--format-output", "query,target,evalue,bits"])
        if os.path.getsize(res_tsv) == 0:
            return np.zeros((len(q_entries), len(LABELS)), dtype=np.float32)
        res = pd.read_csv(res_tsv, sep="\t", header=None, names=["query", "target", "evalue", "bits"])
        return _vote(res, {e: i for i, e in enumerate(q_entries)}, len(q_entries), dict(zip(train_entries, train_y)))


def _foldseek_db(entries, struct_dir, dbpath, tmp):
    """Symlink the needed {Entry}.pdb into a folder and foldseek-createdb it."""
    sdir = os.path.join(tmp, os.path.basename(dbpath) + "_pdbs")
    os.makedirs(sdir, exist_ok=True)
    n = 0
    for e in entries:
        src = os.path.join(struct_dir, f"{e}.pdb")
        if os.path.exists(src):
            os.symlink(os.path.abspath(src), os.path.join(sdir, f"{e}.pdb"))
            n += 1
    sh([foldseek_bin(), "createdb", sdir, dbpath])
    return n


def _name_to_entry(name):
    """Foldseek query/target name (e.g. 'P12345.pdb_A') -> UniProt Entry (files are {Entry}.pdb; accessions have no '_')."""
    return str(name).split()[0].replace(".pdb", "").split("_")[0]


def foldseek_scores(train_entries, train_y, q_entries, struct_dir):
    with tempfile.TemporaryDirectory() as tmp:
        trainDB, qDB, resDB = (os.path.join(tmp, x) for x in ("trainDB", "qDB", "resDB"))
        n_tr, n_q = _foldseek_db(train_entries, struct_dir, trainDB, tmp), _foldseek_db(q_entries, struct_dir, qDB, tmp)
        print(f"      structures found: train {n_tr}/{len(train_entries)}, query {n_q}/{len(q_entries)}", flush=True)
        res_tsv = os.path.join(tmp, "res.tsv")
        sh([foldseek_bin(), "search", qDB, trainDB, resDB, os.path.join(tmp, "t1"), "-e", "100", "--max-seqs", "10", "--threads", str(THREADS)])
        sh([foldseek_bin(), "convertalis", qDB, trainDB, resDB, res_tsv, "--format-output", "query,target,evalue,bits"])
        if os.path.getsize(res_tsv) == 0:
            return np.zeros((len(q_entries), len(LABELS)), dtype=np.float32)
        res = pd.read_csv(res_tsv, sep="\t", header=None, names=["query", "target", "evalue", "bits"])
        res["query"], res["target"] = res["query"].map(_name_to_entry), res["target"].map(_name_to_entry)
        return _vote(res, {e: i for i, e in enumerate(q_entries)}, len(q_entries), dict(zip(train_entries, train_y)))


def hmm_scores(train_entries, train_seqs, train_y, q_entries, q_seqs):
    """Best hmmsearch bit score of every query against the profile HMM of each label, built from <= 500 train positives (-1000 = no hit)."""
    scores = np.full((len(q_entries), len(LABELS)), -1000.0)
    with tempfile.TemporaryDirectory() as tmp:
        q_fa = os.path.join(tmp, "query.fasta")
        _fasta(q_fa, q_entries, q_seqs)
        q_idx = {e: i for i, e in enumerate(q_entries)}
        for ci, label in enumerate(LABELS):
            pos = [i for i in range(len(train_entries)) if train_y[i, ci] == 1]
            if len(pos) < 2:
                continue
            if len(pos) > HMM_MAX_SEQS:
                pos = [pos[i] for i in random.sample(range(len(pos)), HMM_MAX_SEQS)]
            class_fa, aligned, hmm, tbl = (os.path.join(tmp, f"{label}{s}") for s in (".fasta", "_aln.fasta", ".hmm", ".tbl"))
            _fasta(class_fa, [train_entries[i] for i in pos], [train_seqs[i] for i in pos])
            if subprocess.run(f"mafft --auto --thread {MAFFT_THREADS} --quiet {class_fa} > {aligned}", shell=True, capture_output=True).returncode != 0:
                print(f"      [warn] mafft failed for {label}; class skipped", flush=True)
                continue
            if subprocess.run(["hmmbuild", "--amino", "--cpu", str(THREADS), hmm, aligned], capture_output=True).returncode != 0:
                print(f"      [warn] hmmbuild failed for {label}; class skipped", flush=True)
                continue
            sh(["hmmsearch", "--tblout", tbl, "-E", "1000", "--cpu", str(THREADS), hmm, q_fa])
            with open(tbl) as f:
                for line in f:
                    if not line.startswith("#"):
                        parts = line.split()
                        if parts[0] in q_idx:
                            scores[q_idx[parts[0]], ci] = max(scores[q_idx[parts[0]], ci], float(parts[5]))
    return scores


# ============================================================================
# scoring one method
# ============================================================================

def score_method(name, s_va, s_te, y_va, y_te, overall, per_label, scores, hmm=False):
    """Primary: the fixed THRESHOLD (0.5) for probabilities / vote fractions; for the HMM (bit scores) the MCC-optimal threshold on the validation scores.
    mcc_valtuned: the validation-tuned thresholds for every method (for the HMM it is the primary threshold itself)."""
    if hmm:
        valid = s_va > -999
        thr_val = tune_thresholds(y_va, s_va, valid)
        thr_val[~valid.any(0)] = 0.0                                          # no scored validation protein for a label: nothing passes
        sd = np.array([np.std(s_va[valid[:, j], j]) if valid[:, j].sum() > 1 else 0.0 for j in range(len(LABELS))])
        scale = np.where(sd > 0, 2 * sd, 50.0)
        p_te = (1 / (1 + np.exp(-(s_te - thr_val) / scale))).astype(np.float32)   # p >= 0.5 <=> score >= the validation threshold
        prim_thr, probs, tuned_probs, tuned_thr = 0.5, p_te, p_te, 0.5
    else:
        thr_val = tune_thresholds(y_va, s_va)
        prim_thr, probs, tuned_probs, tuned_thr = THRESHOLD, s_te, s_te, thr_val
    overall[name] = all_metrics(y_te, probs, prim_thr)
    per_label[name] = per_label_metrics(y_te, probs, prim_thr)
    overall[name]["mcc_valtuned"] = all_metrics(y_te, tuned_probs, tuned_thr)["mcc"]
    tuned_pl = per_label_metrics(y_te, tuned_probs, tuned_thr)
    for lab in LABELS:
        per_label[name][lab]["mcc_valtuned"] = tuned_pl[lab]["mcc"]
    scores.update({f"{name}__val": s_va.astype(np.float32), f"{name}__test": s_te.astype(np.float32), f"{name}__thr_valtuned": np.asarray(thr_val, dtype=np.float64)})


# ============================================================================
# one iteration
# ============================================================================

def run_iteration(it, smoke, force):
    rel = release_dir()
    out_root = Path(os.environ.get("BAKEOFF_OUT", rel / "benchmark")) / ("_smoke" if smoke else "")
    out_dir = out_root / "iterations"
    out_dir.mkdir(parents=True, exist_ok=True)
    f_det, f_pl = out_dir / f"iter{it:02d}_detailed.csv", out_dir / f"iter{it:02d}_per_label.csv"
    if f_det.exists() and not force:
        print(f"iteration {it}: {f_det.name} exists, skipping (use --force to redo)")
        return
    methods = [m.strip() for m in os.environ.get("BAKEOFF_METHODS", ",".join(ALL_METHODS)).split(",")]
    assert all(m in ALL_METHODS for m in methods), f"unknown method in {methods}; choose from {ALL_METHODS}"
    check_tools(methods, rel)
    sp_cfg = CFG["splits"]
    splits_json = Path(os.environ.get("TVT_SPLITS_JSON", rel / "clusters" / f"splits_tvt_sam_benchmark_id30_t{round(sp_cfg['test_frac'] * 100)}_v{round(sp_cfg['val_frac'] * 100)}.json"))
    if not splits_json.exists():
        sys.exit(f"{splits_json} not found: run 1_clustering/tvt_splits.py (make splits) first")
    (out_root / "config_used.json").write_text(json.dumps({"refresh_config": CFG, "release": str(rel), "splits_file": str(splits_json), "epochs": EPOCHS, "patience": PATIENCE}, indent=1))

    df, mats = load_data(rel)
    entries, seqs = df["Entry"].values, df["Sequence"].values
    y = df[LABELS].values.astype(np.float32)
    sp = json.loads(splits_json.read_text())["iterations"][it - 1]
    pos = {e: i for i, e in enumerate(entries)}
    idx = {part: np.array([pos[e] for e in sp[part] if e in pos]) for part in ("train", "validation", "test")}
    missing = {part: sum(e not in pos for e in sp[part]) for part in idx}
    tr, va, te = idx["train"], idx["validation"], idx["test"]
    assert not (set(tr) & set(va)) and not (set(tr) & set(te)) and not (set(va) & set(te)), "partitions overlap"
    epochs, patience = EPOCHS, PATIENCE
    if smoke:
        rng = np.random.RandomState(0)
        tr, va, te = (rng.choice(a, min(n, len(a)), replace=False) for a, n in ((tr, 600), (va, 150), (te, 150)))
        epochs, patience = 20, 5
    y_tr, y_va, y_te = y[tr], y[va], y[te]
    print(f"iteration {it}: train {len(tr)}, validation {len(va)}, test {len(te)} (proteins of the split without data: {missing}) | positives M/P/R/D "
          f"train {y_tr.sum(0).astype(int).tolist()} validation {y_va.sum(0).astype(int).tolist()} test {y_te.sum(0).astype(int).tolist()} | device {DEVICE} | threads {THREADS} | "
          f"epochs {epochs}, patience {patience}, threshold {THRESHOLD} | smoke={smoke}", flush=True)

    overall, per_label, scores = {}, {}, {"y_val": y_va, "y_test": y_te, "val_entries": entries[va], "test_entries": entries[te]}
    curves, training = [], []
    for k, (name, srcs) in enumerate(MLP_SOURCES.items()):
        if name not in methods:
            continue
        t0 = time.time()
        Xtr, Xva, Xte = (np.hstack([mats[s][a] for s in srcs]) for a in (tr, va, te))
        model, hist, stop = train_early_stop(Xtr, y_tr, Xva, y_va, epochs, patience, LEARNING_RATE, SEED + 100 * (it - 1) + k)
        score_method(name, predict(model, Xva), predict(model, Xte), y_va, y_te, overall, per_label, scores)
        emb, ep = name[len("MLP_"):], np.arange(1, epochs + 1)
        curves.append(pd.DataFrame({"iteration": it, "embedding": emb, "epoch": np.repeat(ep, len(LABELS)), "label": np.tile(LABELS, epochs),
                                    "train_loss": hist["train"].reshape(-1), "val_loss": hist["val"].reshape(-1)}))
        training.append({"iteration": it, "embedding": emb, "best_epoch": stop.best_epoch, "stop_epoch": stop.stop_epoch if stop.stop_epoch is not None else np.nan,
                         "best_val_total": stop.best, "val_total_at_end": float(hist["val"][-1].sum()), "epochs_trained": epochs, "n_train": len(tr), "n_val": len(va)})
        print(f"  {name:20s} MCC {overall[name]['mcc']:.3f} (validation-tuned thresholds: {overall[name]['mcc_valtuned']:.3f}) | kept epoch {stop.best_epoch}, "
              f"early stopping fired at {stop.stop_epoch} | {time.time() - t0:.0f}s", flush=True)

    random.seed(SEED + 1000 + it)                                    # profile-HMM class subsampling
    np.random.seed(SEED + it)
    q = np.concatenate([va, te])
    nv = len(va)
    q_entries, q_seqs = entries[q], seqs[q]
    for name, fn, is_hmm in (("MMseqs2_1NN", lambda: mmseqs_scores(entries[tr], seqs[tr], y_tr, q_entries, q_seqs), False),
                             ("Profile_HMM", lambda: hmm_scores(entries[tr], seqs[tr], y_tr, q_entries, q_seqs), True),
                             ("Foldseek_1NN", lambda: foldseek_scores(entries[tr], y_tr, q_entries, str(rel / "structures")), False)):
        if name not in methods:
            continue
        t0 = time.time()
        s = fn()
        score_method(name, s[:nv], s[nv:], y_va, y_te, overall, per_label, scores, hmm=is_hmm)
        print(f"  {name:20s} MCC {overall[name]['mcc']:.3f}  ({time.time() - t0:.0f}s)", flush=True)

    meta = {"split_mode": "cluster", "dataset": "full", "iteration": it}
    pd.DataFrame([{**meta, "classifier": c, "n_train": len(tr), "n_val": len(va), "n_test": len(te), **m} for c, m in overall.items()]).to_csv(f_det, index=False)
    pd.DataFrame([{**meta, "classifier": c, "label": lab, **lm} for c, d in per_label.items() for lab, lm in d.items()]).to_csv(f_pl, index=False)
    if curves:
        pd.concat(curves).to_csv(out_dir / f"iter{it:02d}_curves.csv.gz", index=False)
        pd.DataFrame(training).to_csv(out_dir / f"iter{it:02d}_training.csv", index=False)
    np.savez_compressed(out_dir / f"iter{it:02d}_scores.npz", **scores)
    print(f"wrote {f_det}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--iteration", type=int, required=True, help="1..16 (SLURM_ARRAY_TASK_ID)")
    ap.add_argument("--smoke", action="store_true", help="20 epochs, 600 / 150 / 150 proteins, writes to _smoke/")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    torch.set_num_threads(THREADS)
    run_iteration(a.iteration, a.smoke, a.force)


if __name__ == "__main__":
    main()

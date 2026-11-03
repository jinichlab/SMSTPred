#!/usr/bin/env python3
"""Figure 8 -- Horizyn reaction profiles of the four dark clusters where CLEAN and Horizyn name the same chemistry.

Each panel is one dark cluster (c551, c667, c110, c714). Grey histogram: the cluster-mean Horizyn score of every one of
the 36,552 Rhea reaction entries (mean over the cluster's members of the member-vs-reaction score; log y). Orange rug: the
SAM-methylation reactions (both directional entries of the 677 reactions flagged as SAM methylation = 1,714 entries).
Solid vertical line: the highest-scoring reaction. The substrate -> product reaction drawn in the panel is the reaction the
largest share of members voted for as their single best match (c714: the furocoumarin reaction CLEAN names, which is its
second most-voted reaction and third by score; the highest-scoring reaction, solid line, is a different reaction, so both are marked).
Rhea numbers a reaction as master n with directional entries n+1 (left->right) and n+2 (right->left); Horizyn scores the two
directions almost identically, so reactions are collapsed over direction and labelled by their master ID (as in the text).

Inputs: horizyn_dark_scores.npz (3,478 dark proteins x 36,552 reactions; Zenodo record, in data/ or $SMSTPRED_DATA; it can be
regenerated with 6_dark_cluster_annotation/horizyn/run_horizyn_predict_full.py) and
6_dark_cluster_annotation/horizyn/rhea_reaction_clusters.csv (in this repository).
Output: Fig8_reaction_profiles.png

  python make_fig8_reaction_profiles.py                    the archived figure (c551, c667, c110, c714)
  python make_fig8_reaction_profiles.py --release <dir>    the refreshed figure: the clusters chosen by select_panels.py (CLEAN reliable + MTase and its EC inside Horizyn's dominant
                                                           family; c317, c353, c512, c729 on release 2026_03; titles from panel_titles.csv), scores from
                                                           <dir>/dark_annotation/horizyn_run/horizyn_full_scores.npz, output <dir>/fig8/Fig8_reaction_profiles.png
"""
import argparse
from pathlib import Path
import io
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DATA = Path(os.environ.get("SMSTPRED_DATA", REPO / "data"))
ST6 = REPO / "6_dark_cluster_annotation"

# (cluster, panel letter, title, Rhea master ID of the drawn reaction)
PANELS = [(551, "A", "Tricetin → 3′,5′-di-O-methyltricetin", 32347),
          (667, "B", "7-methylxanthine → theobromine", 24604),
          (110, "C", "Guanidinoacetate → creatine", 10656),
          (714, "D", "5-hydroxy → 5-methoxyfurocoumarin", 18861)]
COL = {551: "#4578aa", 667: "#618d43", 110: "#dbb92c", 714: "#c14445"}      # release mode adds its own (select_panels)
RUG = "#ed8c44"
SAM = Chem.CanonSmiles("C[S+](CC[C@H]([NH3+])C(=O)[O-])C[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1O")
SAH = Chem.CanonSmiles("NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](CSCC[C@H]([NH3+])C(=O)[O-])[C@@H](O)[C@H]1O")
FONT = ImageFont.truetype(str(Path(matplotlib.get_data_path()) / "fonts/ttf/DejaVuSans.ttf"), 30)


def canon(f):
    try:
        return Chem.CanonSmiles(f)
    except Exception:
        return f


# ------------------------------------------------------------------ data
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--release", default=None, help="release folder (e.g. uniprot_refresh/release_2026_03); default = the archived figure")
args = ap.parse_args()
if args.release:
    REL = Path(args.release)
    PANELS = None      # chosen below by select_panels.select (the rule), titles from panel_titles.csv
    SCORES = REL / "dark_annotation/horizyn_run/horizyn_full_scores.npz"
    OUTDIR = REL / "fig8"; OUTDIR.mkdir(exist_ok=True)
else:
    SCORES, OUTDIR = DATA / "horizyn_dark_scores.npz", HERE
if args.release:
    import sys; sys.path.insert(0, str(HERE))
    from select_panels import select, master_map
    MASTER = master_map(pd.read_csv(REL / "reference_tables/rhea-directions.tsv", sep="\t"))
    ann = REL / "dark_annotation"
    chosen = select(pd.read_csv(ann / "clean/clean_summary.csv"), pd.read_csv(ann / "horizyn/protein_clusters_family_final.csv"),
                    pd.read_csv(ST6 / "horizyn/rhea_reaction_clusters.csv"), pd.read_csv(REL / "reference_tables/rhea2ec.tsv", sep="\t", dtype=str), MASTER)
    if not chosen:
        raise SystemExit("no dark cluster meets the Fig 8 rule (CLEAN reliable + MTase, EC inside Horizyn's dominant family): nothing to draw")
    PANELS = [(c, "ABCDEF"[i], None, None) for i, c in enumerate(chosen)]
    CLEAN_NAME = pd.read_csv(ann / "clean/clean_summary.csv").set_index("cluster")["enzyme"].to_dict()
    TITLES = pd.read_csv(HERE / "panel_titles.csv").set_index("rhea_master")["title"].to_dict()      # optional hand-written overrides (empty by default)
    sys.path.insert(0, str(ST6))
    from rhea_text import load_equations, reaction_title
    EQ = load_equations(REL / "reference_tables" / "rhea_equations.tsv")                          # the titles are written from Rhea's equation text
    PALETTE = ["#4578aa", "#618d43", "#dbb92c", "#c14445", "#7b5ea7", "#3b9aa3"]
    for i, c in enumerate(chosen):
        COL[c] = PALETTE[i]
    print("Fig 8 panels:", chosen, flush=True)
    pd.DataFrame({"cluster": chosen}).to_csv(OUTDIR / "fig8_clusters.csv", index=False)      # Fig 7 outlines these (build_cluster_grids.py --release)
z = np.load(SCORES, allow_pickle=True)
S, clus, rid, is_sam = z["scores"], z["clusters"], z["rhea_ids"], z["is_sam_methylation"]
N = S.shape[1]
rc = pd.read_csv(ST6 / "horizyn/rhea_reaction_clusters.csv").set_index("rhea_lr_id").loc[rid]
side = lambda x: ".".join(sorted(x.split(".")))
keys = [tuple(sorted(side(q) for q in smi.split(">>"))) for smi in rc.reaction_smiles]
kid = {k: i for i, k in enumerate(dict.fromkeys(keys))}
rxn = np.array([kid[k] for k in keys])                       # one id per direction-collapsed reaction
sam_rxn = np.isin(rxn, np.unique(rxn[is_sam]))
top1 = S.argmax(1)


def forward_pair(group):
    """Substrate / product SMILES and entry ID of the direction in which SAM is consumed."""
    for e in np.where(rxn == group)[0]:
        l, r = rc.reaction_smiles.iloc[e].split(">>")
        L = [f for f in l.split(".") if f != "[H+]"]; R = [f for f in r.split(".") if f != "[H+]"]
        if any(canon(f) == SAM for f in L):
            sub = [f for f in L if canon(f) != SAM]; prod = [f for f in R if canon(f) != SAH]
            if sub and prod:
                return sub[0], prod[0], int(rid[e])
    raise ValueError(group)


def mol_img(smiles, w, h):
    d = rdMolDraw2D.MolDraw2DCairo(w, h); o = d.drawOptions()
    o.bondLineWidth = 3; o.minFontSize = 26; o.maxFontSize = 34
    d.DrawMolecule(Chem.MolFromSmiles(smiles)); d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText())).convert("RGB")


def reaction_img(sub, prod, w=520, h=330, gap=200):
    cv = Image.new("RGB", (2 * w + gap, h), "white")
    cv.paste(mol_img(sub, w, h), (0, 0)); cv.paste(mol_img(prod, w, h), (w + gap, 0))
    d = ImageDraw.Draw(cv); y = h // 2
    d.line([(w + 28, y), (w + gap - 32, y)], fill="black", width=5)
    d.polygon([(w + gap - 32, y - 13), (w + gap - 6, y), (w + gap - 32, y + 13)], fill="black")
    d.text((w + 44, y + 16), "SAM → SAH", fill="#444444", font=FONT)
    return cv


# ------------------------------------------------------------------ figure
plt.rcParams.update({"font.size": 12})
nrow = (len(PANELS) + 1) // 2
fig, axes = plt.subplots(nrow, 2, figsize=(14, 4.9 * nrow), squeeze=False)
for ax in axes.ravel()[len(PANELS):]:
    ax.axis("off")
SHOWN = []
for ax, (c, letter, title, master) in zip(axes.ravel(), PANELS):
    idx = np.where(clus == c)[0]; mu = S[idx].mean(0); top = int(mu.argmax())
    votes = pd.Series(rxn[top1[idx]]).value_counts()
    ax.hist(mu, bins=140, color="#c8c8c8", edgecolor="none"); ax.set_yscale("log")
    ymin = ax.get_ylim()[0]; ymax = 10 ** 5.3                 # headroom above the histogram for the drawing
    ax.set_ylim(ymin, ymax); ax.set_yticks([1, 10, 100, 1000]); ax.set_yticklabels(["10⁰", "10¹", "10²", "10³"])
    ax.vlines(mu[sam_rxn], ymin, ymin * 2.2, color=RUG, lw=0.5, alpha=0.7)
    ax.axvline(mu[top], color=COL[c], lw=2.6)
    # drawn reaction: the most-voted one (c714: the furocoumarin reaction, a lower vote rank)
    g_show = votes.index[0] if c != 714 or args.release else rxn[np.where(rid == master + 1)[0][0]]
    sub, prod, entry = forward_pair(g_show)
    if master is None:                                        # release mode: the drawn reaction is whatever the members voted for
        master = MASTER[entry]
        title = TITLES.get(master) or reaction_title(EQ.get(master, "")) or f"{CLEAN_NAME.get(c, '')} (RHEA:{master})"
    assert (MASTER[entry] if args.release else entry - 1) == master, (c, entry, master)   # the drawn entry belongs to the master reaction quoted in the text (archived figure: its left->right entry is master + 1)
    SHOWN.append((c, master))
    ents = np.where(rxn == g_show)[0]; e_ = ents[np.argmax(mu[ents])]
    same = rxn[top] == g_show
    if not same:
        ax.axvline(mu[e_], color=COL[c], lw=2.0, ls="--")
    ins = ax.inset_axes([0.20, 0.70, 0.66, 0.29])
    ins.imshow(np.asarray(reaction_img(sub, prod))); ins.axis("off")
    ymid = 10 ** (np.log10(ymin) + 0.85 * (np.log10(ymax) - np.log10(ymin)))
    ax.annotate("", xy=(mu[e_], ymid), xycoords="data", xytext=(0.87, 0.85), textcoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=COL[c], lw=1.8, shrinkA=0, shrinkB=0))
    share = votes[g_show] / len(idx) * 100; z_e = (mu[e_] - mu.mean()) / mu.std()
    lines = [f"RHEA:{master} · voted by {share:.0f}% of members", f"score {mu[e_]:.2f}  (z = {z_e:.1f})"]
    if not same:
        lines += ["solid line: top-scoring reaction"]
    ax.text(0.86, 0.655, "\n".join(lines), transform=ax.transAxes, ha="right", va="top", fontsize=12.5,
            color=COL[c], fontweight="bold", linespacing=1.25)
    ax.set_title(f"c{c} · {title}  (n = {len(idx)})", loc="left", fontweight="bold", fontsize=14, pad=10)
    ax.text(-0.115, 1.075, letter, transform=ax.transAxes, fontsize=26, fontweight="bold", va="bottom")
    ax.set_xlabel("cluster-mean Horizyn score", fontsize=13); ax.set_ylabel("no. of reactions (log)", fontsize=13)
    ax.tick_params(labelsize=12)
fig.legend(handles=[Patch(fc="#c8c8c8", label=f"all {N:,} Rhea reaction entries"),
                    plt.Line2D([], [], color=RUG, lw=2.2, label="SAM-methylation reactions (rug)")],
           loc="lower center", ncol=2, frameon=False, fontsize=13.5)
fig.tight_layout(rect=[0.0, 0.045, 1, 1])
out = OUTDIR / "Fig8_reaction_profiles.png"
fig.savefig(out, dpi=300, facecolor="white"); plt.close(fig)
print("wrote", out)
for c, m in SHOWN:
    v = pd.Series(rxn[top1[np.where(clus == c)[0]]]).value_counts()
    print(f"  c{c}: master RHEA:{m}; distinct top hits {len(v)}; top-voted share {v.iloc[0] / v.sum():.0%}")
print(f"  SAM-methylation: {int(is_sam.sum())} flagged entries -> {len(np.unique(rxn[is_sam]))} reactions, {int(sam_rxn.sum())} directional entries")

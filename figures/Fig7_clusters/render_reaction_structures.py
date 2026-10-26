#!/usr/bin/env python3
"""Render the actual 2D reaction schemes (substrate + SAM -> product + SAH)
for the four clusters where CLEAN and Horizyn independently name the same (or
related) specific chemistry, for Figure 7 panel E.

Rhea reactions used (the single most-voted individual reaction among each
cluster's members, from Horizyn's per-protein rank-1 votes):
  c551  RHEA:32347  tricetin -> 3',5'-di-O-methyltricetin            (84% of members)
  c667  RHEA:24604  7-methylxanthine -> theobromine (caffeine pathway) (65% of members)
  c110  RHEA:10656  guanidinoacetate -> creatine                      (84% of members)
  c714  RHEA:18861  5-hydroxyfurocoumarin -> 5-methoxyfurocoumarin    (20% of members; cluster is mixed)

Output: E_reaction_matches.png
"""
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont
import matplotlib

IO = Path(__file__).parent.resolve()

SAM = "C[S+](CC[C@H]([NH3+])C(=O)[O-])C[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1O"
SAH = "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](CSCC[C@H]([NH3+])C(=O)[O-])[C@@H](O)[C@H]1O"
SAM_C, SAH_C = Chem.CanonSmiles(SAM), Chem.CanonSmiles(SAH)

REACTIONS = [
    dict(num="1", cluster="c551", rhea="RHEA:32347", share="84%",
         rxn_smiles="O=C1C=C(C2=CC(O)=C(O)C(O)=C2)OC2=CC([O-])=CC(O)=C12."
                     "C[S+](CC[C@H]([NH3+])C(=O)[O-])C[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1O."
                     "C[S+](CC[C@H]([NH3+])C(=O)[O-])C[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1O"
                     ">>COC1=CC(C2=CC(=O)C3=C(O)C=C([O-])C=C3O2)=CC(OC)=C1O."
                     "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](CSCC[C@H]([NH3+])C(=O)[O-])[C@@H](O)[C@H]1O."
                     "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](CSCC[C@H]([NH3+])C(=O)[O-])[C@@H](O)[C@H]1O.[H+].[H+]",
         label="Tricetin → 3',5'-di-O-methyltricetin", clean="tricin synthase"),
    dict(num="2", cluster="c667", rhea="RHEA:24604", share="65%",
         rxn_smiles="CN1C=NC2=C1C(=O)NC(=O)N2."
                     "C[S+](CC[C@H]([NH3+])C(=O)[O-])C[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1O"
                     ">>CN1C=NC2=C1C(=O)NC(=O)N2C."
                     "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](CSCC[C@H]([NH3+])C(=O)[O-])[C@@H](O)[C@H]1O.[H+]",
         label="7-methylxanthine → theobromine\n(caffeine-pathway step)", clean="caffeine synthase"),
    dict(num="3", cluster="c110", rhea="RHEA:10656", share="84%",
         rxn_smiles="NC(=[NH2+])NCC(=O)[O-]."
                     "C[S+](CC[C@H]([NH3+])C(=O)[O-])C[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1O"
                     ">>CN(CC(=O)[O-])C(N)=[NH2+]."
                     "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](CSCC[C@H]([NH3+])C(=O)[O-])[C@@H](O)[C@H]1O.[H+]",
         label="Guanidinoacetate → creatine", clean="guanidinoacetate N-methyltransferase"),
    dict(num="4", cluster="c714", rhea="RHEA:18861", share="20% (mixed cluster)",
         rxn_smiles="[1*]/C1=C([2*])/C2=C(\\O)C3=C(OC(=O)C([4*])=C3[3*])/C([5*])=C\\2O1."
                     "C[S+](CC[C@H]([NH3+])C(=O)[O-])C[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1O"
                     ">>*C1=C(*)C2=C(OC1=O)/C(*)=C1/O/C(*)=C(*)\\C1=C\\2OC."
                     "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](CSCC[C@H]([NH3+])C(=O)[O-])[C@@H](O)[C@H]1O.[H+]",
         label="5-hydroxyfurocoumarin → 5-methoxyfurocoumarin",
         clean="5-hydroxyfuranocoumarin 5-O-methyltransferase"),
]

MOL_PX = 560


def split_reaction(rxn_smiles):
    """Return (substrate, SAM, product, SAH) SMILES, dropping bare protons."""
    lhs, rhs = rxn_smiles.split(">>")
    reactants = [f for f in lhs.split(".") if f.strip() != "[H+]"]
    products = [f for f in rhs.split(".") if f.strip() != "[H+]"]

    def is_cofactor(frag, canon):
        try:
            return Chem.CanonSmiles(frag) == canon
        except Exception:
            return False

    sub = next(f for f in reactants if not is_cofactor(f, SAM_C))
    prod = next(f for f in products if not is_cofactor(f, SAH_C))
    return sub, prod


_mol_cache = {}


def draw_mol(smiles, px):
    key = (smiles, px)
    if key in _mol_cache:
        return _mol_cache[key]
    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(px, px)
    d.drawOptions().padding = 0.12
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    png = d.GetDrawingText()
    tmp = IO / "_tmp_mol.png"
    with open(tmp, "wb") as f:
        f.write(png)
    img = Image.open(tmp).convert("RGB")
    tmp.unlink()
    _mol_cache[key] = img
    return img


NCOLS = 2


def draw_reaction_cell(canvas, d, x0, y0, r, fonts):
    f_bold, f_reg, f_small = fonts
    sub_smi, prod_smi = split_reaction(r["rxn_smiles"])
    sub_img = draw_mol(sub_smi, MOL_PX)
    prod_img = draw_mol(prod_smi, MOL_PX)

    bd = 48
    d.ellipse((x0, y0, x0 + bd, y0 + bd), fill="#1a1a1a", outline="white", width=2)
    bb = d.textbbox((0, 0), r["num"], font=f_bold)
    d.text((x0 + bd / 2 - (bb[2] - bb[0]) / 2, y0 + bd / 2 - (bb[3] - bb[1]) / 2 - bb[1]),
           r["num"], font=f_bold, fill="white")

    tx = x0 + bd + 16
    d.text((tx, y0 - 4), r["cluster"], font=f_bold, fill="#1a1a1a")
    d.text((tx, y0 + 40), r["label"], font=f_reg, fill="#1a1a1a")
    n_label_lines = r["label"].count("\n") + 1
    ny = y0 + 40 + 40 * n_label_lines
    d.text((x0, ny), f"CLEAN: {r['clean']}", font=f_small, fill="#555555")
    d.text((x0, ny + 36), f"{r['rhea']}  ·  {r['share']} of members' top vote", font=f_small, fill="#777777")

    mol_y = ny + 95
    arrow_w = 160
    row_top = mol_y
    ay = mol_y + MOL_PX // 2

    x = x0
    canvas.paste(sub_img, (x, row_top)); x += MOL_PX

    d.line((x + 18, ay, x + arrow_w - 18, ay), fill="#1a1a1a", width=4)
    d.polygon([(x + arrow_w - 18, ay - 10), (x + arrow_w - 18, ay + 10), (x + arrow_w, ay)], fill="#1a1a1a")
    d.text((x + arrow_w / 2, ay + 16), "SAM → SAH", font=f_small, fill="#777777", anchor="ma")
    x += arrow_w

    canvas.paste(prod_img, (x, row_top)); x += MOL_PX

    return x, mol_y + MOL_PX + 24


def main():
    fdir = matplotlib.get_data_path() + "/fonts/ttf"
    f_bold = ImageFont.truetype(f"{fdir}/DejaVuSans-Bold.ttf", 38)
    f_reg = ImageFont.truetype(f"{fdir}/DejaVuSans.ttf", 32)
    f_small = ImageFont.truetype(f"{fdir}/DejaVuSans.ttf", 27)
    fonts = (f_bold, f_reg, f_small)

    TOP_PAD = 10
    COL_GAP, ROW_GAP = 90, 70

    cell_w = 34 + 12 + MOL_PX + 160 + MOL_PX + 30
    cell_h = 260 + MOL_PX
    nrows = (len(REACTIONS) + NCOLS - 1) // NCOLS
    W = cell_w * NCOLS + COL_GAP * (NCOLS - 1) + 20
    H = TOP_PAD + 10 + cell_h * nrows + ROW_GAP * (nrows - 1)

    canvas = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(canvas)

    for i, r in enumerate(REACTIONS):
        col, row = i % NCOLS, i // NCOLS
        x0 = 10 + col * (cell_w + COL_GAP)
        y0 = TOP_PAD + 10 + row * (cell_h + ROW_GAP)
        draw_reaction_cell(canvas, d, x0, y0, r, fonts)

    canvas.save(IO / "E_reaction_matches.png", dpi=(200, 200))
    print("wrote E_reaction_matches.png", canvas.size)


if __name__ == "__main__":
    main()

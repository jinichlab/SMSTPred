# Render the plain grey PRMT1 (1OR8 chain A) + faint SAM/SAH panel used in Fig 1 (needs PyMOL).
# Run from the repository root:   SMSTPRED_REPO=$PWD pymol -cq figures/Fig1_overview_schematic/render_prmt1_sam.pml
python
import os
from pymol import cmd

REPO = os.environ.get("SMSTPRED_REPO", ".")
SRC = f"{REPO}/7_motif_analysis/structures/_carrier_src"
OUT = f"{REPO}/figures/Fig1_overview_schematic/biorender_components"

cmd.set("ray_opaque_background", 0)
cmd.set("antialias", 2)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("ray_shadows", 0)
cmd.set("ambient", 0.35)
cmd.bg_color("white")

cmd.load(f"{SRC}/1OR8.cif", "m")
COF = "resn SAM+SAH+SFG"
# keep chain-A enzyme + the cofactor copy sitting in its pocket
cmd.remove(f"not (polymer.protein or {COF})")
cmd.remove("polymer.protein and not chain A")
cmd.remove(f"{COF} and not byres ({COF} within 6 of (chain A and polymer.protein))")

cmd.hide("everything")
cmd.show("cartoon", "polymer.protein")
cmd.color("grey70", "polymer.protein")

# faint SAM: thin muted sticks, pale-wheat carbons, standard heteroatom cues, slightly sheer
cmd.show("sticks", COF)
cmd.set("stick_radius", 0.16, COF)
cmd.color("wheat", f"{COF} and elem C")
cmd.color("skyblue", f"{COF} and elem N")
cmd.color("salmon", f"{COF} and elem O")
cmd.color("paleyellow", f"{COF} and elem S")
cmd.set("stick_transparency", 0.12, COF)

cmd.orient("polymer.protein")
cmd.zoom("all", 4)
cmd.ray(1500, 1500)
cmd.png(f"{OUT}/struct_1OR8_SAM.png", dpi=300)
print("cofactor atoms kept:", cmd.count_atoms(COF))
python end

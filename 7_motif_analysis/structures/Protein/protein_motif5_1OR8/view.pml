# Protein-5 THWKQTVFYM - Protein arginine N-methyltransferase (1OR8; X-ray; SAH + partner protein) - near substrate.  Open:  pymol view.pml
load protein_motif5_1OR8.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper, 1
set cartoon_ring_mode, 3
set cartoon_ring_finder, 1
set ray_shadows, 0
set valence, 1
show cartoon, polymer.protein
color 0xb0b0b0, polymer.protein
select motif, chain A and resi 292-301
color 0x238b45, motif
show sticks, motif
util.cnc motif
select partner, chain B
show cartoon, partner
color wheat, partner
select cof, (hetatm and not solvent)
show sticks, cof
color grey60, cof
util.cnc cof
set_color cofcarbon, [0.710, 0.396, 0.114]
select samcof, (cof and resn 0XU+NEA+QVR+SAH+SAM+SFG)
color cofcarbon, samcof
util.cnc samcof
distance contacts, motif, (partner or cof), 3.6, mode=2
hide labels, contacts
color black, contacts
set_title m, 1, "Protein-5 THWKQTVFYM - Protein arginine N-methyltransferase (1OR8; X-ray; SAH + partner protein) - near substrate"
orient
deselect
# >>> rotate to taste, then:  ray 2400, 1800;  png protein_motif5_1OR8.png, dpi=300

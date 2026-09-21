# Protein-3 MADCIYYE - Protein N-lysine methyltransferase M (7OAT; X-ray; SAH + partner protein) - near substrate/SAM.  Open:  pymol view.pml
load protein_motif3_7OAT.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper, 1
set cartoon_ring_mode, 3
set cartoon_ring_finder, 1
set ray_shadows, 0
set valence, 1
show cartoon, polymer.protein
color 0xb0b0b0, polymer.protein
select motif, chain C and resi 142-149
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
set_title m, 1, "Protein-3 MADCIYYE - Protein N-lysine methyltransferase M (7OAT; X-ray; SAH + partner protein) - near substrate/SAM"
orient
deselect
# >>> rotate to taste, then:  ray 2400, 1800;  png protein_motif3_7OAT.png, dpi=300

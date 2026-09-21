# Other-8 GVNFDLPHVIS - Scoulerine-9-O-methyltransferase 1 ( (6I6K; X-ray; SAH + SLX) - near SAM.  Open:  pymol view.pml
load other_motif8_6I6K.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper, 1
set cartoon_ring_mode, 3
set cartoon_ring_finder, 1
set ray_shadows, 0
set valence, 1
show cartoon, polymer.protein
color 0xb0b0b0, polymer.protein
select motif, chain A and resi 254-264
color 0xe6550d, motif
show sticks, motif
util.cnc motif
select cof, (hetatm and not solvent)
show sticks, cof
color grey60, cof
util.cnc cof
set_color cofcarbon, [0.710, 0.396, 0.114]
select samcof, (cof and resn 0XU+NEA+QVR+SAH+SAM+SFG)
color cofcarbon, samcof
util.cnc samcof
color yellow, resn SLX
util.cnc resn SLX
distance contacts, motif, (cof), 3.6, mode=2
hide labels, contacts
color black, contacts
set_title m, 1, "Other-8 GVNFDLPHVIS - Scoulerine-9-O-methyltransferase 1 ( (6I6K; X-ray; SAH + SLX) - near SAM"
orient
deselect
# >>> rotate to taste, then:  ray 2400, 1800;  png other_motif8_6I6K.png, dpi=300

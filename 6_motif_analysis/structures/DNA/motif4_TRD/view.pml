# LTPREAARLQGFPDDYTF - cofactor-proximal.  Open:  pymol motif4_TRD/view.pml
load motif4_TRD.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper, 1
set cartoon_ring_mode, 3
set cartoon_ring_finder, 1
set ray_shadows, 0
set valence, 1
# protein motif region as grey ribbon
show cartoon, chain A
color grey70, chain A
# the motif: colour + sticks + CPK heteroatoms
select motif, chain A and resi 274-291
color 0x6a51a3, motif
show sticks, motif
util.cnc motif
# DNA duplex as tan base-slabs
select dna, chain B+C
show cartoon, dna
color wheat, dna
# cofactor + target base as sticks
select cof, resn SAH+C36
show sticks, cof
color grey60, cof
util.cnc cof
# polar contacts (H-bonds) motif <-> DNA/cofactor
distance contacts, motif, (dna or cof), 3.6, mode=2
hide labels, contacts
color black, contacts
set_title m, 1, "motif4 LTPREAARLQGFPDDYTF - cofactor-proximal"
orient motif
deselect
# >>> rotate to taste, then:   ray 2400, 1800;   png motif4.png, dpi=300

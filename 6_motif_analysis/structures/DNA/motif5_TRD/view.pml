# KQIGNAVPPPVARA - DNA-cofactor junction.  Open:  pymol motif5_TRD/view.pml
load motif5_TRD.pdb, m
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
select motif, chain A and resi 300-315
color 0xce1256, motif
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
set_title m, 1, "motif5 KQIGNAVPPPVARA - DNA-cofactor junction"
orient motif
deselect
# >>> rotate to taste, then:   ray 2400, 1800;   png motif5.png, dpi=300

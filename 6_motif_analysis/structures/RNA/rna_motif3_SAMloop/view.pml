# RNA Trm5 SAM-binding DxGxGxG loop.  Open: pymol view.pml
load rna_motif3_SAMloop.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper, 1
set cartoon_ring_mode, 3
set cartoon_ring_finder, 1
set ray_shadows, 0
show cartoon, chain A
color grey70, chain A
select motif, chain A and resi 198-212
color 0x2171b5, motif
show sticks, motif
util.cnc motif
label motif and name CA, "%s%s" % (one_letter[resn], resi)
set label_size, 12
select rna, chain B
show cartoon, rna
color wheat, rna
select cof, resn SAM
show sticks, cof
color grey60, cof
util.cnc cof
distance contacts, motif, (rna or cof), 3.6, mode=2
hide labels, contacts
color black, contacts
orient motif
deselect
# >>> rotate, then:  ray 2400,1800;  png rna_motif3_SAMloop.png, dpi=300

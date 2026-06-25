# RNA Trm5 N-terminal region.  Open: pymol view.pml
load rna_motif6.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper, 1
set cartoon_ring_mode, 3
set cartoon_ring_finder, 1
set ray_shadows, 0
show cartoon, chain A
color grey70, chain A
select motif, chain A and resi 45-55
color 0x6a51a3, motif
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
# >>> rotate, then:  ray 2400,1800;  png rna_motif6.png, dpi=300

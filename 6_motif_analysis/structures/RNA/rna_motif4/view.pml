# RNA RsmF region (apo).  Open: pymol view.pml
load rna_motif4.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set ray_shadows,0
show cartoon, chain A
color grey80, chain A
select motif, chain A and resi 85-101
color 0x6a51a3, motif
show sticks, motif
util.cnc motif
label motif and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,12
orient chain A
deselect
# >>> rotate, then: ray 2400,1800; png rna_motif4.png, dpi=300

# RNA RsmF m5C catalytic PC (apo).  Open: pymol view.pml
load rna_motif1_m5Ccat.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set ray_shadows,0
show cartoon, chain A
color grey80, chain A
select motif, chain A and resi 172-185
color 0xcb181d, motif
show sticks, motif
util.cnc motif
label motif and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,12
orient chain A
deselect
# >>> rotate, then: ray 2400,1800; png rna_motif1_m5Ccat.png, dpi=300

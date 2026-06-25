# Protein HemK/PrmC region (AlphaFold).  Open: pymol view.pml
load protein_motif5_HemK.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set ray_shadows,0
show cartoon, chain A
color grey80, chain A
select motif, chain A and resi 10-19
color 0xe6550d, motif
show sticks, motif
util.cnc motif
label motif and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,12
orient chain A
deselect
# >>> rotate, then: ray 2400,1800; png protein_motif5_HemK.png, dpi=300

# Protein HemK catalytic NPPY (AlphaFold).  Open: pymol view.pml
load protein_motif8_HemK_NPPY.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set ray_shadows,0
show cartoon, chain A
color grey80, chain A
select motif, chain A and resi 188-195
color 0x8c564b, motif
show sticks, motif
util.cnc motif
label motif and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,12
orient chain A
deselect
# >>> rotate, then: ray 2400,1800; png protein_motif8_HemK_NPPY.png, dpi=300

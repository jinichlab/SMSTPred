# Other salicylate MTase homolog (SAH+salicylate) [HOMOLOG].  Open: pymol view.pml
load other_motif4.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set ray_shadows,0
show cartoon, chain X
color grey80, chain X
select motif, chain X and resi 128-145
color 0xcb181d, motif
show sticks, motif
util.cnc motif
label motif and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,12
show sticks, resn SAH
color grey50, resn SAH
util.cnc resn SAH
color yellow, resn SAL
util.cnc resn SAL
distance contacts, motif, (resn SAH or resn SAL), 3.6, mode=2
hide labels, contacts
color black, contacts
orient chain X
deselect
# >>> rotate, then: ray 2400,1800; png other_motif4.png, dpi=300

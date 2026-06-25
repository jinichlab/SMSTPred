# Other caffeic-acid OMT homolog (SAH+ferulate) [HOMOLOG].  Open: pymol view.pml
load other_motif1.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set ray_shadows,0
show cartoon, chain E
color grey80, chain E
select motif, chain E and resi 285-295
color 0x2171b5, motif
show sticks, motif
util.cnc motif
label motif and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,12
show sticks, resn SAH
color grey50, resn SAH
util.cnc resn SAH
color yellow, resn FER
util.cnc resn FER
distance contacts, motif, (resn SAH or resn FER), 3.6, mode=2
hide labels, contacts
color black, contacts
orient chain E
deselect
# >>> rotate, then: ray 2400,1800; png other_motif1.png, dpi=300

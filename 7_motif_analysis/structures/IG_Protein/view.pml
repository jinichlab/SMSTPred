# Q54S83 IG on apo 5ZTZ.  Open: pymol view.pml
load Q54S83_IG.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set ray_shadows,0
show cartoon, chain A
color grey80, chain A
select igtop, chain A and resi 384+281+209+382+374+383+385+279+375+191
color 0xcb181d, igtop
show sticks, igtop
util.cnc igtop
label igtop and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,14
orient chain A
deselect
# >>> rotate, then: ray 2400,1800; png Q54S83_IG.png, dpi=300

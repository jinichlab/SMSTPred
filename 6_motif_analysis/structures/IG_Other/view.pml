# P0DPA9 IG on 8PB3 (SAH+substrate).  Open: pymol view.pml
load P0DPA9_IG.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set ray_shadows,0
show cartoon, chain A
color grey80, chain A
select igtop, chain A and resi 220+86+78+127+222+125+22+213+11+13
color 0xcb181d, igtop
show sticks, igtop
util.cnc igtop
label igtop and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,14
show sticks, resn SAH+XP6
util.cnc resn SAH+XP6
color yellow, resn XP6 and elem C
orient chain A
deselect
# >>> rotate, then: ray 2400,1800; png P0DPA9_IG.png, dpi=300

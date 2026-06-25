# Q9V2G1 IG on tRNA-bound 5WT1.  Open: pymol view.pml
load Q9V2G1_IG.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set cartoon_ring_mode,3
set ray_shadows,0
show cartoon, chain A
color grey80, chain A
select igtop, chain A and resi 196+139+189+194+258+270+318+188+250+111
color 0xcb181d, igtop
show sticks, igtop
util.cnc igtop
label igtop and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,14
select rna, chain C+F
show cartoon, rna
color wheat, rna
select cof, resn SAH
show sticks, cof
color grey60, cof
util.cnc cof
orient chain A
deselect
# >>> rotate, then: ray 2400,1800; png Q9V2G1_IG.png, dpi=300

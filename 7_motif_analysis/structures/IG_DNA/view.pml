# T4 Dam (P04392) IG on DNA-bound 1YF3.  Open: pymol view.pml
load P04392_IG.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper, 1
set cartoon_ring_mode, 3
set ray_shadows, 0
# whole protein one colour (grey)
show cartoon, chain A
color grey80, chain A
# top-IG residues highlighted as red sticks (CPK heteroatoms) + labels
select igtop, chain A and resi 127+117+114+126+124+123+128+122+77+10
color 0xcb181d, igtop
show sticks, igtop
util.cnc igtop
label igtop and name CA, "%s%s" % (one_letter[resn], resi)
set label_size, 14
# DNA duplex (tan slabs)
show cartoon, chain C+D
color wheat, chain C+D
# cofactor AdoHcy
select cof, resn SAH
show sticks, cof
color grey50, cof
util.cnc cof
orient chain A
deselect
# >>> rotate, then:  ray 2400,1800;  png P04392_IG.png, dpi=300

# Protein PRMT5 region (+MEP50 partner).  Open: pymol view.pml
load protein_motif6_PRMT5.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set ray_shadows,0
show cartoon, chain A
color grey80, chain A
select motif, chain A and resi 45-53
color 0xce1256, motif
show sticks, motif
util.cnc motif
label motif and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,12
select partner, chain B+C
show cartoon, partner
color wheat, partner
select cof, resn 0XU
show sticks, cof
color grey50, cof
util.cnc cof
distance contacts, motif, (chain B+C or resn 0XU), 3.6, mode=2
hide labels, contacts
color black, contacts
orient chain A
deselect
# >>> rotate, then: ray 2400,1800; png protein_motif6_PRMT5.png, dpi=300

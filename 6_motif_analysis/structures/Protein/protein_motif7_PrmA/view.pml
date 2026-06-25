# Protein PrmA region (+L11 substrate).  Open: pymol view.pml
load protein_motif7_PrmA.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set ray_shadows,0
show cartoon, chain A
color grey80, chain A
select motif, chain A and resi 98-114
color 0x00868b, motif
show sticks, motif
util.cnc motif
label motif and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,12
select partner, chain B
show cartoon, partner
color wheat, partner
distance contacts, motif, (chain B), 3.6, mode=2
hide labels, contacts
color black, contacts
orient chain A
deselect
# >>> rotate, then: ray 2400,1800; png protein_motif7_PrmA.png, dpi=300

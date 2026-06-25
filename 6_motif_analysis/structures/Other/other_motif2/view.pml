# Other PFOMT homolog (SAH, no substrate) [HOMOLOG].  Open: pymol view.pml
load other_motif2.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper,1
set ray_shadows,0
show cartoon, chain A
color grey80, chain A
select motif, chain A and resi 99-111
color 0x6a51a3, motif
show sticks, motif
util.cnc motif
label motif and name CA, "%s%s"%(one_letter[resn],resi)
set label_size,12
show sticks, resn SAH
color grey50, resn SAH
util.cnc resn SAH
distance contacts, motif, (resn SAH), 3.6, mode=2
hide labels, contacts
color black, contacts
orient chain A
deselect
# >>> rotate, then: ray 2400,1800; png other_motif2.png, dpi=300

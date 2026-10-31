# Metabolite-5 GGKERTEEEWKKL - (RS)-norcoclaurine 6-O-methyltransfe (5ICE; X-ray; SAH + 2H4) - distal.  Open:  pymol view.pml
load metabolite_motif5_5ICE.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper, 1
set cartoon_ring_mode, 3
set cartoon_ring_finder, 1
set ray_shadows, 0
set valence, 1
show cartoon, polymer.protein
color 0xb0b0b0, polymer.protein
select motif, chain A and resi 312-324
color 0xe6550d, motif
show sticks, motif
util.cnc motif
select cof, (hetatm and not solvent)
show sticks, cof
color grey60, cof
util.cnc cof
set_color cofcarbon, [0.710, 0.396, 0.114]
select samcof, (cof and resn 0XU+NEA+QVR+SAH+SAM+SFG)
color cofcarbon, samcof
util.cnc samcof
color yellow, resn 2H4
util.cnc resn 2H4
distance contacts, motif, (cof), 3.6, mode=2
hide labels, contacts
color black, contacts
set_title m, 1, "Metabolite-5 GGKERTEEEWKKL - (RS)-norcoclaurine 6-O-methyltransfe (5ICE; X-ray; SAH + 2H4) - distal"
orient
deselect
# >>> rotate to taste, then:  ray 2400, 1800;  png metabolite_motif5_5ICE.png, dpi=300

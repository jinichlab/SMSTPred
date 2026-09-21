# DNA-1 CAGFPCQAFSISG - Type II methyltransferase M.HhaI (M. (10MH; X-ray; SAH + DNA) - near substrate/SAM.  Open:  pymol view.pml
load dna_motif1_10MH.pdb, m
bg_color white
hide everything
set cartoon_side_chain_helper, 1
set cartoon_ring_mode, 3
set cartoon_ring_finder, 1
set ray_shadows, 0
set valence, 1
show cartoon, polymer.protein
color 0xb0b0b0, polymer.protein
select motif, chain A and resi 76-88
color 0x6a51a3, motif
show sticks, motif
util.cnc motif
set_color dnatan, [0.612,0.478,0.357]
select nuc, chain C
show cartoon, nuc
hide sticks, nuc
color dnatan, nuc
select cof, (hetatm and not solvent and not chain C)
show sticks, cof
color grey60, cof
util.cnc cof
set_color cofcarbon, [0.710, 0.396, 0.114]
select samcof, (cof and resn 0XU+NEA+QVR+SAH+SAM+SFG)
color cofcarbon, samcof
util.cnc samcof
distance contacts, motif, (cof), 3.6, mode=2
hide labels, contacts
color black, contacts
set_title m, 1, "DNA-1 CAGFPCQAFSISG - Type II methyltransferase M.HhaI (M. (10MH; X-ray; SAH + DNA) - near substrate/SAM"
orient
deselect
# >>> rotate to taste, then:  ray 2400, 1800;  png dna_motif1_10MH.png, dpi=300

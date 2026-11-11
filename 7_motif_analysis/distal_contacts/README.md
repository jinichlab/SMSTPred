# What do the distal motifs touch?

For every PDB entry of every carrier of a distal motif (`../structure_search/motif_structure_candidates.csv`), `distal_motif_contacts.py` locates the motif on the
carrier chain by its matched sequence and lists its heavy-atom contacts (<= 4.5 A) with other polymer chains (other entities, copies of the carrier, nucleic acids),
ligands and ions, and the surface buried by partners (SASA of the motif residues alone minus in the assembly). It downloads missing mmCIF files from RCSB.
`run_distal_contacts.sh` runs it for the distal motifs of the rebuilt release (DNA-4, Metabolite-7, Protein-6, RNA-5, RNA-7) and of the archived one
(DNA-4, Metabolite-3, Protein-6, RNA-6, RNA-8); set `REL_NEW` and `REL_OLD`. Outputs: `out_*/contacts_per_structure.csv`, `out_*/summary.txt`.

Caveats: contacts are read from the asymmetric unit (a crystal contact is not told from a biological interface); carriers include enzymes of other classes
(FIMO hits at p < 1e-4); "motif_not_resolved" entries lack the motif in the model.

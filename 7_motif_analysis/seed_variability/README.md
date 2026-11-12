# STREME seed variability

STREME (`--seed 0..3`) was rerun on the inputs of two releases (archived: seeds 0-3; rebuilt: seeds 0-2), all other settings as in `../run_streme.sh`,
to ask how much the motif set depends on the seed. `run_streme_array.sh` (SLURM array, one class x variant per task), then `run_scan.sh` (FIMO cross scan of every
variant's motifs, `summarise.py`). `variants.txt` lists the variants as `name source-motif_analysis-folder seed` (`$ARCHIVED_MOTIF_DIR` and `$REBUILT_MOTIF_DIR` are the two `motif_analysis` folders; export them before submitting).

Results: `variability_counts.csv` (motifs with ratio_vs_bg >= 2 per variant and class), `variability_motifs.csv` (every motif's consensus),
`variability_overlap.csv` (motifs of a variant that reappear in the seed-0 run, consensus up to 3 mismatches; a strict, position-wise criterion).
Summary: 22-25 of 32 motifs are class-enriched (>= 2-fold) depending on seed and release; which motifs are found changes more (12-21 of 32 reappear).

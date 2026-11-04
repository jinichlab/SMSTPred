# Fig 5 as power-law fits (release version)

Adapted from the `scaling-laws-fig5` branch (A. Jinich) so it reads a release: `scaling_loglog.py --dir <release>/fig5`, run by `run_fig5.sh`
after `build_proteome_vs_sam_table.py`. Reads the release's per-proteome table (proteome size = proteins in the reference proteome, `n_sam` =
IPR029063 proteins, `n_Metabolite` = ensemble calls at 0.5) and lineage from `data/selected_1259_organisms.csv`. The robustness analyses of that
branch (`analysis_scaling_laws/`) were not ported; they still read the old release.

## One setting, fixed
- log bins of 0.15 decades from 10^3; weighted non-linear least squares on bin means (n = exp(a)(g/1e4)^alpha)
- bins with fewer than **5** proteomes dropped (10, the 2017 NAR protocol, leaves Viridiplantae a single bin; the strict version is `strict_protocol` in the numbers file)
- clade bootstrap, 2,000 replicates, clades = lineage prefix of depth 7 (bacteria, family level) / 8 (eukaryotes); depth 6 reported as a sensitivity
- archaea excluded (6 proteomes)

## Release 2026_03 result (fig5/fig5_scaling_numbers.csv)
| panel | alpha lineage | alpha rest | delta alpha [95% CI] | level at equal size [95% CI] | old OLS ratio |
|---|---|---|---|---|---|
| A Actinomycetota total | 1.42 | 0.69 | +0.73 [+0.27, +0.97] | x1.28 [1.12, 1.49] at 5,987 | 2.21 |
| B Actinomycetota Metabolite | 1.89 | 1.26 | +0.63 [-0.25, +1.66] | x1.85 [1.48, 2.48] | 2.98 |
| C Viridiplantae total | 0.69 | 0.75 | -0.06 [-0.43, +0.34] | x1.64 [0.86, 2.13] at 34,300 | 1.84 |
| D Viridiplantae Metabolite | 0.75 | 0.64 | +0.10 [-0.57, +0.94] | x3.25 [0.90, 6.95] | 8.45 |

Only A shows a different exponent. B is a level difference (exponents indistinguishable), C and D are level differences whose intervals reach 1 for
plants. The old "x-fold faster" numbers are ratios of OLS slopes of raw counts and are not rates.

## Caveats (from the branch's README, still true)
A single power law fits poorly (clades sit at different offsets; reduced chi-square 4 to 40), binned and unbinned estimators can disagree,
the proteomes were chosen by rule and are not independent, and no phylogenetic regression was done.

## Text (as used in the draft)
**Figure 5 caption.** Proteome size and SAM-methyltransferase counts as power laws n = lambda g^alpha. Predicted SAM-methyltransferases (A, C all; B, D with a Metabolite call) versus proteome size on log axes for (A, B) Actinomycetota against other bacteria and (C, D) Viridiplantae against other eukaryotes. Small markers are proteomes, open markers the mean (+/- s.e.) of logarithmic size bins (0.15 decades, at least five proteomes), lines weighted non-linear least-squares fits. Boxes give each exponent alpha, the difference delta alpha, and the level difference at equal proteome size (lineage over rest at the lineage's median size), with 95% intervals from 2,000 bootstrap replicates that resample whole clades. Only the Actinomycetota total count has a steeper exponent (A); in B to D the exponents cannot be told apart and the lineage sits at a higher level.

**Results.** Actinomycetota gain SAM-methyltransferases more steeply with proteome size than other bacteria (exponent 1.42 vs 0.69; difference 0.73, 95% CI 0.27 to 0.97), over a range of only about 1,900 to 11,000 proteins. For metabolite-substrate enzymes the exponents do not differ (1.89 vs 1.26; difference 0.63, CI -0.25 to 1.66), but Actinomycetota carry 1.9 times as many as other bacteria of equal proteome size (CI 1.5 to 2.5). Viridiplantae carry 1.6 times as many SAM-methyltransferases as other eukaryotes of equal proteome size (CI 0.9 to 2.1) and 3.2 times as many metabolite-substrate enzymes (CI 0.9 to 6.9), with no detectable difference in exponent (-0.06 and +0.10), so the plant excess is a level difference whose size is not resolved by these data.

**Discussion line.** Both lineages carry more small-molecule methyltransferases than other organisms of the same proteome size, and in Actinomycetota the total count also rises more steeply with proteome size (Figure 5); for plants the excess is a higher level at equal size rather than a different exponent, and it is not precisely resolved.

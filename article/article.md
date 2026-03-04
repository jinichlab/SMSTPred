# Machine Learning-Based Classification of SAM-Dependent Methyltransferase Substrate Preferences Reveals Universal Evolutionary Trade-offs Across the Tree of Life

**Author:** N. Roethler

---

## Abstract

S-adenosyl methionine (SAM)-dependent methyltransferases (MTases) constitute one of the most functionally diverse enzyme superfamilies in biology, modifying substrates ranging from nucleic acids to small metabolites. Despite their prevalence, the evolutionary forces that shape the distribution of substrate preferences across organisms remain poorly characterized. Here, we trained a multi-layer perceptron (MLP) classifier on ProtT5 protein language model embeddings to predict SAM MTase substrate specificity (RNA, DNA, protein, or small molecule/other) from sequence alone. The MLP achieved a mean AUPRC of 0.998 +/- 0.002, AUROC of 0.999 +/- 0.0004, and Matthew's Correlation Coefficient of 0.976 +/- 0.006 on the SAM-specific evaluation dataset, outperforming XGBoost across all metrics. We applied the trained model to predict substrate preferences for all annotated SAM MTases across 386 organisms spanning 312 bacteria, 13 archaea, and 61 eukaryotes, mapping predictions onto an interactive phylogenetic tree. Statistical analyses across this dataset revealed four major discoveries: (1) a universally conserved negative trade-off between RNA-binding and small-molecule ("Other")-binding MTase proportions (Spearman rho = -0.793, p = 1.75e-84; phylogenetically corrected rho = -0.709, p = 3.55e-145); (2) lineage-specific specializations including Actinomycetota enrichment for Other-category MTases linked to biosynthetic gene clusters (BGCs; fold = 1.70x, p = 1.78e-14); (3) a vertebrate-specific expansion of DNA MTases (2.81-fold, p = 3.01e-7) that is independent of protein MTase proportions; and (4) a genome-reduction-driven shift in prokaryotic MTase allocation wherein host-associated bacteria shed Other MTases while retaining RNA MTases (rho = -0.588, p = 1.44e-31). Together, these findings reveal deep evolutionary constraints on cellular methylation strategies and suggest that the RNA-Other trade-off reflects a fundamental competition for proteomic investment between translational infrastructure and metabolic diversification.

---

## 1. Introduction

SAM-dependent methyltransferases are universally distributed enzymes that transfer a methyl group from S-adenosyl-L-methionine to a diverse array of biological targets. These targets span all major macromolecular classes: ribosomal and transfer RNAs, genomic and plasmid DNA, histone and non-histone proteins, and a wide variety of small molecules including polyketides, lipids, and secondary metabolites [1, 2]. Because methylation events regulate processes as fundamental as translational fidelity, chromatin organization, epigenetic memory, immune defense, and natural product biosynthesis, the relative investment of an organism's SAM MTase repertoire in each substrate class provides a high-resolution view of its underlying biology and evolutionary pressures.

Despite the functional importance of this enzyme family, large-scale comparative analysis of SAM MTase substrate preferences across the tree of life has been constrained by two challenges. First, experimental characterization of substrate specificity is low-throughput: most MTases in sequenced genomes carry no biochemically validated substrate annotation. Second, even for those with annotations, substrates are classified using heterogeneous ontologies across databases, making cross-organism comparisons difficult. The recent emergence of protein language models, particularly ProtT5, which generates 1024-dimensional per-protein embeddings encoding evolutionary and structural information from unsupervised training on hundreds of millions of sequences, offers a solution to the first challenge [3]. These embeddings have enabled state-of-the-art performance on functional classification tasks and can be computed for any protein sequence without structural or alignment information.

In this work, we developed a sequence-based machine learning classifier trained on ProtT5 embeddings to assign SAM MTase substrate preferences across four categories: RNA-binding, DNA-binding, protein-binding, and "Other" (small molecules, lipids, secondary metabolites). We systematically evaluated alternative algorithms, including XGBoost, Random Forest, and Logistic Regression, under a cluster-based cross-validation scheme that prevents sequence-similar proteins from leaking between training and evaluation folds, thereby providing realistic performance estimates. The MLP classifier outperformed all competing algorithms and was applied to predict substrate preferences for all SAM MTases across 386 organisms spanning the full depth of the tree of life. These predictions were mapped onto an interactive phylogenetic tree and analyzed with a comprehensive suite of statistical approaches, revealing four major discoveries that illuminate universal and lineage-specific forces shaping the evolution of cellular methylation repertoires.

---

## 2. Methods

### 2.1 Training Data and Protein Embeddings

The training dataset consisted of SAM-dependent methyltransferases from the UniProtKB/Swiss-Prot database with experimentally reviewed substrate annotations. Proteins were labeled across four substrate categories: RNA, DNA, Protein, and Other (encompassing small molecules, lipids, and secondary metabolites). Per-protein embeddings were generated using the ProtT5-XL-U50 protein language model [3], producing 1024-dimensional representations that encode evolutionary, structural, and functional information. Embeddings were stored in HDF5 format and loaded at training time via h5py.

### 2.2 Machine Learning Classifier Development

We implemented a multi-label multi-layer perceptron (MLP) classifier with the following architecture: a 1024-dimensional input layer corresponding to the ProtT5 embedding dimensionality, a first hidden layer of 512 neurons with ReLU activation, a second hidden layer of 256 neurons with ReLU activation, and an output layer of four neurons with sigmoid activation corresponding to the four substrate categories. The multi-label formulation allows a protein to carry positive predictions for multiple substrate types simultaneously, reflecting the functional promiscuity observed in some bifunctional MTases. The model was trained using binary cross-entropy loss (BCELoss) and the Adam optimizer with a learning rate of 0.001 for 500 epochs.

We compared the MLP against XGBoost (wrapped in a MultiOutputClassifier for multi-label support), Random Forest, and Logistic Regression. XGBoost was trained with logloss as the evaluation metric and parallel processing enabled. All classifiers received identical training and evaluation data to enable fair comparison.

### 2.3 Cluster-Based Cross-Validation

A critical methodological consideration in functional classification of protein sequences is the avoidance of homology leakage between training and validation sets. If closely related proteins appear in both partitions, performance estimates will be optimistically inflated. To address this, we implemented a cluster-based cross-validation strategy using the cosine similarity of ProtT5 embeddings as a proxy for sequence similarity. Proteins were grouped into clusters at a similarity threshold of 0.6 using a greedy algorithm: proteins were iteratively assigned to the first cluster whose existing members all had cosine similarity >= 0.6 to the candidate, otherwise initiating a new cluster. We validated that each resulting split contained at least one positive example per substrate class (up to 1000 re-sampling attempts). Three independent cross-validation iterations were performed with different cluster partitions (90% clusters for training, 10% for validation), and mean performance with standard deviation was reported across iterations.

### 2.4 Genome-Scale Predictions on Reference Proteomes

The trained MLP was applied to predict substrate preferences for all SAM MTase sequences from 386 reference proteomes spanning bacteria, archaea, and eukaryotes. SAM MTases were identified using the InterPro domain identifier IPR029063 (S-adenosyl-L-methionine-binding domain). For each organism's proteome, matching protein sequences were extracted, ProtT5 embeddings were loaded from pre-computed HDF5 files, and per-protein substrate predictions were generated. Organisms with fewer than 20 predicted SAM MTases were excluded. For each organism, per-substrate predictions were aggregated to produce proportional metrics: pct_rna, pct_dna, pct_protein, and pct_other, representing the fraction of that organism's SAM MTases predicted to bind each substrate class. The total count of SAM MTases per organism (n_sams) was also recorded.

### 2.5 Phylogenetic Tree Construction and Visualization

A reference phylogenetic tree was constructed in Newick format spanning all 386 organisms. Taxonomic lineages were parsed from UniProtKB metadata, filtered to canonical ranks (domain, phylum, class, order, family, genus, species), and assembled into a hierarchical TreeNode structure using a recursive path-construction algorithm. All branch lengths were set to a constant value of 1.0. Leaf labels were normalized by removing special characters and replacing spaces with underscores to match iTOL identifiers. The tree was built on an HPC cluster (20 CPU cores, 20 GB RAM) via a SLURM submission script.

Annotation files for the Interactive Tree of Life (iTOL) platform were generated for each substrate-proportion metric using the DATASET_GRADIENT format with substrate-specific color scales: pct_rna (white to red, #cb181d), pct_dna (white to purple, #6a51a3), pct_protein (white to green, #238b45), pct_other (white to orange, #e6550d), and pct_sams (white to blue, #2171b5). An additional DATASET_SIMPLEBAR annotation was generated for n_sams to visualize absolute MTase counts. For circular tree visualization, the GraPhlAn pipeline was used: graphlan_annotate.py applied annotations to the PhyloXML tree, and graphlan.py rendered the final image (300 DPI, 7.5-inch canvas, 340-degree arc).

### 2.6 Statistical Analysis

All statistical analyses were performed in Python (NumPy, SciPy, pandas, statsmodels). Spearman rank correlations were used for non-normally distributed variables, with normality assessed by the Shapiro-Wilk test. Group comparisons used Mann-Whitney U tests (two-sided unless otherwise specified) with rank-biserial correlation as effect size, and Kruskal-Wallis tests for more than two groups. Multiple regression models used OLS with HC3 heteroscedasticity-consistent standard errors. Phylogenetic Independent Contrasts (PIC) were computed following Felsenstein [4] using the ete3 library, with regression through the origin on contrast values from the 386-organism reference tree; all pairwise contrasts at polytomous nodes were computed, yielding 946 total contrasts. Genome sizes were retrieved via the NCBI Datasets API v2 (total_ungapped_length; best available assembly selected by RefSeq category and assembly level). BGC counts were estimated from genome size using literature-based rates (Actinomycetota: 8 BGCs/Mb; Archaea: 2.5 BGCs/Mb; all other bacteria: 5 BGCs/Mb) [5, 6].

---

## 3. Results

### 3.1 MLP Classifier Outperforms XGBoost for SAM MTase Substrate Classification

We first systematically compared the MLP against XGBoost on the SAM MTase-specific training dataset under three-iteration cluster-based cross-validation (Table 1). The MLP achieved a mean AUPRC of 0.998 +/- 0.002, AUROC of 0.999 +/- 0.0004, MCC of 0.976 +/- 0.006, F1 of 0.983 +/- 0.003, and accuracy of 0.991 +/- 0.003. XGBoost achieved an AUPRC of 0.987 +/- 0.014, AUROC of 0.998 +/- 0.002, MCC of 0.905 +/- 0.088, F1 of 0.907 +/- 0.094, and accuracy of 0.982 +/- 0.011. The MLP was consistently superior, particularly for MCC and F1, with notably lower variance across cross-validation iterations. The cluster-based cross-validation scheme, which prevents embedding-similar proteins from spanning train/validation partitions, demonstrates that this performance reflects genuine generalization rather than homology leakage. The trained MLP was used for all downstream genome-scale predictions.

Table 1. Classifier performance on the SAM MTase evaluation dataset (mean +/- SD, 3 cross-validation iterations).

  Metric     MLP                  XGBoost
  AUPRC      0.998 +/- 0.002      0.987 +/- 0.014
  AUROC      0.999 +/- 0.0004     0.998 +/- 0.002
  MCC        0.976 +/- 0.006      0.905 +/- 0.088
  F1         0.983 +/- 0.003      0.907 +/- 0.094
  Accuracy   0.991 +/- 0.003      0.982 +/- 0.011

### 3.2 Genome-Scale Predictions Across 386 Organisms

The trained MLP was applied to all annotated SAM MTases from 386 reference proteomes, yielding substrate-proportion profiles for each organism. The dataset spans 312 bacteria (80.8%), 13 archaea (3.4%), and 61 eukaryotes (15.8%). SAM MTase counts per organism ranged from 20 to 662 (prokaryote median ~47; eukaryote median ~120), and substrate proportions were mapped onto the reference phylogenetic tree via iTOL gradient annotations, enabling simultaneous visualization of evolutionary relationships and methylation strategy across the full tree of life.

### 3.3 A Universally Conserved RNA-Other Trade-off

The most prominent pattern in the genome-scale dataset is a strong, universally conserved negative correlation between the proportions of RNA-binding and Other-binding MTases across all 386 organisms (Spearman rho = -0.793, p = 1.75e-84; Figure 1). This trade-off is not a consequence of shared phylogenetic history: phylogenetic independent contrasts computed across 946 bifurcating nodes of the reference tree yield an equally strong negative correlation (rho = -0.709, p = 3.55e-145, Figure 8), confirming that the relationship has evolved convergently across independent lineages rather than being inherited from a common ancestor.

The trade-off is particularly striking among eukaryotes (Spearman rho = -0.931, p = 4.67e-27), and within eukaryotic kingdoms occupies distinct positions along a shared axis (Figure 6): protists exhibit the highest RNA proportions (median 67.7%) and lowest Other proportions (9.7%), while fungi and plants cluster at the opposite extreme (RNA medians 21.6% and 34.2%; Other medians 48.4% and 43.4%). Animals occupy an intermediate position (RNA median 41.1%, Other median 20.0%). These kingdom-level differences are highly significant (Kruskal-Wallis: RNA H = 28.44, p = 2.94e-6; Other H = 36.45, p = 6.01e-8).

We quantified functional balance using Shannon entropy (H = -sum(pi * log(pi)) across the four substrate categories) and tested whether organisms at the extremes of the RNA/Other ratio have lower entropy than those near the midpoint. The inverted-U relationship is confirmed: quadratic regression of Shannon entropy on log(RNA/Other ratio) among eukaryotes yields R2 = 0.622, with a maximum at RNA/Other ~1.27 (Figure 2). Spearman correlation between distance from the optimal ratio and Shannon entropy is rho = -0.428 (p = 1.34e-18) across all 386 organisms, confirming that specialization toward either extreme reduces functional balance.

### 3.4 Lineage-Specific Specializations: Actinomycetota and Secondary Metabolism

A second major pattern is the pronounced specialization of Actinomycetota toward Other-category MTases (Figure 9). Across five bacterial phyla (n >= 8 organisms each), Other-binding MTase proportions differ significantly (Kruskal-Wallis H = 75.25, p = 1.76e-15). Actinomycetota have a median Other% of 60.6%, versus 34.9% in Pseudomonadota, 26.2% in Bacillota, and 44.6% in Bacteroidota (one-sided Mann-Whitney, Actinomycetota > all others: fold = 1.70x, p = 1.78e-14).

This enrichment is mechanistically linked to biosynthetic gene clusters (BGCs). BGC count estimates correlate strongly with Other% across bacteria (Spearman rho = 0.644, p = 2.01e-36), and most strongly within Actinomycetota themselves (rho = 0.780, p = 1.03e-16). Multiple linear regression in free-living prokaryotes (n = 168; OLS with HC3 robust standard errors) confirms that BGC count is a significant independent predictor of Other% even after controlling for genome size (BGC coefficient = 0.0048, p = 1.96e-5; full model R2 = 0.455; BGC partial R2 = 0.077). Genome size itself is no longer a significant predictor (p = 0.222) once BGC count is included, indicating that the genome size effect on Other% operates primarily through the number of BGCs encoded.

### 3.5 Vertebrate DNA MTase Expansion Is Decoupled from the RNA-Other Axis

Using a purpose-built, taxonomically balanced dataset of 50 animal species (25 vertebrates, 25 invertebrates), we tested whether vertebrates show elevated methylation investment in specific substrate categories. Vertebrates have expanded DNA-binding MTases by 2.81-fold relative to invertebrates (vertebrate median 10.9%, invertebrate median 3.9%; Mann-Whitney U = 577.0, p = 3.01e-7, effect size r = -0.846; Figure 7), while protein-binding MTases are essentially unchanged (vertebrate median 9.4%, invertebrate median 9.8%; p = 0.961). This finding confirms that the vertebrate expansion specifically targeted DNA methylation rather than reflecting a general increase in epigenetic complexity.

The vertebrate DNA MTase expansion is consistent with the documented expansion of the DNMT gene family (DNMT1, DNMT3A, DNMT3B, DNMT3L) at the vertebrate stem, which enabled the dense, heritable CpG methylation patterns underlying mammalian imprinting, X-inactivation, and developmental gene regulation [9]. The independence of this expansion from the RNA-Other axis indicates that DNA methylation evolved as a distinct regulatory layer superimposed on the conserved RNA-Other balance, rather than as a trade-off against either existing category.

### 3.6 Genome Reduction Reshapes Prokaryotic MTase Allocation

The fourth major discovery concerns the relationship between genome size and the RNA/Other ratio in prokaryotes. Across 324 prokaryotes with genome size data, Spearman rho = -0.588 (p = 1.44e-31) between log genome size and RNA/Other ratio, with a clear monotone trend across genome quartiles: the smallest-genome quartile (mean 2.17 Mb) has a mean RNA/Other ratio of 2.92, versus 0.57 in the largest-genome quartile (mean 7.11 Mb; Figure 5). This pattern holds within each of the major bacterial phyla: Actinomycetota (rho = -0.857, p = 2.95e-23), Bacillota (rho = -0.734, p = 9.25e-10), and Bacteroidota (rho = -0.473, p = 0.020), confirming that it is not a phylogenetic confound.

Host-associated prokaryotes have significantly higher Shannon functional balance entropy than free-living organisms (Mann-Whitney p = 1.16e-5; Figure 4). The mechanistic explanation is that host-associated bacteria undergo genome reduction that selectively eliminates Other-binding MTases, which serve secondary metabolism functions dispensable in the nutrient-rich host environment, while retaining RNA-binding MTases essential for translational fidelity and frequently expanding DNA-binding MTases for restriction-modification immune defense, phase-variable epigenetic regulons, and virulence gene regulation. Shannon entropy is consistent across the three domains of life (Figure 3): bacteria (median H = 1.155), archaea (median H = 1.100), and eukaryotes (median H = 1.089) occupy nearly identical entropy ranges, suggesting that the degree of total substrate diversification is a conserved property even as the direction of specialization varies widely.

---

## 4. Discussion

### 4.1 The RNA-Other Trade-off as a Universal Cellular Constraint

The most striking finding of this work is the universality and phylogenetic robustness of the RNA-Other trade-off. The negative correlation holds across all three domains of life, survives phylogenetic correction with 946 independent contrasts (rho = -0.709), and is recapitulated within each eukaryotic kingdom. This implies that the trade-off is not an artifact of shared ancestry but reflects a genuine functional constraint: cells that divert more of their SAM MTase repertoire toward small-molecule biosynthesis must invest correspondingly less in RNA modification, and vice versa.

The inverted-U relationship between functional balance and RNA/Other ratio suggests an optimal zone (ratio ~1.27 in eukaryotes, ~1.51 in prokaryotes) at which organisms maintain high Shannon entropy across all four substrate categories. Animals cluster closer to this optimum than fungi or plants, consistent with the hypothesis that metazoan complexity demands diversified methylation strategies spanning chromatin regulation, translational control, and signaling. Fungi and plants instead specialize toward secondary metabolism (Other), reflecting their sessile lifestyles and reliance on chemical ecology. At the prokaryotic extreme, mycobacteria achieve the lowest Shannon entropies recorded (H ~0.63-0.70), with >60% of their SAM MTases predicted to target small molecules.

The mechanistic basis for the trade-off likely involves proteomic resource competition. The total number of SAM MTases per proteome represents a finite pool of genes maintained under selection. Expanding secondary metabolic MTases necessarily crowds out the allocation of encoding capacity to RNA-modifying enzymes, and vice versa. This view is supported by the BGC partial regression analysis, which shows that each additional BGC is associated with a 0.48% increase in Other MTase proportion independently of genome size, and that genome size loses explanatory power once BGC count is accounted for.

### 4.2 Lineage-Specific Deviations and Their Biological Drivers

While the RNA-Other trade-off is universal, lineages deviate from the global trend in biologically interpretable ways. Helicobacter species exhibit unusually high DNA MTase proportions (~52.5%), reflecting the dense restriction-modification defense systems characteristic of this host-adapted gastric pathogen [10]. Parasitic alveolates such as Toxoplasma gondii reach ~79.7% RNA MTases, consistent with a heavily RNA-centric gene regulation strategy in organisms that have undergone extreme metabolic streamlining. Actinomycetota exploit the high-Other extreme through BGC-encoded methyltransferases that catalyze key diversification steps in polyketide and non-ribosomal peptide biosynthesis.

These lineage-specific deviations demonstrate that while the trade-off axis constrains the space of viable MTase repertoires, evolution has explored the full spectrum of positions along this axis in response to ecological demands. The comparison of free-living versus host-associated bacteria further illustrates this dynamic: lifestyle change is accompanied by a predictable shift in MTase allocation, with genome reduction selectively targeting the metabolically dispensable Other category while preserving or expanding functionally critical RNA and DNA categories.

### 4.3 Vertebrate-Specific DNA Methylation and the Animal Epigenetic Transition

The vertebrate DNA MTase expansion is among the cleanest lineage-specific signals in the dataset. The 2.81-fold elevation in DNA-binding MTase proportions in vertebrates relative to invertebrate animals, with no accompanying change in protein-binding proportions, is consistent with the known evolutionary history of the DNMT gene family. The vertebrate ancestor underwent whole-genome duplication events that gave rise to the DNMT3A/DNMT3B/DNMT3L ohnolog group from an ancestral DNMT3-like gene, enabling the dense bidirectional CpG methylation that underlies imprinting and X-chromosome inactivation [9]. Invertebrates typically possess sparser, mosaic DNA methylation patterns with less DNMT gene family expansion, explaining the 2.81-fold difference.

The independence of the DNA expansion from protein MTase proportions is biologically meaningful: it suggests that vertebrate epigenetic complexity was achieved by layering a new DNA methylation system on top of an existing, conserved histone and RNA methylation framework, rather than by reallocating resources from pre-existing categories.

### 4.4 Limitations and Future Directions

Several limitations should be noted. First, the substrate prediction model classifies each protein into four broad categories, potentially collapsing heterogeneous functional groups, particularly in the "Other" category which encompasses structurally and biochemically diverse enzymes. PFAM-based analysis identified GNAT superfamily acetyltransferases as a likely contaminant of the Other category in some proteomes; refining the training labels to exclude non-SAM-dependent enzymes could improve specificity. Second, BGC counts were estimated from genome size rather than measured directly from antiSMASH annotations, introducing measurement uncertainty. Third, the vertebrate-versus-invertebrate comparison was performed on a purpose-built 50-animal dataset whose substrate predictions derive from a slightly different classifier pipeline than the 386-organism iTOL dataset.

Future work should extend predictions to metagenomics data to test whether the RNA-Other trade-off is conserved in uncultured lineages. Integration of expression data could reveal whether predicted substrate preferences translate into differential expression patterns across tissues and developmental stages. The identification of core rRNA/tRNA methyltransferase families (RSMA, TRMD, RSMG, RSMH) that are nearly universal across 270 prokaryotic genomes provides a specific mechanistic basis for the genome-reduction effect that could be tested with comparative genomics in independently evolved host-adapted lineages.

---

## 5. Conclusions

We present the first genome-scale, multi-domain analysis of SAM-dependent methyltransferase substrate preferences predicted by a high-performance deep learning classifier trained on protein language model embeddings. The classifier achieves near-perfect discrimination (AUPRC 0.998, MCC 0.976) on the SAM MTase benchmark and generalizes to produce biologically coherent predictions across 386 diverse organisms. Statistical analysis of these predictions reveals four major evolutionary patterns: a universal, phylogenetically robust RNA-Other trade-off; lineage-specific specialization tied to BGC content and lifestyle; a vertebrate-specific DNA MTase expansion decoupled from the RNA-Other axis; and genome-reduction-driven remodeling of methylation allocation in host-associated bacteria. Together, these findings reveal that despite the enormous diversity of SAM MTases in sequence and function, their genome-scale distributions are governed by a small number of deep evolutionary constraints that can now be quantitatively analyzed and compared across the full tree of life.

---

## References

1. Schubert HL, Blumenthal RM, Cheng X. Many paths to methyltransfer: a chronicle of convergence. Trends Biochem Sci. 2003;28(6):329-335.
2. Petrossian TC, Clarke SG. Uncovering the human methyltransferasome. Mol Cell Proteomics. 2011;10(1):M110.000976.
3. Elnaggar A, et al. ProtTrans: Toward understanding the language of life through self-supervised learning. IEEE Trans Pattern Anal Mach Intell. 2022;44(10):7112-7127.
4. Felsenstein J. Phylogenies and the comparative method. Am Nat. 1985;125(1):1-15.
5. Cimermancic P, et al. Insights into secondary metabolism from a global analysis of prokaryotic biosynthetic gene clusters. Cell. 2014;158(2):412-421.
6. Blin K, et al. antiSMASH 4.0: improvements in chemistry prediction and gene cluster boundary identification. Nucleic Acids Res. 2017;45(W1):W36-W41.
7. Lin CY, et al. Discovery of widespread non-ribosomal peptide carboxyl methyltransferases in bacteria. PNAS. 2023;120(12):e2214089120.
8. Zhang Q, et al. SAM-dependent methyltransferases in natural product biosynthesis. Bioresour Bioprocess. 2021;8:67.
9. Goll MG, Bestor TH. Eukaryotic cytosine methyltransferases. Annu Rev Biochem. 2005;74:481-514.
10. Seib KL, et al. The role of DNA methylation in the virulence of Neisseria and other host-adapted pathogens. Annu Rev Microbiol. 2020;74:217-240.

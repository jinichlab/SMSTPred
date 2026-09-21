# Interpretability roadmap: ProstT5 SAE + per-residue attribution

Scoping notes for opening the black box of the ProstT5 → MLP substrate classifier.
Two complementary tracks:

- **Track 1 — Sparse Autoencoder (SAE):** *global, feature-level.* "What internal
  ProstT5 features distinguish DNA/RNA/Protein/Other MTases, and which residues fire
  them?" Adapts InterProt (Adams et al. 2025, bioRxiv 2025.02.06.636901).
- **Track 2 — per-residue attribution:** *local, per-protein, model-faithful.* "For
  *this* protein, which residues drove the classifier's call?" Needs **no SAE** — runs
  on the existing classifier today.

Planned trigger: revisit when STREME leaves a set of sequences with **no discriminative
linear motif**. Attribution shows *where* the classifier looks on those; the SAE says
*what feature* it is.

---

## Track 1 — Training a ProstT5 SAE

### Key realization: the SAE core is model-agnostic
InterProt's `sae_model.py` (`SparseAutoencoder` + `loss_fn`) takes a
`(batch, length, d_model)` activation tensor and knows nothing about ESM. **Reusable
verbatim** with `d_model=1024`. Per-token LayerNorm, TopK activation, auxk dead-latent
revival, and weight/grad normalization (Gao et al. 2024 recipe) all transfer unchanged.

Everything ESM-specific is in two places only:
1. `esm_wrapper.ESM2Model.get_layer_activations()` — produces activations → replace with
   a ProstT5 extractor.
2. The `diff_cross_entropy` validation metric — splices reconstructed activations back
   through ESM's remaining layers + MLM head. ProstT5 is encoder–decoder with no clean
   MLM head → **drop it; use FVU / explained-variance + dead-latent count** (standard
   SAE-quality metrics).

So this is a front-end swap + validation simplification, not a from-scratch build.

### Recommended configuration
| Knob | Value | Why |
|---|---|---|
| pLM | `Rostlab/ProstT5` encoder, fp16 | our model |
| `d_model` | **1024** | ProstT5/ProtT5 hidden dim (ESM-650M was 1280) |
| layer | **last encoder hidden state** (primary) | the representation the MLP classifier actually consumes → latents explain *our* model |
| `d_hidden` | **4096** (4× expansion); option 8192 | InterProt's released, well-behaved size; our corpus is small-ish, start conservative |
| `k` (TopK) | **128** default; try 64 | their k-sweep: lower k → more family-specific latents |
| `auxk` / dead-thresh | 256 / 2000 | unchanged — revives dead latents |
| optim | AdamW, lr 2e-4, grad-clip 1.0, 1 epoch | InterProt recipe, unchanged |
| validation | FVU + dead-count | no MLM-head diff-CE for ProstT5 |

**Layer choice:** InterProt used layer 24/33 (~mid-late) because *family-specific*
features peak there. We have two valid targets, both cheap since extraction dominates
cost: (a) the layer the **classifier** reads → latents directly explain DNA/RNA/Protein/
Other decisions [**primary**]; (b) a mid encoder layer (~12/24) for richest family
features. **TODO: confirm the classifier uses last-layer mean-pooled ProstT5** (the
standard usage) so the SAE matches it.

### Corpus — MTase-only is the right call
InterProt trained on 1M generic UniRef50 sequences. We want the opposite: train on **our
MTase corpus** (the 106,608 single-label sequences ≈ **32M residue tokens**). Every
sequence shares the SAM-binding Rossmann fold, so the SAE must spend its dictionary on
**within-MTase variation** — exactly the substrate-discriminative features we want.
Cheaper and more on-target than a generic SAE; 32M tokens is enough for a 4096-latent
dictionary under Gao et al. scaling. (Optionally extend to the full 118,816.)

### Compute & storage
InterProt runs ESM *inside* the training step (fine for 650M). ProstT5 is **~3B params** —
too heavy to recompute per step, and we want to sweep k / d_hidden. So **decouple**:

- **Phase A — extract once (~30–60 min GPU):** run ProstT5 encoder over the 106k seqs,
  `output_hidden_states=True`, strip special tokens, write **per-residue fp16 shards** to
  scratch. ⚠️ **Storage ≈ 64 GB** (32M × 1024 × 2 bytes) — the one real resource cost.
  Reuses the existing ProstT5 embedder, minus mean-pooling, saving full L×1024.
- **Phase B — train SAE from cache (minutes–1 hr/config):** a single 1024×4096 linear
  enc/dec; trivially cheap. Sweep k∈{64,128}, d_hidden∈{4096,8192} from one cache for free.

End-to-end: a few GPU-hours, dominated by extraction + disk, not SAE training.

### Code work (~1–2 focused days)
1. `prostt5_extract.py` — adapt the ProstT5 embedder: prefix `<AA2fold>`,
   `output_hidden_states`, save residue shards + index. (~80 lines, mostly existing)
2. `sae_dataset.py` — stream residue vectors from shards into token-batches (replaces
   InterProt's `SequenceDataModule`, ~60 lines)
3. `train_prostt5_sae.py` — InterProt's `training.py`/`sae_module.py` with ESM removed,
   MSE+auxk kept, FVU/dead-count logging (~100 lines)
4. `sae_model.py` — **copy verbatim**, `d_model=1024`
5. Downstream (the science, separate step): per-latent top-activating sequences,
   substrate-specificity scoring (does a latent fire only on DNA-MTases?), linear-probe-
   on-latents to find substrate-discriminative features (analogs of their
   `make_viz_files`/`oned_probe`), plus our own residue-highlight plots (interprot.com is
   ESM-only). Cross-check vs STREME: do top DNA latents land on the same DPPY catalytic
   residues?

### Reference snippets (from the InterProt clone)
- `SparseAutoencoder.__init__(d_model, d_hidden, k=128, auxk=256, batch_size=256,
  dead_steps_threshold=2000)`
- Loss: `mse + (1/32)*auxk_mse`; per-token LN normalize→reconstruct→rescale; w_dec
  unit-norm columns; auxk revives latents dead > 2000/batch_size steps.
- Their defaults: layer 24, d_hidden 16384, k 128, lr 2e-4, batch 48 seqs, 1 epoch.
- Released usable checkpoint: `esm2_plm1280_l24_sae4096` (ESM2-650M, layer 24, 4096, TopK)
  on HF `liambai/InterProt-ESM2-SAEs`.

---

## Track 2 — Per-residue attribution (no SAE needed)

Answers "which residues contribute to the feature the model attends to" directly on the
existing ProstT5 → MLP classifier.

### The mean-pooling subtlety (important)
The classifier consumes a **mean-pooled** embedding: `z = (1/L) Σ e_i`, `logit = f(z)`.
Because pooling is linear, `∂logit/∂e_i = (1/L)·∂f/∂z` is **identical for every residue** —
so *naive gradient on the pooled vector gives uniform, useless attribution*. The fix is
methods that use each residue's embedding `e_i` (which *do* differ):

### Methods, best first
1. **Integrated Gradients / Gradient×Input on residue embeddings** *(recommended first cut)*
   - Per-residue contribution `c_i = (e_i − baseline_i) · ∫ ∂f/∂z dα` (IG), or the cheap
     linear version `c_i = (1/L)·(∂f/∂z)·e_i` (input×gradient). Residues differ via `e_i`,
     so this *does* localize. Σ c_i ≈ logit → attributions are conservation-checkable.
   - Needs only the existing classifier + per-residue embeddings (transient, no 64 GB
     cache). Cheap, batchable, genome-scale.
2. **Occlusion / in-silico mutagenesis** *(intuitive validation on showcase proteins)*
   - Mask or mutate each residue (or sliding window), re-embed through ProstT5, re-predict,
     measure Δlogit for the class. Model-agnostic and very interpretable ("remove this
     position → DNA call collapses"), but O(L) ProstT5 passes/protein → only for a handful
     of showcase enzymes, not genome-scale.
3. **ProstT5 self-attention maps** *(weakest evidence, complementary)*
   - Shows which residues the *pLM* attends to internally, but attention ≠ task importance
     ("attention is not explanation"), and it reflects the MLM objective, not our
     substrate head. A "what the model wires together" view, not "what drives the call."

### How it pairs with everything else
- **STREME-negative sequences:** run IG → paint per-residue importance on sequence +
  AF/ProstT5 structure. If the signal is a *dispersed/structural* region rather than a
  conserved linear k-mer, that explains why STREME found nothing — and is itself a finding.
- **vs SAE:** attribution = local, per-protein, model-faithful ("which residues *here*");
  SAE = global, feature-level, interpretable ("what feature, across the dataset"). Use
  attribution to find *where*, SAE to name *what*.

### Effort
IG/Grad×Input: ~half a day on the existing classifier (reuse the per-residue extraction
from Phase A, but only transiently). Occlusion: another ~half day for showcase proteins.
No training, no new model.

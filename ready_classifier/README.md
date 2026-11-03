# SMSTPred ready classifier

Give it protein sequences, get back the **substrate type** of each SAM-dependent methyltransferase: **Metabolite** (small molecule), **Protein**, **RNA** or **DNA**
(any combination; a protein can be called for more than one). It is the trained model of the paper, ready to run on your own sequences: no training, no data download except the
protein language model it uses to read the sequences.

It always uses the **newest trained model**: every time the pipeline retrains on a new UniProt release it replaces `model/<release>/` and `model/LATEST.json`, so `git pull` is all
that is needed to update. `python get_model.py` shows which model you have.

## 1. Install (once)

```bash
conda env create -f environment.yml
conda activate smstpred-classify
```
A GPU makes it fast (about 10 sequences per second); without one it still works, slowly (use `--device cpu`). The first run downloads the protein language model
`Rostlab/prot_t5_xl_half_uniref50-enc` (about 1.2 GB) from Hugging Face into your cache; later runs reuse it.

## 2. Prepare your input

A plain **FASTA** file of protein sequences:

```
>my_protein_1 any description you like
MIEIKDKQLTGLRFIDLFAGLGGFRLALESCGAECVYSNEWDKYAQEVYEMNFGEKPEG...
>my_protein_2
MSRLPPA...
```
* One header line per protein starting with `>`; the **first word** after `>` is the identifier that comes back in the results (repeated identifiers get `_2`, `_3` appended).
* The sequence can be on one line or wrapped over several. Capital or lower case, `-` and `*` are ignored.
* Standard amino acids; `U`, `Z`, `O`, `B`, `X` are accepted. Anything else makes the protein "not scored" (it stays in the output with the reason).
* **Length 30 to 2,000 amino acids.** Shorter or longer sequences are listed as "not scored", never dropped silently.
* The model was trained on proteins that carry the **SAM-dependent methyltransferase fold** (InterPro IPR029063). Run it on such proteins. On other proteins it will still give an answer,
  but the answer means nothing: the `in_domain` column below is the warning.

## 3. Run

```bash
python classify.py my_proteins.fasta                  # writes my_proteins_predictions.csv
python classify.py my_proteins.fasta -o results.csv --device cpu
bash run_example.sh                                   # the six example sequences in example/
```

## 4. Read the output

One row per input sequence, in input order (`example/example_predictions.csv` is a real output):

| Column | Meaning |
|---|---|
| `Entry` | your identifier |
| `length` | amino acids |
| `status` | `scored`, or why the sequence was not scored |
| `score_Metabolite`, `score_Protein`, `score_RNA`, `score_DNA` | the model's probability (0 to 1) that the enzyme acts on that substrate type: the mean of six ensemble members |
| `call_Metabolite` ... `call_DNA` | 1 if the score is at least 0.5, else 0. The four labels are decided independently, so a protein can have none, one or several calls |
| `predicted_substrate` | the calls in words, e.g. `RNA`, `Metabolite+RNA`, or `none` |
| `disagreement` | how much the six members disagree (average spread of their scores). Near 0 = they agree; above about 0.05 the call is shaky |
| `ad_maxsim` | the highest cosine similarity of your protein to any training protein |
| `in_domain` | 1 if `ad_maxsim` is at least the model's applicability-domain cutoff (about 0.81, the 1st percentile of training proteins' similarity to the rest of the training set); **0 means your protein is unlike anything the model was trained on: treat its calls with caution** |

How to use them: trust a call most when `in_domain` is 1 and `disagreement` is small. A high score means the model resembles the enzymes it was trained on that act on that substrate type; it is **not**
a prediction of the exact substrate molecule, and it is not experimental evidence. A protein with `none` is one the model could not place in any of the four types.

## 5. How good is it

On held-out proteins from protein families the model had never seen (cluster-disjoint cross-validation, honest out-of-fold predictions) the current model reaches the MCC listed in
`model/LATEST.json` (`out_of_fold_metrics`, per label as well). Performance on your proteins will be lower the further they are from the training families (see `in_domain`).

## 6. Where the model comes from

`model/LATEST.json` names the release (UniProt release, training set size, labels, checksums of every file). It is written by `workflow/publish_classifier.py`, the last step of the training stage of the
pipeline in this repository (`workflow/run_refresh.sh`), so the model here is always the one of the newest completed run. If a model file is damaged or missing, `classify.py` says so and (once the data deposit is
public and its address is in `LATEST.json`) fetches it again. Use your own trained ensemble with `--model-dir <folder>` or `SMSTPRED_MODEL_DIR`.

Files: `classify.py` (the command), `get_model.py` (finds and checks the newest model), `model/` (the ensemble), `example/` (six sequences and the real output for them), `environment.yml`, `run_example.sh`.

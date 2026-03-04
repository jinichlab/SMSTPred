import os
import torch
import pandas as pd
import h5py
from pathlib import Path
from tqdm import tqdm
from transformers import T5Tokenizer, T5EncoderModel
import numpy as np

# ========= CONFIG =========
DATA_DIR = Path("../substrate_type_predictor_testing/reference_files/reference_proteomes")
EMBED_DIR = Path("../substrate_type_predictor_testing/reference_files/embeddings")
EMBED_DIR.mkdir(exist_ok=True)
OUT_SUFFIX = "_embeddings.h5"
IPR_TARGET = "IPR029063"
BATCH_SIZE = 2  # Adjust based on GPU memory

# ========= Load ProtT5 =========
print("🧠 Loading ProtT5 model and tokenizer...")
tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False)
model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
model = model.eval()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = model.to(device)
print(f"✅ Model loaded on: {device}")

# ========= Helper =========
def preprocess_sequence(sequence: str) -> str:
    return ' '.join(list(sequence.replace(" ", "").replace("\n", "")))

def batch_embed(sequences: list[str]) -> list[torch.Tensor]:
    tokenized = tokenizer(sequences, return_tensors='pt', padding=True, truncation=True)
    input_ids = tokenized.input_ids.to(device)
    attention_mask = tokenized.attention_mask.to(device)

    with torch.no_grad():
        outputs = model(input_ids=input_ids, attention_mask=attention_mask).last_hidden_state
        mask = attention_mask.unsqueeze(-1).expand(outputs.size())
        summed = torch.sum(outputs * mask, dim=1)
        counts = torch.clamp(mask.sum(dim=1), min=1e-9)
        mean_embs = summed / counts
    return mean_embs.cpu()

# ========= Run per proteome =========
for csv_file in tqdm(sorted(DATA_DIR.glob("*.csv"))):
    df = pd.read_csv(csv_file)
    if "InterPro" not in df.columns or "Sequence" not in df.columns:
        print(f"⚠️ Skipping {csv_file.name} (missing columns)")
        continue

    total_proteins = len(df)
    df_sams = df[df['InterPro'].fillna('').str.contains(IPR_TARGET)]
    total_sams = len(df_sams)
    df_sams = df_sams[df_sams['Length'] <= 1000]
    print(f"🔬 Filtered SAMs ≤ 1000 aa: {len(df_sams)} remaining")

    if df_sams.empty:
        print(f"⚠️ No SAM MTases in {csv_file.name}")
        continue

    out_h5_path = EMBED_DIR / f"{csv_file.stem}{OUT_SUFFIX}"
    print(f"\n📂 Processing {csv_file.name}")
    print(f"   → Total proteins: {total_proteins}")
    print(f"   → SAM MTases found: {total_sams}")

    sequences = []
    ids = []
    for _, row in df_sams.iterrows():
        seq_id = row["Entry"]
        sequence = row["Sequence"]
        if pd.isna(sequence) or not isinstance(sequence, str) or len(sequence.strip()) < 10:
            continue
        sequences.append(preprocess_sequence(sequence))
        ids.append(seq_id)

    print(f"   → Sequences valid for embedding: {len(sequences)}")

    embeddings_written = 0
    with h5py.File(out_h5_path, "w") as h5file:
        for i in range(0, len(sequences), BATCH_SIZE):
            batch_seqs = sequences[i:i+BATCH_SIZE]
            batch_ids = ids[i:i+BATCH_SIZE]
            try:
                batch_embs = batch_embed(batch_seqs)
                for j, emb in enumerate(batch_embs):
                    # ---------------- Robustness checks ----------------
                    if emb is None:
                        print(f"⚠️ {batch_ids[j]}: embedding is None, skipping")
                        continue
                    if not isinstance(emb, torch.Tensor):
                        print(f"⚠️ {batch_ids[j]}: embedding not a torch.Tensor, converting")
                        emb = torch.tensor(emb, dtype=torch.float32)
                    try:
                        np_emb = emb.numpy()
                        if np_emb.ndim != 1:
                            print(f"❌ {batch_ids[j]}: non-1D embedding: shape={np_emb.shape}")
                            continue
                        if np_emb.shape[0] != 1024:
                            print(f"❌ {batch_ids[j]}: unexpected embedding size: {np_emb.shape}")
                            continue
                        if np.isnan(np_emb).any() or np.isinf(np_emb).any():
                            print(f"❌ {batch_ids[j]}: contains NaN or Inf")
                            continue

                        h5file.create_dataset(batch_ids[j], data=np_emb.astype(np.float32))
                        embeddings_written += 1  # ← Add this

                    except Exception as e:
                        print(f"❌ Error saving {batch_ids[j]}: {e}")

            except Exception as batch_error:
                print(f"❌ Error embedding batch starting with {batch_ids[0]}: {batch_error}")

    print(f"✅ {out_h5_path.name} written with {embeddings_written} embeddings")
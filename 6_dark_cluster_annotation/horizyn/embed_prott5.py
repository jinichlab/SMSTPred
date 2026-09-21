"""ProtT5-XL-half embeddings, matching the checkpoint Horizyn was trained on
(Rostlab/prot_t5_xl_half_uniref50-enc; see data/README.md line 93). Standard
Rostlab/ProtTrans mean-pooling recipe: replace rare amino acids with X,
space-separate residues (T5 tokenizer requirement), run the encoder, mean-pool
per-residue states over the true sequence length (excludes padding and the
trailing </s> special token).

NOT part of the horizyn package itself -- horizyn only *consumes* protein
embeddings (via `target_encoder`), it doesn't generate them. This is the
missing "generate prots_t5.h5 for your own proteins" step the repo's data
README leaves to the user.
"""
import re
import torch
from transformers import T5EncoderModel, T5Tokenizer

CHECKPOINT = "Rostlab/prot_t5_xl_half_uniref50-enc"


def load_prott5(device):
    tokenizer = T5Tokenizer.from_pretrained(CHECKPOINT, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(CHECKPOINT)
    model = model.to(device)
    model = model.half() if device.type == "cuda" else model.float()  # half ckpt needs fp32 cast on CPU
    model.eval()
    return tokenizer, model


def embed_sequences(sequences, tokenizer, model, device, batch_size=16):
    """Returns a (N, 1024) float32 tensor, one row per input sequence, in the
    SAME order as the input `sequences` list. Internally batches in length-sorted
    order (much less padding waste when lengths vary widely, e.g. 70-900+ aa)
    then un-permutes back to the caller's original order."""
    clean = [re.sub(r"[UZOB]", "X", s.strip().upper()) for s in sequences]
    order = sorted(range(len(clean)), key=lambda i: len(clean[i]))  # short -> long
    spaced = [" ".join(list(clean[i])) for i in order]
    out_rows = [None] * len(sequences)

    for i in range(0, len(spaced), batch_size):
        idx_batch = order[i:i + batch_size]
        batch = spaced[i:i + batch_size]
        true_lens = [len(clean[k]) for k in idx_batch]
        enc = tokenizer.batch_encode_plus(
            batch, add_special_tokens=True, padding="longest", return_tensors="pt"
        )
        input_ids = enc["input_ids"].to(device)
        attention_mask = enc["attention_mask"].to(device)
        with torch.no_grad():
            hidden = model(input_ids=input_ids, attention_mask=attention_mask).last_hidden_state
        for j, (orig_idx, L) in enumerate(zip(idx_batch, true_lens)):
            out_rows[orig_idx] = hidden[j, :L].mean(dim=0).float().cpu()
        print(f"  embedded {min(i + batch_size, len(spaced))}/{len(spaced)}", flush=True)

    return torch.stack(out_rows)

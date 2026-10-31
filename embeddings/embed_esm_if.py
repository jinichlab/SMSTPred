#!/usr/bin/env python3
"""
ESM-IF1 structure embeddings (512-D, mean over residues of the encoder output, chain A) for AlphaFold DB models.

  python embed_esm_if.py --csv labels_SAM_full.csv --out embeddings_SAM_ESM-IF.h5

Pipelined download -> ESM-IF embed -> delete for AFDB structures (model version v6).

CPU threads stream PDBs from AlphaFold DB into a bounded queue; the
GPU main thread pulls each PDB, runs ESM-IF, writes a 512-D mean-
pooled embedding into an HDF5 file keyed by UniProt Entry, then
deletes the PDB. Resumable: entries already in the H5 are skipped.

Pass --limit N to embed only the first N Entries (test mode).
"""

import argparse
import queue
import tempfile
import threading
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import requests
import torch

import esm
import esm.inverse_folding
from esm.inverse_folding.util import CoordBatchConverter

META_CSV    = 'uniprot_metadata.csv'      # only for --reviewed-only
AFDB_URL    = 'https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.pdb'

N_DOWNLOAD_WORKERS = 6
QUEUE_MAXSIZE      = 16
DOWNLOAD_TIMEOUT   = 60


def download(entry, tmpdir):
    """Fetch one AFDB PDB; return Path or None on miss / error."""
    url = AFDB_URL.format(acc=entry)
    try:
        r = requests.get(url, timeout=DOWNLOAD_TIMEOUT)
        if r.status_code == 404:
            return None
        r.raise_for_status()
    except requests.RequestException as e:
        print(f'  ! download {entry}: {e}', flush=True)
        return None
    p = Path(tmpdir) / f'{entry}.pdb'
    p.write_bytes(r.content)
    return p


def producer(entries, work_q, tmpdir):
    for entry in entries:
        pdb = download(entry, tmpdir)
        work_q.put((entry, pdb))


def encode_structure(model, alphabet, coords, device):
    """Workaround for fair-esm 2.0 bug where CoordBatchConverter doesn't
    move every returned tensor to the requested device. We move them
    all explicitly before the encoder forward pass."""
    batch_converter = CoordBatchConverter(alphabet)
    coords_b, confidence, _, _, padding_mask = batch_converter(
        [(coords, None, None)], device=device)
    coords_b    = coords_b.to(device)
    confidence  = confidence.to(device)
    padding_mask = padding_mask.to(device)
    with torch.no_grad():
        out = model.encoder.forward(coords_b, padding_mask, confidence,
                                    return_all_hiddens=False)
    # out['encoder_out'][0]: shape (L+2, batch=1, 512); strip <cath>/<af2>
    return out['encoder_out'][0][1:-1, 0, :]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int, default=None,
                    help='embed only the first N Entries (test mode)')
    ap.add_argument('--reviewed-only', action='store_true',
                    help='restrict to Swiss-Prot reviewed entries '
                         '(uses uniprot_metadata.csv; guaranteed in AFDB)')
    ap.add_argument('--device', choices=['cuda', 'cpu'], default='cuda',
                    help='force device. cpu is ~30x slower but bypasses a '
                         'fair-esm 2.0 + GPU device-mismatch bug')
    ap.add_argument('--csv', required=True, help='CSV with a column of UniProt accessions (default column: Entry)')
    ap.add_argument('--out', required=True, help='output H5 path')
    ap.add_argument('--entry-col', default='Entry', help='accession column name')
    ap.add_argument('--hitrate', default=None,
                    help='optional AFDB hit-rate CSV (columns Entry, status) used to pre-filter known 404s')
    args = ap.parse_args()
    out_h5 = args.out

    entries = (pd.read_csv(args.csv, usecols=[args.entry_col])[args.entry_col]
               .dropna().astype(str).drop_duplicates().tolist())

    # Pre-filter known AFDB misses to avoid wasted 404 round-trips.
    if args.hitrate and Path(args.hitrate).exists():
        hr = pd.read_csv(args.hitrate)
        in_afdb = set(hr.loc[hr['status'] == 200, 'Entry'])
        before = len(entries)
        entries = [e for e in entries if e in in_afdb]
        print(f'AFDB-filtered: {len(entries):,} / {before:,} '
              f'(dropped {before - len(entries):,} known 404s)', flush=True)

    if args.reviewed_only:
        m = pd.read_csv(META_CSV, dtype=str, low_memory=False)
        key = 'From' if 'From' in m.columns else 'Entry'
        rev = set(m.loc[m['Reviewed'] == 'reviewed', key])
        entries = [e for e in entries if e in rev]
        print(f'filtered to Swiss-Prot reviewed: {len(entries):,}', flush=True)
    if args.limit:
        entries = entries[:args.limit]
    print(f'queued: {len(entries):,} Entries', flush=True)

    # Resume: skip entries already in the H5.
    if Path(out_h5).exists():
        with h5py.File(out_h5, 'r') as f:
            done = set(f.keys())
        entries = [e for e in entries if e not in done]
        print(f'  {len(done):,} already in {out_h5}; '
              f'{len(entries):,} remaining', flush=True)
    if not entries:
        print('nothing to do', flush=True)
        return

    print('loading ESM-IF1...', flush=True)
    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.eval()
    use_cuda = args.device == 'cuda' and torch.cuda.is_available()
    if use_cuda:
        model = model.cuda()
        print(f'  cuda: {torch.cuda.get_device_name(0)}', flush=True)
    else:
        print(f'  device: cpu (requested={args.device}, '
              f'cuda_available={torch.cuda.is_available()})', flush=True)

    n_done = n_miss = n_fail = 0
    t0 = time.time()
    with tempfile.TemporaryDirectory() as tmpdir:
        work_q = queue.Queue(maxsize=QUEUE_MAXSIZE)
        chunks = [entries[i::N_DOWNLOAD_WORKERS]
                  for i in range(N_DOWNLOAD_WORKERS)]
        threads = []
        for chunk in chunks:
            t = threading.Thread(target=producer,
                                 args=(chunk, work_q, tmpdir), daemon=True)
            t.start()
            threads.append(t)

        with h5py.File(out_h5, 'a') as h5:
            remaining = len(entries)
            while remaining > 0:
                entry, pdb = work_q.get()
                remaining -= 1
                seen = len(entries) - remaining
                if pdb is None:
                    n_miss += 1
                    print(f'  [{seen}/{len(entries)}] {entry}: no AFDB',
                          flush=True)
                    continue
                try:
                    coords, _ = esm.inverse_folding.util.load_coords(
                        str(pdb), 'A')
                    device = next(model.parameters()).device
                    rep = encode_structure(model, alphabet, coords, device)
                    emb = rep.mean(dim=0).detach().cpu().numpy().astype(
                        np.float32)
                    h5.create_dataset(entry, data=emb, compression='gzip')
                    n_done += 1
                    rate = n_done / max(time.time() - t0, 1e-6)
                    print(f'  [{seen}/{len(entries)}] {entry}: '
                          f'L={len(coords)} dim={emb.shape[0]}  '
                          f'({rate:.2f} prot/s)', flush=True)
                except Exception as e:
                    n_fail += 1
                    print(f'  ! embed {entry}: {e}', flush=True)
                finally:
                    try:
                        pdb.unlink()
                    except OSError:
                        pass

        for t in threads:
            t.join()

    dt = time.time() - t0
    print(f'\nfinished in {dt:.1f}s', flush=True)
    print(f'  embedded: {n_done}', flush=True)
    print(f'  no AFDB:  {n_miss}', flush=True)
    print(f'  failed:   {n_fail}', flush=True)
    with h5py.File(out_h5, 'r') as f:
        print(f'  {out_h5}: {len(f.keys())} embeddings total', flush=True)


if __name__ == '__main__':
    main()

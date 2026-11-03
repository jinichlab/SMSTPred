#!/usr/bin/env python3
"""
Small UniProt REST client for the refresh pull (Stage 0): polite rate limit, retries with backoff, release tracking,
cursor-paged TSV download for the big queries.

Every response carries the UniProt release it was served from (X-UniProt-Release); the client counts them so the pull can
refuse to mix two releases. Nothing here is specific to the SAM project; pull_uniprot.py holds the queries.
"""
import gzip
import hashlib
import os
import threading
import time
from collections import Counter
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE = "https://rest.uniprot.org"
USER_AGENT = "SMSTPred-uniprot-refresh/1.0 (+https://github.com/jinichlab/SMSTPred)"


class Client:
    def __init__(self, min_interval=0.2, timeout=300):
        self.s = requests.Session()
        retry = Retry(total=6, connect=6, read=6, backoff_factor=2.0, status_forcelist=(429, 500, 502, 503, 504),
                      allowed_methods=frozenset(["GET"]), respect_retry_after_header=True)
        self.s.mount("https://", HTTPAdapter(max_retries=retry, pool_connections=8, pool_maxsize=8))
        self.s.headers.update({"User-Agent": USER_AGENT})
        self.timeout = timeout
        self.min_interval = min_interval          # seconds between the starts of two requests, across all threads
        self._lock = threading.Lock()
        self._last = 0.0
        self.releases = Counter()                 # (release, release date) -> number of responses

    def _throttle(self):
        with self._lock:
            wait = self._last + self.min_interval - time.monotonic()
            if wait > 0:
                time.sleep(wait)
            self._last = time.monotonic()

    def get(self, url, params=None):
        self._throttle()
        r = self.s.get(url, params=params, timeout=self.timeout)
        r.raise_for_status()
        rel = r.headers.get("X-UniProt-Release")
        if rel:
            with self._lock:
                self.releases[(rel, r.headers.get("X-UniProt-Release-Date", ""))] += 1
        return r

    # ------------------------------------------------------------------ queries
    def count(self, query):
        """(number of entries, release, release date) for a UniProtKB query."""
        r = self.get(f"{BASE}/uniprotkb/search", {"query": query, "size": 1, "fields": "accession"})
        return int(r.headers["X-Total-Results"]), r.headers.get("X-UniProt-Release"), r.headers.get("X-UniProt-Release-Date")

    def proteome_meta(self, upid):
        """A few fields of the proteome record (type, protein count, last modified, redundant-to), or {} if it does not exist."""
        try:
            d = self.get(f"{BASE}/proteomes/{upid}").json()
        except requests.HTTPError as e:
            return {"proteome_http_status": e.response.status_code}
        red = d.get("redundantTo")
        return {"proteome_type": d.get("proteomeType"), "proteome_protein_count": d.get("proteinCount"),
                "proteome_modified": d.get("modified"), "proteome_gene_count": d.get("geneCount"),
                "proteome_redundant_to": (red.get("id") if isinstance(red, dict) else red)}

    def stream_tsv(self, query, fields):
        """Whole result of a (small) query as TSV text through the stream endpoint, and the release that served it."""
        r = self.get(f"{BASE}/uniprotkb/stream", {"query": query, "fields": fields, "format": "tsv"})
        return r.text, r.headers.get("X-UniProt-Release")

    def stream_tsv_gz(self, query, fields, out_gz, retries=3, log=print):
        """Download a big query in one streamed request (the endpoint UniProt recommends for large downloads) into a gzipped TSV,
        atomically. The row count must equal the count the API announces for the query; otherwise (or after a dropped connection)
        the whole download is repeated, up to `retries` times. Same return value as paged_tsv."""
        out_gz = Path(out_gz)
        tmp = out_gz.with_suffix(out_gz.suffix + ".part")
        expected, _, _ = self.count(query)
        last = None
        for attempt in range(1, retries + 1):
            t0 = time.time()
            try:
                self._throttle()
                with self.s.get(f"{BASE}/uniprotkb/stream", params={"query": query, "fields": fields, "format": "tsv"},
                                stream=True, timeout=self.timeout) as r:
                    r.raise_for_status()
                    release, rdate = r.headers.get("X-UniProt-Release"), r.headers.get("X-UniProt-Release-Date")
                    newlines, tail, mb = 0, b"", 0
                    with gzip.open(tmp, "wb") as fh:
                        for chunk in r.iter_content(1 << 20):
                            fh.write(chunk)
                            newlines += chunk.count(b"\n")
                            tail = chunk[-1:] or tail
                            mb += 1
                            if mb % 100 == 0:
                                log(f"    {mb} MB received")
                rows = newlines - 1 + (0 if tail == b"\n" else 1)          # minus the header line
                with self._lock:
                    self.releases[(release, rdate)] += 1
                if rows != expected:
                    raise RuntimeError(f"expected {expected} rows but received {rows}")
                os.replace(tmp, out_gz)
                return {"rows": rows, "expected": expected, "release": release, "release_date": rdate,
                        "seconds": round(time.time() - t0, 1)}
            except Exception as e:                      # dropped connection, truncated stream, release switch ...
                last = e
                log(f"    attempt {attempt}/{retries} failed: {type(e).__name__}: {str(e)[:150]}")
                time.sleep(10 * attempt)
        raise RuntimeError(f"streamed download of {query!r} failed {retries} times: {last}")

    def paged_tsv(self, query, fields, out_gz, page_size=500, limit=None, log=print):
        """Download a big query page by page (cursor pagination) into a gzipped TSV, atomically.

        Returns dict(rows, expected, release, release_date, seconds). Raises if the row count differs from the count the
        API announced for the query (a release switch or a dropped page would show up here)."""
        out_gz = Path(out_gz)
        tmp = out_gz.with_suffix(out_gz.suffix + ".part")
        t0 = time.time()
        params = {"query": query, "fields": fields, "format": "tsv", "size": page_size}
        url, rows, expected, release, rdate, page = f"{BASE}/uniprotkb/search", 0, None, None, None, 0
        with gzip.open(tmp, "wt", newline="\n") as fh:
            while url:
                r = self.get(url, params if page == 0 else None)
                if page == 0:
                    expected = int(r.headers["X-Total-Results"])
                    release, rdate = r.headers.get("X-UniProt-Release"), r.headers.get("X-UniProt-Release-Date")
                lines = r.text.split("\n")
                if page == 0:
                    fh.write(lines[0] + "\n")
                body = [ln for ln in lines[1:] if ln]
                fh.write("\n".join(body) + ("\n" if body else ""))
                rows += len(body)
                page += 1
                if page % 20 == 0:
                    log(f"    page {page}: {rows}/{expected} rows")
                url = r.links.get("next", {}).get("url")
                if limit and rows >= limit:
                    break
        if not limit and rows != expected:
            raise RuntimeError(f"expected {expected} rows for {query!r} but received {rows}")
        os.replace(tmp, out_gz)
        return {"rows": rows, "expected": expected, "release": release, "release_date": rdate, "seconds": round(time.time() - t0, 1)}


def sha256(path, block=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(block), b""):
            h.update(b)
    return h.hexdigest()

#!/usr/bin/env python3
"""Print the release folder the chain works on: $RELEASE_DIR if set, else the newest release in $UNIPROT_PULL_DIR ($RELEASE picks one by its name, e.g. 2026_03).
Used by the stage scripts that run after the pull has made the release folder (so the chain need not know its name when it is submitted)."""
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "1_clustering"))
from mmseqs_clusters import newest_release  # noqa: E402

if os.environ.get("RELEASE_DIR"):
    print(os.environ["RELEASE_DIR"])
else:
    print(newest_release(Path(os.environ.get("UNIPROT_PULL_DIR", REPO / "data" / "uniprot_pull")), os.environ.get("RELEASE")))

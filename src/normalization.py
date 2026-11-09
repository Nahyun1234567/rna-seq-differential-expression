"""Loading, sample-metadata matching, filtering and normalization helpers."""

from pathlib import Path

import numpy as np
import pandas as pd

CONDITION_CODES = {"C": "control", "W": "wt", "L": "variant"}


def load_counts(path):
    """Read the raw count matrix (genes x samples)."""
    return pd.read_csv(path, sep="\t", index_col=0)


def load_metadata(path, counts=None):
    """Read sample metadata and, if ``counts`` is given, align it to its columns.

    Raises if any count column has no metadata row, so a mislabelled sample
    can never slip through silently.
    """
    meta = pd.read_csv(path).set_index("sample")
    if counts is not None:
        missing = set(counts.columns) - set(meta.index)
        if missing:
            raise ValueError(f"Samples missing from metadata: {sorted(missing)}")
        meta = meta.loc[counts.columns]
        meta.index.name = "sample"
    return meta


def condition_from_sample_name(name):
    """Infer condition from the sample code in the column name, e.g. '...-HepG2-W3_S12' -> 'wt'."""
    code = name.split("-")[-1][0]
    return CONDITION_CODES[code]


def cpm(counts):
    """Counts per million."""
    return counts / counts.sum(axis=0) * 1e6


def log2_cpm(counts, pseudocount=1):
    return np.log2(cpm(counts) + pseudocount)


def filter_low_counts(counts, min_cpm=1.0, min_samples=3):
    """Keep genes with CPM >= ``min_cpm`` in at least ``min_samples`` samples.

    ``min_samples`` is usually the size of the smallest group (3 here), so a
    gene expressed in only one condition is still kept.
    """
    keep = (cpm(counts) >= min_cpm).sum(axis=1) >= min_samples
    return counts.loc[keep]


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)
    return Path(path)

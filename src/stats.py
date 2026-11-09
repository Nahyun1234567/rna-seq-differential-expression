"""Differential expression: DESeq2 (pydeseq2) and a simple Welch t-test baseline."""

import numpy as np
import pandas as pd
from scipy import stats

# (numerator, denominator): log2FC > 0 means higher in the numerator group
CONTRASTS = [
    ("wt", "control"),
    ("variant", "control"),
    ("variant", "wt"),
]


def bh_fdr(pvalues):
    """Benjamini-Hochberg adjusted p-values (NaNs are kept as NaN)."""
    p = np.asarray(pvalues, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


def ttest_de(log_expr, meta, numerator, denominator, condition_col="condition"):
    """Welch t-test per gene on log-scale expression (e.g. log2 CPM)."""
    a = log_expr.loc[:, meta[condition_col] == numerator]
    b = log_expr.loc[:, meta[condition_col] == denominator]
    t, p = stats.ttest_ind(a, b, axis=1, equal_var=False)
    res = pd.DataFrame(
        {
            "baseMean_log2": log_expr.mean(axis=1),
            "log2FoldChange": a.mean(axis=1) - b.mean(axis=1),
            "stat": t,
            "pvalue": p,
        },
        index=log_expr.index,
    )
    res["padj"] = bh_fdr(res["pvalue"])
    return res.sort_values("padj")


def run_deseq2(counts, meta, condition_col="condition", ref_level="control", n_cpus=1):
    """Fit one DESeq2 model on all samples. Returns the fitted DeseqDataSet."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference

    dds = DeseqDataSet(
        counts=counts.T.astype(int),  # pydeseq2 wants samples x genes
        metadata=meta[[condition_col]],
        design=f"~{condition_col}",
        refit_cooks=True,
        inference=DefaultInference(n_cpus=n_cpus),
        quiet=True,
    )
    dds.deseq2()
    return dds


def deseq2_contrast(dds, numerator, denominator, condition_col="condition", n_cpus=1):
    """Wald test for one contrast on a fitted DeseqDataSet."""
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    ds = DeseqStats(
        dds,
        contrast=[condition_col, numerator, denominator],
        inference=DefaultInference(n_cpus=n_cpus),
        quiet=True,
    )
    ds.summary()
    return ds.results_df.sort_values("padj")


def call_degs(res, padj=0.05, lfc=1.0):
    """Label each gene as up / down / ns."""
    sig = (res["padj"] < padj) & (res["log2FoldChange"].abs() > lfc)
    direction = np.where(res["log2FoldChange"] > 0, "up", "down")
    return pd.Series(np.where(sig, direction, "ns"), index=res.index, name="regulation")


def deg_summary(results, padj=0.05, lfc=1.0):
    """Count up/down DEGs for a dict {contrast_name: results_df}."""
    rows = {}
    for name, res in results.items():
        reg = call_degs(res, padj, lfc).value_counts()
        rows[name] = {"up": reg.get("up", 0), "down": reg.get("down", 0)}
    out = pd.DataFrame(rows).T
    out["total"] = out["up"] + out["down"]
    return out

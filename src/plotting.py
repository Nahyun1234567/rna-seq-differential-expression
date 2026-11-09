"""Reusable plots for QC and differential expression."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

CONDITION_COLORS = {"control": "#7f7f7f", "wt": "#1f77b4", "variant": "#d62728"}
REGULATION_COLORS = {"up": "#d62728", "down": "#1f77b4", "ns": "#c7c7c7"}


def pca_plot(log_expr, meta, condition_col="condition", ax=None, n_top=500):
    """PCA on the ``n_top`` most variable genes (samples as points)."""
    top = log_expr.var(axis=1).sort_values(ascending=False).index[:n_top]
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(log_expr.loc[top].T)
    ax = ax or plt.subplots(figsize=(7, 6))[1]
    for cond in meta[condition_col].unique():
        m = (meta[condition_col] == cond).to_numpy()
        ax.scatter(pcs[m, 0], pcs[m, 1], s=100, label=cond,
                   color=CONDITION_COLORS.get(cond))
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    ax.set_title(f"PCA (top {n_top} variable genes)")
    ax.legend()
    return ax


def volcano_plot(res, title="", padj=0.05, lfc=1.0, ax=None, label_top=10):
    """Volcano plot of log2FC vs -log10(padj)."""
    df = res.dropna(subset=["padj"]).copy()
    df["neglog10"] = -np.log10(df["padj"].clip(lower=1e-300))
    sig = (df["padj"] < padj) & (df["log2FoldChange"].abs() > lfc)
    df["reg"] = np.where(sig, np.where(df["log2FoldChange"] > 0, "up", "down"), "ns")

    ax = ax or plt.subplots(figsize=(7, 6))[1]
    for reg in ["ns", "down", "up"]:
        sub = df[df["reg"] == reg]
        ax.scatter(sub["log2FoldChange"], sub["neglog10"], s=6, alpha=0.6,
                   color=REGULATION_COLORS[reg], label=f"{reg} ({len(sub)})")
    ax.axhline(-np.log10(padj), ls="--", lw=0.8, color="k")
    for x in (-lfc, lfc):
        ax.axvline(x, ls="--", lw=0.8, color="k")
    for gene, row in df[df["reg"] != "ns"].nsmallest(label_top, "padj").iterrows():
        ax.annotate(gene, (row["log2FoldChange"], row["neglog10"]), fontsize=7)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(padj)")
    ax.set_title(title)
    ax.legend(loc="upper left", fontsize=8)
    return ax


def ma_plot(res, title="", padj=0.05, ax=None):
    df = res.dropna(subset=["padj"])
    sig = df["padj"] < padj
    ax = ax or plt.subplots(figsize=(7, 5))[1]
    ax.scatter(np.log10(df.loc[~sig, "baseMean"] + 1), df.loc[~sig, "log2FoldChange"],
               s=4, color=REGULATION_COLORS["ns"])
    ax.scatter(np.log10(df.loc[sig, "baseMean"] + 1), df.loc[sig, "log2FoldChange"],
               s=4, color=REGULATION_COLORS["up"])
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("log10(mean normalized count + 1)")
    ax.set_ylabel("log2 fold change")
    ax.set_title(title)
    return ax


def zscore_rows(df):
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1).replace(0, np.nan), axis=0)


def deg_heatmap(log_expr, genes, meta, condition_col="condition", title=""):
    """Clustered heatmap (row z-scores) of the given genes."""
    import seaborn as sns

    order = meta.sort_values(condition_col).index
    data = zscore_rows(log_expr.loc[genes, order]).dropna()
    col_colors = meta.loc[order, condition_col].map(CONDITION_COLORS)
    g = sns.clustermap(data, cmap="RdBu_r", center=0, col_cluster=False,
                       col_colors=col_colors, yticklabels=len(data) <= 60,
                       figsize=(8, 10))
    g.figure.suptitle(title, y=1.02)
    return g

#!/usr/bin/env python3
"""
CLUSTER_QC.py - isoseq cluster2 QC (step 04)
Author: Talha Lone (Parvez lab)

Reads results/qc/cluster_stats_<DATASET>.tsv from 04b_EXTRACT_CLUSTER_STATS.sh
for ALL and for any group you clustered on its own.
Writes results/qc/Cluster_01..07 PNGs.

Keep in mind: if groups were sequenced to different depths, their isoform
counts and singleton rates will differ for that reason alone. Treat the
per-group panels as QC, not as biology. Comparisons between conditions should
be done on the pooled ALL isoform set.

Run:  cd 04_CLUSTERING && python CLUSTER_QC.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import qc_utils as qc  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

QC_DIR = "results/qc"
KEYS = ("flnc_input", "total_isoforms", "support_1", "support_2_5", "support_6_10",
        "support_11_50", "support_50plus", "mean_length", "min_length", "max_length",
        "len_n", "len_0_500", "len_501_1000", "len_1001_2000", "len_2001_3000",
        "len_3001_5000", "len_5001_10000", "len_10000plus")


def load(ds):
    d = qc.read_kv(f"{QC_DIR}/cluster_stats_{ds}.tsv")
    if not d:
        return None
    for k in KEYS:
        if not isinstance(d.get(k), (int, float)):
            d[k] = 0
    return d


DATA = {ds: d for ds in qc.DATASETS for d in [load(ds)] if d}
DS = list(DATA)


def total(g): return DATA[g]["total_isoforms"]
def singletons(g): return DATA[g]["support_1"]
def multi(g): return total(g) - singletons(g)
def singleton_pct(g): return 100 * singletons(g) / total(g) if total(g) else 0
def flnc(g): return DATA[g]["flnc_input"]
def compression(g): return flnc(g) / total(g) if total(g) else 0


def colors():
    return [qc.group_color(g) for g in DS]


def labels():
    return [qc.dataset_label(g) for g in DS]


def width():
    return max(9, 3.2 * len(DS))


def fig1(path):
    head = ["Dataset", "FLNC in", "Isoforms", "Multi-read", "Singletons", "Singleton %",
            "Reads / isoform", "Mean len (bp)"]
    rows = [head]
    for g in DS:
        rows.append([g, f"{flnc(g):,}" if flnc(g) else "n/a", f"{total(g):,}", f"{multi(g):,}",
                     f"{singletons(g):,}", f"{singleton_pct(g):.1f}%",
                     f"{compression(g):.1f}x" if flnc(g) else "n/a",
                     f"{DATA[g]['mean_length']:.0f}"])
    qc.table_fig(rows, "isoseq cluster2 summary", path, row_groups=DS)


def fig2(path):
    vals = [total(g) / 1e6 for g in DS]
    fig, ax = plt.subplots(figsize=(width(), 6))
    top = max(vals) * 1.2
    ax.set_ylim(0, top)
    ax.bar(range(len(DS)), vals, color=colors(), edgecolor="black", linewidth=2, width=0.55)
    for i, v in enumerate(vals):
        ax.text(i, v + top * 0.02, f"{v:.2f}M", ha="center", va="bottom", fontweight="bold", fontsize=12)
    for i, g in enumerate(DS):
        if flnc(g):
            ax.text(i, top * 0.06, f"{compression(g):.1f} reads/isoform", ha="center",
                    fontweight="bold", color="white",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor=qc.group_color(g), alpha=0.85))
    ax.set_xticks(range(len(DS)))
    ax.set_xticklabels(labels())
    ax.set_ylabel("Clustered isoforms (millions)", fontweight="bold")
    ax.set_title("Clustered isoforms per dataset", fontweight="bold")
    qc.save(fig, path)


def fig3(path):
    x = np.arange(len(DS))
    w = 0.35
    fin = [flnc(g) / 1e6 for g in DS]
    iso = [total(g) / 1e6 for g in DS]
    fig, ax = plt.subplots(figsize=(width(), 6))
    top = max(fin + iso) * 1.22
    ax.set_ylim(0, top)
    ax.bar(x - w / 2, fin, w, label="FLNC reads in", color=[qc.group_light(g) for g in DS], edgecolor="black")
    ax.bar(x + w / 2, iso, w, label="isoforms out", color=colors(), edgecolor="black")
    for i in range(len(DS)):
        ax.text(x[i] - w / 2, fin[i] + top * 0.01, f"{fin[i]:.1f}M", ha="center", va="bottom",
                fontsize=9, fontweight="bold")
        ax.text(x[i] + w / 2, iso[i] + top * 0.01, f"{iso[i]:.2f}M", ha="center", va="bottom",
                fontsize=9, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels())
    ax.set_ylabel("reads / isoforms (millions)", fontweight="bold")
    ax.set_title("FLNC reads in vs isoforms out", fontweight="bold")
    ax.legend()
    qc.save(fig, path)


def fig4(path):
    m = [multi(g) / 1e6 for g in DS]
    s = [singletons(g) / 1e6 for g in DS]
    fig, ax = plt.subplots(figsize=(width(), 6))
    top = max(total(g) for g in DS) / 1e6 * 1.15
    ax.set_ylim(0, top)
    ax.bar(range(len(DS)), m, color=colors(), edgecolor="black", linewidth=2, width=0.55,
           label="multi-read")
    ax.bar(range(len(DS)), s, bottom=m, color=[qc.group_light(g) for g in DS], edgecolor="black",
           linewidth=2, width=0.55, label="singletons")
    for i, g in enumerate(DS):
        ax.text(i, m[i] / 2, f"{m[i]:.2f}M\n({100 - singleton_pct(g):.1f}%)", ha="center",
                va="center", color="white", fontweight="bold", fontsize=9)
        ax.text(i, m[i] + s[i] / 2, f"{s[i]:.2f}M\n({singleton_pct(g):.1f}%)", ha="center",
                va="center", fontweight="bold", fontsize=9)
        ax.text(i, total(g) / 1e6 + top * 0.02, f"total {total(g) / 1e6:.2f}M", ha="center",
                va="bottom", fontweight="bold", fontsize=10)
    ax.set_xticks(range(len(DS)))
    ax.set_xticklabels(labels())
    ax.set_ylabel("isoforms (millions)", fontweight="bold")
    ax.set_title("Singletons vs multi-read isoforms", fontweight="bold")
    ax.legend(loc="upper right")
    qc.save(fig, path)


def _grouped(path, bins, keys, xlabel, title, note=None, legend_loc="upper left"):
    x = np.arange(len(bins))
    w = 0.8 / len(DS)
    fig, ax = plt.subplots(figsize=(13, 6))
    for i, g in enumerate(DS):
        ax.bar(x + (i - (len(DS) - 1) / 2) * w, [DATA[g][k] / 1e6 for k in keys], w,
               label=g, color=qc.group_color(g), edgecolor="black")
    ax.set_xticks(x)
    ax.set_xticklabels(bins)
    ax.set_ylabel("isoforms (millions)", fontweight="bold")
    ax.set_xlabel(xlabel, fontweight="bold")
    ax.set_title(title, fontweight="bold")
    ax.legend(loc=legend_loc)
    if note:
        ax.text(0.98, 0.97, note, transform=ax.transAxes, ha="right", va="top", fontsize=10,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.95, edgecolor="gray"))
    qc.save(fig, path)


def fig5(path):
    _grouped(path, ["1\n(singleton)", "2-5", "6-10", "11-50", ">50"],
             ["support_1", "support_2_5", "support_6_10", "support_11_50", "support_50plus"],
             "FLNC reads per isoform", "Read support per isoform", legend_loc="upper right")


def fig6(path):
    note = "\n".join(f"mean {g}: {DATA[g]['mean_length']:.0f} bp" for g in DS)
    _grouped(path, ["<=500", "501-1k", "1-2k", "2-3k", "3-5k", "5-10k", ">10k"],
             ["len_0_500", "len_501_1000", "len_1001_2000", "len_2001_3000",
              "len_3001_5000", "len_5001_10000", "len_10000plus"],
             "isoform length (bp)", "Isoform length", note)


def fig7(path):
    fig, axes = plt.subplots(1, 3, figsize=(max(13, 4.5 * len(DS)), 6))
    panels = [
        ([total(g) / 1e6 for g in DS], "Isoforms", "millions", lambda v: f"{v:.2f}M"),
        ([singleton_pct(g) for g in DS], "Singleton %", "%", lambda v: f"{v:.1f}%"),
        ([DATA[g]["mean_length"] for g in DS], "Mean isoform length", "bp", lambda v: f"{v:.0f} bp"),
    ]
    for ax, (vals, title, ylab, f) in zip(axes, panels):
        top = 100 if ylab == "%" else max(vals) * 1.2
        ax.set_ylim(0, top)
        ax.bar(range(len(DS)), vals, color=colors(), edgecolor="black", width=0.55)
        for i, v in enumerate(vals):
            ax.text(i, v + top * 0.02, f(v), ha="center", va="bottom", fontweight="bold", fontsize=9)
        ax.set_title(title, fontweight="bold")
        ax.set_xticks(range(len(DS)))
        ax.set_xticklabels(DS, rotation=30 if len(DS) > 4 else 0)
        ax.set_ylabel(ylab, fontweight="bold")
    fig.suptitle("Dataset comparison (QC only; differences can just be sequencing depth)",
                 fontsize=14, fontweight="bold")
    qc.save(fig, path)


def main():
    if not DS:
        raise SystemExit(f"no cluster_stats_*.tsv in {QC_DIR} (run 04b_EXTRACT_CLUSTER_STATS.sh)")
    os.makedirs(QC_DIR, exist_ok=True)
    qc.use_style()
    for g in DS:
        print(f"{g}: isoforms={total(g):,} singletons={singletons(g):,} ({singleton_pct(g):.1f}%) "
              f"reads/isoform={compression(g):.1f} mean_len={DATA[g]['mean_length']:.0f}")
        if not flnc(g):
            print(f"  WARNING: flnc_input missing for {g}, reads-in figures will be empty")
    fig1(f"{QC_DIR}/Cluster_01_summary.png")
    fig2(f"{QC_DIR}/Cluster_02_total_isoforms.png")
    fig3(f"{QC_DIR}/Cluster_03_compression.png")
    fig4(f"{QC_DIR}/Cluster_04_singleton_composition.png")
    fig5(f"{QC_DIR}/Cluster_05_support_distribution.png")
    fig6(f"{QC_DIR}/Cluster_06_length_distribution.png")
    fig7(f"{QC_DIR}/Cluster_07_dataset_comparison.png")
    print(f"figures in {QC_DIR}/")


if __name__ == "__main__":
    main()

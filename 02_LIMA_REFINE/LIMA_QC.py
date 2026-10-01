#!/usr/bin/env python3
"""
LIMA_QC.py - demultiplexing QC (step 02a)
Author: Talha Lone (Parvez lab)

Reads results/demux/<cell>.lima.summary and <cell>.lima.counts for every cell
and adds counts up per sample using samples.tsv. Read lengths are taken from
the first LIMA_LENGTH_READS reads of each demuxed BAM (all BAMs of a sample
pooled), so it runs fine on a login node. Needs samtools on PATH for the
length figures; the rest works without it.
Writes results/qc/Lima_01..10 PNGs.

Run:  cd 02_LIMA_REFINE && python LIMA_QC.py
      LIMA_LENGTH_READS=50000 python LIMA_QC.py     (quicker)
"""
import glob
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import qc_utils as qc  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

DEMUX_DIR = "results/demux"
OUT_DIR = "results/qc"
LENGTH_READS = int(os.environ.get("LIMA_LENGTH_READS", 200_000))   # per sample
MIN_ASSIGN_PCT = 80.0
MIN_FL_READS = 1_000_000
SAMTOOLS = shutil.which("samtools")


# ---------------------------------------------------------------- parsing
def parse_summaries():
    """lima.summary numbers, summed over cells."""
    out = {}
    for f in sorted(glob.glob(f"{DEMUX_DIR}/*.lima.summary")):
        for line in open(f):
            m = re.match(r"\s*(.+?)\s*:\s*([\d,]+)", line)
            if m:
                k = m.group(1).strip()
                out[k] = out.get(k, 0) + int(m.group(2).replace(",", ""))
    return out


def parse_counts():
    """reads per sample from lima.counts (summed over cells)"""
    counts, unknown = {}, {}
    for f in sorted(glob.glob(f"{DEMUX_DIR}/*.lima.counts")):
        cell = os.path.basename(f)[:-len(".lima.counts")]
        for line in open(f):
            p = line.rstrip("\n").split("\t")
            if len(p) < 5 or p[0] == "IdxFirst":
                continue
            pair = f"{p[2]}--{p[3]}"
            s = qc.sample_for(cell, pair)
            if s is None:
                unknown[f"{cell}:{pair}"] = int(p[4])
                continue
            counts[s] = counts.get(s, 0) + int(p[4])
    for k, v in unknown.items():
        print(f"  not in samples.tsv: {k} ({v:,} reads)")
    return counts


def sample_bams():
    out = {}
    for bam in sorted(glob.glob(f"{DEMUX_DIR}/*.bam")):
        if "unbarcoded" in bam:
            continue
        s = qc.sample_for(*qc.split_name(bam))
        if s:
            out.setdefault(s, []).append(bam)
    return out


def read_lengths(bams, n):
    if not SAMTOOLS or not bams:
        return np.array([])
    per = max(1, n // len(bams))
    L = []
    for bam in bams:
        r = subprocess.run(f'"{SAMTOOLS}" view "{bam}" 2>/dev/null | head -n {per} '
                           f"| awk '{{print length($10)}}'",
                           shell=True, capture_output=True, text=True)
        L += [int(x) for x in r.stdout.split()]
    return np.array(L)


def sample_ticks(ax, ss):
    """sample names under the bars, slanted if they are long"""
    long_names = max(len(s) for s in ss) > 10
    ax.set_xticks(range(len(ss)))
    ax.set_xticklabels([qc.sample_label(s) for s in ss], fontsize=9,
                       rotation=30 if long_names else 0, ha="right" if long_names else "center")


def present(d):
    return [s for s in qc.SAMPLES_BY_GROUP if d.get(s) is not None and
            (len(d[s]) > 0 if hasattr(d[s], "__len__") else d[s])]


def total_in(lima):
    return next((lima[k] for k in lima if k.lower().startswith("reads input")), None)


def total_assigned(lima):
    return next((lima[k] for k in lima if "above all thresholds" in k.lower()), None)


def total_unassigned(lima):
    return next((lima[k] for k in lima if "below any threshold" in k.lower()), None)


def group_legend(ax, extra=None):
    h = [Patch(color=qc.group_color(g), label=g) for g in qc.GROUPS]
    ax.legend(handles=h + (extra or []))


# ---------------------------------------------------------------- figures
def fig1_summary(lima, counts, path):
    inp, asg, una = total_in(lima), total_assigned(lima), total_unassigned(lima)
    rows = [["Metric", "Value", "Expected", "Status"]]
    st = lambda ok: "PASS" if ok else "REVIEW"
    ncell = len(glob.glob(f"{DEMUX_DIR}/*.lima.summary"))
    rows.append(["Cells", f"{ncell}", "-", "PASS"])
    if inp:
        rows.append(["Input S-reads", f"{inp:,}", "-", "PASS"])
    if asg and inp:
        rows.append(["Assigned", f"{asg:,}", "-", "PASS"])
        rows.append(["Assignment rate", f"{100 * asg / inp:.2f}%", f">{MIN_ASSIGN_PCT:.0f}%",
                     st(100 * asg / inp > MIN_ASSIGN_PCT)])
    if una and inp:
        rows.append(["Unassigned", f"{una:,}", "-", "PASS"])
        rows.append(["Unassigned rate", f"{100 * una / inp:.2f}%", f"<{100 - MIN_ASSIGN_PCT:.0f}%",
                     st(100 * una / inp < 100 - MIN_ASSIGN_PCT)])
    nexp = len(qc.SAMPLES)
    rows.append(["Samples recovered", f"{len(counts)}", f"{nexp}", st(len(counts) == nexp)])
    if counts:
        rows.append(["Assigned to sheet samples", f"{sum(counts.values()):,}", "~ assigned", "PASS"])
    qc.table_fig(rows, "Lima demultiplexing summary", path, status_col=3,
                 col_widths=[0.34, 0.22, 0.2, 0.12], fontsize=11)


def fig2_distribution(counts, path):
    ss = present(counts)
    if not ss:
        return
    vals = np.array([counts[s] for s in ss]) / 1e6
    tot = vals.sum()
    fig, ax = plt.subplots(figsize=(max(9, 1.4 * len(ss)), 6))
    top = vals.max() * 1.25
    ax.set_ylim(0, top)
    ax.bar(range(len(ss)), vals, color=[qc.sample_color(s) for s in ss],
           edgecolor="black", linewidth=1.5)
    for i, v in enumerate(vals):
        ax.text(i, v + top * 0.012, f"{v:.2f}M\n({100 * v / tot:.1f}%)", ha="center",
                va="bottom", fontweight="bold", fontsize=9)
    ax.axhline(vals.mean(), color="black", ls="--", lw=2)
    qc.group_blocks(ax, ss, top)
    sample_ticks(ax, ss)
    ax.set_ylabel("Reads assigned (millions)", fontweight="bold")
    ax.set_title("Reads per sample after lima", fontweight="bold")
    group_legend(ax, [plt.Line2D([], [], color="black", ls="--", label=f"mean {vals.mean():.2f}M")])
    qc.save(fig, path)


def fig3_failures(lima, path):
    keys = [k for k in lima if lima[k] > 0 and
            k.strip().startswith(("Below ", "Undesired", "Without", "Wrong"))]
    if not keys:
        return
    keys = sorted(keys, key=lambda k: lima[k], reverse=True)
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.bar(range(len(keys)), [lima[k] / 1e3 for k in keys], color="#7A7A7A", edgecolor="black")
    for i, k in enumerate(keys):
        ax.text(i, lima[k] / 1e3, f"{lima[k] / 1e3:.0f}K", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([k.replace("Below ", "").replace("Undesired ", "") for k in keys],
                       rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("reads (thousands)", fontweight="bold")
    ax.set_title("Unassigned reads by reason (all cells)", fontweight="bold")
    qc.save(fig, path)


def fig4_group_balance(counts, path):
    groups = [g for g in qc.GROUPS if any(qc.GROUP_OF[s] == g for s in counts)]
    if not groups:
        return
    data = [[counts[s] / 1e6 for s in qc.SAMPLES_BY_GROUP if qc.GROUP_OF[s] == g and s in counts]
            for g in groups]
    fig, ax = plt.subplots(figsize=(max(7, 2.2 * len(groups)), 6))
    bp = ax.boxplot(data, patch_artist=True, widths=0.5,
                    medianprops=dict(color="black", linewidth=2))
    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels([f"{g} (n={len(d)})" for g, d in zip(groups, data)])
    for patch, g in zip(bp["boxes"], groups):
        patch.set_facecolor(qc.group_color(g))
        patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for x, g in enumerate(groups, 1):
        ss = [s for s in qc.SAMPLES_BY_GROUP if qc.GROUP_OF[s] == g and s in counts]
        for s in ss:
            ax.scatter(x + rng.uniform(-0.05, 0.05), counts[s] / 1e6, color=qc.sample_color(s),
                       edgecolor="black", zorder=5, s=60)
    allc = [v for d in data for v in d]
    cv = 100 * np.std(allc) / np.mean(allc) if allc else 0
    ax.set_ylabel("Reads per sample (millions)", fontweight="bold")
    ax.set_title(f"Depth per group (CV across all samples {cv:.1f}%)", fontweight="bold")
    qc.save(fig, path)


def fig5_efficiency(lima, path):
    inp, asg = total_in(lima), total_assigned(lima)
    if not (inp and asg):
        return
    rate = 100 * asg / inp
    ok = rate > MIN_ASSIGN_PCT
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.axis("off")
    ax.text(0.5, 0.75, "Overall assignment rate", ha="center", fontsize=18,
            fontweight="bold", transform=ax.transAxes)
    ax.text(0.5, 0.45, f"{rate:.2f}% assigned   ({asg:,} / {inp:,})", ha="center",
            fontsize=18, transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=1", facecolor="lightgreen" if ok else "#F4C84D",
                      alpha=0.6, edgecolor="black", linewidth=2))
    ax.text(0.5, 0.15, f"PASS (above {MIN_ASSIGN_PCT:.0f}%)" if ok else f"REVIEW (below {MIN_ASSIGN_PCT:.0f}%)",
            ha="center", fontsize=13, fontweight="bold",
            color="green" if ok else "#A65900", transform=ax.transAxes)
    qc.save(fig, path)


def fig6_length_table(lens, path):
    ss = present(lens)
    if not ss:
        return
    rows = [["Sample", "Group", "n (subsample)", "Mean (bp)", "Median", "SD", "IQR"]]
    for s in ss:
        L = lens[s]
        rows.append([s, qc.GROUP_OF[s], f"{len(L):,}", f"{L.mean():.0f}", f"{np.median(L):.0f}",
                     f"{L.std():.0f}", f"{np.percentile(L, 25):.0f}-{np.percentile(L, 75):.0f}"])
    qc.table_fig(rows, "Read length by sample (subsampled)", path,
                 row_groups=[qc.GROUP_OF[s] for s in ss])


def fig7_length_grid(lens, path):
    ss = present(lens)
    if not ss:
        return
    ncol = min(3, len(ss))
    nrow = int(np.ceil(len(ss) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6 * ncol, 4.5 * nrow), squeeze=False)
    axes = axes.flatten()
    for i, s in enumerate(ss):
        ax, L = axes[i], lens[s]
        ax.hist(L, bins=50, color=qc.sample_color(s), alpha=0.8, edgecolor="black")
        ax.axvline(L.mean(), color="red", ls="--", lw=2, label=f"mean {L.mean():.0f}")
        ax.set_title(f"{s} ({qc.GROUP_OF[s]})", fontsize=11, fontweight="bold")
        ax.set_xlabel("read length (bp)")
        ax.set_ylabel("reads")
        ax.legend(fontsize=8)
        ax.set_xlim(0, min(10000, L.max()))
    for j in range(len(ss), len(axes)):
        axes[j].axis("off")
    fig.suptitle("Read length distributions", fontsize=16, fontweight="bold")
    qc.save(fig, path)


def fig8_length_box(lens, path):
    ss = present(lens)
    if not ss:
        return
    fig, ax = plt.subplots(figsize=(max(10, 1.4 * len(ss)), 6))
    bp = ax.boxplot([lens[s] for s in ss], positions=range(len(ss)),
                    patch_artist=True, widths=0.6, showfliers=False,
                    medianprops=dict(color="black", linewidth=2))
    for patch, s in zip(bp["boxes"], ss):
        patch.set_facecolor(qc.sample_color(s))
        patch.set_alpha(0.7)
    sample_ticks(ax, ss)
    ymax = ax.get_ylim()[1] * 1.1
    ax.set_ylim(0, ymax)
    qc.group_blocks(ax, ss, ymax)
    ax.set_ylabel("read length (bp)", fontweight="bold")
    ax.set_title("Read length by sample", fontweight="bold")
    qc.save(fig, path)


def fig9_mean_length(lens, path):
    ss = present(lens)
    if not ss:
        return
    means = [lens[s].mean() for s in ss]
    sds = [lens[s].std() for s in ss]
    fig, ax = plt.subplots(figsize=(max(10, 1.4 * len(ss)), 6))
    top = max(m + d for m, d in zip(means, sds)) * 1.2
    ax.set_ylim(0, top)
    ax.bar(range(len(ss)), means, yerr=sds, capsize=5, color=[qc.sample_color(s) for s in ss],
           edgecolor="black", error_kw={"ecolor": "black"})
    for i, (m, d) in enumerate(zip(means, sds)):
        ax.text(i, m + d + top * 0.01, f"{m:.0f}", ha="center", va="bottom",
                fontweight="bold", fontsize=9)
    ax.axhline(np.mean(means), color="black", ls="--", lw=2,
               label=f"overall mean {np.mean(means):.0f} bp")
    qc.group_blocks(ax, ss, top)
    sample_ticks(ax, ss)
    ax.set_ylabel("mean read length (bp)", fontweight="bold")
    ax.set_title("Mean read length by sample (+/- SD)", fontweight="bold")
    ax.legend(loc="lower right")
    qc.save(fig, path)


def fig10_low_samples(counts, path):
    ss = present(counts)
    if not ss:
        return
    vals = np.array([counts[s] for s in ss], float)
    mean = vals.mean()
    thr = mean - vals.std()
    fig, ax = plt.subplots(figsize=(max(10, 1.4 * len(ss)), 6))
    top = vals.max() / 1e6 * 1.2
    ax.set_ylim(0, top)
    bars = ax.bar(range(len(ss)), vals / 1e6, color=[qc.sample_color(s) for s in ss],
                  edgecolor="black", linewidth=2)
    for i, v in enumerate(vals):
        if v < thr:
            bars[i].set_edgecolor("#CC0000")
            bars[i].set_linewidth(3.5)
        ax.text(i, v / 1e6 + top * 0.012, f"{v / 1e6:.2f}M", ha="center", va="bottom",
                fontweight="bold", fontsize=9)
    ax.axhline(mean / 1e6, color="red", ls="--", lw=2.5, label=f"mean {mean / 1e6:.2f}M")
    ax.axhline(thr / 1e6, color="orange", ls=":", lw=2, label=f"mean - 1 SD {thr / 1e6:.2f}M")
    ax.axhline(MIN_FL_READS / 1e6, color="green", ls="--", lw=2,
               label=f"floor {MIN_FL_READS / 1e6:.0f}M")
    qc.group_blocks(ax, ss, top)
    sample_ticks(ax, ss)
    ax.set_ylabel("reads (millions)", fontweight="bold")
    ax.set_title("Low-depth check (red outline = below mean - 1 SD)", fontweight="bold")
    ax.legend(loc="lower right")
    qc.save(fig, path)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    qc.use_style()
    lima = parse_summaries()
    if not lima:
        raise SystemExit(f"no *.lima.summary in {DEMUX_DIR} (run 02a first)")
    counts = parse_counts()
    print("reads per sample:")
    for s in qc.SAMPLES_BY_GROUP:
        print(f"  {s:<25} {qc.GROUP_OF[s]:<10} {counts.get(s, 0):>12,}")
        if s not in counts:
            print(f"    WARNING: no reads for {s}")

    fig1_summary(lima, counts, f"{OUT_DIR}/Lima_01_summary.png")
    fig2_distribution(counts, f"{OUT_DIR}/Lima_02_distribution.png")
    fig3_failures(lima, f"{OUT_DIR}/Lima_03_failure_modes.png")
    fig4_group_balance(counts, f"{OUT_DIR}/Lima_04_group_balance.png")
    fig5_efficiency(lima, f"{OUT_DIR}/Lima_05_assignment_efficiency.png")

    if not SAMTOOLS:
        print("samtools not on PATH, skipping the read length figures (6-9)")
        lens = {}
    else:
        print(f"reading lengths ({LENGTH_READS:,} reads per sample)...")
        lens = {s: read_lengths(b, LENGTH_READS) for s, b in sample_bams().items()}
    fig6_length_table(lens, f"{OUT_DIR}/Lima_06_length_table.png")
    fig7_length_grid(lens, f"{OUT_DIR}/Lima_07_length_distributions.png")
    fig8_length_box(lens, f"{OUT_DIR}/Lima_08_length_boxplot.png")
    fig9_mean_length(lens, f"{OUT_DIR}/Lima_09_mean_length.png")
    fig10_low_samples(counts, f"{OUT_DIR}/Lima_10_low_depth_check.png")
    print(f"figures in {OUT_DIR}/")


if __name__ == "__main__":
    main()

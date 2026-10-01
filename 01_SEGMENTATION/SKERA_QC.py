#!/usr/bin/env python3
"""
SKERA_QC.py - library-level QC for skera split (step 01)
Author: Talha Lone (Parvez lab)

Reads, for every cell in results/segmented/:
  <cell>.segmented.summary.json    read counts, mean array size, % full arrays
  <cell>.segmented.ligations.csv   adapter ligation table
  results/qc/<cell>_array_hist.tsv exact segments-per-read table from 01b (optional)
Writes results/qc/<cell>_QC1..QC4 PNGs.

If the 01b table is not there, the segments-per-read distribution is estimated
from the first few million S-reads (needs samtools on PATH). That part is
cheap and fine on a login node; it does not read the whole BAM.

Per-sample numbers come after lima (step 02), not here.

Run:  cd 01_SEGMENTATION && python SKERA_QC.py
"""
import glob
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import qc_utils as qc  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

SEG_DIR = "results/segmented"
OUT_DIR = "results/qc"

EXPECTED_ARRAY = qc.CFG["EXPECTED_ARRAY_SIZE"]
ARRAY_SIZE_MIN = EXPECTED_ARRAY - 1      # mean array size below this = look closer
FULL_ARRAY_MIN = 90.0                    # % reads with a full array
SAMTOOLS_THREADS = 4
SUBSAMPLE_NREADS = 3_000_000             # S-reads read for the quick estimate

DARK = qc.ALL_COLOR    # before demux everything is one pooled library -> purple
LIGHT = qc.ALL_LIGHT
MAS8_LABELS = ["A", "B", "C", "D", "E", "F", "G", "H", "Q"]

# run-length count of S-reads per ZMW (ZMW = 2nd field of the read name)
AWK_HIST = (r'{split($1,a,"/"); z=a[2]; '
            r'if(NR>1 && z!=p){h[c]++; c=0} c++; p=z} '
            r'END{if(c>0)h[c]++; for(k in h) print k"\t"h[k]}')


def load_metrics(path):
    with open(path) as fh:
        s = json.load(fh)
    return {a["id"]: a["value"] for a in s["attributes"]}


def array_histogram(bam, cache):
    """(hist, exact). Uses the 01b table if present, otherwise a subsample."""
    if os.path.exists(cache):
        df = pd.read_csv(cache, sep="\t")
        return dict(zip(df["segments"], df["num_reads"])), True
    if shutil.which("samtools") is None:
        print("  samtools not on PATH, skipping the distribution plot")
        return None, None
    print(f"  estimating segments per read from the first {SUBSAMPLE_NREADS:,} S-reads")
    cmd = (f'samtools view -@ {SAMTOOLS_THREADS} "{bam}" 2>/dev/null '
           f'| head -n {SUBSAMPLE_NREADS} '
           f'| awk -F"\\t" \'{AWK_HIST}\'')
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if not res.stdout.strip():
        print("  subsample failed:", (res.stderr or "no output")[:300])
        return None, None
    hist = {}
    for line in res.stdout.strip().splitlines():
        k, v = line.split("\t")
        hist[int(k)] = int(v)
    return hist, False


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    qc.use_style()
    summaries = sorted(glob.glob(os.path.join(SEG_DIR, "*.segmented.summary.json")))
    if not summaries:
        raise SystemExit(f"no *.segmented.summary.json in {SEG_DIR} (run 01a first)")

    for summary_path in summaries:
        m = load_metrics(summary_path)
        base = os.path.basename(summary_path).replace(".segmented.summary.json", "")

        reads, sreads = m["reads"], m["s_reads"]
        mean_arr, pfa, slen = m["mean_array_size"], m["percent_full_array"], m["mean_len_s_reads"]
        amp = sreads / reads if reads else 0
        arr_ok, fa_ok = mean_arr >= ARRAY_SIZE_MIN, pfa >= FULL_ARRAY_MIN
        overall = "PASS" if (arr_ok and fa_ok) else "REVIEW"

        print("=" * 72)
        print(f"skera QC  {base}")
        print("=" * 72)
        print(f"  HiFi reads in          : {reads:>14,}")
        print(f"  S-reads out            : {sreads:>14,}")
        print(f"  fold                   : {amp:>13.2f}x   (expect ~{EXPECTED_ARRAY})")
        print(f"  mean array size        : {mean_arr:>13.2f}    ({'ok' if arr_ok else 'LOW'})")
        print(f"  % reads w/ full array  : {pfa:>13.2f}%   (want >= {FULL_ARRAY_MIN:.0f}%, "
              f"{'ok' if fa_ok else 'LOW'})  <- the one to look at")
        print(f"  mean S-read length     : {slen:>13,} bp")
        print(f"  -> {overall}")

        # QC1: table
        rows = [
            ["Metric", "Value", "Expected", "Status"],
            ["Input HiFi reads", f"{reads:,}", "-", "OK"],
            ["Output S-reads", f"{sreads:,}", f"~{EXPECTED_ARRAY}x input", "OK"],
            ["Fold increase", f"{amp:.2f}x", f"~{EXPECTED_ARRAY}x", "OK" if arr_ok else "LOW"],
            ["Mean array size", f"{mean_arr:.2f}", f"~{EXPECTED_ARRAY}", "OK" if arr_ok else "LOW"],
            ["% reads w/ full array", f"{pfa:.2f}%", f">={FULL_ARRAY_MIN:.0f}%", "OK" if fa_ok else "LOW"],
            ["Mean S-read length", f"{slen:,} bp", "1,000-3,000 bp", "OK"],
        ]
        fig, ax = plt.subplots(figsize=(11, 4.5))
        ax.axis("off")
        tbl = ax.table(cellText=rows, cellLoc="center", loc="center",
                       colWidths=[0.34, 0.22, 0.24, 0.12])
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(11)
        tbl.scale(1, 2.4)
        for j in range(4):
            c = tbl[(0, j)]
            c.set_facecolor(DARK)
            c.set_text_props(weight="bold", color="white")
        for i in range(1, len(rows)):
            for j in range(4):
                c = tbl[(i, j)]
                if j == 3:
                    c.set_facecolor("#90EE90" if rows[i][3] == "OK" else "#F4A6A6")
                    c.set_text_props(weight="bold")
                else:
                    c.set_facecolor("#F2F0F7" if i % 2 == 0 else "white")
        ax.set_title(f"Skera segmentation QC ({overall})\n{base}", fontsize=15,
                     fontweight="bold", pad=16)
        qc.save(fig, os.path.join(OUT_DIR, f"{base}_QC1_summary.png"))

        # QC2: reads in vs out
        fig, ax = plt.subplots(figsize=(7, 6))
        bars = ax.bar(["Input\nHiFi reads", "Output\nS-reads"], [reads / 1e6, sreads / 1e6],
                      color=[LIGHT, DARK], edgecolor="black", linewidth=1.5, width=0.6)
        top = sreads / 1e6 * 1.22
        ax.set_ylim(0, top)
        for b in bars:
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + top * 0.012,
                    f"{b.get_height():.1f}M", ha="center", va="bottom", fontweight="bold")
        ax.text(0.5, sreads / 1e6 * 1.08, f"{amp:.2f}x", ha="center", va="center",
                fontsize=13, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.5", facecolor="gold", alpha=0.7))
        ax.set_ylabel("Reads (millions)", fontweight="bold")
        ax.set_title("Segmentation fold increase", fontsize=14, fontweight="bold", pad=16)
        qc.save(fig, os.path.join(OUT_DIR, f"{base}_QC2_amplification.png"))

        # QC3: mean array size + S-read length
        fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 5))
        axL.bar(["Mean array size"], [mean_arr], color=DARK, edgecolor="black", width=0.5)
        axL.axhline(EXPECTED_ARRAY, color="red", ls="--", lw=1.5, label=f"expected {EXPECTED_ARRAY}")
        axL.set_ylim(0, max(EXPECTED_ARRAY, mean_arr) * 1.2)
        axL.text(0, mean_arr + 0.15, f"{mean_arr:.2f}", ha="center", va="bottom",
                 fontweight="bold", fontsize=13)
        axL.set_ylabel("Segments per HiFi read", fontweight="bold")
        axL.set_title("Mean segments per read", fontweight="bold")
        axL.legend(loc="lower right")
        axR.bar(["Mean S-read length"], [slen], color=LIGHT, edgecolor="black", width=0.5)
        axR.set_ylim(0, max(3000, slen * 1.3))
        axR.text(0, slen * 1.02, f"{slen:,} bp", ha="center", va="bottom",
                 fontweight="bold", fontsize=13)
        axR.set_ylabel("Length (bp)", fontweight="bold")
        axR.set_title("Mean S-read length", fontweight="bold")
        qc.save(fig, os.path.join(OUT_DIR, f"{base}_QC3_means.png"))

        # QC4: segments-per-read distribution + ligation heatmap
        bam = summary_path.replace(".segmented.summary.json", ".segmented.bam")
        hist, exact = array_histogram(bam, os.path.join(OUT_DIR, f"{base}_array_hist.tsv"))

        fig, (axD, axH) = plt.subplots(1, 2, figsize=(16, 7))
        if hist:
            ks = sorted(hist)
            cnt = np.array([hist[k] for k in ks], dtype=float)
            pct = 100 * cnt / cnt.sum()
            colors = [DARK if k == EXPECTED_ARRAY else LIGHT for k in ks]
            axD.bar(ks, pct, color=colors, edgecolor="black", width=0.8)
            axD.axvline(mean_arr, color="red", ls="--", lw=2.5, label=f"mean {mean_arr:.2f}")
            axD.set_xlabel("Segments per HiFi read", fontweight="bold")
            axD.set_ylabel("% of HiFi reads", fontweight="bold")
            src = "exact" if exact else f"estimated from {int(cnt.sum()):,} reads"
            axD.set_title(f"Segments per read ({src})", fontweight="bold")
            axD.set_xticks(ks)
            axD.legend()
        else:
            axD.text(0.5, 0.5, "segments-per-read distribution not available\n"
                               "(load samtools, or run 01b_ARRAY_HIST.sh)",
                     ha="center", va="center", transform=axD.transAxes)
            axD.set_axis_off()

        lig_path = summary_path.replace(".summary.json", ".ligations.csv")
        try:
            lig = pd.read_csv(lig_path)
            cols = {c.lower(): c for c in lig.columns}
            if "ligations" in cols and sum("adapter" in c for c in cols) >= 2:
                ac = [cols[c] for c in cols if "adapter" in c]
                a1 = cols.get("adapter_1", ac[0])
                a2 = cols.get("adapter_2", ac[1])
                mat = lig.pivot_table(index=a1, columns=a2, values=cols["ligations"],
                                      aggfunc="sum", fill_value=0)
            else:
                mat = lig.set_index(lig.columns[0])
            vals = mat.values.astype(float)
            rs = vals.sum(axis=1, keepdims=True)
            frac = np.divide(vals, rs, out=np.zeros_like(vals), where=rs != 0)
            if frac.shape == (len(MAS8_LABELS), len(MAS8_LABELS)):
                xlab = ylab = MAS8_LABELS
            else:
                xlab = [str(c) for c in mat.columns]
                ylab = [str(r) for r in mat.index]
            im = axH.imshow(frac, cmap="Purples", vmin=0, vmax=1, aspect="auto")
            axH.set_xticks(range(frac.shape[1]))
            axH.set_xticklabels(xlab, fontsize=10)
            axH.set_yticks(range(frac.shape[0]))
            axH.set_yticklabels(ylab, fontsize=10)
            axH.set_xlabel("Following adapter", fontweight="bold")
            axH.set_ylabel("Leading adapter", fontweight="bold")
            axH.set_title("Adapter ligations (row-normalised)\na clean run is one bright diagonal band",
                          fontweight="bold")
            for i in range(frac.shape[0]):
                for j in range(frac.shape[1]):
                    if frac[i, j] >= 0.15:
                        axH.text(j, i, f"{frac[i, j]:.2f}", ha="center", va="center",
                                 color="white" if frac[i, j] >= 0.6 else "black", fontsize=8)
            fig.colorbar(im, ax=axH, label="fraction of the leading adapter's ligations")
        except Exception as e:
            axH.text(0.5, 0.5, f"ligation heatmap not available\n({type(e).__name__})",
                     ha="center", va="center", transform=axH.transAxes)
            axH.set_axis_off()
            print(f"  ligation heatmap skipped ({type(e).__name__}: {e})")
        qc.save(fig, os.path.join(OUT_DIR, f"{base}_QC4_distribution_ligation.png"))

    print("\ndone. For the exact distribution run 01b_ARRAY_HIST.sh, then rerun this.")


if __name__ == "__main__":
    main()

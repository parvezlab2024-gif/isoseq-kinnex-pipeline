#!/usr/bin/env python3
"""
REFINE_QC.py - isoseq refine QC (step 02b, plus 02c if you ran it)
Author: Talha Lone (Parvez lab)

Reads the refine reports, adds them up per sample (a sample can be on
more than one cell):
  results/refined/*.flnc.filter_summary.report.json          (--require-polya)
  results/refined_nopolya/*.flnc.filter_summary.report.json  (optional, 02c)
Only the JSON reports are read, not the BAMs, so it is quick.
Writes results/qc/Refine_01..08 PNGs. Figures 7-8 need the 02c run.

Run:  cd 02_LIMA_REFINE && python REFINE_QC.py
"""
import glob
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import qc_utils as qc  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REFINED_DIR = "results/refined"
NOPOLYA_DIR = "results/refined_nopolya"
OUT_DIR = "results/qc"
MIN_FLNC = 1_000_000     # per sample, below this the sample is flagged LOW

KEYS = {"num_reads_fl": "fl", "num_reads_flnc": "flnc", "num_reads_flnc_polya": "polya"}


def read_reports(d):
    """{sample: {fl, flnc, polya}} summed over every report for that sample"""
    out = {}
    for f in sorted(glob.glob(f"{d}/*.flnc.filter_summary.report.json")):
        bam = f.replace(".filter_summary.report.json", ".bam")
        s = qc.sample_for(*qc.split_name(bam))
        if s is None:
            print(f"  not in samples.tsv, ignored: {os.path.basename(f)}")
            continue
        try:
            attrs = {a["id"]: a["value"] for a in json.load(open(f))["attributes"]}
        except Exception as e:
            print(f"  could not read {f}: {e}")
            continue
        row = out.setdefault(s, {})
        for k, name in KEYS.items():
            if k in attrs:
                row[name] = row.get(name, 0) + int(attrs[k])
    return out


def build():
    rows = read_reports(REFINED_DIR)
    nopolya = read_reports(NOPOLYA_DIR) if os.path.isdir(NOPOLYA_DIR) else {}
    for s, d in nopolya.items():
        if s in rows and d.get("flnc"):
            rows[s]["flnc_nopolya"] = d["flnc"]
    return rows


def sample_ticks(ax, ss):
    """sample names under the bars, slanted if they are long"""
    long_names = max(len(s) for s in ss) > 10
    ax.set_xticks(range(len(ss)))
    ax.set_xticklabels([qc.sample_label(s) for s in ss], fontsize=9,
                       rotation=30 if long_names else 0, ha="right" if long_names else "center")


def present(rows):
    return [s for s in qc.SAMPLES_BY_GROUP if s in rows]


def f1_table(rows, path):
    ss = present(rows)
    head = ["Sample", "Group", "FL", "FLNC", "FLNC+polyA", "FL->final %",
            "Chimera %", "polyA %", "Status"]
    tab = [head]
    for s in ss:
        d = rows[s]
        fl, flnc, pa = d.get("fl"), d.get("flnc"), d.get("polya")
        tab.append([s, qc.GROUP_OF[s],
                    f"{fl:,}" if fl else "n/a",
                    f"{flnc:,}" if flnc else "n/a",
                    f"{pa:,}" if pa else "n/a",
                    f"{100 * pa / fl:.2f}%" if (fl and pa) else "n/a",
                    f"{100 * (fl - flnc) / fl:.3f}%" if (fl and flnc) else "n/a",
                    f"{100 * pa / flnc:.2f}%" if (flnc and pa) else "n/a",
                    "PASS" if (pa or 0) >= MIN_FLNC else "LOW"])
    qc.table_fig(tab, "isoseq refine summary (--require-polya)", path,
                 status_col=len(head) - 1, row_groups=[qc.GROUP_OF[s] for s in ss])


def _bar(rows, valfn, ylab, title, path, pct=False, fmt="{:.2f}", height=6):
    ss = [s for s in present(rows) if valfn(rows[s]) is not None]
    if not ss:
        print(f"  skipped '{title}': no data")
        return
    vals = np.array([valfn(rows[s]) for s in ss], float)
    unit = "%" if pct else ""
    fig, ax = plt.subplots(figsize=(max(10, 1.4 * len(ss)), height))
    if pct:
        # zoom in, these are usually all within a percent of each other
        lo = max(0, vals.min() - max(0.5, 0.1 * (vals.max() - vals.min())))
        hi = min(100.0, vals.max()) + (vals.max() - lo) * 0.25 + 0.01
        ax.set_ylim(lo, hi)
    else:
        lo, hi = 0, vals.max() * 1.25
        ax.set_ylim(lo, hi)
    ax.bar(range(len(ss)), vals, color=[qc.sample_color(s) for s in ss],
           edgecolor="black", linewidth=1.5)
    for i, v in enumerate(vals):
        ax.text(i, v + (hi - lo) * 0.01, (fmt + unit).format(v), ha="center",
                va="bottom", fontweight="bold", fontsize=9)
    ax.axhline(vals.mean(), color="red", ls="--", lw=2,
               label=("mean " + fmt + unit).format(vals.mean()))
    qc.group_blocks(ax, ss, hi)
    sample_ticks(ax, ss)
    ax.set_ylabel(ylab, fontweight="bold")
    ax.set_title(title, fontweight="bold")
    ax.legend(loc="lower right")
    qc.save(fig, path)


def f2_yield(rows, path):
    _bar(rows, lambda d: d["polya"] / 1e6 if d.get("polya") else None,
         "FLNC reads (millions)", "FLNC reads per sample (--require-polya)", path)


def f3_conv(rows, path):
    _bar(rows, lambda d: 100 * d["polya"] / d["fl"] if (d.get("fl") and d.get("polya")) else None,
         "FL -> final FLNC (%)", "FL to FLNC+polyA conversion", path, pct=True, height=5)


def f4_polya(rows, path):
    _bar(rows, lambda d: 100 * d["polya"] / d["flnc"] if (d.get("flnc") and d.get("polya")) else None,
         "reads with poly(A) (%)", "Poly(A) tail detection", path, pct=True, fmt="{:.3f}", height=5)


def f5_chimera(rows, path):
    _bar(rows, lambda d: 100 * (d["fl"] - d["flnc"]) / d["fl"] if (d.get("fl") and d.get("flnc")) else None,
         "chimeric reads (%)", "Chimeric read rate (worry if >10%)", path, fmt="{:.3f}", height=5)


def f6_waterfall(rows, path):
    ss = [s for s in present(rows) if rows[s].get("fl") and rows[s].get("polya")]
    if not ss:
        print("  skipped waterfall: no refine reports")
        return
    x = np.arange(len(ss))
    w = 0.27
    fl = np.array([rows[s]["fl"] for s in ss]) / 1e6
    flnc = np.array([rows[s].get("flnc", rows[s]["polya"]) for s in ss]) / 1e6
    pa = np.array([rows[s]["polya"] for s in ss]) / 1e6
    fig, ax = plt.subplots(figsize=(max(11, 1.8 * len(ss)), 6))
    ax.bar(x - w, fl, w, label="FL", color="#B0B0B0", edgecolor="black")
    ax.bar(x, flnc, w, label="FLNC (chimeras removed)", color="#9B86C2", edgecolor="black")
    ax.bar(x + w, pa, w, label="FLNC + polyA (used downstream)", color=qc.ALL_COLOR, edgecolor="black")
    top = fl.max() * 1.2
    ax.set_ylim(0, top)
    qc.group_blocks(ax, ss, top)
    sample_ticks(ax, ss)
    ax.set_ylabel("reads (millions)", fontweight="bold")
    ax.set_title("FL -> FLNC -> FLNC+polyA", fontweight="bold")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, frameon=False)
    qc.save(fig, path)


def f7_polya_table(rows, path):
    ss = [s for s in present(rows) if rows[s].get("polya") and rows[s].get("flnc_nopolya")]
    if not ss:
        print("  skipped poly(A) comparison table: no 02c output")
        return
    tab = [["Sample", "Group", "FLNC with --require-polya", "FLNC without", "Lost", "Lost %"]]
    for s in ss:
        w, wo = rows[s]["polya"], rows[s]["flnc_nopolya"]
        tab.append([s, qc.GROUP_OF[s], f"{w:,}", f"{wo:,}", f"{wo - w:,}", f"{100 * (wo - w) / wo:.3f}%"])
    qc.table_fig(tab, "Reads lost to the poly(A) requirement", path,
                 row_groups=[qc.GROUP_OF[s] for s in ss])


def f8_polya_chart(rows, path):
    ss = [s for s in present(rows) if rows[s].get("polya") and rows[s].get("flnc_nopolya")]
    if not ss:
        print("  skipped poly(A) comparison chart: no 02c output")
        return
    x = np.arange(len(ss))
    w = 0.38
    wv = np.array([rows[s]["polya"] for s in ss]) / 1e6
    wo = np.array([rows[s]["flnc_nopolya"] for s in ss]) / 1e6
    fig, ax = plt.subplots(figsize=(max(11, 1.8 * len(ss)), 6))
    top = wo.max() * 1.2
    ax.set_ylim(0, top)
    ax.bar(x - w / 2, wv, w, label="with --require-polya", color=qc.ALL_COLOR, edgecolor="black")
    ax.bar(x + w / 2, wo, w, label="without", color=qc.ALL_LIGHT, edgecolor="black")
    for i in range(len(ss)):
        diff = int(round((wo[i] - wv[i]) * 1e6))
        ax.text(i, max(wv[i], wo[i]) + top * 0.01, f"+{diff:,}", ha="center", va="bottom",
                fontsize=8, color="darkred", fontweight="bold")
    sample_ticks(ax, ss)
    ax.set_ylabel("FLNC reads (millions)", fontweight="bold")
    ax.set_title("FLNC reads with vs without --require-polya", fontweight="bold")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, frameon=False)
    qc.save(fig, path)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    qc.use_style()
    rows = build()
    if not rows:
        raise SystemExit(f"no refine reports in {REFINED_DIR} (run 02b first)")
    for s in qc.SAMPLES_BY_GROUP:
        d = rows.get(s)
        if d is None:
            print(f"  WARNING: no refine output for {s}")
            continue
        print(f"  {s:<25} {qc.GROUP_OF[s]:<10} FL={d.get('fl')} FLNC={d.get('flnc')} "
              f"polyA={d.get('polya')} no-polyA-run={d.get('flnc_nopolya')}")
    f1_table(rows, f"{OUT_DIR}/Refine_01_summary.png")
    f2_yield(rows, f"{OUT_DIR}/Refine_02_flnc_yield.png")
    f3_conv(rows, f"{OUT_DIR}/Refine_03_fl_to_flnc.png")
    f4_polya(rows, f"{OUT_DIR}/Refine_04_polya_detection.png")
    f5_chimera(rows, f"{OUT_DIR}/Refine_05_chimera_rate.png")
    f6_waterfall(rows, f"{OUT_DIR}/Refine_06_waterfall.png")
    f7_polya_table(rows, f"{OUT_DIR}/Refine_07_polya_comparison_table.png")
    f8_polya_chart(rows, f"{OUT_DIR}/Refine_08_polya_comparison_chart.png")
    print(f"figures in {OUT_DIR}/")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
MAPPING_QC.py -- pbmm2 mapping QC plots (step 05)
Talha Lone (Parvez lab)

Reads  results/qc/<BUILD_LABEL>/mapstats_<DATASET>.tsv   (from 05d_EXTRACT_MAP_STATS.sh)
       DATASET = ALL, each group, and each sample (whatever exists)
Writes results/qc/<BUILD_LABEL>/MAPPING_*.png

Run:   cd 05_MAPPING; python MAPPING_QC.py
Clustered datasets (ALL + groups) get the full figure set. Per-sample BAMs
from 05c, if present, show up in the summary table and mapping-rate plot.
"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import qc_utils as qc  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

qc.use_style()
plt.rcParams['font.family'] = 'DejaVu Sans'

LABEL = qc.CFG["BUILD_LABEL"]
QC_DIR = os.path.join('results', 'qc', LABEL)
MAP_RATE_FLOOR = 95.0   # expected minimum primary mapping rate
HEADER = '#4B2E83'


def g(d, k, default=0):
    v = d.get(k, default)
    return v if isinstance(v, (int, float)) else default


DATA = {}
for ds in qc.DATASETS + qc.SAMPLES_BY_GROUP:
    p = os.path.join(QC_DIR, f'mapstats_{ds}.tsv')
    if os.path.exists(p):
        DATA[ds] = qc.read_kv(p)
if not DATA:
    raise SystemExit(f'No mapstats_*.tsv in {QC_DIR}. Run 05d_EXTRACT_MAP_STATS.sh first.')

CLUSTERED = [d for d in qc.DATASETS if d in DATA]
INDIV = [s for s in qc.SAMPLES_BY_GROUP if s in DATA]

print('=' * 80)
print(f'PBMM2 MAPPING QC  ({LABEL})')
print('clustered :', ', '.join(CLUSTERED) or '-')
print('samples   :', ', '.join(INDIV) or '-')
print('=' * 80)


def mapq_cats(d):
    """(total, unique MAPQ60, intermediate 1-59, multi MAPQ0)"""
    total = sum(v for k, v in d.items() if k.startswith('mapq_') and isinstance(v, (int, float)))
    q60 = g(d, 'mapq_60')
    q0 = g(d, 'mapq_0')
    return total, q60, total - q60 - q0, q0


def ds_label(ds):
    return qc.dataset_label(ds)


def offsets(n, width):
    return [(i - (n - 1) / 2) * width for i in range(n)]


def out_path(name):
    return os.path.join(QC_DIR, name)


# ---------------------------------------------------------------------------
# 1. summary table
# ---------------------------------------------------------------------------
def fig_summary():
    rows = CLUSTERED + INDIV
    fig, ax = plt.subplots(figsize=(16, 0.7 * len(rows) + 2))
    ax.axis('off')
    header = ['Dataset', 'Primary', 'Primary\nMapped', 'Primary\nUnmapped',
              'Map Rate\n(%)', 'Supplementary', 'Supp %\n(multi-map)', 'Status']
    tdata = [header]
    status = []
    for ds in rows:
        d = DATA[ds]
        P = g(d, 'flagstat_primary')
        SU = g(d, 'flagstat_supplementary')
        rate = g(d, 'flagstat_map_rate')
        st = 'PASS' if rate >= MAP_RATE_FLOOR else 'CHECK'
        status.append(st)
        tdata.append([ds, f'{P:,}', f"{g(d, 'flagstat_primary_mapped'):,}",
                      f"{g(d, 'flagstat_primary_unmapped'):,}", f'{rate:.2f}%',
                      f'{SU:,}', f'{(SU / P * 100 if P else 0):.2f}%', st])
    t = ax.table(cellText=tdata, cellLoc='center', loc='center',
                 colWidths=[0.20, 0.12, 0.12, 0.12, 0.10, 0.13, 0.12, 0.09])
    t.auto_set_font_size(False)
    t.set_fontsize(10)
    t.scale(1, 2.2)
    for j in range(len(header)):
        c = t[(0, j)]
        c.set_facecolor(HEADER)
        c.set_text_props(weight='bold', color='white')
    for i, ds in enumerate(rows, start=1):
        col = qc.dataset_color(ds)
        for j in range(len(header)):
            c = t[(i, j)]
            if j == len(header) - 1:
                c.set_facecolor('#90EE90' if status[i - 1] == 'PASS' else '#F4C76E')
                c.set_text_props(weight='bold')
            else:
                c.set_facecolor('#FFFFFF')
            c.set_edgecolor(col)
            if j == 0:
                c.set_text_props(weight='bold', color=col)
    ax.set_title(f'pbmm2 Mapping Summary  ({LABEL})', fontsize=16, fontweight='bold', pad=20)
    qc.save(fig, out_path('MAPPING_01_summary.png'), dpi=200)


# ---------------------------------------------------------------------------
# 2. mapping rate per dataset
# ---------------------------------------------------------------------------
def fig_maprate():
    rows = CLUSTERED + INDIV
    rates = [g(DATA[d], 'flagstat_map_rate') for d in rows]
    mapped = [g(DATA[d], 'flagstat_primary_mapped') for d in rows]
    fig, ax = plt.subplots(figsize=(max(10, 1.5 * len(rows) + 3), 5.5))
    x = np.arange(len(rows))
    bars = ax.bar(x, rates, color=[qc.dataset_color(d) for d in rows], alpha=0.9,
                  edgecolor='black', linewidth=1.5, width=0.6)
    # keep the 95% line in view
    lo = max(50.0, min(MAP_RATE_FLOOR - 1, min(rates) - 2))
    ax.set_ylim(lo, 100.0)
    for b, r, m in zip(bars, rates, mapped):
        ax.text(b.get_x() + b.get_width() / 2., lo + (100 - lo) * 0.04,
                f'{r:.2f}%\n({m / 1e6:.2f}M)', ha='center', va='bottom',
                fontweight='bold', fontsize=8.5, color='white')
    ax.axhline(MAP_RATE_FLOOR, color='red', linestyle='--', linewidth=2,
               label=f'Min expected ({MAP_RATE_FLOOR:.0f}%)', zorder=10)
    ax.axhline(np.mean(rates), color='orange', linestyle='--', linewidth=2,
               label=f'Mean {np.mean(rates):.2f}%', zorder=10)
    if CLUSTERED and INDIV:
        ax.axvline(len(CLUSTERED) - 0.5, color='black', linestyle=':', linewidth=1.5, alpha=0.6)
    ax.set_ylabel('Primary Mapping Rate (%)', fontsize=12, fontweight='bold')
    ax.set_title(f'Primary Mapping Rate per Dataset  ({LABEL})', fontsize=14, fontweight='bold', pad=10)
    ax.set_xticks(x)
    ax.set_xticklabels([ds_label(d) for d in rows], fontsize=9)
    ax.legend(fontsize=10, loc='upper left', bbox_to_anchor=(1.01, 1.0))
    ax.grid(True, alpha=0.3, axis='y')
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, p: f'{v:.1f}%'))
    qc.save(fig, out_path('MAPPING_02_maprate.png'), dpi=200)


# ---------------------------------------------------------------------------
# 3. alignment breakdown, clustered
# ---------------------------------------------------------------------------
def fig_breakdown():
    rows = CLUSTERED
    if not rows:
        return
    pm = np.array([g(DATA[d], 'flagstat_primary_mapped') for d in rows], float) / 1e6
    su = np.array([g(DATA[d], 'flagstat_supplementary') for d in rows], float) / 1e6
    un = np.array([g(DATA[d], 'flagstat_primary_unmapped') for d in rows], float) / 1e6
    tot = pm + su + un
    x = np.arange(len(rows))
    width = 0.55
    fig, ax = plt.subplots(figsize=(max(9, 2.2 * len(rows) + 3), 7))
    ax.bar(x, pm, width, color=[qc.dataset_color(d) for d in rows], alpha=0.9,
           edgecolor='black', linewidth=1.5, zorder=3)
    ax.bar(x, su, width, bottom=pm, color='#B6A0D3', alpha=0.9, edgecolor='black',
           linewidth=1.5, zorder=3, label='Supplementary (multi-map)')
    ax.bar(x, un, width, bottom=pm + su, color='#D73027', alpha=0.9, edgecolor='black',
           linewidth=1.5, zorder=3, label='Unmapped')
    top = tot.max() * 1.25
    ax.set_ylim(0, top)
    for i, ds in enumerate(rows):
        s = DATA[ds]
        P = g(s, 'flagstat_primary')
        ax.text(i, pm[i] / 2, f"{pm[i]:.2f}M\n({g(s, 'flagstat_map_rate'):.2f}%)",
                ha='center', va='center', fontweight='bold', fontsize=9, color='white')
        supp_pct = g(s, 'flagstat_supplementary') / P * 100 if P else 0
        un_pct = g(s, 'flagstat_primary_unmapped') / P * 100 if P else 0
        ax.text(i, tot[i] + top * 0.015,
                f"supp: {g(s, 'flagstat_supplementary') / 1e3:.0f}K ({supp_pct:.2f}%)\n"
                f"unmapped: {int(g(s, 'flagstat_primary_unmapped')):,} ({un_pct:.3f}%)",
                ha='center', va='bottom', fontsize=8.5, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='black', alpha=0.85))
    ax.set_xticks(x)
    ax.set_xticklabels([ds_label(d) for d in rows], fontsize=10)
    ax.set_ylabel('Isoforms (Millions)', fontsize=12, fontweight='bold')
    ax.set_title(f'Alignment Category Breakdown - Clustered  ({LABEL})', fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    qc.save(fig, out_path('MAPPING_03_breakdown.png'), dpi=200)


# ---------------------------------------------------------------------------
# 4. MAPQ composition + table, clustered
# ---------------------------------------------------------------------------
def fig_mapq():
    rows = CLUSTERED
    if not rows:
        return
    cats = ['Unique (MAPQ 60)', 'Intermediate (1-59)', 'Multi-map (MAPQ 0)']
    fig, (axb, axt) = plt.subplots(2, 1, figsize=(11, 7 + 0.4 * len(rows)),
                                   gridspec_kw={'height_ratios': [3, 2]})
    x = np.arange(len(cats))
    width = 0.8 / len(rows)
    for off, ds in zip(offsets(len(rows), width), rows):
        total, q60, mid, q0 = mapq_cats(DATA[ds])
        pct = [v / total * 100 if total else 0 for v in (q60, mid, q0)]
        bars = axb.bar(x + off, pct, width, color=qc.dataset_color(ds), alpha=0.9,
                       edgecolor='black', linewidth=1, label=ds)
        for b, p in zip(bars, pct):
            axb.text(b.get_x() + b.get_width() / 2., p + 1, f'{p:.1f}%',
                     ha='center', va='bottom', fontsize=8, fontweight='bold')
    axb.set_ylabel('% of primary alignments', fontsize=12, fontweight='bold')
    axb.set_xticks(x)
    axb.set_xticklabels(cats, fontsize=10)
    axb.set_ylim(0, 105)
    axb.set_title(f'MAPQ Composition - Clustered  ({LABEL})', fontsize=14, fontweight='bold', pad=10)
    axb.legend(fontsize=10)
    axb.grid(True, alpha=0.3, axis='y')

    axt.axis('off')
    header = ['Dataset', 'Total', 'Unique (60)', 'Unique %', 'Intermed.', 'Multi (0)', 'Multi %']
    tdata = [header]
    for ds in rows:
        total, q60, mid, q0 = mapq_cats(DATA[ds])
        tdata.append([ds, f'{total:,}', f'{q60:,}', f'{(q60 / total * 100 if total else 0):.2f}%',
                      f'{mid:,}', f'{q0:,}', f'{(q0 / total * 100 if total else 0):.2f}%'])
    t = axt.table(cellText=tdata, cellLoc='center', loc='center',
                  colWidths=[0.13, 0.16, 0.16, 0.12, 0.15, 0.15, 0.12])
    t.auto_set_font_size(False)
    t.set_fontsize(9.5)
    t.scale(1, 2.0)
    for j in range(len(header)):
        c = t[(0, j)]
        c.set_facecolor(HEADER)
        c.set_text_props(weight='bold', color='white')
    for i, ds in enumerate(rows, start=1):
        for j in range(len(header)):
            c = t[(i, j)]
            c.set_facecolor('#FFFFFF')
            c.set_edgecolor(qc.dataset_color(ds))
            if j == 0:
                c.set_text_props(weight='bold', color=qc.dataset_color(ds))
    qc.save(fig, out_path('MAPPING_04_mapq.png'), dpi=200)


# ---------------------------------------------------------------------------
# 5. reads per chromosome, clustered
# ---------------------------------------------------------------------------
SEX_ORDER = {'X': 0, 'Y': 1, 'Z': 2, 'W': 3}


def chrom_key(name):
    """natural order: 1..N, then X Y Z W, then M/MT, then anything else, 'other' last"""
    if name == 'other':
        return (4, 0, name)
    core = re.sub(r'^chr', '', name, flags=re.I)
    if core.isdigit():
        return (0, int(core), name)
    if core.upper() in SEX_ORDER:
        return (1, SEX_ORDER[core.upper()], name)
    if core.upper() in ('M', 'MT'):
        return (2, 0, name)
    return (3, 0, name)


def fig_chrom():
    rows = CLUSTERED
    if not rows:
        return
    chroms = set()
    for ds in rows:
        chroms.update(k[4:] for k in DATA[ds] if k.startswith('chr_'))
    chroms = sorted(chroms, key=chrom_key)
    if not chroms:
        return
    fig, ax = plt.subplots(figsize=(min(26, max(12, 0.35 * len(chroms) * len(rows) + 4)), 7))
    x = np.arange(len(chroms))
    width = 0.8 / len(rows)
    for off, ds in zip(offsets(len(rows), width), rows):
        vals = np.array([g(DATA[ds], f'chr_{c}') for c in chroms], float) / 1e3
        ax.bar(x + off, vals, width, color=qc.dataset_color(ds), alpha=0.88,
               edgecolor='black', linewidth=0.5, label=ds_label(ds).replace('\n', ' '))
    # divider before the non-nuclear / unplaced bins
    n_nuc = sum(1 for c in chroms if chrom_key(c)[0] < 2)
    if 0 < n_nuc < len(chroms):
        ax.axvline(n_nuc - 0.5, color='black', linestyle='--', linewidth=1.5, alpha=0.6)
    ax.set_ylabel('Mapped Isoforms (Thousands)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Chromosome', fontsize=12, fontweight='bold')
    ax.set_title(f'Chromosome Distribution of Mapped Isoforms  ({LABEL})',
                 fontsize=14, fontweight='bold', pad=10)
    ax.set_xticks(x)
    ax.set_xticklabels(chroms, fontsize=9 if len(chroms) <= 30 else 7, rotation=45, ha='right')
    ax.set_xlim(-0.6, len(chroms) - 0.4)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    qc.save(fig, out_path('MAPPING_05_chromosome.png'), dpi=200)


# ---------------------------------------------------------------------------
# 6. rate / unmapped / supplementary, clustered
# ---------------------------------------------------------------------------
def fig_triplot():
    rows = CLUSTERED
    if not rows:
        return
    colors = [qc.dataset_color(d) for d in rows]
    prim = [g(DATA[d], 'flagstat_primary') for d in rows]
    rates = [g(DATA[d], 'flagstat_map_rate') for d in rows]
    unmap = [g(DATA[d], 'flagstat_primary_unmapped') for d in rows]
    supp = [g(DATA[d], 'flagstat_supplementary') for d in rows]
    unmap_pct = [u / p * 100 if p else 0 for u, p in zip(unmap, prim)]
    supp_pct = [s / p * 100 if p else 0 for s, p in zip(supp, prim)]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(1, 3, figsize=(max(15, 5 * len(rows)), 5.2))

    ax = axes[0]
    bars = ax.bar(x, rates, color=colors, alpha=0.9, edgecolor='black', linewidth=1.5, width=0.62)
    lo = max(80.0, min(99.0, min(rates) - 0.10))
    ax.set_ylim(lo, 100.0)
    for b, r in zip(bars, rates):
        ax.text(b.get_x() + b.get_width() / 2, lo + (100 - lo) * 0.05, f'{r:.2f}%',
                ha='center', va='bottom', fontweight='bold', color='white', fontsize=11)
    mr = float(np.mean(rates))
    ax.axhline(mr, color='orange', linestyle='--', linewidth=2, label=f'mean {mr:.2f}%')
    ax.set_title('Primary Mapping Rate', fontsize=13, fontweight='bold', color=HEADER)
    ax.set_ylabel('Mapping Rate (%)', fontsize=11, fontweight='bold')
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, p: f'{v:.2f}%'))
    ax.legend(fontsize=9, loc='upper right')

    ax = axes[1]
    top = max(unmap) or 1
    bars = ax.bar(x, unmap, color=colors, alpha=0.9, edgecolor='black', linewidth=1.5, width=0.62)
    ax.set_ylim(0, top * 1.20)
    for b, u, pc in zip(bars, unmap, unmap_pct):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + top * 0.02,
                f'{u:,}\n({pc:.3f}%)', ha='center', va='bottom', fontweight='bold', fontsize=10)
    ax.set_title('Unmapped Isoforms', fontsize=13, fontweight='bold', color=HEADER)
    ax.set_ylabel('Count', fontsize=11, fontweight='bold')

    ax = axes[2]
    supp_k = [s / 1e3 for s in supp]
    top = max(supp_k) or 1
    bars = ax.bar(x, supp_k, color=colors, alpha=0.9, edgecolor='black', linewidth=1.5, width=0.62)
    ax.set_ylim(0, top * 1.20)
    for b, s, pc in zip(bars, supp_k, supp_pct):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + top * 0.02,
                f'{s:.1f}K\n({pc:.2f}%)', ha='center', va='bottom', fontweight='bold', fontsize=10)
    ax.set_title('Supplementary Alignments', fontsize=13, fontweight='bold', color=HEADER)
    ax.set_ylabel('Count (Thousands)', fontsize=11, fontweight='bold')

    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(rows, fontsize=11)
        ax.grid(True, alpha=0.3, axis='y')
    fig.suptitle(f'Mapping QC Summary - Clustered  ({LABEL})', fontsize=15, fontweight='bold')
    qc.save(fig, out_path('MAPPING_06_triplot.png'), dpi=200)


fig_summary()
fig_maprate()
fig_breakdown()
fig_mapq()
fig_chrom()
fig_triplot()

print('\n' + '=' * 80)
print('MAPPING QC SUMMARY')
print('=' * 80)
for ds in CLUSTERED + INDIV:
    d = DATA[ds]
    t, q60, _, _ = mapq_cats(d)
    flag = '' if g(d, 'flagstat_map_rate') >= MAP_RATE_FLOOR else '   <-- below expected'
    print(f"  {ds:<24} rate={g(d, 'flagstat_map_rate'):6.2f}%  "
          f"unmapped={int(g(d, 'flagstat_primary_unmapped')):>9,}  "
          f"MAPQ60={(q60 / t * 100 if t else 0):5.1f}%{flag}")
print('\nDone.')

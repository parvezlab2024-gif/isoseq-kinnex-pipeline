#!/usr/bin/env python3
"""
COLLAPSE_QC.py -- isoseq collapse QC plots (step 06)
Talha Lone (Parvez lab)

Reads  results/qc/<BUILD_LABEL>/collapse_stats_<GROUP>.tsv  (from 06b_EXTRACT_COLLAPSE_STATS.sh)
       GROUP = ALL plus any group that was collapsed on its own
Writes results/qc/<BUILD_LABEL>/COLLAPSE_*.png

Run:   cd 06_COLLAPSE; python COLLAPSE_QC.py
Per-sample numbers come from the ALL collapse (flnc_count.txt columns).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import qc_utils as qc  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

qc.use_style()
plt.rcParams['font.family'] = 'DejaVu Sans'

LABEL = qc.CFG["BUILD_LABEL"]
QC_DIR = os.path.join('results', 'qc', LABEL)
HEADER = '#4B2E83'

SPAN_KEYS = ['span_le200', 'span_201_500', 'span_501_1000', 'span_1001_2000',
             'span_2001_5000', 'span_5001_10k', 'span_10k_15k', 'span_15kplus']
SPAN_LABELS = ['<=200', '201-500', '501-1000', '1001-2000',
               '2001-5000', '5001-10k', '10k-15k', '>15k']
IPL_KEYS = ['ipl_1', 'ipl_2', 'ipl_3', 'ipl_4_5', 'ipl_6_10', 'ipl_11_20', 'ipl_20plus']
IPL_LABELS = ['1', '2', '3', '4-5', '6-10', '11-20', '>20']
SUPP_KEYS = ['supp_1', 'supp_2_5', 'supp_6_10', 'supp_11_50', 'supp_51_100', 'supp_100plus']
SUPP_LABELS = ['1', '2-5', '6-10', '11-50', '51-100', '100+']


def g(d, k, default=0):
    v = d.get(k, default)
    return v if isinstance(v, (int, float)) else default


DATA = {}
for ds in qc.DATASETS:
    p = os.path.join(QC_DIR, f'collapse_stats_{ds}.tsv')
    if os.path.exists(p):
        DATA[ds] = qc.read_kv(p)
if not DATA:
    raise SystemExit(f'No collapse_stats_*.tsv in {QC_DIR}. Run 06b_EXTRACT_COLLAPSE_STATS.sh first.')
ROWS = [d for d in qc.DATASETS if d in DATA]

print('=' * 80)
print(f'ISOSEQ COLLAPSE QC  ({LABEL})')
print('datasets:', ', '.join(ROWS))
print('=' * 80)


def singleton_pct(d):
    post = g(d, 'post_isoforms')
    return g(d, 'supp_1') / post * 100 if post else 0


def mono_pct(d):
    tot = g(d, 'exon_mono') + g(d, 'exon_multi')
    return g(d, 'exon_mono') / tot * 100 if tot else 0


def offsets(n, width):
    return [(i - (n - 1) / 2) * width for i in range(n)]


def out_path(name):
    return os.path.join(QC_DIR, name)


# ---------------------------------------------------------------------------
# 1. summary table
# ---------------------------------------------------------------------------
def fig_summary():
    fig, ax = plt.subplots(figsize=(18, 0.9 * len(ROWS) + 2))
    ax.axis('off')
    header = ['Dataset', 'Pre-Collapse\n(mapped)', 'Post-Collapse\n(isoforms)',
              'Redundancy\nRemoved', 'Reduction\n(%)', 'Gene\nLoci',
              'Mean Iso\n/ Locus', 'Mono-exon\n(%)', 'Singleton\n(%)']
    tdata = [header]
    for ds in ROWS:
        d = DATA[ds]
        pre = g(d, 'pre_isoforms')
        post = g(d, 'post_isoforms')
        loci = g(d, 'loci')
        red = pre - post
        tdata.append([ds, f'{pre:,}', f'{post:,}', f'{red:,}',
                      f'{(red / pre * 100 if pre else 0):.1f}%', f'{loci:,}',
                      f'{(post / loci if loci else 0):.1f}', f'{mono_pct(d):.1f}%',
                      f'{singleton_pct(d):.1f}%'])
    t = ax.table(cellText=tdata, cellLoc='center', loc='center',
                 colWidths=[0.10, 0.13, 0.13, 0.12, 0.09, 0.10, 0.09, 0.10, 0.09])
    t.auto_set_font_size(False)
    t.set_fontsize(10.5)
    t.scale(1, 2.6)
    for j in range(len(header)):
        c = t[(0, j)]
        c.set_facecolor(HEADER)
        c.set_text_props(weight='bold', color='white')
    for i, ds in enumerate(ROWS, start=1):
        col = qc.dataset_color(ds)
        for j in range(len(header)):
            c = t[(i, j)]
            c.set_facecolor('#FFFFFF')
            c.set_edgecolor(col)
            if j == 0:
                c.set_text_props(weight='bold', color=col)
    ax.set_title(f'IsoSeq Collapse Summary  ({LABEL})', fontsize=16, fontweight='bold', pad=20)
    qc.save(fig, out_path('COLLAPSE_01_summary.png'), dpi=200)


# ---------------------------------------------------------------------------
# 2. pre vs post collapse
# ---------------------------------------------------------------------------
def fig_pre_post():
    pre = [g(DATA[d], 'pre_isoforms') / 1e6 for d in ROWS]
    post = [g(DATA[d], 'post_isoforms') / 1e6 for d in ROWS]
    x = np.arange(len(ROWS))
    width = 0.36
    fig, ax = plt.subplots(figsize=(max(9, 2.5 * len(ROWS) + 3), 6))
    b1 = ax.bar(x - width / 2, pre, width, label='Pre-collapse (clustered, mapped)',
                color=[qc.group_light(d) for d in ROWS], alpha=0.9, edgecolor='black', linewidth=1.8)
    b2 = ax.bar(x + width / 2, post, width, label='Post-collapse (unique isoforms)',
                color=[qc.group_color(d) for d in ROWS], alpha=0.9, edgecolor='black', linewidth=1.8)
    ytop = (max(pre + post) or 1) * 1.22
    ax.set_ylim(0, ytop)
    for b, v in zip(b1, pre):
        ax.text(b.get_x() + b.get_width() / 2., b.get_height() + ytop * 0.015,
                f'{v:.2f}M', ha='center', va='bottom', fontweight='bold', fontsize=9)
    for b, v, ds in zip(b2, post, ROWS):
        p0 = g(DATA[ds], 'pre_isoforms')
        q0 = g(DATA[ds], 'post_isoforms')
        rp = (p0 - q0) / p0 * 100 if p0 else 0
        cx = b.get_x() + b.get_width() / 2.
        ax.text(cx, b.get_height() + ytop * 0.055, f'{v:.2f}M',
                ha='center', va='bottom', fontweight='bold', fontsize=9)
        ax.text(cx, b.get_height() + ytop * 0.013, f'-{rp:.0f}%',
                ha='center', va='bottom', fontweight='bold', fontsize=8.5, color='#B23B3B')
    ax.set_ylabel('Isoforms (Millions)', fontsize=12, fontweight='bold')
    ax.set_title(f'Pre- vs Post-Collapse Isoforms  ({LABEL})', fontsize=14, fontweight='bold', pad=10)
    ax.set_xticks(x)
    ax.set_xticklabels([qc.dataset_label(d) for d in ROWS], fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    qc.save(fig, out_path('COLLAPSE_02_pre_vs_post.png'), dpi=200)


# ---------------------------------------------------------------------------
# 3. FLNC support per isoform + singletons
# ---------------------------------------------------------------------------
def fig_flnc_support():
    if not any(k in DATA[d] for d in ROWS for k in SUPP_KEYS):
        print('  (support figure skipped: no flnc_count stats)')
        return
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={'width_ratios': [2.3, 1]})
    x = np.arange(len(SUPP_KEYS))
    width = 0.8 / len(ROWS)
    sing, cnt = [], []
    for off, ds in zip(offsets(len(ROWS), width), ROWS):
        d = DATA[ds]
        tot = sum(g(d, k) for k in SUPP_KEYS) or 1
        axl.bar(x + off, [g(d, k) / tot * 100 for k in SUPP_KEYS], width,
                color=qc.group_color(ds), alpha=0.9, edgecolor='black', linewidth=1,
                label=f'{ds} ({qc.n_samples(ds)} samples)')
        sing.append(g(d, 'supp_1') / tot * 100)
        cnt.append(g(d, 'supp_1'))
    axl.set_ylim(0, 100)
    axl.set_ylabel('% of Collapsed Isoforms', fontsize=12, fontweight='bold')
    axl.set_xlabel('FLNC Read Support per Isoform', fontsize=12, fontweight='bold')
    axl.set_title('Support Distribution', fontsize=13, fontweight='bold', pad=10)
    axl.set_xticks(x)
    axl.set_xticklabels(SUPP_LABELS, fontsize=10)
    axl.legend(fontsize=10)
    axl.grid(True, alpha=0.3, axis='y')

    bars = axr.bar(range(len(ROWS)), sing, color=[qc.group_color(d) for d in ROWS],
                   alpha=0.9, edgecolor='black', linewidth=1.5, width=0.6)
    lo = max(0, min(sing) - 1.5)
    hi = max(sing) + 2.5
    axr.set_ylim(lo, hi)
    for b, s, c in zip(bars, sing, cnt):
        axr.text(b.get_x() + b.get_width() / 2., s + (hi - lo) * 0.02,
                 f'{s:.1f}%\n({c / 1e6:.2f}M)', ha='center', va='bottom', fontsize=9, fontweight='bold',
                 bbox=dict(boxstyle='round', fc='white', ec='gray', alpha=0.9))
    axr.set_ylabel('% of Singleton Isoforms', fontsize=12, fontweight='bold')
    axr.set_title('Singletons (support = 1)', fontsize=13, fontweight='bold', pad=10)
    axr.set_xticks(range(len(ROWS)))
    axr.set_xticklabels(ROWS, fontsize=11)
    axr.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, p: f'{v:.1f}%'))
    axr.grid(True, alpha=0.3, axis='y')
    fig.suptitle(f'FLNC Read Support per Collapsed Isoform  ({LABEL})', fontsize=15, fontweight='bold')
    qc.save(fig, out_path('COLLAPSE_03_flnc_support.png'), dpi=200)


# ---------------------------------------------------------------------------
# 4. isoforms per gene locus
# ---------------------------------------------------------------------------
def fig_isoforms_per_locus():
    x = np.arange(len(IPL_KEYS))
    width = 0.8 / len(ROWS)
    fig, ax = plt.subplots(figsize=(12, 6))
    for off, ds in zip(offsets(len(ROWS), width), ROWS):
        ax.bar(x + off, [g(DATA[ds], k) / 1e3 for k in IPL_KEYS], width,
               color=qc.group_color(ds), alpha=0.9, edgecolor='black', linewidth=1, label=ds)
    ax.set_ylabel('Gene Loci (Thousands)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Isoforms per Gene Locus', fontsize=12, fontweight='bold')
    ax.set_title(f'Isoform Diversity per Gene Locus  ({LABEL})', fontsize=14, fontweight='bold', pad=10)
    ax.set_xticks(x)
    ax.set_xticklabels(IPL_LABELS, fontsize=10)
    txt = '\n'.join(f"{ds}: {g(DATA[ds], 'post_isoforms') / g(DATA[ds], 'loci'):.1f}"
                    for ds in ROWS if g(DATA[ds], 'loci'))
    ax.text(0.98, 0.95, f'mean iso/locus\n{txt}', transform=ax.transAxes, ha='right', va='top',
            fontsize=9, bbox=dict(boxstyle='round', fc='white', ec='gray', alpha=0.9))
    ax.legend(fontsize=10, loc='upper center')
    ax.grid(True, alpha=0.3, axis='y')
    qc.save(fig, out_path('COLLAPSE_04_isoforms_per_locus.png'), dpi=200)


# ---------------------------------------------------------------------------
# 5. per-sample isoforms / reads (from the ALL collapse)
# ---------------------------------------------------------------------------
def fig_cross_sample():
    if 'ALL' not in DATA:
        return
    d = DATA['ALL']
    samples = [s for s in qc.SAMPLES_BY_GROUP if f'ps_{s}_iso' in d]
    if not samples:
        print('  (per-sample figure skipped: no ps_<sample>_iso keys in ALL; '
              'check that flnc_count.txt column names match samples.tsv)')
        return
    colors = [qc.sample_color(s) for s in samples]
    iso = [g(d, f'ps_{s}_iso') / 1e3 for s in samples]
    flnc = [g(d, f'ps_{s}_flnc') / 1e6 for s in samples]
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(max(14, 1.3 * len(samples) * 2), 5.8))
    for ax, vals, ylab, unit, title in (
            (axa, iso, 'Isoforms Detected (Thousands)', 'K', 'Isoforms Detected per Sample'),
            (axb, flnc, 'FLNC Assigned (Millions)', 'M', 'FLNC Reads Assigned per Sample')):
        fmt = (lambda v: f'{v:.0f}K') if unit == 'K' else (lambda v: f'{v:.2f}M')
        bars = ax.bar(range(len(samples)), vals, color=colors, alpha=0.9,
                      edgecolor='black', linewidth=1.5, width=0.7)
        ytop = (max(vals) or 1) * 1.25
        ax.set_ylim(0, ytop)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2., b.get_height() + ytop * 0.02,
                    fmt(v), ha='center', va='bottom', fontweight='bold', fontsize=9)
        ax.axhline(np.mean(vals), color='gray', linestyle='--', linewidth=1.8,
                   label=f'mean {fmt(np.mean(vals))}')
        qc.group_blocks(ax, samples, ytop)
        ax.set_ylabel(ylab, fontsize=11, fontweight='bold')
        ax.set_title(title, fontsize=13, fontweight='bold')
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(samples, fontsize=8, rotation=30, ha='right')
        ax.legend(fontsize=9, loc='lower right')
        ax.grid(True, alpha=0.3, axis='y')
    fig.suptitle(f'Cross-Sample Consistency (from pooled ALL collapse)  ({LABEL})',
                 fontsize=15, fontweight='bold')
    qc.save(fig, out_path('COLLAPSE_05_cross_sample.png'), dpi=200)


# ---------------------------------------------------------------------------
# 6. ALL vs each group
# ---------------------------------------------------------------------------
def fig_comparison():
    n = len(ROWS)
    fig, axes = plt.subplots(1, 3, figsize=(max(15, 1.6 * n * 3), 5))

    def panel(ax, vals, title, ylab, fmt):
        bars = ax.bar(range(n), vals, color=[qc.group_color(d) for d in ROWS],
                      alpha=0.9, edgecolor='black', linewidth=1.8, width=0.6)
        ytop = (max(vals) or 1) * 1.22
        ax.set_ylim(0, ytop)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2., b.get_height() + ytop * 0.02,
                    fmt(v), ha='center', va='bottom', fontweight='bold', fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_ylabel(ylab, fontsize=11, fontweight='bold')
        ax.set_xticks(range(n))
        ax.set_xticklabels(ROWS, fontsize=11)
        ax.grid(True, alpha=0.3, axis='y')

    panel(axes[0], [g(DATA[d], 'post_isoforms') / 1e6 for d in ROWS],
          'Collapsed Isoforms', 'Isoforms (Millions)', lambda v: f'{v:.2f}M')
    panel(axes[1], [g(DATA[d], 'loci') / 1e3 for d in ROWS],
          'Unique Gene Loci', 'Gene Loci (Thousands)', lambda v: f'{v:.1f}K')
    panel(axes[2], [singleton_pct(DATA[d]) for d in ROWS],
          'Singleton Isoforms (%)', 'Singleton %', lambda v: f'{v:.1f}%')
    fig.suptitle(f"Collapse Comparison: {' / '.join(ROWS)}  ({LABEL})", fontsize=15, fontweight='bold')
    qc.save(fig, out_path('COLLAPSE_06_comparison.png'), dpi=200)


# ---------------------------------------------------------------------------
# 7. genomic span
# ---------------------------------------------------------------------------
def fig_span():
    fig, (axl, axm) = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={'width_ratios': [2.5, 1]})
    x = np.arange(len(SPAN_KEYS))
    width = 0.8 / len(ROWS)
    for off, ds in zip(offsets(len(ROWS), width), ROWS):
        d = DATA[ds]
        post = g(d, 'post_isoforms')
        axl.bar(x + off, [g(d, k) / post * 100 if post else 0 for k in SPAN_KEYS], width,
                color=qc.group_color(ds), alpha=0.9, edgecolor='black', linewidth=1.2, label=ds)
    axl.set_ylabel('% of Collapsed Isoforms', fontsize=12, fontweight='bold')
    axl.set_xlabel('Genomic Span (bp, includes introns)', fontsize=12, fontweight='bold')
    axl.set_title('Transcript Genomic Span Distribution', fontsize=13, fontweight='bold', pad=10)
    axl.set_xticks(x)
    axl.set_xticklabels(SPAN_LABELS, fontsize=9, rotation=20, ha='right')
    axl.legend(fontsize=10)
    axl.grid(True, alpha=0.3, axis='y')

    means = [g(DATA[d], 'span_mean') / 1e3 for d in ROWS]
    bars = axm.bar(range(len(ROWS)), means, color=[qc.group_color(d) for d in ROWS],
                   alpha=0.9, edgecolor='black', linewidth=1.8, width=0.55)
    ytop = (max(means) or 1) * 1.22
    axm.set_ylim(0, ytop)
    for b, v in zip(bars, means):
        axm.text(b.get_x() + b.get_width() / 2., b.get_height() + ytop * 0.02,
                 f'{v:.1f} kb', ha='center', va='bottom', fontweight='bold', fontsize=10)
    axm.set_xticks(range(len(ROWS)))
    axm.set_xticklabels(ROWS, fontsize=11)
    axm.set_ylabel('Mean Genomic Span (kb)', fontsize=12, fontweight='bold')
    axm.set_title('Mean Span', fontsize=13, fontweight='bold', pad=10)
    axm.grid(True, alpha=0.3, axis='y')
    fig.suptitle(f'Transcript Genomic Span (Post-Collapse)  ({LABEL})', fontsize=15, fontweight='bold')
    qc.save(fig, out_path('COLLAPSE_07_span.png'), dpi=200)


fig_summary()
fig_pre_post()
fig_flnc_support()
fig_isoforms_per_locus()
fig_cross_sample()
fig_comparison()
fig_span()

print('\n' + '=' * 80)
print('COLLAPSE QC SUMMARY')
print('=' * 80)
for ds in ROWS:
    d = DATA[ds]
    pre = g(d, 'pre_isoforms')
    post = g(d, 'post_isoforms')
    loci = g(d, 'loci')
    print(f"  {ds:<10} pre={pre:>10,}  post={post:>10,}  "
          f"removed={(pre - post) / pre * 100 if pre else 0:4.1f}%  loci={loci:>7,}  "
          f"iso/locus={(post / loci if loci else 0):4.1f}  singleton={singleton_pct(d):4.1f}%  "
          f"mean_span={g(d, 'span_mean') / 1e3:.1f}kb")
print('\nDone.')

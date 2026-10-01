#!/usr/bin/env python3
"""
CLASSIFY_QC.py -- pigeon classify / filter QC plots (step 07)
Talha Lone (Parvez lab)

Reads  results/qc/<BUILD_LABEL>/classify_stats_<GROUP>.tsv   (from 07e_EXTRACT_CLASSIFY_STATS.sh)
       results/<BUILD_LABEL>/<GROUP>/clustered_<GROUP>.saturation.txt  (from 07d, optional)
Writes results/qc/<BUILD_LABEL>/CLASSIFY_*.png

Run:   cd 07_CLASSIFY; python CLASSIFY_QC.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import qc_utils as qc  # noqa: E402
import matplotlib  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

qc.use_style()
plt.rcParams['font.family'] = 'DejaVu Sans'

LABEL = qc.CFG["BUILD_LABEL"]
QC_DIR = os.path.join('results', 'qc', LABEL)
HEADER = '#4B2E83'

# structural categories. GenicIntron is tiny, so it is folded into Other.
CATS = ['FSM', 'ISM', 'NIC', 'NNC', 'Antisense', 'GenicGenomic', 'Intergenic', 'Other']
CAT_SHORT = {'FSM': 'FSM', 'ISM': 'ISM', 'NIC': 'NIC', 'NNC': 'NNC', 'Antisense': 'AS',
             'GenicGenomic': 'Genic', 'Intergenic': 'IG', 'Other': 'Other'}
CAT_COL = {'FSM': '#27AE60', 'ISM': '#2980B9', 'NIC': '#8E44AD', 'NNC': '#E67E22',
           'Antisense': '#E74C3C', 'GenicGenomic': '#F39C12', 'Intergenic': '#95A5A6', 'Other': '#BDC3C7'}

FILT_ORDER = ['post_passed', 'reason_intrapriming', 'reason_lowcov', 'reason_rtswitch']
FILT_LABEL = {'post_passed': 'Passed', 'reason_intrapriming': 'Intrapriming',
              'reason_lowcov': 'Low cov/NC', 'reason_rtswitch': 'RT switching'}


def g(d, k, default=0):
    v = d.get(k, default)
    return v if isinstance(v, (int, float)) else default


DATA = {}
for ds in qc.DATASETS:
    p = os.path.join(QC_DIR, f'classify_stats_{ds}.tsv')
    if os.path.exists(p):
        DATA[ds] = qc.read_kv(p)
if not DATA:
    raise SystemExit(f'No classify_stats_*.tsv in {QC_DIR}. Run 07e_EXTRACT_CLASSIFY_STATS.sh first.')
ROWS = [d for d in qc.DATASETS if d in DATA]

print('=' * 80)
print(f'PIGEON CLASSIFY & FILTER QC  ({LABEL})')
print('datasets:', ', '.join(ROWS))
print('=' * 80)


def ramp(ds, n=4):
    """dark -> light shades of the dataset colour (pie slices, stacked bars)"""
    if ds == 'ALL':
        cmap = matplotlib.colormaps['Purples'] if hasattr(matplotlib, 'colormaps') \
            else plt.get_cmap('Purples')
    else:
        name = qc.GROUP_CMAPS[qc.GROUPS.index(ds) % len(qc.GROUP_CMAPS)]
        cmap = matplotlib.colormaps[name] if hasattr(matplotlib, 'colormaps') else plt.get_cmap(name)
    return [matplotlib.colors.to_hex(cmap(x)) for x in np.linspace(0.9, 0.3, n)]


def cat_val(d, prefix, cat):
    v = g(d, f'{prefix}_{cat}')
    if cat == 'Other':
        v += g(d, f'{prefix}_GenicIntron')
    return v


def pct(n, tot):
    return f'{n:,} ({n / tot * 100:.1f}%)' if tot else f'{n:,}'


def out_path(name):
    return os.path.join(QC_DIR, name)


# ---------------------------------------------------------------------------
# 1. summary table
# ---------------------------------------------------------------------------
def fig_summary():
    fig, ax = plt.subplots(figsize=(4 + 3.2 * len(ROWS), 9))
    ax.axis('off')
    body = [['Metric'] + [f'{d} ({qc.n_samples(d)} smp)' for d in ROWS],
            ['-- CLASSIFY (pre-filter) --'] + [''] * len(ROWS)]
    body.append(['Input isoforms'] + [f"{g(DATA[d], 'pre_input'):,}" for d in ROWS])
    for c in ['FSM', 'ISM', 'NIC', 'NNC', 'Intergenic']:
        body.append(['FSM (full-splice)' if c == 'FSM' else CAT_SHORT[c]] +
                    [pct(cat_val(DATA[d], 'pre', c), g(DATA[d], 'pre_input')) for d in ROWS])
    body.append(['Known genes (pre)'] + [f"{g(DATA[d], 'pre_genes_known'):,}" for d in ROWS])
    body.append(['-- FILTER --'] + [''] * len(ROWS))
    body.append(['Passed filter'] +
                [pct(g(DATA[d], 'post_passed'), g(DATA[d], 'post_input')) for d in ROWS])
    for rk, rl in [('reason_intrapriming', 'Removed: Intrapriming'),
                   ('reason_lowcov', 'Removed: Low cov/NC'),
                   ('reason_rtswitch', 'Removed: RT switch')]:
        body.append([rl] + [pct(g(DATA[d], rk), g(DATA[d], 'post_removed')) for d in ROWS])
    genes_post = {d: g(DATA[d], 'post_genes_known') + g(DATA[d], 'post_genes_novel') for d in ROWS}
    body.append(['Known genes (post)'] +
                [pct(g(DATA[d], 'post_genes_known'), genes_post[d]) for d in ROWS])
    body.append(['Novel loci (post)'] +
                [pct(g(DATA[d], 'post_genes_novel'), genes_post[d]) for d in ROWS])

    sec = {i for i, r in enumerate(body) if r[0].startswith('--')}
    t = ax.table(cellText=body, cellLoc='center', loc='center', bbox=[0, 0, 1, 1])
    t.auto_set_font_size(False)
    t.set_fontsize(9)
    for (r, c), cell in t.get_celld().items():
        cell.set_edgecolor('#CCCCCC')
        cell.set_linewidth(0.5)
        if r == 0:
            cell.set_facecolor(HEADER)
            cell.set_text_props(color='white', fontweight='bold')
        elif r in sec:
            cell.set_facecolor('#EDE7F6')
            cell.set_text_props(fontweight='bold', color=HEADER)
        elif r % 2 == 0:
            cell.set_facecolor('#F5F5F5')
        else:
            cell.set_facecolor('white')
        if c == 0:
            cell.set_text_props(ha='left')
    ax.set_title(f'Pigeon Classify & Filter - QC Summary  ({LABEL})',
                 fontsize=14, fontweight='bold', pad=14, color=HEADER)
    qc.save(fig, out_path('CLASSIFY_01_summary.png'), dpi=200)


# ---------------------------------------------------------------------------
# 2. structural categories, pre vs post filter
# ---------------------------------------------------------------------------
def fig_categories():
    fig, axes = plt.subplots(1, len(ROWS), figsize=(5.3 * len(ROWS), 6.2), squeeze=False)
    axes = axes[0]
    fig.suptitle(f'Isoform Structural Categories: Pre- vs Post-Filter  ({LABEL})',
                 fontsize=15, fontweight='bold', y=1.02)
    x = np.arange(len(CATS))
    w = 0.38
    cc = [CAT_COL[c] for c in CATS]
    for ax, ds in zip(axes, ROWS):
        d = DATA[ds]
        pv = [cat_val(d, 'pre', c) for c in CATS]
        qv = [cat_val(d, 'post', c) for c in CATS]
        pp = [v / (sum(pv) or 1) * 100 for v in pv]
        qp = [v / (sum(qv) or 1) * 100 for v in qv]
        ax.bar(x - w / 2, pp, w, color=cc, alpha=0.35, edgecolor='black', linewidth=0.8)
        ax.bar(x + w / 2, qp, w, color=cc, edgecolor='black', linewidth=0.8)
        for xi, q2 in enumerate(qp):
            if q2 > 1.5:
                ax.text(xi + w / 2, q2 + 0.6, f'{q2:.0f}', ha='center', va='bottom',
                        fontsize=7, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels([CAT_SHORT[c] for c in CATS], fontsize=8.5, rotation=30, ha='right')
        ax.set_ylim(0, (max(pp + qp) or 1) * 1.2)
        ax.set_title(f"{ds}\n{int(g(d, 'pre_input')):,} -> {int(g(d, 'post_passed')):,}",
                     fontsize=11, fontweight='bold')
        ax.set_ylabel('% of isoforms', fontsize=10)
        kept = g(d, 'post_passed') / g(d, 'pre_input') * 100 if g(d, 'pre_input') else 0
        ax.text(0.97, 0.97, f'Retained\n{kept:.1f}%', transform=ax.transAxes, ha='right', va='top',
                fontsize=9, fontweight='bold', color='#1B5E20',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='#E8F5E9', edgecolor='#27AE60',
                          linewidth=1.5, alpha=0.9))
    handles = [mpatches.Patch(color=CAT_COL[c], label=CAT_SHORT[c]) for c in CATS]
    handles += [mpatches.Patch(facecolor='grey', alpha=0.35, edgecolor='black', label='Pre-filter'),
                mpatches.Patch(facecolor='grey', edgecolor='black', label='Post-filter')]
    fig.legend(handles=handles, loc='lower center', ncol=10, fontsize=8.5,
               bbox_to_anchor=(0.5, -0.07), frameon=True)
    qc.save(fig, out_path('CLASSIFY_02_categories.png'), dpi=200)


# ---------------------------------------------------------------------------
# 3. filter removal reasons
# ---------------------------------------------------------------------------
def fig_filter_reasons():
    fig, axes = plt.subplots(1, len(ROWS), figsize=(5.3 * len(ROWS), 5.6), squeeze=False)
    axes = axes[0]
    fig.suptitle(f'Pigeon Filter: Removal Reasons  ({LABEL})', fontsize=15, fontweight='bold', y=1.06)
    for ax, ds in zip(axes, ROWS):
        d = DATA[ds]
        sizes = [g(d, k) for k in FILT_ORDER]
        tot = sum(sizes) or 1
        ax.pie(sizes, colors=ramp(ds), startangle=90,
               wedgeprops=dict(edgecolor='white', linewidth=2), radius=1.0)
        ax.set_title(ds, fontsize=12, fontweight='bold', pad=14)
        ax.text(0, -1.28, f"Input: {int(g(d, 'post_input')):,}", ha='center', va='top',
                fontsize=9.5, fontweight='bold')
        leg = [f'{FILT_LABEL[k]}: {int(g(d, k)):,} ({g(d, k) / tot * 100:.1f}%)' for k in FILT_ORDER]
        ax.legend(leg, loc='upper center', bbox_to_anchor=(0.5, -0.16), fontsize=8.5, frameon=True)
    qc.save(fig, out_path('CLASSIFY_03_filter_reasons.png'), dpi=200)


# ---------------------------------------------------------------------------
# 4. known vs novel gene loci
# ---------------------------------------------------------------------------
def fig_known_novel():
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(max(12, 2.4 * len(ROWS) + 6), 6))
    x = np.arange(len(ROWS))
    for ax, stage, title in ((axl, 'pre', 'Pre-Filter'), (axr, 'post', 'Post-Filter')):
        kn = [g(DATA[d], f'{stage}_genes_known') / 1e3 for d in ROWS]
        nv = [g(DATA[d], f'{stage}_genes_novel') / 1e3 for d in ROWS]
        ax.bar(x, kn, 0.6, color=[ramp(d)[0] for d in ROWS], edgecolor='black', linewidth=1)
        ax.bar(x, nv, 0.6, bottom=kn, color=[ramp(d)[2] for d in ROWS], edgecolor='black', linewidth=1)
        top = max(k + n for k, n in zip(kn, nv)) or 1
        for i, (k, n) in enumerate(zip(kn, nv)):
            if k + n > 0:
                ax.text(i, k + n + top * 0.02, f'{k / (k + n) * 100:.0f}% known',
                        ha='center', va='bottom', fontsize=8.5, fontweight='bold')
        ax.set_ylim(0, top * 1.12)
        ax.set_xticks(x)
        ax.set_xticklabels(ROWS, fontsize=11)
        ax.set_ylabel('Gene loci (thousands)', fontsize=11, fontweight='bold')
        ax.set_title(title, fontsize=13, fontweight='bold')
        ax.legend(handles=[mpatches.Patch(facecolor='#555555', edgecolor='black', label='Known (darker)'),
                           mpatches.Patch(facecolor='#C8C8C8', edgecolor='black', label='Novel (lighter)')],
                  fontsize=9)
        ax.grid(True, alpha=0.3, axis='y')
    fig.suptitle(f'Known vs Novel Gene Loci  ({LABEL})', fontsize=15, fontweight='bold')
    qc.save(fig, out_path('CLASSIFY_04_known_novel.png'), dpi=200)


# ---------------------------------------------------------------------------
# 5. per-sample passed isoforms / FL reads (from the ALL filtered classification)
# ---------------------------------------------------------------------------
def fig_cross_sample():
    if 'ALL' not in DATA:
        return
    d = DATA['ALL']
    samples = [s for s in qc.SAMPLES_BY_GROUP if f'ps_{s}_isopass' in d]
    if not samples:
        print('  (per-sample figure skipped: no FL.<sample> columns matched samples.tsv)')
        return
    colors = [qc.sample_color(s) for s in samples]
    iso = [g(d, f'ps_{s}_isopass') / 1e3 for s in samples]
    flr = [g(d, f'ps_{s}_flpass') / 1e6 for s in samples]
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(max(14, 1.3 * len(samples) * 2), 5.8))
    for ax, vals, ylab, unit, title in (
            (axa, iso, 'Passed Isoforms (Thousands)', 'K', 'Passed Isoforms Detected per Sample'),
            (axb, flr, 'FL Reads Passed (Millions)', 'M', 'FL Reads in Passed Isoforms per Sample')):
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
    fig.suptitle(f'Cross-Sample Consistency (post-filter, from pooled ALL)  ({LABEL})',
                 fontsize=15, fontweight='bold')
    qc.save(fig, out_path('CLASSIFY_05_cross_sample.png'), dpi=200)


# ---------------------------------------------------------------------------
# 6. categories by isoform vs by read (post-filter)
# ---------------------------------------------------------------------------
def fig_iso_vs_read(stage='post'):
    fig, axes = plt.subplots(2, len(ROWS), figsize=(max(7, 4.7 * len(ROWS)), 7), squeeze=False)
    for ci, ds in enumerate(ROWS):
        d = DATA[ds]
        for ri, (pref, rlab, unit) in enumerate([(stage, 'by isoform', 'isoforms'),
                                                 (stage + 'read', 'by read', 'reads')]):
            ax = axes[ri][ci]
            vals = [cat_val(d, pref, c) for c in CATS]
            pcts = [v / (sum(vals) or 1) * 100 for v in vals]
            ax.bar(range(len(CATS)), pcts, color=[CAT_COL[c] for c in CATS],
                   edgecolor='black', linewidth=0.8)
            for xi, p in enumerate(pcts):
                if p > 1.5:
                    ax.text(xi, p + max(pcts) * 0.02, f'{p:.0f}', ha='center', va='bottom',
                            fontsize=7.5, fontweight='bold')
            ax.set_xticks(range(len(CATS)))
            ax.set_xticklabels([CAT_SHORT[c] for c in CATS], fontsize=8, rotation=30, ha='right')
            ax.set_ylim(0, (max(pcts) or 1) * 1.18)
            if ci == 0:
                ax.set_ylabel(f'{rlab}\n% of {unit}', fontsize=10.5, fontweight='bold')
            if ri == 0:
                ax.set_title(f'{ds} ({qc.n_samples(ds)} samples)', fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.3, axis='y')
    stage_disp = 'Post-filter' if stage == 'post' else 'Pre-filter'
    fig.suptitle(f'Structural Categories: by Isoform vs by Read - {stage_disp}  ({LABEL})',
                 fontsize=15, fontweight='bold', y=1.0)
    qc.save(fig, out_path('CLASSIFY_06_iso_vs_read.png'), dpi=200)


# ---------------------------------------------------------------------------
# 7. saturation curves (pigeon report, 07d)
# ---------------------------------------------------------------------------
def read_saturation(ds):
    p = os.path.join('results', LABEL, ds, f'clustered_{ds}.saturation.txt')
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        head = fh.readline().split()
        rows = [list(map(float, line.split())) for line in fh if line.strip()]
    if not rows:
        return None
    arr = np.array(rows)
    return {h: arr[:, i] for i, h in enumerate(head)}


def fig_saturation():
    sat = {ds: read_saturation(ds) for ds in qc.DATASETS}
    sat = {k: v for k, v in sat.items() if v is not None}
    if not sat:
        print('  (saturation figure skipped: no *.saturation.txt, run 07d)')
        return
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for ax, what, title in ((axes[0], 'genes', 'Genes'), (axes[1], 'isoforms', 'Isoforms')):
        for ds, s in sat.items():
            col = qc.dataset_color(ds)
            x = s['reads'] / 1e6
            ax.plot(x, s[f'unique_{what}'] / 1e3, color=col, lw=2.2, label=f'{ds} all')
            if f'unique_{what}_known' in s:
                ax.plot(x, s[f'unique_{what}_known'] / 1e3, color=col, lw=1.6, ls='--',
                        label=f'{ds} known')
        ax.set_xlabel('FL reads subsampled (millions)', fontsize=11, fontweight='bold')
        ax.set_ylabel(f'Unique {what} (thousands)', fontsize=11, fontweight='bold')
        ax.set_title(f'{title} discovered vs depth', fontsize=13, fontweight='bold')
        ax.set_ylim(bottom=0)
        ax.set_xlim(left=0)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9, loc='lower right')
    fig.suptitle(f'Saturation (pigeon report, singletons excluded)  ({LABEL})',
                 fontsize=15, fontweight='bold')
    qc.save(fig, out_path('CLASSIFY_07_saturation.png'), dpi=200)


fig_summary()
fig_categories()
fig_filter_reasons()
fig_known_novel()
fig_cross_sample()
fig_iso_vs_read('post')
fig_saturation()

print('\n' + '=' * 80)
print('CLASSIFY QC SUMMARY')
print('=' * 80)
for ds in ROWS:
    d = DATA[ds]
    pin = g(d, 'pre_input')
    pas = g(d, 'post_passed')
    fsm = cat_val(d, 'post', 'FSM') / (pas or 1) * 100
    print(f"  {ds:<10} input={pin:>9,}  passed={pas:>9,} ({(pas / pin * 100 if pin else 0):4.1f}%)  "
          f"FSM_post={fsm:4.1f}%  known_genes_post={g(d, 'post_genes_known'):>6,}")
print('\nDone.')

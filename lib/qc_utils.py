"""
qc_utils.py
Shared bits for the *_QC.py plotting scripts: reads config.sh and the sample
sheet, and hands out consistent colours per group/sample.

The QC scripts add lib/ to sys.path and do `import qc_utils as qc`.
No need to edit this; change config.sh / samples.tsv instead.
"""
import os
import re
import subprocess

import matplotlib
matplotlib.use("Agg")  # no display on the cluster
import matplotlib.colors
import matplotlib.pyplot as plt
from matplotlib import cm

LIB_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(LIB_DIR)

ALL_COLOR = "#4E2A84"   # pooled data (ALL) is always purple
ALL_LIGHT = "#B6A6D9"
# one colormap per group, in sample-sheet order. Purples is kept for ALL.
GROUP_CMAPS = ["Blues", "Reds", "Greens", "Oranges", "RdPu", "YlOrBr", "GnBu", "Greys"]

_KEYS = ["PROJECT_DIR", "SAMPLE_SHEET", "BUILD_LABEL", "EXPECTED_ARRAY_SIZE",
         "MAIN_CHROM_REGEX", "GENOME_FA", "SORTED_GTF"]


def load_config():
    """Source lib/functions.sh in bash and pull out the settings we need."""
    script = 'source "%s/functions.sh" >/dev/null || exit 1\n' % LIB_DIR
    for k in _KEYS:
        script += 'printf "%s=%%s\\n" "${%s}"\n' % (k, k)
    script += 'printf "HIFI_BAMS=%s\\n" "$(IFS=$\'\\x1f\'; echo "${HIFI_BAMS[*]}")"\n'
    res = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
    if res.returncode != 0:
        raise SystemExit("could not read config.sh:\n" + res.stderr)
    cfg = {}
    for line in res.stdout.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            cfg[k] = v
    cfg["HIFI_BAMS"] = [x for x in cfg.get("HIFI_BAMS", "").split("\x1f") if x]
    try:
        cfg["EXPECTED_ARRAY_SIZE"] = int(cfg.get("EXPECTED_ARRAY_SIZE") or 8)
    except ValueError:
        cfg["EXPECTED_ARRAY_SIZE"] = 8
    return cfg


CFG = load_config()


def load_sheet(path=None):
    """Rows of the sample sheet as dicts (barcode, sample, group, cell)."""
    path = path or CFG["SAMPLE_SHEET"]
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 3 or p[0] == "barcode":
                continue
            cell = p[3] if len(p) > 3 and p[3] else "all"
            rows.append(dict(barcode=p[0], sample=p[1], group=p[2], cell=cell))
    if not rows:
        raise SystemExit(f"no samples found in {path}")
    return rows


SHEET = load_sheet()
SAMPLES = list(dict.fromkeys(r["sample"] for r in SHEET))      # sheet order
GROUPS = list(dict.fromkeys(r["group"] for r in SHEET))
GROUP_OF = {r["sample"]: r["group"] for r in SHEET}
# samples ordered by group, so group blocks sit together on plots
SAMPLES_BY_GROUP = [s for g in GROUPS for s in SAMPLES if GROUP_OF[s] == g]
DATASETS = ["ALL"] + GROUPS


def sample_for(cell, pair):
    for r in SHEET:
        if r["barcode"] == pair and (r["cell"] == "all" or r["cell"] in cell):
            return r["sample"]
    return None


def split_name(path):
    """<cell>.<pair>[.flnc].bam -> (cell, pair)"""
    b = os.path.basename(path)
    b = re.sub(r"\.bam$", "", b)
    b = re.sub(r"\.flnc$", "", b)
    cell, _, pair = b.rpartition(".")
    return cell, pair


def _cmap(group):
    name = GROUP_CMAPS[(GROUPS.index(group) if group in GROUPS else 0) % len(GROUP_CMAPS)]
    try:
        return matplotlib.colormaps[name]
    except AttributeError:          # matplotlib < 3.5
        return cm.get_cmap(name)


def group_color(group):
    if group == "ALL":
        return ALL_COLOR
    return matplotlib.colors.to_hex(_cmap(group)(0.70))


def group_light(group):
    if group == "ALL":
        return ALL_LIGHT
    return matplotlib.colors.to_hex(_cmap(group)(0.30))


def sample_color(sample):
    g = GROUP_OF.get(sample, GROUPS[0])
    members = [s for s in SAMPLES_BY_GROUP if GROUP_OF[s] == g]
    n = len(members)
    k = members.index(sample) if sample in members else 0
    # spread shades from dark to mid so none disappear on white
    x = 0.85 - (0.45 * k / (n - 1)) if n > 1 else 0.70
    return matplotlib.colors.to_hex(_cmap(g)(x))


def dataset_color(ds):
    if ds == "ALL" or ds in GROUPS:
        return group_color(ds)
    return sample_color(ds)


def n_samples(ds):
    if ds == "ALL":
        return len(SAMPLES)
    return sum(1 for s in SAMPLES if GROUP_OF[s] == ds)


def dataset_label(ds):
    if ds == "ALL" or ds in GROUPS:
        return f"{ds}\n({n_samples(ds)} samples)"
    return f"{ds}\n{GROUP_OF.get(ds, '')}"


def sample_label(sample):
    return f"{sample}\n{GROUP_OF.get(sample, '')}"


def group_blocks(ax, samples, ymax):
    """Dotted lines between group blocks + a group tag at the top of each block.
    `samples` is the x order used on the plot (x = 0..n-1). ymax is not used
    any more (kept so old calls still work)."""
    if len(GROUPS) < 2:
        return
    start = 0
    for i in range(1, len(samples) + 1):
        if i == len(samples) or GROUP_OF[samples[i]] != GROUP_OF[samples[start]]:
            g = GROUP_OF[samples[start]]
            # y in axes fraction so it also works when the y axis does not start at 0
            ax.text((start + i - 1) / 2, 0.96, g, ha="center", va="top",
                    transform=ax.get_xaxis_transform(), fontsize=11, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor=group_color(g),
                              alpha=0.3, edgecolor="black"))
            if i < len(samples):
                ax.axvline(i - 0.5, color="black", ls=":", lw=1.5, alpha=0.5)
            start = i


def table_fig(rows, title, path, status_col=None, row_groups=None, header_color="#333333",
              col_widths=None, fontsize=10):
    """Draw a list-of-lists as a matplotlib table. rows[0] is the header.
    status_col: index of a PASS/REVIEW column to colour green/red.
    row_groups: optional list (len rows-1) of group names to tint rows."""
    ncol = len(rows[0])
    fig, ax = plt.subplots(figsize=(max(8, 1.6 * ncol), 0.8 + 0.5 * len(rows)))
    ax.axis("off")
    t = ax.table(cellText=rows, cellLoc="center", loc="center", colWidths=col_widths)
    t.auto_set_font_size(False)
    t.set_fontsize(fontsize)
    t.scale(1, 2.2)
    for j in range(ncol):
        c = t[(0, j)]
        c.set_facecolor(header_color)
        c.set_text_props(weight="bold", color="white")
    for i in range(1, len(rows)):
        for j in range(ncol):
            c = t[(i, j)]
            if status_col is not None and j == status_col:
                ok = str(rows[i][j]).upper() in ("PASS", "OK")
                c.set_facecolor("#90EE90" if ok else "#F4A6A6")
                c.set_text_props(weight="bold")
            elif row_groups is not None:
                c.set_facecolor(group_light(row_groups[i - 1]))
            else:
                c.set_facecolor("#F2F2F2" if i % 2 else "white")
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    save(fig, path)


def save(fig, path, dpi=300):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}")


def use_style():
    for s in ("seaborn-v0_8-darkgrid", "seaborn-darkgrid"):
        if s in plt.style.available:
            plt.style.use(s)
            return


def read_kv(path):
    """key<TAB>value stat files written by the *_EXTRACT_*.sh scripts."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            try:
                out[p[0]] = float(p[1]) if "." in p[1] else int(p[1])
            except ValueError:
                out[p[0]] = p[1]
    return out

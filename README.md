# PacBio Iso-Seq (Kinnex) long-read pipeline

Talha Lone, Parvez lab. Last updated 2026-09-17.

These are the scripts we use to go from PacBio Revio Kinnex full-length RNA
data (HiFi BAMs) to a classified, filtered set of isoforms with per-sample
read counts. There is a QC step at every stage. We most recently ran this on
6 mouse B-cell samples (3 control, 3 knockout) on one Revio SMRT cell, and
some of those numbers are given below so you know roughly what to expect.

Everything runs on SLURM (written for Northwestern Quest) with a single
conda environment.

```
HiFi BAM
  01  skera split           concatemers -> S-reads
  02  lima                  demultiplex by cDNA barcode, name the samples
      isoseq refine         trim polyA, remove chimeras -> FLNC reads
  03  file lists            which FLNC BAMs belong to which sample/group
  04  isoseq cluster2       pool all samples, cluster -> transcripts
  05  pbmm2                 map transcripts to the genome
  06  isoseq collapse       unique isoforms + per-sample FLNC counts
  07  pigeon                classify (FSM/ISM/NIC/NNC...), filter, saturation
```

## Folder layout

```
forShibo/
  README.md          this file
  config.sh          <- EDIT: paths, reference, options
  samples.tsv        <- EDIT: barcode -> sample -> group
  lib/               shared helpers (no need to touch)
  SETUP/             conda environment
  RAW/               your HiFi BAMs (empty)
  REFERENCE/         genome + annotation (empty, read REFERENCE/README.md)
  01_SEGMENTATION/   01a, 01b, SKERA_QC.py, adapters/
  02_LIMA_REFINE/    02a, 02b, 02c, LIMA_QC.py, REFINE_QC.py, primers/
  03_MERGE_FOFN/     03_CREATE_FOFN.sh
  04_CLUSTERING/     04a, 04b, CLUSTER_QC.py
  05_MAPPING/        05a-05d, MAPPING_QC.py
  06_COLLAPSE/       06a, 06b, COLLAPSE_QC.py
  07_CLASSIFY/       07a-07e, CLASSIFY_QC.py
```

Each step folder gets a `logs/` (SLURM output) and a `results/` folder when
you run it. QC tables and figures go to `results/qc/`.

## Before you start

### 1. Environment

```
bash SETUP/setup_environment.sh
```

See `SETUP/README.md`. The jobs activate the env on their own.

### 2. Reference

Put a genome FASTA and a GTF in `REFERENCE/`. **Read `REFERENCE/README.md`
first.** pigeon only accepts GTFs in a particular format; GENCODE works
as is, NCBI RefSeq needs a fix (`07a`).

### 3. Edit `config.sh`

Everything below is in `config.sh`. Relative paths are relative to this
folder.

| setting | what to put there |
|---|---|
| `CONDA_SH` | path to `etc/profile.d/conda.sh` of your mamba/conda. The default is the Quest one |
| `CONDA_ENV` | env name from step 1 (`pacbio_workflow`) |
| `MODULES` | modules to `module load` in every job. Leave empty if samtools is in the env (it is) |
| `HIFI_BAMS` | one line per HiFi BAM (per SMRT cell / Kinnex index). Don't include `*.unassigned.bam` |
| `SAMPLE_SHEET` | normally `samples.tsv` |
| `SKERA_ADAPTERS` | `mas8_primers.fasta` (included) for Kinnex full-length RNA. MAS-Seq 16-fold libraries need mas16 from PacBio |
| `EXPECTED_ARRAY_SIZE` | 8 for Kinnex, only used for the QC plot |
| `LIMA_PRIMERS` | `IsoSeq_v2_primers_12.fasta` (included, bc01-bc12). Change only if your core used a different primer set |
| `MIN_POLYA_LENGTH` | 20 (isoseq default) |
| `BUILD_LABEL` | short name for genome + annotation, e.g. `GRCh38_v47`. Used in folder names |
| `GENOME_FA` | genome FASTA, uncompressed |
| `ANNOTATION_GTF` | pigeon-compatible GTF (the unsorted one) |
| `CAGE_BED`, `POLYA_MOTIFS` | optional, leave `""` to skip. See REFERENCE/README.md |
| `MAIN_CHROM_REGEX` | which chromosomes get their own bar in the mapping QC. Default covers `chr1`/`1`...`X`,`Y`,`M`,`MT`. For NCBI RefSeq genomes (`NC_...` names) use `'^NC_'` |
| `CLUSTER_OPTS` | `--singletons` (we keep them and drop them later) |
| `COLLAPSE_OPTS` | `--do-not-collapse-extra-5exons` (see note under step 06) |
| `THREADS_*` | threads per tool. If you change these, change `--cpus-per-task` in the script too |

### 4. Edit `samples.tsv`

Tab separated, one row per sample:

```
barcode                       sample   group   cell
IsoSeqX_bc01_5p--IsoSeqX_3p   CTRL_1   CTRL    all
IsoSeqX_bc02_5p--IsoSeqX_3p   CTRL_2   CTRL    all
IsoSeqX_bc04_5p--IsoSeqX_3p   KO_1     KO      all
```

- `barcode` is the lima pair name: `IsoSeqX_bcNN_5p--IsoSeqX_3p`. Get the
  barcode for each sample from the core's sample sheet.
- `sample` ends up as a column name in the count matrix. No spaces or commas.
- `group` is your condition. It's used for the optional per-group runs and
  for colours in the QC plots. `ALL` is reserved.
- `cell` only matters if the same barcode was used for different samples on
  different SMRT cells. Then put part of the movie name there
  (e.g. `m84241_260523_024948_s1`). Otherwise `all`.

A sample that was sequenced on more than one cell just gets one row (with
`cell` = `all`); its reads are combined automatically.

### 5. SLURM account

All scripts have `#SBATCH --account=b1042` and `--partition=genomics`
(our Quest allocation). Change them all at once with:

```
sed -i 's/--account=b1042/--account=YOUR_ALLOC/; s/--partition=genomics/--partition=YOUR_PARTITION/' 0*/*.sh
```

## Running it

**Always `cd` into the step folder before `sbatch`.** The scripts find
`config.sh` and write `logs/` relative to where you submit from. If you submit
from somewhere else the job dies right away, sometimes without writing a log.

Array jobs need `--array=1-N`, where N is the number of inputs. The scripts
whose N isn't obvious will print it for you with `--count`.

### 01 Segmentation

```
cd 01_SEGMENTATION
sbatch --array=1-<number of HIFI_BAMS> 01a_RUN_SKERA.sh
sbatch --array=1-<number of HIFI_BAMS> 01b_ARRAY_HIST.sh     # optional, for the QC plot
python SKERA_QC.py
```

Check `percent_full_array` in `results/segmented/*.summary.json`. It should be
above ~90%, and the mean array size should be close to 8. If full arrays are
near 0% you have the wrong adapter set. Ours: 9.68 M HiFi reads -> 75.3 M
S-reads (7.78x), 94.5% full arrays.

### 02 Demultiplex + refine

```
cd 02_LIMA_REFINE
sbatch --array=1-<number of HIFI_BAMS> 02a_RUN_LIMA.sh
# when lima is done:
bash 02b_RUN_REFINE.sh --count                  # prints N
sbatch --array=1-N 02b_RUN_REFINE.sh
sbatch --array=1-N 02c_RUN_REFINE_NOPOLYA.sh    # optional, QC only
python LIMA_QC.py
python REFINE_QC.py
```

`02a` writes your sample names into the BAMs. Don't skip or work around it:
that is where collapse gets the per-sample column names later.

Check:
- lima `*.lima.summary`: most reads should pass ("above all thresholds").
  Ours was 97.6% (73.4 of 75.3 M).
- Every sample in `samples.tsv` should get a BAM. A missing one means a
  wrong barcode in the sheet.
- Refine: FL -> FLNC loss (chimeras) should be tiny (ours was <0.1%), and
  the reads per sample shouldn't be wildly uneven. We had 10.6-15.1 M FLNC
  per sample.

### 03 File lists

```
cd 03_MERGE_FOFN
bash 03_CREATE_FOFN.sh
```

Runs in a few seconds, no sbatch needed. It stops if any sample in the sheet has no
refined BAM.

### 04 Clustering

```
cd 04_CLUSTERING
sbatch 04a_RUN_CLUSTER.sh ALL
sbatch 04b_EXTRACT_CLUSTER_STATS.sh     # after 04a
python CLUSTER_QC.py
```

**Cluster everything together (`ALL`).** That gives one set of isoforms,
and collapse then counts reads from each sample against that same set. That
is what you want for any differential analysis. `sbatch 04a_RUN_CLUSTER.sh
<GROUP>` clusters one group on its own. We only used that for exploring;
isoforms from separate runs can't be matched up one-to-one.

This is the heavy step. For ~73 M FLNC reads we gave it 256 GB and it took
12 h. Exit code 137 means it ran out of memory.

### 05 Mapping

```
cd 05_MAPPING
sbatch 05a_INDEX_GENOME.sh               # once; can run while 04 is going
sbatch 05b_MAP_CLUSTERED.sh ALL          # after 04a
bash 05c_MAP_INDIVIDUAL.sh --count       # optional per-sample BAMs for IGV
sbatch --array=1-N 05c_MAP_INDIVIDUAL.sh
sbatch 05d_EXTRACT_MAP_STATS.sh          # after 05b (and 05c if you ran it)
python MAPPING_QC.py
```

Mapping rate for the clustered transcripts should be very high. Ours was
99.96%. Anything under ~95% probably means the wrong genome or contamination.

### 06 Collapse

```
cd 06_COLLAPSE
sbatch 06a_COLLAPSE.sh ALL
sbatch 06b_EXTRACT_COLLAPSE_STATS.sh     # needs 05d done
python COLLAPSE_QC.py
```

Output that matters: `results/<BUILD_LABEL>/ALL/clustered_ALL.<BUILD_LABEL>.collapsed.flnc_count.txt`,
the isoform x sample count matrix. Its header should list your sample
names. If it just says `BioSample_1`, the names were lost at lima (the
script warns about this).

`--do-not-collapse-extra-5exons` keeps 5' truncated versions as separate
isoforms instead of merging them into the longer model. We did this because
we cared about alternative first exons. Without the flag, isoseq merges them
and you end up with fewer isoforms. pigeon filter and later count filters
deal with most of the truncated ones either way.

### 07 Classification (pigeon)

```
cd 07_CLASSIFY
bash 07a_FIX_NCBI_GTF.sh REFERENCE/<ncbi>.gtf   # ONLY for NCBI RefSeq GTFs, then update config.sh
sbatch 07b_PIGEON_PREPARE_REF.sh                  # once per annotation
sbatch 07c_PIGEON_CLASSIFY.sh ALL                 # after 06a and 07b
sbatch 07d_PIGEON_REPORT.sh ALL                   # saturation curve
bash 07e_EXTRACT_CLASSIFY_STATS.sh
python CLASSIFY_QC.py
```

Main outputs in `results/<BUILD_LABEL>/ALL/`:

| file | what it is |
|---|---|
| `clustered_ALL_classification.txt` | every isoform, with category, gene, per-sample FL counts |
| `clustered_ALL_classification.filtered_lite_classification.txt` | isoforms that passed the filter |
| `clustered_ALL.sorted.filtered_lite.gff` | the same isoforms as a GFF (load into IGV) |
| `clustered_ALL_classification.filtered_lite_reasons.txt` | why each removed isoform was removed |
| `clustered_ALL.summary.txt`, `*_classification.filtered.summary.txt` | counts before/after filtering |
| `clustered_ALL.saturation.txt` | genes/isoforms vs. reads |

pigeon filter removes likely intra-priming, RT switching and low-coverage /
non-canonical junction isoforms. In our data 45.7% of 2.56 M collapsed
isoforms passed. By read, ~68% of reads were full splice matches (FSM) before
filtering and ~77% after.

Expect a lot of ISM (incomplete splice match). Many of them are 5' truncated
reads rather than real alternative starts. If you have CAGE peaks for your
species, set `CAGE_BED` so you can tell the two apart (REFERENCE/README.md).

## Starting from segmented or demultiplexed reads

Revio can run skera and lima on the instrument. Check what the core sent.

- **Segmented reads** (`*.segmented.bam` / "S-reads"): skip 01. Link or copy
  each file to `01_SEGMENTATION/results/segmented/<cell>.segmented.bam`, where
  `<cell>` is the HiFi BAM name without `.hifi_reads` and `.bam`
  (`m84xxx_s1.hifi_reads.bcM0001.bam` -> `m84xxx_s1.bcM0001`). Keep
  `HIFI_BAMS` filled in; it's used to name things.
- **Already demultiplexed but not refined** (per-barcode BAMs): make sure the
  sample names are set (`samtools view -H file.bam | grep '^@RG'` should show
  `SM:` your sample name). If they aren't, rerun lima yourself from the
  segmented reads. It's cheaper than finding out at step 06. Otherwise put
  them in `02_LIMA_REFINE/results/demux/` as `<cell>.<barcode pair>.bam` and
  start at 02b.
- **FLNC reads**: put them in `02_LIMA_REFINE/results/refined/` as
  `<cell>.<barcode pair>.flnc.bam` and start at 03.

## Resources we used

Wall times from our run (one Revio cell, 6 samples, ~73 M FLNC reads).
These are the memory amounts we requested, not what the jobs actually used.

| step | cpus | memory | time |
|---|---|---|---|
| skera | 4 | 16G | 3.7 h |
| lima | 8 | 16G | 2.5-3 h |
| refine (per sample) | 4 | 16G | ~1 h |
| cluster2 ALL | 16 | 256G | 12 h |
| pbmm2 clustered ALL | 16 | 64G | 50 min |
| pbmm2 per sample | 16 | 48G | ~2 h each |
| collapse ALL | 8 | 64G | 8 min |
| pigeon classify + filter ALL | 8 | 32G | 7 min |

The `#SBATCH` time limits have plenty of slack. More cells or samples mostly
make cluster2 slower and hungrier.

## Things to watch out for

Most of these don't produce an error. The results just come out wrong.

- Wrong skera adapter set: ~0% full arrays (we checked mas16 against our
  mas8 library).
- Hand-edited primer FASTA: an old edited copy of ours had a truncated primer
  and a typo in one barcode. Use the PacBio file as shipped.
- Lima without biosample names: collapse puts every read into a single
  `BioSample_1` column.
- CAGE peaks from the wrong assembly: almost no FSMs land in peaks (1.5% vs ~37%).
- Submitting from the wrong folder: SLURM can't create the log file, so the
  job fails without leaving a log.
- Conda env on a scratch filesystem that gets purged: tools stop starting.

## What's not included

This stops at a classified, filtered isoform set plus a count matrix.
We did differential transcript usage (DRIMSeq, IsoformSwitchAnalyzeR)
and gene-level DESeq2 afterwards in R, starting from
`flnc_count.txt` / the filtered classification and `03_MERGE_FOFN/fofn_files/design.tsv`.
We can share those scripts too, but they depend more on the experiment.

## Links

- Iso-Seq / pigeon docs: https://isoseq.how
- Kinnex primer and adapter files: https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/

Questions: Talha Lone (Parvez lab).

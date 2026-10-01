# REFERENCE

This folder is empty on purpose. Put your genome and annotation here (or
anywhere else) and point `GENOME_FA` and `ANNOTATION_GTF` in `config.sh` at
them.

The scripts write index files next to these files, so the folder has to be
writable:

- `05a` writes `<genome>.fai` and `<genome>.ISOSEQ.mmi`
- `07b` writes `<annotation>.sorted.gtf` and `<annotation>.sorted.gtf.pgi`

## What goes in here

| file | needed? | notes |
|---|---|---|
| genome FASTA | yes | uncompressed. Primary assembly if there is one (see below) |
| annotation GTF | yes | same assembly and same chromosome names as the FASTA |
| CAGE / TSS peaks BED | optional | set `CAGE_BED`. Fills `within_CAGE_peak` in the classification |
| polyA motif list | optional | set `POLYA_MOTIFS`. Fills `polyA_motif` |

**Genome.** Use the primary assembly (no alt haplotypes or patches) if your
source offers one. Alt contigs make a lot of reads map equally well to two
places and end up with MAPQ 0. Some examples:

- GENCODE human: `GRCh38.primary_assembly.genome.fa.gz`
- GENCODE mouse: `GRCm39.primary_assembly.genome.fa.gz` (we used `GRCm39.genome.fa`,
  which only has the reference chromosomes and scaffolds, and that was also fine)
- Ensembl: `*.dna.primary_assembly.fa.gz`, or `*.dna.toplevel.fa.gz` for species
  that don't have a primary assembly file

Unzip it first (`gunzip -k`). pigeon and samtools can read bgzip files but not
plain gzip, and it isn't worth the trouble.

## Making the annotation work with pigeon

pigeon is fussy about the GTF. Most failures in `07b`/`07c` come from this.
What it needs:

1. **GTF, not GFF3.** If you only have GFF3, convert it:
   `gffread annotation.gff3 -T -o annotation.gtf` (gffread is in the conda env).
   Then check that `gene_name` survived (step 3).
2. **Chromosome names must match the FASTA exactly.** `chr1` vs `1` is the
   usual mismatch (UCSC-style vs Ensembl-style). Compare the two:
   ```
   cut -f1 genome.fa.fai | head
   grep -v '^#' annotation.gtf | cut -f1 | sort -u | head
   ```
   NCBI RefSeq uses accessions (`NC_000001.11`) in both files, which is fine
   as long as both came from the same NCBI download. `07b` stops if the first
   GTF chromosome is not in the genome.
3. **Every line needs `gene_id`, `gene_name` and (apart from gene lines)
   `transcript_id`.** Check it:
   ```
   grep -v '^#' annotation.gtf | grep -vc 'gene_name "'     # want 0
   ```
4. **There have to be `exon` lines.** `transcript` and `gene` lines are good to have.
5. **Uncompressed, with no blank lines, and nothing that isn't 9 tab-separated columns.**

What that means for the common sources:

- **GENCODE**: works as is. We used `gencode.vM38.annotation.gtf` directly.
  The `##` header lines did not cause problems.
- **Ensembl**: gene_id/gene_name/transcript_id are all there. Genes that
  have no symbol are missing `gene_name`. We haven't run an Ensembl GTF
  through this exact setup, so check step 3 above. If a lot of lines are
  missing it, you can copy gene_id into gene_name with a one-line awk, like
  `07a` does.
- **NCBI RefSeq**: does NOT work as downloaded. It writes `gene "xyz"` instead
  of `gene_name "xyz"`, and it puts an empty `transcript_id ""` on gene lines.
  Run the fix script, then point `ANNOTATION_GTF` at the `.pigeon.gtf` it
  writes:
  ```
  bash 07_CLASSIFY/07a_FIX_NCBI_GTF.sh REFERENCE/GCF_xxxxx_genomic.gtf
  ```
  We had to do this for zebrafish (GRCz12tu is only available from NCBI).
- **Custom / lab annotations** (StringTie, liftoff, etc.): go through the
  checklist above. StringTie output usually has no `gene_name`.

Point `ANNOTATION_GTF` at the unsorted `.gtf`. `07b` runs `pigeon prepare`,
which writes `<name>.sorted.gtf` + `.pgi`, and classify uses those.

If `pigeon prepare` still complains, look at the first line it mentions.
It's usually one odd record, e.g. a feature with start > end or a gene with no
exons. Removing it is fine.

## Optional files

**CAGE peaks** (`CAGE_BED`). These let pigeon say whether a transcript 5' end
falls in a known TSS peak. That is the best evidence for telling a real
alternative TSS apart from a 5' truncated read (a lot of ISMs are the
second kind). Human and mouse peaks are on refTSS (https://reftss.riken.jp) for
the current assemblies. Things to watch:

- It has to be the **same assembly** as your genome. We once used an old mm10
  file with mm39 data by mistake. Nothing errored; only 1.5% of FSM
  transcripts landed in a peak, against ~37% with the right mm39 file. If your
  FSM-in-peak number is tiny, check this first.
- Same chromosome naming as the genome.
- Plain BED, no header or `track` line. Sort it (`sort -k1,1 -k2,2n`).
- `07b` runs `pigeon prepare` on it (makes `<name>.sorted.bed` + `.pgi`).

**polyA motifs** (`POLYA_MOTIFS`). A text file with one motif per line,
strongest first. The one SQANTI3 ships with is a good default:

```
AATAAA
ATTAAA
AGTAAA
TATAAA
CATAAA
GATAAA
AATATA
AATACA
AATAGA
AAAAAG
ACTAAA
AAGAAA
AATGAA
TTTAAA
AAAACA
GGGGCT
```

PacBio also has ready-made pigeon reference sets for human and mouse, with the
genome, GTF, CAGE and polyA files already in the right format. See the pigeon
page on https://isoseq.how. They're the easiest option if one matches the
annotation version you want.

## Quick check before running 05a / 07b

```
ls -lh REFERENCE/
head -2 REFERENCE/<genome>.fa
grep -v '^#' REFERENCE/<annotation>.gtf | head -3
grep -v '^#' REFERENCE/<annotation>.gtf | awk -F'\t' 'NF!=9' | wc -l   # want 0
```

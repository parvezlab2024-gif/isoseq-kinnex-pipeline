# SETUP

Build the conda environment once before running anything.

```
cd /path/to/forShibo
bash SETUP/setup_environment.sh            # makes an env called pacbio_workflow
```

On Quest, mamba is at `/hpc/software/mamba/23.1.0` (`module load mamba` or
source its `etc/profile.d/conda.sh`). The build takes 10-20 minutes and pulls
~370 MB.

If you give the env a different name (`bash SETUP/setup_environment.sh myenv`)
or install conda somewhere else, update `CONDA_ENV` and `CONDA_SH` in
`config.sh`. Every job activates the env itself; you don't have to activate it
before `sbatch`.

What's in it (`isoseq_environment.yml`):

| tool | version | used in |
|---|---|---|
| skera (pbskera) | 1.4.0 | 01 |
| lima | 2.13.0 | 02 |
| isoseq | 4.3.0 | 02 refine, 04 cluster2, 06 collapse |
| pbmm2 | 1.17.0 | 05 |
| pigeon (pbpigeon) | 1.4.0 | 07 |
| samtools / htslib | 1.22.1 | everywhere (counts, indexes) |
| gffread | 0.12.7 | only if you need to turn a GFF3 into a GTF |
| python 3.11 + numpy, pandas, matplotlib | | the `*_QC.py` plots |

These are the versions we ran. Newer PacBio releases usually work, but
option names and output file names do change between releases (pigeon
especially), so if you upgrade, check the `--help` of anything that fails.

Quick check after building:

```
conda activate pacbio_workflow
isoseq --version; pigeon --version; samtools --version | head -1
```

Problems we've run into:

- `conda activate` does nothing in a batch job: the script needs to source
  `conda.sh` first. The pipeline does that using `CONDA_SH` from config.sh.
- Python plots fail with import errors that mention `~/.local/lib/python...`:
  packages installed with `pip install --user` are overriding the env. Run
  `export PYTHONNOUSERSITE=1` before the QC scripts.
- Env stored on a scratch filesystem that gets purged: the env breaks without
  an obvious error (binaries get zeroed). Keep it in home or project space.

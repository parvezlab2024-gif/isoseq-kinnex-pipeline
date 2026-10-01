# RAW

Put the sequencing core's delivery here (or leave it where it is and use
full paths in `config.sh`).

For a Revio Kinnex run you get one folder per SMRT cell, with a `hifi_reads/`
folder inside. The pipeline starts from the HiFi BAMs:

```
RAW/1_A01/hifi_reads/m84xxxx_xxxxxx_xxxxxx_s1.hifi_reads.bcM0001.bam
```

List each one under `HIFI_BAMS` in `config.sh`. Leave out the
`*.unassigned.bam` files; those are reads the instrument couldn't assign to a
Kinnex index.

If the core already ran skera and/or lima on the instrument (look for
`segmented_reads/` or per-barcode `fl*.bam` files), you can skip steps. See
"Starting from segmented or demultiplexed reads" in the main README.

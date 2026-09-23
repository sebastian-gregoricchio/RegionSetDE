# Androgen receptor example data

Androgen receptor (AR) ChIP-seq in LNCaP cells, used by the peak based workflow of RegionSetDE:
the sample sheet, the peak calls and the alignments the package counts over.

## Where it comes from

| | |
| --- | --- |
| Accession | [GSE284522](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284522) |
| Publication | Eickhoff *et al.*, *Communications Biology* **8**, 1043 (2025), doi:10.1038/s42003-025-08449-2 |
| Cells | LNCaP, hormone deprived for 72 h in stripped serum, then DMSO or 5 nM R1881 |
| Design | DMSO 4 h, R1881 4 h, R1881 24 h, three biological replicates each, plus one input |
| Condition | The non-targeting arm of the knockdown design of the study, hence the `NT` in the original file names |
| Assembly | GRCh38 (Ensembl 102), chromosomes named without the `chr` prefix |
| Reads | Paired-end, 54 bp |

| Sample | GEO |
| --- | --- |
| `AR_DMSO_r1`, `_r2`, `_r3` | GSM8685837, GSM8685838, GSM8685839 |
| `AR_R1881_4h_r1`, `_r2`, `_r3` | GSM8685840, GSM8685841, GSM8685842 |
| `AR_R1881_24h_r1`, `_r2`, `_r3` | GSM8685843, GSM8685844, GSM8685845 |
| `input` | GSM8685835 |

## How it was processed

Alignment, filtering and peak calling were run with
[SPACCa](https://github.com/sebastian-gregoricchio/SPACCa), which aligns with `bwa-mem`, keeps
reads at `MAPQ >= 20`, marks duplicates, and calls peaks with MACS2 in `BAMPE` mode with default
parameters. The `.narrowPeak` files here are those calls, made on the complete libraries.

## How it was cut down to a shippable size

`inst/scripts/make-example-bams.sh` holds the commands. In short:

1. The window `19:46,000,000-58,000,000` was cut out of each sorted BAM.
2. Proper pairs only (`-f 2`, the same set MACS2 counted in `BAMPE` mode), subsampled to 12 % of
   the read pairs with `samtools view -s 0.12 --subsample-seed 42`. Subsampling acts on whole
   templates, so no mate is ever orphaned.
3. The sequence, the base qualities, the alignment tags and the read names were thrown away, and
   the header was reduced to the sort order and the one contig the reads sit on. Counting reads a
   region needs the position and the CIGAR, nothing else, so the counts are unaffected by any of
   this beyond the subsampling.
4. Indexed with `samtools index -c`, which writes a CSI. On a window this far along a chromosome a
   BAI carries an entry per 16 kb from the start of the chromosome and comes out nine times larger.
5. The peak files were filtered to the same window. They were **not** recalled on the slice, so
   they carry the depth and the statistics the caller actually had.

What this means when reading numbers off the example: counts over the regions are the counts of
the real libraries scaled by 0.12, and ratios between samples and between regions are preserved,
because every library was subsampled at the same rate. Library sizes are those of the window, not
of the genome. The peak scores belong to libraries roughly eight times deeper than the alignments
shipped next to them, which is why a region can hold few fragments and still be a confident call.

The window was chosen for size, not for biology. Anything read off these files describes the
example, not the androgen response of LNCaP cells, and the paper is the place to look for that.

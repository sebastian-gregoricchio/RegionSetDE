# Region lists shipped with RegionSetDE

The files under `blacklist/` and `greenlist/` are redistributed as they were published, converted to
gzipped BED and nothing else. `regionLists.tsv` indexes them and is what `loadBlacklist()`,
`loadGreenlist()` and `availableRegionLists()` read. Adding a list means dropping a gzipped BED file
in one of the two folders and writing one row in that table.

## blacklist/*_encode.v2.bed.gz

The ENCODE blacklist, version 2, for hg19, hg38, mm10, dm3, dm6, ce10 and ce11. Regions of anomalous
coverage found across many experiments, from the pipeline of the Boyle laboratory.

* Source: https://github.com/Boyle-Lab/Blacklist, `lists/`, downloaded on 2026-09-22.
* Licence: GPL-3, the licence of that repository, which is the licence of this package as well.
* Reference: Amemiya HM, Kundaje A, Boyle AP (2019). The ENCODE blacklist: identification of
  problematic regions of the genome. *Scientific Reports* 9:9354. doi:10.1038/s41598-019-45839-z

## blacklist/hs1_excluderanges.v1.bed.gz

The exclusion set of T2T-CHM13v2.0, the assembly UCSC calls hs1. ENCODE never published one, so this
is the set the excluderanges authors produced by running the same Boyle laboratory software on the
T2T assembly: the high signal and low mappability regions, labelled as in the ENCODE files.

* Source: the `T2T.excluderanges` record of the excluderanges package, exported from AnnotationHub as
  BED on 2026-09-22.
* Licence: MIT, the licence of that package.
* Reference: Dozmorov MG et al. (2023). excluderanges: exclusion sets for T2T-CHM13, GRCm39, and
  other genome assemblies. *Bioinformatics* 39(4):btad198. doi:10.1093/bioinformatics/btad198
* The file covers chr1 to chr22 and chrX. CHM13 comes from a hydatidiform mole and carries no Y
  chromosome of its own, so the chrY added to v2.0 from HG002 has no exclusion regions here.

## blacklist/{hg38,mm39}_{cutrun,cuttag}.v1.bed.gz and greenlist/

The greenlists and the matched high signal regions of the CUT&RUN greenlist paper, one pair per
assay. The greenlist holds the bins whose background is consistent across experiments, found by
Shannon entropy over hundreds of public libraries and kept at least 5 kb away from genes; the
blacklist holds the regions of high and inconsistent signal in the same compendium.

* Source: https://github.com/fndemello/CUT-RUN_greenlist, main folder, downloaded on 2026-09-22.
* Licence: the repository carries no licence file. The paper is open access and the files are its
  supporting data, redistributed here unchanged and with the citation attached. Ask the authors
  before relying on this for anything beyond academic use.
* Reference: de Mello FN, Tahira AC, Berzoti-Coelho MG, Verjovski-Almeida S (2024). The CUT&RUN
  greenlist: genomic regions of consistent noise are effective normalizing factors for quantitative
  epigenome mapping. *Briefings in Bioinformatics* 25(2):bbad538. doi:10.1093/bib/bbad538

## What is not here

There is no greenlist for mm10 and none for mm39 CUT&Tag, because neither has been published. An
mm10 analysis needs the mm39 list lifted over, and a lifted file belongs here only with the chain
file and the date written in this document, marked as derived rather than original.

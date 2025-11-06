## Mobile element sequences
File `AluY_L1_SVA_cons_noPA.fa` contains consensus sequences of three major
human active mobile element families (AluY, L1HS and SVA) excluding the polyA
tails. It can be used with the `-T` option to add mobile element insertion
(MEI) information in the INFO field of somatic/mosaic variant calls.

## Non-centromeric regions

File `chm13v2.reg.nocen.bed` *excludes* the approximate locations of centromeric
satellite repeats and acrocentric short arms in CHM13 v2.0, which was *manually*
constructed based on the [official satellite][cen-sat] annotation, the
[DNA-BRNN][dna-brnn] satellite annotation and the [minigraph pangenome
graph][HPRC-mg] from the HPRC year-1 data.

File `hs38.reg.nocen.bed` was constructed similarly from DNA-BRNN annotation and
excludes regions where minigraph alignment faded.

File `hs37.reg.nocen.bed` was constructed by running DNA-BRNN on GRCh37 (hs37d5.fa).

[cen-sat]: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.0.bed
[dna-brnn]: https://github.com/lh3/dna-nn
[HPRC-mg]: https://zenodo.org/records/10693675
[zenodo]: https://zenodo.org/records/10963019

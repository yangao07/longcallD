* Tumour purity $\alpha$: est. -> max. -> est. ... -> converge
* Genotyping & Phasing & Haplotyping with EM
  * Initial haplotype (0|1, 1|0, 0|0, 1|1)
  * Assign each read to one allele (0 or 1)
  * Phasing based on assigned alleles
  * Startover based on refined haplotypes

* SV EM
  * Initial genotype: high-qual/confidence SVs
  * Assign nearby SV with undetermined boundaries/low-qual/confidence to high-qual/confidence SVs
  * New determined SV if not assigned to any existing one and pass filters
  * Update genotype and startover

* Multiple sequence alignment in variant calling?

* VCF + variant.fa: no long insertion seq. in VCF

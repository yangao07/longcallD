TODO:
- [ ] Rescue reads without hap/Remove reads with ambiguous hap
- [ ] Rescue SNPs in noisy region/Remove SNPs with low-confidence read evidence (T,T 1|2)
- [x] Phase set
- [ ] Stitch phase sets/bam chunks
- [ ] Link supp and primary
- [x] 2 alt bases (no ref base)
- [x] I/D bases during hap-assignment
- [ ] Local optima
- [ ] Hom. SNPs
- [ ] Germline indels

- [ ] Only phase bam with/without input VCF

---

* Tumour purity $\alpha$: est. -> max. -> est. ... -> converge
* Genotyping & Phasing & Haplotyping with EM
  * Initial haplotype (0|1, 1|0, 0|0, 1|1)
  * Assign each read to one allele (0 or 1)
  * Phasing based on assigned alleles
  * Startover based on refined haplotypes

* SV EM/iteration
  * Initial genotype: high-qual/confidence SVs
  * Assign nearby SV with undetermined boundaries/low-qual/confidence to high-qual/confidence SVs
  * New determined SV if not assigned to any existing one and pass filters
  * Update genotype and startover

* Multiple sequence alignment in variant calling?

* VCF + variant.fa: no long insertion seq. in VCF

* Joint alignment phasing and variant calling
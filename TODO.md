TODO:
- [ ] Multiple input BAMs (independent samples)
- [ ] Rescue reads without hap/Remove reads with ambiguous hap
- [ ] Rescue SNPs in noisy region/Remove SNPs with low-confidence read evidence (T,T 1|2)
- [x] Phase set
- [ ] Phase quality (FORMAT:PQ)
- [ ] Stitch phase sets/bam chunks
- [ ] Link supp and primary
- [x] 2 alt bases (no ref base)
- [x] I/D bases during hap-assignment
- [ ] Local optima
- [ ] Hom. SNPs
- [ ] Germline indels
- [ ] MSA for noisy regions
- [ ] Phasing using both SNPs and indels
- [ ] Rescue SNPs near gaps, they were considered as low-qual and ignored in the first round
- [ ] Large insertions -> re-alignment for clipping sequence: >30%(>2) Large INS around a candidate long INS/DEL, triger re-alignment (including clipping sequence)
- [ ] Only phase bam with/without input VCF
- [ ] Output assembly?
- [ ] LV instead of edlib for short seq.?
- [ ] Edlib: return path only if ED > 0?
- [ ] No candidate SVs (count > 2): do MSA (abPOA), use phase-based consensus as candidate
- [ ] Around SV: check left/right for longer match, extend/re-align the other side
- [ ] Annotation of inserted sequence? L1/Alu, etc.
- [ ] Weighted score during co-calling & phasing, variant-wise weight: base/mapping Q, noisy density, around large clipping
  - Low-weight variant contribute less than high-weight variants during phasing
- [ ] Coverage-based analysis to filter out wrong-mapping reads
  - With higher coverage
  - Large number of SNPs/indels or long insertions/deletions
  - Not in a repeat region
  - Non-repeat/homopolyer XID density
- [ ] BAM and Ref.FA not match
- [ ] Learn sequencing error profile from BAM
- [ ] Methylation
- [ ] Merge overlapping variants, i.e., I followed by D, caused by make_var_from_MSA_two_cons
- [ ] Split whole genome into chunks, use BAM index to load each chunk within each thread, instead of pre-loading N (i.e.,64) chunks
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



* DEBUG sites
  * 18925699 PhaseSet
  * 18923692 No HAP
  * 18923837 No HAP
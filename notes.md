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

```pesudocode
MQ_per read: skip or not

For each read (in a block, e.g., 100kb)
  For each mismatch sites (X in cigar)
    skip if BQ < thres
    If site not exist in site_to_base:
      Insert site (bin_tree_ins)
    Update site_to_base
    // site_to_base: bases, base_cnts // 0:ref, 1-n:alt
    //               n_bases, m_bases
    Update read_to_site
    // read_to_site: 0/1/2..: ref/alt1/alt2..
    //               1 bit if only max. 1 alt allele.
    //               2 or more bits when necessarg
Filtering site_to_base: site_used_for_phasing
For each read in read_to_site:
  If first read
    Assign with H1
    Pre_read = read
    Pre_hap = H1
    continue
  Compare with Pre_read
  If similarity > Thres:
    Assign with Pre_hap
  Else
    Assign with 3-Pre_hap
```

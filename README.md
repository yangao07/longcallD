<!-- # LongcallD: local-haplotagging-based small and structural variant calling -->
<!-- [![Latest Release](https://img.shields.io/github/v/tag/yangao07/longcalld?label=Release)](https://github.com/yangao07/longcallD/releases/latest) -->
[![Latest Release](https://img.shields.io/github/release/yangao07/longcallD.svg?label=Release)](https://github.com/yangao07/longcallD/releases/latest)
[![Github All Releases](https://img.shields.io/github/downloads/yangao07/longcallD/total.svg?label=Download)](https://github.com/yangao07/longcallD/releases)
[![Bioconda Version](https://img.shields.io/conda/vn/bioconda/longcallD.svg?style=flag&label=Bioconda)](https://anaconda.org/bioconda/longcalld)
[![Bioconda Install](https://img.shields.io/conda/dn/bioconda/longcallD.svg?style=flag&label=Conda-install)](https://anaconda.org/bioconda/longcalld)
[![C/C++ CI](https://github.com/yangao07/longcallD/actions/workflows/linux-CI.yml/badge.svg)](https://github.com/yangao07/longcallD/actions/workflows/linux-CI.yml)
[![C/C++ CI](https://github.com/yangao07/longcallD/actions/workflows/macos-CI.yml/badge.svg)](https://github.com/yangao07/longcallD/actions/workflows/macos-CI.yml)
[![License](https://img.shields.io/badge/License-MIT-black.svg)](https://github.com/yangao07/longcallD/blob/main/LICENSE)
<!-- [![Published in Bioinformatics](https://img.shields.io/badge/Published%20in-Bioinformatics-blue.svg)](https://dx.doi.org/10.1093/bioinformatics/btaa963) -->
<!-- [![GitHub Issues](https://img.shields.io/github/issues/yangao07/longcallD.svg?label=Issues)](https://github.com/yangao07/longcallD/issues) -->

## Updates (release v0.0.7)

* sort reads internally before processing to fix potential inconsistency when multiple input BAM/CRAM files are provided
<!-- * Fix corrupted VCF output in v0.0.5 -->
<!-- * Fix missing MEI header in VCF output -->
<!-- * Improved run time and memory usage (especially when mosaic variant calling enabled) -->
<!-- * Add `--input-is-list` and `-X` to support multiple input BAM/CRAM files of the same sample for variant calling -->


## Getting Started
```sh
# Download pre-built executables and test data (recommended)
# Linux-x64
wget https://github.com/yangao07/longcallD/releases/download/v0.0.7/longcallD-v0.0.7_x64-linux.tar.gz
tar -zxvf longcallD-v0.0.7_x64-linux.tar.gz && cd longcallD-v0.0.7_x64-linux
# MacOS-arm64
wget https://github.com/yangao07/longcallD/releases/download/v0.0.7/longcallD-v0.0.7_arm64-macos.tar.gz
tar -zxvf longcallD-v0.0.7_arm64-macos.tar.gz && cd longcallD-v0.0.7_arm64-macos

# PacBio HiFi reads
./longcallD call ./test_data/chr11_2M.fa ./test_data/HG002_chr11_hifi_test.bam --hifi > HG002_hifi_test.vcf
# Oxford Nanopore reads
./longcallD call ./test_data/chr11_2M.fa ./test_data/HG002_chr11_ont_test.bam --ont > HG002_ont_test.vcf
```
<!-- # man page for detailed command line options
man ./longcallD.1
``` -->

## Table of Contents
- [Updates (release v0.0.7)](#updates-release-v007)
- [Getting Started](#getting-started)
- [Table of Contents](#table-of-contents)
- [Introduction](#introduction)
- [Installation](#installation)
  - [Pre-built executables (recommended)](#pre-built-executables-recommended)
  - [Bioconda](#bioconda)
  - [Build from source](#build-from-source)
- [Usage](#usage)
  - [Variant calling with PacBio HiFi/Nanopore long reads](#variant-calling-with-pacbio-hifinanopore-long-reads)
  - [Multiple input BAM/CRAM files of the same sample](#multiple-input-bamcram-files-of-the-same-sample)
  - [Low-allele-fraction mosaic variant calling](#low-allele-fraction-mosaic-variant-calling)
  - [Region-specific variant calling](#region-specific-variant-calling)
  - [Output phased (\& refined) long-read BAM/CRAM](#output-phased--refined-long-read-bamcram)
  - [Variant calling from remote files](#variant-calling-from-remote-files)
- [Memory usage](#memory-usage)
- [Acknowledgements](#acknowledgements)
- [Contact](#contact)


## Introduction
LongcallD is a **local-haplotagging-based variant caller** designed for detecting small variants and structural variants (SVs)
using long-read sequencing data. It supports both **PacBio HiFi** and **Oxford Nanopore** reads.

LongcallD phases long reads into haplotypes using SNPs and small indels before calling SVs. It outputs phased variant calls in VCF format, including SNPs, small indels, and large SVs (currently only supporting insertions and deletions).

LongcallD (≥v0.0.5) can also call low-allele-fraction mosaic variant when `-s/--mosaic` is used.
Currently, only SNVs and large indels are supported, no mosaic small indels will be called.
Specifically, longcallD can sensitively identify mosaic mobile element insertions (MEIs).
Providing the annotation sequence of common mobile elements, i.e., Alu/L1/SVA, using `-T` is highly recommanded, which is included [here](https://github.com/yangao07/longcallD/tree/main/anno).
## Installation

### Pre-built executables (recommended)
**Linux-x64**
```
wget https://github.com/yangao07/longcallD/releases/download/v0.0.7/longcallD-v0.0.7_x64-linux.tar.gz
tar -zxvf longcallD-v0.0.7_x64-linux.tar.gz
```
**MacOS-arm64**
```
wget https://github.com/yangao07/longcallD/releases/download/v0.0.7/longcallD-v0.0.7_arm64-macos.tar.gz
tar -zxvf longcallD-v0.0.7_arm64-macos.tar.gz
```

**Linux-arm64/macOS-x64**

There is no pre-built executable for Linux-arm64 or macOS-x64, please try conda or build from source.

### Bioconda
**For Linux and macOS**
```
conda install -c bioconda longcalld
```

### Build from source
To compile longcallD from source, ensure you have **GCC/clang(9.0+)** and **zlib/libbz2/liblzma/libcurl** (for htslib) installed. 
```
wget https://github.com/yangao07/longcallD/releases/download/v0.0.7/longcallD-v0.0.7.tar.gz
tar -zxvf longcallD-v0.0.7.tar.gz
cd longcallD-v0.0.7; make
```

## Usage
LongcallD requires a **reference genome (FASTA)** and a **long-read BAM/CRAM** file as inputs. It outputs **phased variant calls in VCF format**.
### Variant calling with PacBio HiFi/Nanopore long reads
```
longcallD call -t16 ref.fa hifi.bam > hifi.vcf         # default for PacBio HiFi reads (--hifi)
longcallD call -t16 ref.fa ont.bam --ont > ont.vcf     # for ONT reads
```

### Multiple input BAM/CRAM files of the same sample
You can provide multiple BAM/CRAM files of the same sample for variant calling using `--input-is-list` or `-X`:
```
longcallD call -t16 --input-is-list ref.fa bam_list.txt > sample.vcf
# where bam_list.txt contains:
# sample_part1.bam
# sample_part2.bam
# sample_part3.bam
```
or
```
longcallD call -t16 ref.fa sample_part1.bam -X sample_part2.bam -X sample_part3.bam > sample.vcf
```

### Low-allele-fraction mosaic variant calling
With `-s`, longcallD will detect both germline and low-allele-fraction somatic/mosaic variants.

For each somatic/mosaic variant, a `SOMATIC` tag will be added to the INFO field in the output VCF.
```
longcallD call -s -t16 ref.fa hifi.bam > hifi.vcf
longcallD call -s -t16 ref.fa hifi.bam -T AluY_L1_SVA_cons_noPA.fa > hifi.vcf # add MEI information in INFO field
```

### Region-specific variant calling
LongcallD supports region-based variant calling, similar to `samtools view`.
```
longcallD call -t16 ref.fa hifi.bam chr11:10,229,956-10,256,221 > hifi_reg.vcf
longcallD call -t16 ref.fa hifi.bam chr11:10,229,956-10,256,221 chr12:10,576,356-10,583,438 > hifi_regs.vcf
longcallD call -t16 ref.fa hifi.bam --region-file reg.bed > hifi_regs.vcf
longcallD call -t16 ref.fa hifi.bam --autosome > hifi_autosome.vcf
```

### Output phased (& refined) long-read BAM/CRAM
LongcallD performs read phasing during variant calling and can outputs phased long reads in BAM/CRAM.

With `--refine-aln`, it can further output refined read alignment based on multiple sequence alignment within each haplotype, which is especially useful for low-complexity regions like homopolymers and tandem repeats.
```
longcallD call -t16 ref.fa hifi.bam --hifi -b hifi_phased.bam > hifi.vcf                  # output phased HiFi reads (BAM tag: HP & PS)
longcallD call -t16 ref.fa ont.bam --ont --refine-aln -b ont_phased_refined.bam > ont.vcf # output phased & refined ONT reads (BAM tag: HP & PS)
```
### Variant calling from remote files
```
ref=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
bam=https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam
longcallD call -t16 $ref $bam chr11:10,229,956-10,256,221 chr12:10,576,356-10,583,438 > hifi_regs.vcf
```

## Memory usage
Because longcallD performs multiple sequence alignment and re-alignment, both are memory-intensive tasks, it generally requires more memory than other variant callers. Peak memory consumption depends primarily on the number of threads (`-t/--threads`), sequencing coverage, and read length. For human whole-genome datasets at ~40× coverage, longcallD typically uses about **1 GB** (**HiFi**) or **2 GB** (**ONT R10**) of memory _per thread_ for germline variant calling.

Memory usage and runtime increase further when mosaic variant calling is enabled.

If you encounter memory constraints, you may restrict processing to specific genomic regions using `--region-file`. A region list for the human genome that excludes centromeres is available [here](https://github.com/yangao07/longcallD/blob/main/anno/).

## Acknowledgements
LongcallD is dependent on the following libraries, we are grateful to all the developers/maintainers:

* [htslib](https://github.com/samtools/htslib): read/write BAM/CRAM/VCF
* [abPOA](https://github.com/yangao07/abPOA): consensus calling
* [WFA](https://github.com/smarco/WFA2-lib): pairwise alignment
* [edlib](https://github.com/Martinsos/edlib): fast sequence similarity calculation
* [cgranges](https://github.com/lh3/cgranges): interval operations
* [sdust](https://github.com/lh3/sdust): identify low-complexity regions

## Contact

For any questions or support, please contact:

* Yan Gao yangao@ds.dfci.harvard.edu

* Heng Li hli@jimmy.harvard.edu

<!-- # LongcallD: local-haplotagging-based small and structural variant calling -->

<!-- [![Latest Release](https://img.shields.io/github/release/yangao07/longcallD.svg?label=Release)](https://github.com/yangao07/longcallD/releases/latest) -->
[![Github All Releases](https://img.shields.io/github/downloads/yangao07/longcallD/total.svg?label=Download)](https://github.com/yangao07/longcallD/releases)
[![Bioconda Install](https://img.shields.io/conda/dn/bioconda/longcallD.svg?style=flag&label=Bioconda)](https://anaconda.org/bioconda/longcalld)
[![C/C++ CI](https://github.com/yangao07/longcallD/actions/workflows/linux-CI.yml/badge.svg)](https://github.com/yangao07/longcallD/actions/workflows/linux-CI.yml)
[![C/C++ CI](https://github.com/yangao07/longcallD/actions/workflows/macos-CI.yml/badge.svg)](https://github.com/yangao07/longcallD/actions/workflows/macos-CI.yml)
[![License](https://img.shields.io/badge/License-MIT-black.svg)](https://github.com/yangao07/longcallD/blob/main/LICENSE)
<!-- [![Published in Bioinformatics](https://img.shields.io/badge/Published%20in-Bioinformatics-blue.svg)](https://dx.doi.org/10.1093/bioinformatics/btaa963) -->
<!-- [![GitHub Issues](https://img.shields.io/github/issues/yangao07/longcallD.svg?label=Issues)](https://github.com/yangao07/longcallD/issues) -->

## Updates (pre-release v0.0.5)

* Add -Oz for compressed VCF output
* Add INFO/TSD for Alu/L1/SVA INS/DEL
* Add support for BAM/cs tag
* Fix compiling in macOS-x64

## Getting Started
```sh
# Download pre-built executables and test data (recommended)
# Linux-x64
wget https://github.com/yangao07/longcallD/releases/download/v0.0.4/longcallD-v0.0.4_x64-linux.tar.gz
tar -zxvf longcallD-v0.0.4_x64-linux.tar.gz && cd longcallD-v0.0.4_x64-linux
# MacOS-arm64
wget https://github.com/yangao07/longcallD/releases/download/v0.0.4/longcallD-v0.0.4_arm64-macos.tar.gz
tar -zxvf longcallD-v0.0.4_arm64-macos.tar.gz && cd longcallD-v0.0.4_arm64-macos

# PacBio HiFi reads
./longcallD call ./test_data/chr11_2M.fa ./test_data/HG002_chr11_hifi_test.bam --hifi > HG002_hifi_test.vcf
# Oxford Nanopore reads
./longcallD call ./test_data/chr11_2M.fa ./test_data/HG002_chr11_ont_test.bam --ont > HG002_ont_test.vcf
```
<!-- # man page for detailed command line options
man ./longcallD.1
``` -->

## Table of Contents
- [Updates (pre-release v0.0.5)](#updates-pre-release-v005)
- [Getting Started](#getting-started)
- [Table of Contents](#table-of-contents)
- [Introduction](#introduction)
- [Installation](#installation)
  - [Pre-built executables (recommended)](#pre-built-executables-recommended)
  - [Bioconda](#bioconda)
  - [Build from source](#build-from-source)
- [Usage](#usage)
  - [Variant calling with HiFi/Nanopore long reads](#variant-calling-with-hifinanopore-long-reads)
  - [Region-specific variant calling](#region-specific-variant-calling)
  - [Variant calling and output phased long reads](#variant-calling-and-output-phased-long-reads)
  - [Variant calling from remote files](#variant-calling-from-remote-files)
- [Acknowledgements](#acknowledgements)
- [Contact](#contact)


## Introduction
LongcallD is a **local-haplotagging-based variant caller** designed for detecting small variants and structural variants (SVs)
using long-read sequencing data. It supports both **PacBio HiFi** and **Oxford Nanopore** reads.

LongcallD phases long reads into haplotypes using SNPs and small indels before calling SVs. It outputs phased variant calls in VCF format, including SNPs, small indels, and large SVs (currently only supporting insertions and deletions).
## Installation

### Pre-built executables (recommended)
**Linux-x64**
```
wget https://github.com/yangao07/longcallD/releases/download/v0.0.4/longcallD-v0.0.4_x64-linux.tar.gz
tar -zxvf longcallD-v0.0.4_x64-linux.tar.gz
```
**MacOS-arm64**
```
wget https://github.com/yangao07/longcallD/releases/download/v0.0.4/longcallD-v0.0.4_arm64-macos.tar.gz
tar -zxvf longcallD-v0.0.4_arm64-macos.tar.gz
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
It is recommended to use the [latest release](https://github.com/yangao07/longcallD/releases).
```
wget https://github.com/yangao07/longcallD/releases/download/v0.0.4/longcallD-v0.0.4.tar.gz
tar -zxvf longcallD-v0.0.4.tar.gz
cd longcallD-v0.0.4; make
```

## Usage
LongcallD requires a **reference genome (FASTA)** and a **long-read BAM/CRAM** file as inputs. It outputs **phased variant calls in VCF format**.
### Variant calling with HiFi/Nanopore long reads
```
longcallD call -t16 ref.fa hifi.bam > hifi.vcf         # default for PacBio HiFi reads (--hifi)
longcallD call -t16 ref.fa ont.bam --ont > ont.vcf     # for ONT reads
```

### Region-specific variant calling
LongcallD supports region-based variant calling, similar to `samtools view`.
```
longcallD call -t16 ref.fa hifi.bam chr11:10,229,956-10,256,221 > hifi_reg.vcf
longcallD call -t16 ref.fa hifi.bam chr11:10,229,956-10,256,221 chr12:10,576,356-10,583,438 > hifi_regs.vcf
longcallD call -t16 ref.fa hifi.bam --region-file reg.bed > hifi_regs.vcf
longcallD call -t16 ref.fa hifi.bam --autosome > hifi_autosome.vcf
```

### Variant calling and output phased long reads
```
longcallD call -t16 ref.fa hifi.bam --hifi -b hifi_phased.bam > hifi.vcf  # output phased HiFi reads (BAM tag: HP & PS)
longcallD call -t16 ref.fa ont.bam --ont -b ont_phased.bam > ont.vcf      # output phased ONT reads (BAM tag: HP & PS)
```
### Variant calling from remote files
```
ref=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
bam=https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam
longcallD call -t16 $ref $bam chr11:10,229,956-10,256,221 chr12:10,576,356-10,583,438 > hifi_regs.vcf
```

## Acknowledgements
LongcallD is dependent on the following libraries, we are grateful to all the developers/maintainers:

* [htslib](https://github.com/samtools/htslib): read/write BAM/CRAM/VCF
* [abPOA](https://github.com/yangao07/abPOA): consensus calling
* [WFA](https://github.com/smarco/WFA2-lib): pairwise alignment
* [cgranges](https://github.com/lh3/cgranges): interval operations
* [sdust](https://github.com/lh3/sdust): identify low-complexity regions

## Contact

For any questions or support, please contact:

* Yan Gao yangao@ds.dfci.harvard.edu

* Heng Li hli@jimmy.harvard.edu

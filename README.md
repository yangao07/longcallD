# LongcallD: local-haplotagging-based small and structural variant calling

[![Latest Release](https://img.shields.io/github/release/yangao07/longcallD.svg?label=Release)](https://github.com/yangao07/longcallD/releases/latest)
[![C/C++ CI](https://github.com/yangao07/longcallD/actions/workflows/linux-CI.yml/badge.svg)](https://github.com/yangao07/longcallD/actions/workflows/linux-CI.yml)
[![C/C++ CI](https://github.com/yangao07/longcallD/actions/workflows/macos-CI.yml/badge.svg)](https://github.com/yangao07/longcallD/actions/workflows/macos-CI.yml)
[![License](https://img.shields.io/badge/License-MIT-black.svg)](https://github.com/yangao07/longcallD/blob/main/LICENSE)
<!-- [![Github All Releases](https://img.shields.io/github/downloads/yangao07/longcallD/total.svg?label=Download)](https://github.com/yangao07/longcallD/releases) -->
<!-- [![BioConda Install](https://img.shields.io/conda/dn/bioconda/longcallD.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/longcallD) -->
<!-- [![Published in Bioinformatics](https://img.shields.io/badge/Published%20in-Bioinformatics-blue.svg)](https://dx.doi.org/10.1093/bioinformatics/btaa963) -->
<!-- [![GitHub Issues](https://img.shields.io/github/issues/yangao07/longcallD.svg?label=Issues)](https://github.com/yangao07/longcallD/issues) -->
## Updates (v0.0.1)

* Pre-release v0.0.1

## Getting Started
```sh
git clone --recursive https://github.com/yangao07/longcallD
cd longcallD && make
# PacBio HiFi reads
./bin/longcallD call ref.fa hifi.bam --hifi > hifi.vcf
# Oxford Nanopore reads
./bin/longcallD call ref.fa ont.bam --ont > ont.vcf
```
<!-- # man page for detailed command line options
man ./longcallD.1
``` -->

## Table of Contents
- [LongcallD: local-haplotagging-based small and structural variant calling](#longcalld-local-haplotagging-based-small-and-structural-variant-calling)
  - [Updates (v0.0.1)](#updates-v001)
  - [Getting Started](#getting-started)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Installation](#installation)
    - [Pre-built binary executable file for Linux or MacOS (recommended)](#pre-built-binary-executable-file-for-linux-or-macos-recommended)
    - [Building longcallD from source files](#building-longcalld-from-source-files)
  - [General Usage](#general-usage)
    - [Variant calling from HiFi/ONT long reads](#variant-calling-from-hifiont-long-reads)
    - [Variant calling and phasing long reads](#variant-calling-and-phasing-long-reads)
    - [Variant calling from remote files](#variant-calling-from-remote-files)
  - [Contact](#contact)

## Introduction
LongcallD is a local-haplotagging-based small and structural variant (SV) caller using long reads.
It works with both PacBio HiFi and Oxford Nanopore genomic long reads. LongcallD first uses SNPs
and small indels to phase long reads into two haplotypes, then call SVs from
phased long reads.

LongcallD outputs phased variant calls in VCF file, containing SNPs, small indels, and large SVs. For SVs, it currently only outputs insertions and deletions.

## Installation

### Pre-built binary executable file for Linux or MacOS (recommended)
For Linux:
```
wget https://github.com/yangao07/longcallD/releases/download/v0.0.1/longcallD-v0.0.1_x64-linux.tar.gz
tar -zxvf longcallD-v0.0.1_x64-linux.tar.gz
```
For MacOS:
```
wget https://github.com/yangao07/longcallD/releases/download/v0.0.1/longcallD-v0.0.1_arm64-macos.tar.gz
tar -zxvf longcallD-v0.0.1_arm64-macos.tar.gz
```

### Building longcallD from source files
You can also build longcallD from source files in Linux or MacOS.
Make sure you have gcc/clang and zlib installed before compiling.
It is recommended to download the [latest release](https://github.com/yangao07/longcallD/releases).
```
wget https://github.com/yangao07/longcallD/releases/download/v0.0.1/abPOA-v0.0.1.tar.gz
tar -zxvf longcallD-v0.0.1.tar.gz
cd longcallD-v0.0.1; make
```

## General Usage
LongcallD takes a reference genome FASTA (can be gzip'd) and a long-read BAM as input
and outputs phased variant calls in a VCF file.
LongcallD seamlessly works with variant calling region(s), in the same way as `samtools view`.
### Variant calling from HiFi/ONT long reads
```
longcallD call -t16 ref.fa hifi.bam > hifi.vcf         # for PacBio HiFi reads, --hifi is default
longcallD call -t16 ref.fa ont.bam --ont > ont.vcf     # for ONT reads
longcallD call -t16 ref.fa hifi.bam chr11:10,229,956-10,256,221 > hifi_reg.vcf   # variant calling in a region
longcallD call -t16 ref.fa hifi.bam chr11:10,229,956-10,256,221 chr12:10,576,356-10,583,438 > hifi_regs.vcf  # variant calling in multiple regions
```
### Variant calling and phasing long reads
```
longcallD call -t16 ref.fa hifi.bam --hifi -b hifi_phased.bam > hifi.vcf  # output phased HiFi reads, with BAM tag HP & PS
longcallD call -t16 ref.fa ont.bam --ont -b ont_phased.bam > ont.vcf      # output phased ONT reads, with BAM tag HP & PS 
```
### Variant calling from remote files
```
ref=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
bam=https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam
longcallD call -t16 $ref $bam chr11:10,229,956-10,256,221 chr12:10,576,356-10,583,438 > hifi_regs.vcf
```

## Contact
Yan Gao yangao@ds.dfci.harvard.edu

Heng Li hli@jimmy.harvard.edu

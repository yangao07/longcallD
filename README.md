# LongcallD: local-haplotagging-based small and structural variant calling

[![Latest Release](https://img.shields.io/github/release/yangao07/longcallD.svg?label=Release)](https://github.com/yangao07/longcallD/releases/latest)
<!-- [![Github All Releases](https://img.shields.io/github/downloads/yangao07/longcallD/total.svg?label=Download)](https://github.com/yangao07/longcallD/releases) -->
<!-- [![BioConda Install](https://img.shields.io/conda/dn/bioconda/longcallD.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/longcallD) -->
<!-- [![Published in Bioinformatics](https://img.shields.io/badge/Published%20in-Bioinformatics-blue.svg)](https://dx.doi.org/10.1093/bioinformatics/btaa963) -->
<!-- [![GitHub Issues](https://img.shields.io/github/issues/yangao07/longcallD.svg?label=Issues)](https://github.com/yangao07/longcallD/issues) -->
[![C/C++ CI](https://github.com/yangao07/longcallD/actions/workflows/linux-CI.yml/badge.svg)](https://github.com/yangao07/longcallD/actions/workflows/linux-CI.yml)
[![C/C++ CI](https://github.com/yangao07/longcallD/actions/workflows/macos-CI.yml/badge.svg)](https://github.com/yangao07/longcallD/actions/workflows/macos-CI.yml)
[![License](https://img.shields.io/badge/License-MIT-black.svg)](https://github.com/yangao07/longcallD/blob/main/LICENSE)
## Updates (v0.0.1)

* Pre-release v0.0.1

## <a name="started"></a>Getting Started
```sh
git clone --recursive https://github.com/yangao07/longcallD
cd longcallD && make
# PacBio HiFi reads
./bin/longcallD ref.fa hifi.bam --hifi > hifi.vcf
# Oxford Nanopore reads
./bin/longcallD ref.fa ont.bam --ont > ont.vcf
```
<!-- # man page for detailed command line options
man ./longcallD.1
``` -->

## Table of Contents

## Introduction
LongcallD is a local-haplotagging-based small and structural variant (SV) caller using long reads.
It works with both PacBio HiFi and Oxford Nanopore genomic long reads. LongcallD first uses SNPs
and small indels to phase long reads into two haplotypes, then call SVs from
phased long reads.

LongcallD outputs phased variant calls in VCF file, containing SNPs, small indels, and large SVs. For SVs, it currently only outputs insertions and deletions.

## Installation

### Pre-built binary executable file for Linux/Unix or MacOS (recommended)
For linux:
```
wget https://github.com/yangao07/longcallD/releases/download/v0.0.1/longcallD-v0.0.1_x64-linux.tar.gz
tar -zxvf longcallD-v0.0.1_x64-linux.tar.gz
```
or for macos:
```
wget https://github.com/yangao07/longcallD/releases/download/v0.0.1/longcallD-v0.0.1_arm64-macos.tar.gz
tar -zxvf longcallD-v0.0.1_arm64-macos.tar.gz
```

### Building longcallD from source files
You can also build longcallD from source files. 
Make sure you have gcc (>=6.4.0) and zlib installed before compiling.
It is recommended to download the [latest release](https://github.com/yangao07/longcallD/releases).
```
wget https://github.com/yangao07/longcallD/releases/download/v0.0.1/abPOA-v0.0.1.tar.gz
tar -zxvf longcallD-v0.0.1.tar.gz
cd longcallD-v0.0.1; make
```

## General usage
### Variant calling from HiFi/ONT long reads
```
longcallD call -t16 ref.fa hifi.bam --hifi > hifi.vcf  # for PacBio HiFi reads
longcallD call -t16 ref.fa ont.bam --ont > ont.vcf     # for ONT reads
```
### Variant calling and phasing long reads
```
longcallD call -t16 ref.fa hifi.bam --hifi -b hifi_phased.bam > hifi.vcf  # output phased HiFi reads, with BAM tag HP & PS
longcallD call -t16 ref.fa ont.bam --ont -b ont_phased.bam > ont.vcf      # output phased ONT reads, with BAM tag HP & PS 
```

## Contact
Yan Gao yangao@ds.dfci.harvard.edu

Heng Li hli@jimmy.harvard.edu

-# cnr-ibba/nf-treeseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.3.0](https://github.com/cnr-ibba/nf-treeseq/releases/tag/v0.3.0) - dev

Release of cnr-ibba/nf-treeseq, created with the [nf-core](https://nf-co.re/) template.
Added support for `threads` method for tree sequence inference, which does not require
ancestral allele information.

### `Added`

- use _major allele_ as ancestral allele
- add _CI_ and tests
- add _nf-core_ `TEMPLATE` branch
- add `nf-schema` plugin for parameter validation
- support for _VCF_ input through sample sheet
- define subworkflow to deal with VCF and PLINK input files
- track all software versions used in the pipeline, including tskit related packages
- pack tree output using `tszip`
- add support for _institutional configuration_
- add `threads` as ancestor method with demography file support and configurable
  `mutation_rate`/`query_interval` parameters
- support for `threads` specific options
- add `threads/infer` and `threads/convert` local modules
- add `threads` local subworkflow
- add `plink2/vcf` module to convert VCF into `pgen` format
- implement _recombination rate_ parameter through pipeline
- add metadata to trees inferred with others method (threads)
- provide constant `thread` _Ne_ using parameter
- update _docs_ and README files

### `Fixed`

- move `tskit` related scripts to `tskit` docker image
- change the default `tsdate_method` to `variational_gamma` in `tsinfer` approaches
- increase default memory allocations across process labels
- revise ancient parameters (condense `mutation_rate` and `recombination_rate` into
  single parameters for both `tsinfer` and `threads` approaches)
- annotate trees relying on VCF sample order
- filter out singletons SNPs in threads method

### `Dependencies`

- update `tskit` to `1.0.0b3` relying on custom docker images
- track pipeline dependencies using `wave` and dockerfiles (`threads`)

### `Deprecated`

- remove `nf-validation` plugin
- remove old validation library
- drop _git lfs_ support
- remove deprecated `dumpsoftwareversion` module

## [v0.2.1](https://github.com/cnr-ibba/nf-treeseq/releases/tag/v0.2.1) - 2024-09-05

### `Dependencies`

- update `tskit` based images

## [v0.2.0](https://github.com/cnr-ibba/nf-treeseq/releases/tag/v0.2.0) - 2024-09-03

The pipeline now supports the _compara_ method for determining ancestral alleles, in addition to the reference genome and est-sfs methods. The tsinfer step is now modularized based on the chosen ancestral inference method.

### `Added`

- add `LICENSE` file with MIT license
- add `tsinfer_compara` local module
- add `compara` local subworkflow

### `Fixed`

- rename `pos` to `position` in `estsfs` related modules
- refactor `create_tstree` function in `helper.py` to handle different ancestral
  allele methods
- enforce validation of `ancestor_method` parameter

### `Dependencies`

- update `tskitetude` related scripts
- update `tskit` based images
- update `README.md` to reflect new changes
- update `conf/modules.config` to reflect new module names

## [v0.1.0](https://github.com/cnr-ibba/nf-treeseq/releases/tag/v0.1.0) - 2024-07-24

First release of nf-treeseq pipeline, which comes out from [bunop/TSKITetude](https://github.com/bunop/TSKITetude)
project. Set up the core structure and implemented
tree sequence inference using est-sfs or reference-based approach. Processing
start from a PLINK binary file and a reference genome where focal samples are
selected, which ancenstor alleles are determined using est-sfs or the reference genome.
Data are then converted to VCF format, imputed and phased using beagle, and finally
tree sequences are inferred using tsinfer through custom scripts.

### `Added`

- implement `est-sfs` and `reference` _local_ subworkflows for tree sequence inference
- add custom helper scripts in `bin/` for various tasks
- add `cnr-ibba/est-sfs` module
- add `estsfs_input` and `estsfs_output` local modules
- add `tsinfer_estsfs` and `tsinfer_reference` local modules
- add `plink_subset` local module
- add `bcftools/merge` nf-core module
- add `bcftools/norm` nf-core module
- add `bcftools/reheader` nf-core module
- add `bcftools/split` nf-core module
- add `beagle5/beagle` nf-core module
- add `plink/recode` nf-core module
- add `samtools/faidx` nf-core module
- add `tabix/tabix` nf-core module
- add test files in `tests/` folder

### `Dependencies`

- git lfs is required for large test files

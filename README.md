# nf-treeseq

A Nextflow pipeline for generating Tree Sequences from PLINK and VCF files.

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/cnr-ibba/nf-treeseq)
[![GitHub Actions CI Status](https://github.com/cnr-ibba/nf-treeseq/actions/workflows/nf-test.yml/badge.svg)](https://github.com/cnr-ibba/nf-treeseq/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/cnr-ibba/nf-treeseq/actions/workflows/linting.yml/badge.svg)](https://github.com/cnr-ibba/nf-treeseq/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/cnr-ibba/nf-treeseq)

## Introduction

**cnr-ibba/nf-treeseq** is a bioinformatics pipeline that infers phylogenetic tree
sequences from genotype data in PLINK or VCF format. The pipeline performs quality
control, format conversion, variant normalization, imputation and phasing using
Beagle, and finally generates tree sequences using `tsinfer` from the `tskit` library.
It supports multiple methods for determining ancestral alleles using `tsinfer`,
including _reference-based_, _frequency-based (major allele)_, estimation via
`est-sfs` with _outgroup samples_, or custom _user-provided ancestral states_.
In addition, it is also possible to skip the `tsinfer` inference steps and create
tree sequences using `threads`.

<!-- TODO cnr-ibba: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/guidelines/graphic_design/workflow_diagrams#examples for examples.   -->

**Main pipeline steps:**

1. **Input Processing**: Reads PLINK binary files (.bed/.bim/.fam) and applies
   sample filtering and SNP quality control (est-sfs/tsinfer approach). Start from VCF
   files for all the other methods.
2. **Variant Normalization**: Normalizes ALT/REF alleles and validates chromosome sizes against reference genome using `bcftools`
3. **Imputation & Phasing**: Imputes missing genotypes and phases haplotypes using Beagle
4. **Ancestral Allele Inference**: Determines ancestral alleles via one of 4 methods for the `tsinfer` approach (reference, major allele, `est-sfs`, or custom), or skips ancestral allele inference entirely when using the `threads` approach
5. **Tree Sequence Generation**: Creates tree sequence files using `tsinfer` with inferred ancestral states or `threads` directly from phased haplotypes.

## Background

This pipeline is designed to infer _Tree Sequences_ from population genetic data,
providing a compact and efficient representation of genealogical relationships
across the genome. Tree sequences, as implemented by the `tskit` library can
be inferred by `tsinfer` and `threads`, and offer a powerful framework for
population genomic analyses
by encoding coalescent histories along the genome.

The pipeline accepts PLINK binary genotype files (.bed/.bim/.fam format) or VCF files
as input, where all samples of interest are contained within a single file set.

A key requirement for `tsinfer` is the specification of ancestral alleles at
each variant site. The pipeline offers flexible approaches to meet this requirement,
accommodating different data scenarios and biological questions
(see ["About Ancestral Alleles"](#about-ancestral-alleles) section below).

### About Ancestral Alleles

`tsinfer` requires ancestral alleles to generate tree sequence files. Currently,
the pipeline supports four different methods for determining ancestral alleles:

1. **Using the reference allele**: The REF allele in the VCF file is used as
   the ancestral allele (default method).
2. **Using the major allele**: The most frequent allele in the dataset is used
   as the ancestral allele.
3. **Using `est-sfs`**: This method estimates the site frequency spectrum and
   infers ancestral alleles. It requires the presence of outgroup samples
   (ancestral to the rest of the data) in the dataset.
4. **Using `custom`**: This method requires an additional CSV file containing
   the ancestral alleles.

### About the `threads` method

In addition to `tsinfer`, the pipeline also supports generating tree sequences
using `threads`. This method does not require ancestral alleles and directly
creates tree sequences from phased haplotypes. This can be useful when ancestral
alleles are not available or when a simpler approach is desired.

## Getting the Pipeline

You can obtain this pipeline by cloning the GitHub repository:

```bash
git clone cnr-ibba/nf-treeseq
```

Alternatively, you can use the `nextflow pull` command:

```bash
nextflow pull cnr-ibba/nf-treeseq
```

For more information on installing and running Nextflow pipelines, including
dealing with revisions, refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html).

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation)
> on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

### Quick Start

The pipeline can be run with parameters provided either via command line or
(recommended) via a JSON parameters file. A minimal command to run the pipeline
looks like this:

```bash
nextflow run cnr-ibba/nf-treeseq \
   -profile singularity \
   -params-file params.json \
   --outdir results
```

Where:

- `-profile singularity` specifies the execution environment (alternative: `docker` or/and institutional profiles)
- `-params-file params.json` points to your parameter configuration file
- `--outdir results` sets the output directory for results

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option.
> Custom config files including those provided by the `-c` Nextflow option can
> be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

### Preparing Input Files

The pipeline requires PLINK binary format genotype files as input.
Ensure you have the following files ready:

- **PLINK binary files**: `.bed`, `.bim`, and `.fam` files with the same prefix (e.g., `mydata.bed`, `mydata.bim`, `mydata.fam`) for the `tsinfer/est-sfs` approach.
- **VCF files with indexes**: One or more VCF with index file containing all samples for all other approaches
- **Reference genome**: A FASTA file (optionally compressed) for allele normalization
- **Sample/Population file**: Tab-separated text file **without a header row** containing `FID`, `IID`, like required by PLINK to specify which samples to retain and their population assignment. `FID` and `IID` must match the corresponding IDs in your PLINK/VCF data.
- **Outgroup files** (optional, for `est-sfs` method): One to three TSV files with `FID` and `IID` columns identifying outgroup samples
- **Custom ancestral allele file** (optional, for `custom` method): CSV file with ancestral allele information

#### PLINK samplesheet

Prepare a samplesheet like the following for PLINK data:

```csv
sample,plink_bfile
test,tests/test_dataset
```

Each row represents a PLINK binary file prefix (i.e., without the `.bed/.bim/.fam`
extensions) while the sample is an identifier for that dataset. This file can be
used only for `tsinfer/est-sfs` approach.

#### VCF samplesheet

Prepare a samplesheet like the following for VCF data:

```csv
sample,vcf,index
test,tests/test_dataset.vcf.gz,tests/test_dataset.vcf.gz.tbi
test2,tests/test_dataset2.vcf.gz,tests/test_dataset2.vcf.gz.tbi
```

### Creating a Parameters File

While all parameters can be passed via the command line, it is **strongly recommended**
to use a JSON parameters file for reproducibility and clarity.
Create a file named `params.json` with at minimum the following required parameters:

```json
{
  "input": "path/to/samplesheet.csv",
  "sample2fid": "path/to/sample2fid.tsv",
  "genome": "path/to/reference_genome.fasta.gz",
  "ancestor_method": "reference",
  "outdir": "results"
}
```

#### Core Parameters Explained:

- **`input`** (required): Path to the samplesheet CSV file
  - For `tsinfer/est-sfs` approach: PLINK samplesheet with `sample` and `plink_bfile` columns
  - For all other approaches: VCF samplesheet with `sample`, `vcf`, and `index` columns

- **`sample2fid`** (required): TSV file with `FID` and `IID` columns to filter samples
  and annotate populations

- **`genome`** (required): Reference genome FASTA file for variant normalization
  - Can be compressed (`.fasta.gz`) or uncompressed (`.fasta`)
  - Must match the genome build used for genotyping

- **`ancestor_method`** (required): Method for determining ancestral alleles
  - Options: `"reference"` (default), `"major"`, `"est-sfs"`, `"custom"`, `"threads"`
  - See [Specifying Ancestral Alleles](#specifying-ancestral-alleles) section below

- **`outdir`** (optional): Output directory for results (default: `"results"`)

#### Optional Quality Control Parameters:

- **`plink_species`**: PLINK species-specific options valid for your data
  (`est-sfs/tsinfer` approach only and threaded VCF-to-PLINK conversion)

- **`plink_geno`**: Maximum missing rate per SNP (default: `0.1`)
  ```json
  "plink_geno": 0.05
  ```
  SNPs with missing rate above this threshold will be excluded

### Specifying Ancestral Alleles

The pipeline requires ancestral alleles to generate tree sequences. At least one of the following methods must be used to infer ancestral alleles:

#### 1. Using the Reference Genome

To use the reference genome for inferring ancestral alleles, simply set the `ancestor_method` flag:

```json
{
  "ancestor_method": "reference"
}
```

This is the default method if the `ancestor_method` parameter is not specified.

#### 2 Using the Major Allele

To use the major allele for inferring ancestral alleles, set the `ancestor_method` flag to `major`:

```json
{
  "ancestor_method": "major"
}
```

#### 3. Using `est-sfs` to Infer Ancestral Alleles

To infer ancestral alleles using `est-sfs`, enable the `with_estsfs` flag and specify one or more outgroup sample files (TSV format with `FID` and `IID` columns). You can provide up to three outgroup files:

```json
{
  "ancestor_method": "est-sfs",
  "outgroup1": "<outgroup1 samples file>",
  "outgroup2": "<outgroup2 samples file>",
  "outgroup3": "<outgroup3 samples file>"
}
```

#### 4. Using `custom` to Infer Ancestral Alleles

To use `custom` for inferring ancestral alleles, provide a CSV file with the following format:

```csv
chrom,position,alleles,anc_allele
26,209049,A/G,C
26,268822,A/G,C
26,285471,A/G,G
26,361728,G/T,G
```

After generating the file, specify it using the `ancestor_file` parameter:

```json
{
  "ancestor_method": "custom",
  "ancestor_file": "<custom ancestor file>"
}
```

#### 5. Using `threads` Method

To generate tree sequences using the `threads` method, which does not require ancestral alleles, set the `ancestor_method` flag to `threads`:

```json
{
  "ancestor_method": "threads"
}
```

### Additional Parameters

Additional parameters can be set in the configuration file to control the pipeline or specify the output directory. To see all available options, run:

```bash
nextflow run cnr-ibba/nf-treeseq --help
```

For more advanced options, including hidden parameters:

```bash
nextflow run cnr-ibba/nf-treeseq --help --validationShowHiddenParams
```

## Running the Pipeline

Once your `params.json` file is configured, execute the pipeline:

```bash
nextflow run cnr-ibba/nf-treeseq \
   -profile singularity \
   -params-file params.json
```

To override specific parameters from the command line:

```bash
nextflow run cnr-ibba/nf-treeseq \
   -profile singularity \
   -params-file params.json \
   --plink_geno 0.05 \
   --outdir custom_output
```

To test the pipeline with the included test dataset:

```bash
nextflow run cnr-ibba/nf-treeseq \
   -profile test,singularity \
   --outdir test_results
```

## Credits

cnr-ibba/nf-treeseq was originally written by Paolo Cozzi.

We thank the following people for their extensive assistance in the development of this pipeline:

- Hannes Becher
- Arianna Manunza
- Jana Obšteter
- Gabriela Mafra Fortuna
- Gregor Gorianc
- Filippo Biscarini

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use cnr-ibba/nf-treeseq for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

# nf-treeseq

A Nextflow pipeline for generating Tree Sequences from PLINK and VCF files.

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![Run with Docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![Run with Singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Background

This pipeline is designed to infer _Tree Sequences_ from genotype data. It is currently tailored for PLINK genotype files, where all relevant samples are contained within a single file. The pipeline converts the PLINK file into a VCF file, corrects ALT/REF alleles, and checks chromosome sizes. It then uses Beagle to impute and phase any missing data before running `tsinfer` to create Tree Sequences from the VCF file.

### About Ancestral Alleles

`tsinfer` requires ancestral alleles to generate tree sequence files. Currently, the pipeline supports three methods for determining ancestral alleles:

1. **Using the reference genome**: The REF allele in the VCF file is used as the ancestral allele.
2. **Using `est-sfs`**: This method estimates the site frequency spectrum and infers ancestral alleles. It requires the presence of outgroup samples (ancestral to the rest of the data) in the dataset.
3. **Using `compara`**: This method requires an additional CSV file containing the ancestral alleles.

## Getting the Pipeline

You can obtain this pipeline by cloning the GitHub repository:

```bash
git clone cnr-ibba/nf-treeseq
```

Alternatively, you can use the `nextflow pull` command:

```bash
nextflow pull cnr-ibba/nf-treeseq
```

For more information on installing and running Nextflow pipelines, including dealing with revisions, refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html).

## Usage

While all parameters can be passed via the command line, it is recommended to use a configuration file. The configuration file should be a simple JSON file containing at least the following parameters:

```json
{
  "plink_bfile": "<binary plink prefix>",
  "plink_species": "<plink species options>",
  "plink_keep": "<plink keep file>",
  "plink_geno": 0.1,
  "genome": "<genome file>"
}
```

### Explanation of Parameters:

- **`plink_bfile`**: The binary PLINK file prefix used as the `--bfile` parameter.
- **`plink_species`**: Species-specific options for PLINK (e.g., `--species sheep` or `--chr-set 26 no-xy no-mt --allow-no-sex`).
- **`plink_keep`**: A TSV file with `FID` and `IID` columns indicating the samples to keep.
- **`plink_geno`**: The PLINK `--geno` parameter (default: 0.1), which excludes SNPs with a higher missing rate.
- **`genome`**: The genome file used by `bcftools` for allele normalization (setting ALT/REF alleles) and chromosome size correction.

### Specifying Ancestral Alleles

The pipeline requires ancestral alleles to generate tree sequences. At least one of the following methods must be used to infer ancestral alleles:

#### 1. Using the Reference Genome

To use the reference genome for inferring ancestral alleles, simply set the `reference_ancestor` flag:

```json
{
  "reference_ancestor": true
}
```

#### 2. Using `est-sfs` to Infer Ancestral Alleles

To infer ancestral alleles using `est-sfs`, enable the `with_estsfs` flag and specify one or more outgroup sample files (TSV format with `FID` and `IID` columns). You can provide up to three outgroup files:

```json
{
  "with_estsfs": true,
  "outgroup1": "<outgroup1 samples file>",
  "outgroup2": "<outgroup2 samples file>",
  "outgroup3": "<outgroup3 samples file>"
}
```

#### 3. Using `compara` to Infer Ancestral Alleles

To use `compara` for inferring ancestral alleles, provide a CSV file with the following format:

```csv
chrom,position,alleles,anc_allele
26,209049,A/G,C
26,268822,A/G,C
26,285471,A/G,G
26,361728,G/T,G
```

After generating the file, specify it using the `compara_ancestor` parameter:

```json
{
  "compara_ancestor": "<compara file>"
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

Once your configuration file is set up, run the pipeline with:

```bash
nextflow run cnr-ibba/nf-treeseq -profile <profile> -params-file <config.json>
```

- `<profile>`: The execution environment profile (e.g., `docker` or `singularity`).
- `<config.json>`: The configuration file you created.

You can also override specific parameters directly in the command line:

```bash
nextflow run cnr-ibba/nf-treeseq -profile singularity -params-file config.json --plink_geno 0.2
```

# cnr-ibba/nf-treeseq

[![GitHub Actions CI Status](https://github.com/cnr-ibba/nf-treeseq/actions/workflows/nf-test.yml/badge.svg)](https://github.com/cnr-ibba/nf-treeseq/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/cnr-ibba/nf-treeseq/actions/workflows/linting.yml/badge.svg)](https://github.com/cnr-ibba/nf-treeseq/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/cnr-ibba/nf-treeseq)

## Introduction

**cnr-ibba/nf-treeseq** is a bioinformatics pipeline that ...

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/guidelines/graphic_design/workflow_diagrams#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run cnr-ibba/nf-treeseq \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

cnr-ibba/nf-treeseq was originally written by Paolo Cozzi.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

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

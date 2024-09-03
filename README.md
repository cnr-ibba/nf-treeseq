
# nf-treeseq

A Nextflow pipeline for generating Tree Sequences from PLINK and VCF files.

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![Run with Docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![Run with Singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Background

This pipeline is designed to infer *Tree Sequences* from genotype data. It is currently tailored for PLINK genotype files, where all relevant samples are contained within a single file. The pipeline converts the PLINK file into a VCF file, corrects ALT/REF alleles, and checks chromosome sizes. It then uses Beagle to impute and phase any missing data before running `tsinfer` to create Tree Sequences from the VCF file.

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

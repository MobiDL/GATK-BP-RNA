# Workflow : GATK-best-practices-RNA

This workflow do VariantCalling on Illumina RNA sequencing data.

## Installation

### Clone the repo

```bash
git clone --recursive https://github.com/MobiDL/GATK-BP-RNA.git
```
### Install dependencies

__Conda is recommended:__ [see here to install conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

You will need to create a new environment based on conda.

```bash
conda env create -f environment.yml
conda activate GATK-best-practices-RNA
```

## Configure your inputs

You can create your input file replacing editing the template or creating your own inputs file.

### Minimal

```json
{
	"GATK_best_practices_RNA.fastqR1": "File",
	"GATK_best_practices_RNA.fastqR2": "File? (optional)",
	"GATK_best_practices_RNA.dbsnp": "File",
	"GATK_best_practices_RNA.dbsnpIdx": "File",
	"GATK_best_practices_RNA.knownSites": "Array[File]",
	"GATK_best_practices_RNA.knownSitesIdx": "Array[File]",
	"GATK_best_practices_RNA.genomeDir": "File",
	"GATK_best_practices_RNA.outputPath": "String",
	"GATK_best_practices_RNA.refDict": "File",
	"GATK_best_practices_RNA.refFai": "File",
	"GATK_best_practices_RNA.refFasta": "File",
}
```

### Extended

A full option templates is provided (inputs.json.tpl).

This template is separating in 4 categories (blank line) :
1. Global pipeline inputs (i.e. minimal)
2. Global pipeline options
3. Specific tasks inputs
4. Specific tasks options

## Launch

### Local

```bash
conda activate GATK-best-practices-RNA
cromwell run \
	-Dconfig.file=backends.conf/local.conf \
	--inputs /path/to/inputs.json \
	GATK-BP-RNA.wdl
```

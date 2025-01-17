<img width="685" alt="image" src="https://github.com/gks-anvil/vrs_anvil_toolkit/assets/47808/909db052-972c-4508-a2f4-8a389de03320">


# VRS AnVIL Toolkit

## Project Overview

This Python package is designed to process Variant Call Format (VCF) files and perform lookup operations on Genomic Variation Representation Service (GA4GH VRS) identifiers. The GA4GH VRS identifiers provide a standardized way to represent genomic variations, making it easier to exchange and share genomic information.

In addition, this project facilitates the retrieval of evidence associated with genomic alleles by leveraging the Genomic Data Representation and Knowledge Base (GA4GH MetaKB) service. GA4GH MetaKB provides a comprehensive knowledge base that links genomic variants to relevant clinical variant interpretations.

## Features

1. **VCF File Processing:**
   - Streamlines reading and parsing of VCF files, to extract relevant genomic information.

2. **GA4GH VRS Identifier Lookup:**
   - Utilizes the GA4GH VRS API to perform lookups for each genomic variation mentioned in the VCF file.
   - Retrieves standardized identifiers for the alleles, enhancing interoperability with GA4GH-compliant systems.
   - GA4GH MetaKB Service Integration:  Utilizes the GA4GH MetaKB retrieve evidence associated with specified genomic alleles.

3. **Output Generation:**
   - Generates summary metrics about throughput, errors, evidence, and hits.
   - Presents the retrieved evidence in a structured format, providing access to information about studies, publications, and other relevant details.
   - Create cohort allele frequency objects

4. **Additional Features**
   - Provides configurable options like threading and caching for processing VCFs.
   - Implements robust error handling to address issues like invalid input files, invalid variants, and more.

## Getting Started

### Prerequisites

- Python 3.10 or later
- Internet connectivity for data dependency setup (seqrepo)

### Installation

1. Get the repository either by...
   1. Source code
   ```bash
   git clone https://github.com/gks-anvil/vrs_anvil_toolkit
   cd vrs_anvil_toolkit
   ```
   2. PyPi
   ```bash
   pip install vrs_anvil_toolkit
   ```

2. Install dependencies either...
   1. for local use
   ```bash
   # install postgresql@14 (required for vrs-python)
   brew install postgresql@14
   bash scripts/setup.sh
   ```
   2. for use on Terra
   ```bash
   bash terra/setup.sh
   ```

### Usage
**General**
All usage has the following general steps...

1. Create a manifest to configure your VCF processing run
1. Use the `vrs_bulk` CLI to create a metrics file of related evidence
1. Use the metrics files for downstream analysis

The follow steps are explained in detail below, with some additional info on using vrs-python to directly annotate VCFs with VRS IDs.

**Manifest**

The configuration of each VCF processing run run is controlled by a `manifest.yaml` file. Most importantly, this file specifies the...
- input VCF file(s) to process
- working directories
- performance and strictness configurations

Use this commented [sample manifest](tests/fixtures/manifest.yaml) as a starting point on the specific variables you can specify per run.

**CLI**

Below are a list of command line utilities that may be useful
```bash
# activate the environment
source venv/bin/activate

# run the vrs_bulk command in the foreground
vrs_bulk annotate

# run the vrs_bulk command in parallel, one process per VCF file
vrs_bulk annotate --scatter

# run the vrs_bulk command in parallel in the background
nohup vrs_bulk annotate --scatter & # press enter to continue

# get the status of the processes for the most recent scatter run
vrs_bulk ps
```

The command line utility supports Google Cloud URIs and running commands in the background to interop with Terra out-of-the-box. This is described in the CLI usage above. For an example notebook, see `vrs-anvil-demo.ipynb` on the `vrs-anvil` workspace.

## Cohort Allele Frequency Generation

### Description
Given a variant of interest and an optional phenotype of interest, get its aggregated allele frequency information as a cohort allele frequency object (CAF).

### General Prerequisites
- Variant of interest
- Valid VRS-annotated joint VCF
  - Assumes chr field is prepended with chr (eg `chr1`)
  - genotyping laid out per-sample
- Precomputed VRS-VCF index (created using [vrsix](https://github.com/gks-anvil/vrsix))
- Phenotype of interest, if desired
- Plugins for project-specific transformations (see [here](README.md#plugins-for-unique-data-inputs) for more info)

### Use Cases
1. Given a variant ID and VCF path, get the CAF for the entire cohort
   - Get VCF row corresponding to variant ID using a variant -> VCF row index
   - Get phenotypes corresponding to each participants using the phenotypes by patient table
   - Aggregate counts for participants using their genotypes
   - Create CAF object using counts

2. Given a variant ID, VCF path, **and participant list**, get the CAF for a subset of participants (subcohort)
   - Same as 1 with subcohort defined by the participant list

3. Given a variant ID, VCF path, **and phenotype**, get the CAF for a subcohort with a specified phenotype
   - Same as 1, but for a subcohort of samples with the phenotype


## Plugins for Unique Data Inputs

### Description
Given the broad variety of data representation used by data generators, we want to be able to generate CAFs for any data consortium. One way this is possible is through the use of a plugin architecture.

A plugin architecture allows users to customize the aggregation of allelic data specified to their own data model. This addresses the problem where...
1. a user has a project-specific phenotype data model
   1. Example: sample-level rare disease data is stored in an unaggregated Terra data table
2. a user is interested in a subset of the cohort based on particular filters
   1. Example: samples must have phenotype A and a minimum read depth to be included in the subcohort
3. each sample's genotype must be aggregated uniquely depending on particular traits
   1. Example: sex is not represented within the VCF, so a user needs to integrate sample-level phenotype data to get accurate counts for chrX variants

These map to three different methods in a given `Plugin` class:
1. `__init__`: Given any set of parameters, create a phenotype index that maps each sample to its list of phenotypes.
2. `include_sample`: given a sample's variant-level or phenotype data, determine whether to include the sample in the allele count
   1. This takes a `pysam.VariantRecord` as input to represent a particular variant record in a VCF
   2. For more details, consult the pysam [VariantRecord docs](https://pysam.readthedocs.io/en/latest/api.html#pysam.VariantRecord).
3. `process_sample_genotype`: determine how to sum the alleles of a sample's genotype using variant-level or phenotype data
   1. This also makes use of a `pysam.VariantRecord` as input
   2. An `alt_index` is also passed in as an input, which describes the index for the allele of interest.
   3. The alt index matches the genotype according to [VCF specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf). For instance, a sample with the 2nd alt might have a genotype containing a 2, ie `(2,1)`, `(2,0)`, `(2,2)`, etc.

For the methods signatures and default implementations, take a look at the [`BasePlugin`](plugin_system/plugins/base_plugin.py) class. For an example plugin, take a look at the [`GregorPlugin`](plugin_system/plugins/gregor_plugin.py). To implement your own plugin....

### Getting Started
1. Copy [`gregor_plugin.py`](plugin_system/plugins/gregor_plugin.py) to the same directory
2. Rename the plugin class and name (eg `MyProjectPlugin` and `my_project_plugin.py`)
3. Implement the three methods mentioned above, calling any default implementations as necessary
   1. [`BasePlugin`](plugin_system/base_plugin.py) is by default the parent class, so you can use the `BasePlugin`'s implementations by calling `super().<method_to_invoke>`
   2. [`GregorPlugin`](plugin_system/plugins/gregor_plugin.py) is a worked example of specific real-world implementation, refer to that for alternative ways to customize allele frequency generation.
4. See this [starter code](scripts/plugin_usage_example.py) for a worked example on how to use the plugin. The two main components of using the plugin are...
   1. For your `MyProjectPlugin` plugin, instantiate it with the `PluginManager` and any input parameters specified.
   2. Call `get_cohort_allele_frequency` with `plugin="MyProjectPlugin"` as a parameter.
   3.

# GREGoR-specific Details

### Work in Progress
- For chromosomes with ploidy of 1 (mitochondrial calling or sex chromosomes), focus allele counts (AC) and locus allele counts (AN) can have a maximum value of 1. Focus allele counts are 1 when the genotype has at least a single allele match (0/1, 1/1, or 1) otherwise it is none.


## Processing VCF Files ([vrs-python](https://github.com/ga4gh/vrs-python))###

vrs-python is a GA4GH GKS package centered around creating Variant Representation specification (VRS) IDs: consistent, globally unique identifiers for variation. Some of its functionality includes variant ID translation and VCF annotation. Used as a dependency in `vrs_bulk`, it can also be used as a standalone package.

For Python usage, see [vrs_vcf_annotator.py](scripts/vrs_vcf_annotator.py) for an example.

For CLI usage:
```bash
python3 -m ga4gh.vrs.extras.vcf_annotation --vcf_in tests/fixtures/1kGP.chr1.1000.vcf --vcf_out annotated_output.vcf.gz --vrs_pickle_out allele_dicts.pkl --seqrepo_root_dir ~/seqrepo/latest
```

The above is an example using an example vcf. Replace the `--vcf_out` and `vrs_pickle_out` here with your desired output file path, where the output vcf can be BCF (`vcf.gz`) or VCF (`vcf`)

Also, see the [VRS Annotator](https://dockstore.org/workflows/github.com/gks-anvil/vrs-annotator/VRSAnnotator:main?tab=info) workflow on Dockstore for a way to do this on Terra.

## Contributing

This project is open to contributions from the research community. If you are interested in contributing to the project, please contact the project team.
See the [contributing guide](CONTRIBUTING.md) for more information on how to contribute to the project.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.md) file for details.

<img width="685" alt="image" src="https://github.com/ohsu-comp-bio/vrs-python-testing/assets/47808/909db052-972c-4508-a2f4-8a389de03320">


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
   git clone https://github.com/ohsu-comp-bio/vrs_anvil_toolkit
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
Given a variant of interest, Create a cohort allele frequency object, subsettable by participant list and a phenotype code.

### General Prerequisites
- Variant ID of interest
- VCF path to file
  - chr field is prepended with chr (eg `chr1`)
  - genotyping laid out per-patient (eg a row has column `PATIENT_1` with value `0/1`)
- Precomputed VRS-VCF index (created using [vrsix](https://github.com/gks-anvil/vrsix))
- Access to phenotypes table either through Terra (default) or as a local file (structured according to the [GREGOR data model](https://gregorconsortium.org/data-model))

### Use Cases
1. Given a variant ID and VCF path, get the allele frequency for the entire cohort
   - Get VCF row corresponding to variant ID using a variant -> VCF row index
   - Get phenotypes corresponding to each participants using the phenotypes by patient table
   - Aggregate counts for participants using their genotypes
   - Create CAF object using counts

2. Given a variant ID, VCF path, **and participant list**, get the allele frequency for a subset of participants (subcohort)
   - Same as 1, just subsetted on a participant list

3. Given a variant ID, VCF path, **and phenotype**, get the allele frequency for the cohort conditional on the phenotype
   - Same as 1, but only increase the counts for the variant of interest if a given patient has the specified phenotype

### Arguments
 - `variant_id` (String): variant ID of interest (VRS ID)
 - `vcf_path` (String): path to VCF file
 - `phenotype_table` (String, optional): where to pull phenotype information from. Defaults to None.
 - `participant_list` (List of Strings, optional): Subset of participants to use. Defaults to None.
 - `phenotype` (String, optional): Specific phenotype to subset on. Defaults to None.

### Example Usage on Terra
```python
# imports
from vrs_anvil.evidence import create_patient_phenotype_index, get_cohort_allele_frequency

# get variant of interest
allele_translator = AlleleTranslator(seqrepo_dataproxy)
variant_ids = ["chr3-10172-AC-A"]
allele = allele_translator.translate_from(variant_id)
vrs_id = allele.id

# specify data paths
vcf_path = "/path/to/vcf"
vcf_index_path = "path/to/vcf/index"

# creating an index
pheno_index_path = "/home/jupyter/pheno.json"
create_patient_phenotype_index(as_set=True, save_path=pheno_index_path)

# generating cohort allele frequency
from vrs_anvil.evidence import create_patient_phenotype_index, get_cohort_allele_frequency
get_cohort_allele_frequency(
   variant_id = vrs_id,
   vcf_path = vcf_path,
   vcf_index_path = vcf_index_path,
   phenotype_index_path=pheno_index_path
)
```
### Work in Progress
- For chromosomes with ploidy of 1 (mitochondrial calling or sex chromosomes), focus allele counts (AC) and locus allele counts (AN) can have a maximum value of 1. Focus allele counts are 1 when the genotype has at least a single allele match (0/1, 1/1, or 1) otherwise it is none.


**Processing VCF Files ([vrs-python](https://github.com/ga4gh/vrs-python))**

vrs-python is a GA4GH GKS package centered around creating Variant Representation specification (VRS) IDs: consistent, globally unique identifiers for variation. Some of its functionality includes variant ID translation and VCF annotation. Used as a dependency in `vrs_bulk`, it can also be used as a standalone package.

For Python usage, see [vrs_vcf_annotator.py](scripts/vrs_vcf_annotator.py) for an example.

For CLI usage:
```bash
python3 -m ga4gh.vrs.extras.vcf_annotation --vcf_in tests/fixtures/1kGP.chr1.1000.vcf --vcf_out annotated_output.vcf.gz --vrs_pickle_out allele_dicts.pkl --seqrepo_root_dir ~/seqrepo/latest
```

The above is an example using an example vcf. Replace the `--vcf_out` and `vrs_pickle_out` here with your desired output file path, where the output vcf can be BCF (`vcf.gz`) or VCF (`vcf`)

Also, see the [VRS Annotator](https://dockstore.org/workflows/github.com/ohsu-comp-bio/vrs-annotator/VRSAnnotator:main?tab=info) workflow on Dockstore for a way to do this on Terra.

### Contributing

This project is open to contributions from the research community. If you are interested in contributing to the project, please contact the project team.
See the [contributing guide](CONTRIBUTING.md) for more information on how to contribute to the project.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.md) file for details.

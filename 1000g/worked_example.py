import json
import os
import subprocess

from ga4gh.vrs.extras.translator import AlleleTranslator
from ga4gh.vrs.dataproxy import create_dataproxy
from plugin_system.plugin_manager import PluginManager
from vrs_anvil.evidence import get_cohort_allele_frequency


# run this in 1000g directory
assert os.getcwd().endswith(
    "1000g"
), "to ensure the plugin can be located, please run this in the 1000g directory"

# set varaible for variant data input
variant_id = "chr1-20094-TAA-T"
vcf_path = "../tests/fixtures/1kGP.chr1.1000.vrs.vcf.gz"

# set path to write VCF index to
vcf_index_path = "1000g_chr1_index.db"

# set phenotype-specific inputs
phenotype = "USA"  # to create subcohorts
phenotype_table = "population_descriptor.tsv"  # downloaded from https://anvil.terra.bio/#workspaces/anvil-datastorage/AnVIL_1000G_PRIMED-data-model/data

# create vcf index from vcf at the specified path using vrsix
command = ["vrsix", "load", f"--db-location={vcf_index_path}", vcf_path]
try:
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    print("vrsix command executed successfully!")
except subprocess.CalledProcessError as e:
    print("Error executing vrsix command:", e.stderr)

# # get VRS ID from variant of interest
seqrepo_rest_service_url = "seqrepo+https://services.genomicmedlab.org/seqrepo"
seqrepo_dataproxy = create_dataproxy(uri=seqrepo_rest_service_url)
allele_translator = AlleleTranslator(seqrepo_dataproxy)
allele = allele_translator.translate_from(variant_id)
vrs_id = allele.id

# instantiate 1000G plugin class with phenotype table
plugin = PluginManager().load_plugin("ThousandGenomesPlugin")
simple_plugin = plugin(phenotype_table)

# generating cohort allele frequency using 1000G plugin
caf = get_cohort_allele_frequency(
    variant_id=vrs_id,
    vcf_path=vcf_path,
    vcf_index_path=vcf_index_path,
    plugin=simple_plugin,
    phenotype=phenotype,
)

print(f"CAF:")
print(json.dumps(caf.model_dump(exclude_none=True), indent=2))

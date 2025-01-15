from ga4gh.vrs.extras.translator import AlleleTranslator
from ga4gh.vrs.dataproxy import create_dataproxy
from plugin_system.plugin_manager import PluginManager
from vrs_anvil.evidence import PLUGIN_DIR, get_cohort_allele_frequency

# get VRS ID from variant of interest
variant_id = "chr3-10172-AC-A"

seqrepo_rest_service_url = "seqrepo+https://services.genomicmedlab.org/seqrepo"
seqrepo_dataproxy = create_dataproxy(uri=seqrepo_rest_service_url)
allele_translator = AlleleTranslator(seqrepo_dataproxy)
allele = allele_translator.translate_from(variant_id)
vrs_id = allele.id

# specify data paths
vcf_path = "/path/to/vcf"
vcf_index_path = "path/to/vcf/index"
phenotype_table_path = "path/to/pheno/data"

# get GREGoR plugin class
plugin_manager = PluginManager(PLUGIN_DIR)
plugin = plugin_manager.load_plugin("GregorPlugin")

# instantiate plugin with gregor data
gregor_plugin = plugin(phenotype_table_path=phenotype_table_path) # here you can specify no args on Terra to pull from the "phenotype" Data Table in your workspace

# generating cohort allele frequency using GREGoR plugin
caf = get_cohort_allele_frequency(
   variant_id = vrs_id,
   vcf_path = vcf_path,
   vcf_index_path = vcf_index_path,
   plugin=gregor_plugin
)

print("caf:", caf)

import os
import pytest
import subprocess

from ga4gh.vrs.extras.translator import AlleleTranslator
from ga4gh.vrs.dataproxy import create_dataproxy
from plugin_system.plugin_manager import PluginManager
from vrs_anvil.evidence import PLUGIN_DIR, get_cohort_allele_frequency


@pytest.fixture()
def existing_vcf_index_path():
    return "tests/fixtures/chr1_multi_sample_index.db"


@pytest.fixture()
def focus_allele_count():
    return 3


@pytest.fixture()
def locus_allele_count():
    return 4


@pytest.fixture()
def phenotype():
    return "HP:0001263"


@pytest.fixture()
def phenotype_index():
    return {
        "sample_A": ["HP:0001263", "HP:0000002"],
        "sample_B": ["HP:0000001"],
    }


@pytest.fixture()
def plugin_methods():
    return [
        "get_phenotype_index",
        "include_sample",
        "process_sample_genotype",
    ]


@pytest.fixture()
def plugin_names():
    return ["BasePlugin", "GregorPlugin", "SimplePlugin"]


@pytest.fixture()
def variant_id():
    return "chr1-20094-TAA-T"


@pytest.fixture()
def vcf_index_path():
    # this path should be empty before
    path = "tests/fixtures/chr1_index.db"

    if os.path.exists(path):
        raise Exception(
            f"VCF index path ({path}) should not exist, please delete ebfore continuing"
        )

    yield path

    if os.path.exists(path):
        os.remove(path)


@pytest.fixture()
def vcf_path():
    return "tests/fixtures/chr1_multi_sample_vrs.vcf.gz"


@pytest.fixture()
def vrs_id():
    return "ga4gh:VA.1WRXw4TC5DYsjC2QLyKux9C0xETLN3Yt"


def test_plugin_manager_can_load_plugins(plugin_methods, plugin_names):
    for plugin_name in plugin_names:
        plugin_manager = PluginManager(PLUGIN_DIR)
        plugin = plugin_manager.load_plugin(plugin_name)

        print("plugin:", plugin)

        for method in plugin_methods:
            assert hasattr(
                plugin, method
            ), f"{plugin_name} does not have method '{method}'"


def test_simple_plugin_can_generate_caf_no_phenotype_index(
    vrs_id, vcf_path, existing_vcf_index_path, focus_allele_count, locus_allele_count
):
    # load plugin class
    plugin_manager = PluginManager(PLUGIN_DIR)
    plugin = plugin_manager.load_plugin(
        "SimplePlugin"
    )  # SimplePlugin inherits all of base plugin's methods

    # instantiate plugin with no phenotype index
    simple_plugin_no_pheno = plugin()

    # for generated vcf path and vcf_index_path
    caf = get_cohort_allele_frequency(
        variant_id=vrs_id,
        vcf_path=vcf_path,
        vcf_index_path=existing_vcf_index_path,
        plugin=simple_plugin_no_pheno,
    )

    print(f"CAF no phenotype: {caf}")

    check_allele_counts(caf, vrs_id, focus_allele_count, locus_allele_count)


def test_simple_plugin_can_generate_cafs_with_phenotype_index(
    vrs_id,
    vcf_path,
    existing_vcf_index_path,
    focus_allele_count,
    locus_allele_count,
    phenotype,
    phenotype_index,
):
    # load plugin class
    plugin_manager = PluginManager(PLUGIN_DIR)
    plugin = plugin_manager.load_plugin(
        "SimplePlugin"
    )  # SimplePlugin inherits all of base plugin's methods

    #####
    # 1 #
    #####
    # create caf using plugin with phenotype index but no phenotype specified
    simple_plugin = plugin(phenotype_index=phenotype_index)
    caf = get_cohort_allele_frequency(
        variant_id=vrs_id,
        vcf_path=vcf_path,
        vcf_index_path=existing_vcf_index_path,
        plugin=simple_plugin,
    )

    print(f"CAF #2: {caf}")

    # check if phenotypes are loaded
    all_phenotypes = [phenotype, "HP:0000002", "HP:0000001"]
    for pheno in all_phenotypes:
        assert (
            pheno in caf["ancillaryResults"]["phenotypes"]
        ), f'{pheno} not found in caf["ancillary_results"]["phenotypes"]: {caf["ancillary_results"]["phenotypes"]}'

    # check if allele counts are accurate
    check_allele_counts(caf, vrs_id, focus_allele_count, locus_allele_count)

    #####
    # 2 #
    #####
    # now, create caf with specified phenotype for subcohort
    caf = get_cohort_allele_frequency(
        variant_id=vrs_id,
        vcf_path=vcf_path,
        vcf_index_path=existing_vcf_index_path,
        plugin=simple_plugin,
        phenotype=phenotype,
    )

    print(f"CAF with phenotype: {caf}")

    # check if phenotypes are loaded
    assert [phenotype] == caf["ancillaryResults"][
        "phenotypes"
    ], f'{phenotype} should be the only phenotype in caf["ancillary_results"]["phenotypes"], instead got: {caf["ancillary_results"]["phenotypes"]}'

    # check if allele counts are accurate
    check_allele_counts(
        caf=caf, focus_allele=vrs_id, focus_allele_count=1, locus_allele_count=2
    )


def test_plugin_worked_example(
    variant_id,
    vcf_path: str,
    vcf_index_path: str,
    phenotype_index: dict[str, list[str]],
    phenotype: str,
):
    """given a vcf, variant of interest, and phenotype of interest,
    create a

    Args:
        vcf_path (str): path to the VCF of interest
        vcf_index_path (str): path to unwritten index
        phenotype_index (dict[str, list[str]]): _description_
        phenotype (str): _description_
    """

    # get VRS ID from variant of interest
    seqrepo_rest_service_url = "seqrepo+https://services.genomicmedlab.org/seqrepo"
    seqrepo_dataproxy = create_dataproxy(uri=seqrepo_rest_service_url)
    allele_translator = AlleleTranslator(seqrepo_dataproxy)
    allele = allele_translator.translate_from(variant_id)
    vrs_id = allele.id

    # create vcf index from vcf at the specified path using vrsix
    command = ["vrsix", "load", f"--db-location={vcf_index_path}", vcf_path]
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        print("Command executed successfully!")
        print("Output:", result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error executing command.")
        print("Error message:", e.stderr)

    # get SimplePlugin class
    plugin_manager = PluginManager(PLUGIN_DIR)
    plugin = plugin_manager.load_plugin("SimplePlugin")

    # instantiate simple plugin with phenotype index
    simple_plugin = plugin(phenotype_index)

    # generating cohort allele frequency using GREGoR plugin
    caf = get_cohort_allele_frequency(
        variant_id=vrs_id,
        vcf_path=vcf_path,
        vcf_index_path=vcf_index_path,
        plugin=simple_plugin,
        phenotype=phenotype,
    )

    print("caf:", caf)

    check_allele_counts(
        caf=caf, focus_allele=vrs_id, focus_allele_count=1, locus_allele_count=2
    )


# check if allele counts are accurate
def check_allele_counts(
    caf: dict, focus_allele: str, focus_allele_count: int, locus_allele_count: int
) -> None:
    assert (
        caf["focusAllele"] == focus_allele
    ), f"Incorrect CAF: expected focusAllele {focus_allele} but got {caf["focusAllele"]}"

    assert (
        caf["focusAlleleCount"] == focus_allele_count
    ), f"Incorrect CAF: expected focusAlleleCount {focus_allele_count} but got {caf["focusAlleleCount"]}"
    assert (
        caf["locusAlleleCount"] == locus_allele_count
    ), f"Incorrect CAF: expected focusAlleleCount {locus_allele_count} but got {caf["focusAlleleCount"]}"

    allele_frequency = focus_allele_count * 1.0 / locus_allele_count
    assert (
        caf["alleleFrequency"] == allele_frequency
    ), f"Incorrect CAF: expected alleleFrequency {allele_frequency} but got {caf["alleleFrequency"]}"

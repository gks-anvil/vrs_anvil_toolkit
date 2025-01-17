import pytest

from plugin_system.plugin_manager import PluginManager
from vrs_anvil.evidence import PLUGIN_DIR, get_cohort_allele_frequency


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
def vrs_id():
    return "ga4gh:VA.1WRXw4TC5DYsjC2QLyKux9C0xETLN3Yt"


@pytest.fixture()
def vcf_path():
    return "tests/fixtures/chr1_multi_sample_vrs.vcf.gz"


@pytest.fixture()
def vcf_index_path():
    return "tests/fixtures/chr1_multi_sample_index.db"


@pytest.fixture()
def focus_allele_count():
    return 3


@pytest.fixture()
def locus_allele_count():
    return 4


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
    vrs_id, vcf_path, vcf_index_path, focus_allele_count, locus_allele_count
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
        vcf_index_path=vcf_index_path,
        plugin=simple_plugin_no_pheno,
    )

    print(f"CAF no phenotype: {caf}")

    check_allele_counts(caf, vrs_id, focus_allele_count, locus_allele_count)


def test_simple_plugin_can_generate_cafs_with_phenotype_index(
    vrs_id, vcf_path, vcf_index_path, focus_allele_count, locus_allele_count
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
    simple_plugin = plugin(
        phenotype_index={
            "sample_A": ["HP:0001263", "HP:0000002"],
            "sample_B": ["HP:0000001"],
        }
    )
    caf = get_cohort_allele_frequency(
        variant_id=vrs_id,
        vcf_path=vcf_path,
        vcf_index_path=vcf_index_path,
        plugin=simple_plugin,
    )

    print(f"CAF #2: {caf}")

    # check if phenotypes are loaded
    all_phenotypes = ["HP:0001263", "HP:0000002", "HP:0000001"]
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
    phenotype = all_phenotypes[0]

    caf = get_cohort_allele_frequency(
        variant_id=vrs_id,
        vcf_path=vcf_path,
        vcf_index_path=vcf_index_path,
        plugin=simple_plugin,
        phenotype=phenotype,
    )

    print(f"CAF with phenotype: {caf}")

    # check if phenotypes are loaded
    assert [phenotype] == caf["ancillaryResults"][
        "phenotypes"
    ], f'{phenotype} should be the only phenotype in caf["ancillary_results"]["phenotypes"], instead got: {caf["ancillary_results"]["phenotypes"]}'

    # check if allele counts are accurate
    check_allele_counts(caf, vrs_id, 1, 2)


# check if allele counts are accurate
def check_allele_counts(caf, fa, fac, lac):
    assert (
        caf["focusAllele"] == fa
    ), f"Incorrect CAF: expected focusAllele {fa} but got {caf["focusAllele"]}"
    assert (
        caf["focusAlleleCount"] == fac
    ), f"Incorrect CAF: expected focusAlleleCount {fac} but got {caf["focusAlleleCount"]}"
    assert (
        caf["locusAlleleCount"] == lac
    ), f"Incorrect CAF: expected focusAlleleCount {lac} but got {caf["focusAlleleCount"]}"

    allele_frequency = fac * 1.0 / lac
    assert (
        caf["alleleFrequency"] == allele_frequency
    ), f"Incorrect CAF: expected alleleFrequency {allele_frequency} but got {caf["alleleFrequency"]}"

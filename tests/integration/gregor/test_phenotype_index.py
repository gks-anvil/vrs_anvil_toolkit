import os
from pathlib import Path
from pysam import VariantFile

from plugin_system.plugins.gregor_plugin import GregorPlugin
from plugin_system.utils import load_dict, save_dict


def test_gregor_plugin_creates_correct_phenotype_index(
    chrY_vcf_path: str, gregor_plugin: GregorPlugin
):
    """given all sample IDs in a VCF, phenotypes index should be created with
    the expected number of unique phenotypes"""

    vcf = VariantFile(chrY_vcf_path)
    assert "VRS_Allele_IDs" in vcf.header.info, (
        "no VRS_Allele_IDs key in INFO found, "
        "please ensure that this is an VRS annotated VCF"
    )

    for record in vcf.fetch():
        phenotypes_set = set()
        samples = [sample for sample in record.samples]
        phenotype_index = gregor_plugin.get_phenotype_index()

        for sample_id in samples:
            if sample_id in phenotype_index:
                phenotypes_set.update(phenotype_index[sample_id])
        print("len(phenotypes_set):", len(phenotypes_set))

        expected_num_phenotypes = 980
        assert (
            len(phenotypes_set) == expected_num_phenotypes
        ), f"Expected {expected_num_phenotypes} phenotypes, got {len(phenotypes_set)}"

        break


def test_loading_gregor_phenotype_index_by_path(
    gregor_plugin: GregorPlugin, tmp_path: Path
):
    index = gregor_plugin.get_phenotype_index()
    save_path = tmp_path / "index.json"
    save_dict(index, save_path)

    loaded_index = load_dict(save_path)

    assert (
        index == loaded_index
    ), "saved index does not match loaded index... use -vv flag for a better diff"
    os.remove(save_path)

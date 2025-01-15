import os
from pathlib import Path
from pysam import VariantFile

from plugin_system.plugins.gregor_plugin import GregorPlugin
from plugin_system.utils import load_dict, save_dict


def test_gregor_plugin_creates_correct_phenotype_index(
    chrY_vcf_path: str, gregor_plugin: GregorPlugin
):
    """creates phenotypes index and calculate list of unique phenotypes for the first five rows,
    checking total unique phenotypes the last row"""

    vcf = VariantFile(chrY_vcf_path)
    assert "VRS_Allele_IDs" in vcf.header.info, (
        "no VRS_Allele_IDs key in INFO found, "
        "please ensure that this is an VRS sannotated VCF"
    )

    for i, record in enumerate(vcf.fetch()):
        vrs_ids = record.info["VRS_Allele_IDs"]
        phenotypes_set = set()
        patients = [
            patient
            for patient, genotype in record.samples.items()
            if 1 in genotype.allele_indices
        ]
        print(f"{record.chrom}-{record.pos}-{record.ref}-{record.alts}")

        print("VRS IDs:", vrs_ids)

        phenotype_index = gregor_plugin.get_phenotype_index()

        print("patients:", patients)
        for patient_id in patients:
            if patient_id in phenotype_index:
                phenotypes_set.update(phenotype_index[patient_id])
        print("len(phenotypes_set):", len(phenotypes_set))

        if i == 4:
            expected_num_phenotypes = 528
            assert (
                len(phenotypes_set) == expected_num_phenotypes
            ), f"Expected {expected_num_phenotypes} phenotypes, got {len(phenotypes_set)}"

            break


def test_loading_gregor_phenotype_index_by_path(
    gregor_plugin: GregorPlugin, tmp_path: Path
):
    os.chdir(tmp_path)

    index = gregor_plugin.get_phenotype_index()
    save_path = "index.json"
    save_dict(index, save_path)

    loaded_index = load_dict(save_path)

    assert (
        index == loaded_index
    ), "saved index does not match loaded index... use -vv flag for a better diff"
    os.remove(save_path)

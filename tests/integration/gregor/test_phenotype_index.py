import os
from pysam import VariantFile

from vrs_anvil.evidence import (
    create_patient_phenotype_index,
    get_patient_phenotype_index,
)


def test_get_phenotypes_from_vcf_row(remote_chry_vcf_path, phenotype_table):
    """creates phenotypes index and calculate list of unique phenotypes for the first five rows,
    checking total unique phenotypes the last row"""

    vcf = VariantFile(remote_chry_vcf_path)
    assert "VRS_Allele_IDs" in vcf.header.info, (
        "no VRS_Allele_IDs key in INFO found, "
        "please ensure that this is an VRS annotated VCF"
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

        phenotype_index = get_patient_phenotype_index(phenotype_table)

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


def test_save_and_use_phenotype_index_from_path(tmp_path, phenotype_table):
    """creates and saves a phenotype index from phenotype csv, then loads it back in"""

    # write to index
    os.chdir(tmp_path)
    save_path = "index.json"
    phenotype_index = create_patient_phenotype_index(
        phenotype_table, save_path=save_path
    )

    patient, saved_phenos = list(phenotype_index.items())[0]
    assert isinstance(
        saved_phenos, list
    ), "phenotypes in index are not returned as a list by default"

    # read in from index
    assert os.path.exists(
        save_path
    ), f"phenotype index file not being written to disk to save path {save_path}"
    loaded_pheno_index = get_patient_phenotype_index(cached_dict=save_path, as_set=True)

    # ensure phenotypes returned as a set
    assert patient in loaded_pheno_index
    loaded_phenos = loaded_pheno_index[patient]
    assert isinstance(
        loaded_phenos, set
    ), f"expected phenotype index values to be a set, but are {type(loaded_phenos)} instead despite specifying as_set=True"

    # check values are loaded as expected
    assert set(saved_phenos) == loaded_phenos

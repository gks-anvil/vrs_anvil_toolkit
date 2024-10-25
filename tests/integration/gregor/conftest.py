import os
import pytest

from pysam import VariantFile


@pytest.fixture
def remote_chry_vcf_path():
    """Return a GS URI to a remote VCF. Currently using annotated chrY GREGoR VCF"""

    return "gs://fc-secure-0b96edf7-d9b1-4ee2-94c9-c458fa44ccb1/gregor_joint_callset/gregor_consortium_u06_sorted_chrY_GRU_VRS.vcf.gz"


@pytest.fixture
def chr3_vcf_path():
    return "gs://fc-secure-0b96edf7-d9b1-4ee2-94c9-c458fa44ccb1/gregor_joint_callset/gregor_consortium_u06_sorted_chr3_GRU_VRS.vcf.gz"


@pytest.fixture()
def vrs_id_solo_alt(remote_chry_vcf_path):
    """VRS ID extracted from VCF row with only one alt"""

    for i, record in enumerate(VariantFile(remote_chry_vcf_path)):
        if i == 10:
            return record.info["VRS_Allele_IDs"][1]  # 1 since 0 is a ref


@pytest.fixture()
def phenotype_table():
    assert (
        "PHENOTYPE_TABLE" in os.environ
    ), "no PHENOTYPE_TABLE bash variable defining the path of the phenotype table, make sure to export"

    return os.environ["PHENOTYPE_TABLE"]

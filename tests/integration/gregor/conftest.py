import os
from pathlib import Path
import pytest

from pysam import VariantFile


@pytest.fixture
def chrY_vcf_path():
    """Return a GS URI to a remote VCF. Currently using annotated chrY GREGoR VCF"""

    return "gs://fc-secure-0b96edf7-d9b1-4ee2-94c9-c458fa44ccb1/gregor_joint_callset/gregor_consortium_u06_sorted_chrY_GRU_VRS.vcf.gz"


@pytest.fixture
def chr3_vcf_path():
    return "gs://fc-secure-0b96edf7-d9b1-4ee2-94c9-c458fa44ccb1/gregor_joint_callset/gregor_consortium_u06_sorted_chr3_GRU_VRS.vcf.gz"


@pytest.fixture()
def vrs_id_chr3(
    chr3_vcf_path,
):
    """VRS ID extracted from VCF row with only one alt"""

    chrom = "chr3"
    pos = 10456
    for record in VariantFile(chr3_vcf_path).fetch(chrom, pos - 1, pos):
        if record.ref == "CTT":
            print("record:", record)
            return record.info["VRS_Allele_IDs"][1]  # 1 since 0 is a ref

    raise ("couldn't find record for vrs_id_solo_alt")
    # for i, record in enumerate(VariantFile(remote_chry_vcf_path)):
    #     if i == 10:
    #         return record.info["VRS_Allele_IDs"][1]  # 1 since 0 is a ref


@pytest.fixture
def vrs_vcf_index() -> str:
    tests_dir = Path(os.path.dirname(__file__)).parent.parent
    index_path = os.path.join(
        tests_dir.absolute(), "fixtures/gregor/chr3_chrY_index.db"
    )

    print("vrs_vcf_index path:", index_path)
    if os.path.exists(index_path):
        return str(index_path)

    assert (
        "PHENOTYPE_TABLE" in os.environ
    ), "No VRS VCF index found, see tests/fixtures/gregor/README.md for setup"

    return os.environ["VRS_VCF_INDEX"]


@pytest.fixture()
def phenotype_table() -> str:
    tests_dir = Path(os.path.dirname(__file__)).parent.parent
    phenotype_table_path = os.path.join(
        tests_dir.absolute(), "fixtures/gregor/phenotypes.tsv"
    )

    print("phenotype_table_path:", phenotype_table_path)
    if os.path.exists(phenotype_table_path):
        return str(phenotype_table_path)

    assert (
        "PHENOTYPE_TABLE" in os.environ
    ), "No phenotype table found, see tests/fixtures/gregor/README.md for setup"

    return os.environ["PHENOTYPE_TABLE"]

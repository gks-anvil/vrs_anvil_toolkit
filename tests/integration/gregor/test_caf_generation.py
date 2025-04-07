import json
import os
import pytest

from typing import Generator
from pysam import VariantFile, VariantRecord
from pytest import approx

from vrs_anvil.evidence import get_cohort_allele_frequency

############
# FIXTURES #
############


@pytest.fixture
def chr3():
    """Return chromosome in the VCF."""
    return "chr3"


@pytest.fixture
def start():
    """Return the start range to query."""
    return 3188848


@pytest.fixture
def stop():
    """Return the end range to query."""
    return 3189029


@pytest.fixture
def expected_record_count():
    """Return the expected record count."""
    return 3


@pytest.fixture()
def vrs_id_chrY(chrY_vcf_path):
    """VRS ID extracted from VCF row with multiple alts"""

    for i, record in enumerate(VariantFile(chrY_vcf_path)):
        if i == 4:
            return record.info["VRS_Allele_IDs"][
                2
            ]  # index 2 refers to 2nd alt since 0 is a ref


def participants(record: VariantRecord) -> Generator[str, None, None]:
    """Return the participants that `have` this allele."""
    assert "GT" in record.format, "Genotype (GT) is required"

    for participant, values in record.samples.items():
        assert "GT" in values, "Genotype (GT) is required"
        # see https://samtools.github.io/hts-specs/VCFv4.1.pdf

        if any(values["GT"]):
            yield participant


#########
# TESTS #
#########


def test_remote_vcf(chrY_vcf_path, start, stop, expected_record_count):
    """Read a remote vcf file, query a range of alleles, check that at least 1 participant exists for each allele."""
    assert "GCS_OAUTH_TOKEN" in os.environ, (
        "GCS_OAUTH_TOKEN required: "
        "export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)"
        "see https://github.com/pysam-developers/pysam/issues/592#issuecomment-353693674 https://support.terra.bio/hc/en-us/articles/360042647612-May-4-2020"
    )
    try:
        vcf = VariantFile(chrY_vcf_path)  # auto-detect input format
        # fetch returns pysam.libcbcf.VariantRecord
        records = [_ for _ in vcf.fetch("chrY", start, stop)]
        assert len(records) == expected_record_count

        for variant_record in records:
            my_participants = [_ for _ in participants(variant_record)]
            assert len(my_participants) < len(
                variant_record.samples
            ), "Not all participants have this allele"
            assert len(my_participants) > 0, "No participants have this allele"
    except ValueError as e:
        print("ValueError: has GCS_OAUTH_TOKEN expired?", e)
        raise e


def test_allele_counts_first_5_rows(chr3_vcf_path, vrs_vcf_index, gregor_plugin):
    """test that the calculated allele counts with no phenotype specified matches
    the actual counts stored in the INFO of the first 10 rows. Works for diploid (non-sex) variants
    """

    vcf = VariantFile(chr3_vcf_path)
    has_ref = "REF" in vcf.header.info["VRS_Allele_IDs"].description

    for i, record in enumerate(vcf.fetch()):
        print("~~~~~~~~~ row ", i, "~~~~~~~~~")

        # use only alt VRS IDs, not REFs
        vrs_allele_ids = record.info["VRS_Allele_IDs"]
        print("vrs_allele_ids:", vrs_allele_ids)
        if has_ref:
            vrs_allele_ids = vrs_allele_ids[1:]

        # for each alt ID, ensure stored allele counts match calculated allele counts
        for alt_index, allele_id in enumerate(vrs_allele_ids):
            print("alt id:", allele_id)
            if not allele_id:
                continue

            caf = get_cohort_allele_frequency(
                allele_id,
                chr3_vcf_path,
                vcf_index_path=vrs_vcf_index,
                plugin=gregor_plugin,
            )

            print("alt_index", alt_index)
            print("AC:", record.info["AC"][alt_index], caf.focusAlleleCount)
            print("AN:", record.info["AN"], caf.locusAlleleCount)
            focus_allele_count = record.info["AC"][alt_index]
            locus_allele_count = record.info["AN"]

            check_caf_allele_data(caf, focus_allele_count, locus_allele_count)
        if i == 5:
            break


def test_correct_caf_given_chr3_variant(
    vrs_id_chr3, chr3_vcf_path, vrs_vcf_index, gregor_plugin
):
    """test caf generation with default parameters and no phenotype specified"""

    # get and log caf
    caf = get_cohort_allele_frequency(
        vrs_id_chr3,
        chr3_vcf_path,
        vcf_index_path=vrs_vcf_index,
        plugin=gregor_plugin,
    )
    print_caf(caf)

    # sanity checks
    assert (
        caf.type == "CohortAlleleFrequencyStudyResult"
    ), f"object of type CohortAlleleFrequency not returned, returned {caf.type} instead"
    assert (
        caf.focusAlleleCount <= caf.locusAlleleCount
    ), f"Focus allele count ({caf['focusAlleleCount']}) is larger than locus allele count ({caf['locusAlleleCount']})"

    print("focusAlleleCount:", caf.focusAlleleCount)
    print("locusAlleleCount:", caf.locusAlleleCount)

    # check allele counts and frequency
    check_caf_allele_data(caf, expected_fac=183, expected_lac=896)

    # ensure list of phenotypes stored
    assert (
        "phenotypes" in caf.ancillaryResults
    ), "no phenotype key stored in caf.ancillaryResults"
    assert (
        len(caf.ancillaryResults["phenotypes"]) > 0
    ), 'no phenotypes stored in caf.ancillaryResults["phenotypes"]'


def test_correct_caf_given_chr3_variant_and_pheno(
    vrs_id_chr3, chr3_vcf_path, vrs_vcf_index, gregor_plugin
):
    """test caf generation for diploid variant with a specified phenotype"""

    phenotype = "HP:0001263"

    caf = get_cohort_allele_frequency(
        vrs_id_chr3,
        chr3_vcf_path,
        vcf_index_path=vrs_vcf_index,
        plugin=gregor_plugin,
        phenotype=phenotype,
    )
    print_caf(caf)

    # sanity checks
    assert (
        caf.type == "CohortAlleleFrequencyStudyResult"
    ), f"object of type CohortAlleleFrequency not returned, returned {caf.type} instead"
    assert (
        caf.focusAlleleCount <= caf.locusAlleleCount
    ), f"Focus allele count ({caf['focusAlleleCount']}) is larger than locus allele count ({caf['locusAlleleCount']})"

    print("focusAlleleCount:", caf.focusAlleleCount)
    print("locusAlleleCount:", caf.locusAlleleCount)

    # check focus counts, locus counts, and allele frequency
    check_caf_allele_data(caf=caf, expected_fac=2, expected_lac=26)

    # check phenotype is stored in cohort data
    assert phenotype in caf.cohort.id
    assert phenotype in caf.cohort.name


def test_correct_allele_freq_for_multi_alts_chrY_variant(
    vrs_id_chrY, chrY_vcf_path, vrs_vcf_index, gregor_plugin
):
    """for a vcf row with multiple alts, test caf generation with default parameters and no phenotype specified"""

    # creat caf
    caf = get_cohort_allele_frequency(
        vrs_id_chrY,
        chrY_vcf_path,
        vcf_index_path=vrs_vcf_index,
        plugin=gregor_plugin,
    )

    # logs
    print(f"CAF generated for {caf.focusAllele.root}")
    print("focusAlleleCount:", caf.focusAlleleCount)
    print("locusAlleleCount:", caf.locusAlleleCount)

    # check allele frequency
    expected_allele_freq = 0.0491
    actual_allele_freq = approx(caf.focusAlleleFrequency, abs=1e-4)
    assert (
        actual_allele_freq == expected_allele_freq
    ), f"incorrect allele frequency, expected {expected_allele_freq} got {actual_allele_freq}"

    # ensure list of phenotypes stored
    assert (
        "phenotypes" in caf.ancillaryResults
    ), "no phenotype key stored in caf.ancillaryResults"
    assert (
        len(caf.ancillaryResults["phenotypes"]) > 0
    ), 'no phenotypes stored in caf.ancillaryResults["phenotypes"]'


def test_correct_allele_freq_for_multi_alts_chrY_variant_and_phenotype(
    vrs_id_chrY, chrY_vcf_path, vrs_vcf_index, gregor_plugin
):
    """test caf generation specifying both a variant and a phenotype of interest"""

    phenotype = "HP:0001263"
    caf = get_cohort_allele_frequency(
        vrs_id_chrY,
        chrY_vcf_path,
        vcf_index_path=vrs_vcf_index,
        plugin=gregor_plugin,
        phenotype=phenotype,
    )
    print_caf(caf)

    expected_allele_freq = 0.1034
    actual_allele_freq = approx(caf.focusAlleleFrequency, abs=1e-4)
    assert (
        actual_allele_freq == expected_allele_freq
    ), f"incorrect allele frequency, expected {expected_allele_freq} got {actual_allele_freq}"


###########
# HELPERS #
###########


def check_caf_allele_data(caf, expected_fac, expected_lac):
    # get allele frequencies
    actual_af = approx(caf.focusAlleleFrequency, abs=1e-4)
    expected_af = expected_fac * 1.0 / expected_lac

    # set up values to compare
    actual_expected_pairs_with_name = [
        (caf.focusAlleleCount, expected_fac, "focusAlleleCount"),
        (caf.locusAlleleCount, expected_lac, "locusAlleleCount"),
        (actual_af, expected_af, "focusAlleleFrequency"),
    ]

    # assert calculated values are expected
    for actual, expected, name in actual_expected_pairs_with_name:
        assert (
            actual == expected
        ), f"incorrect {name}, expected {expected_fac} got {caf.focusAlleleFrequency}"


def print_caf(caf):
    print("CAF:")
    print(json.dumps(caf.model_dump(exclude_none=True), indent=2))

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

        # TODO - this test is a bit naive, should we be more robust. consider AD, GQ, RGQ?
        if any(values["GT"]):
            yield participant


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
            print("AC:", record.info["AC"][alt_index], caf["focusAlleleCount"])
            print("AN:", record.info["AN"], caf["locusAlleleCount"])
            ac = record.info["AC"][alt_index]
            an = record.info["AN"]

            assert (
                caf["focusAlleleCount"] == ac
            ), f"row {i} alt {alt_index} has different focus allele counts, expected {ac}, got {caf['focusAlleleCount']}"
            assert (
                caf["locusAlleleCount"] == an
            ), f"row {i} alt {alt_index} has different locus allele counts, expected {an}, got {caf['locusAlleleCount']}"
        if i == 5:
            break


def test_correct_caf_given_chr3_variant(
    vrs_id_chr3, chr3_vcf_path, vrs_vcf_index, gregor_plugin
):
    """test caf generation with default parameters and no phenotype specified"""
    caf = get_cohort_allele_frequency(
        vrs_id_chr3,
        chr3_vcf_path,
        vcf_index_path=vrs_vcf_index,
        plugin=gregor_plugin,
    )
    print(json.dumps(caf))

    # sanity checks
    assert (
        caf["type"] == "CohortAlleleFrequency"
    ), f"object of type CohortAlleleFrequency not returned, returned {caf['type']} instead"
    assert (
        caf["focusAlleleCount"] <= caf["locusAlleleCount"]
    ), f"Focus allele count ({caf['focusAlleleCount']}) is larger than locus allele count ({caf['locusAlleleCount']})"

    print("focusAlleleCount:", caf["focusAlleleCount"])
    print("locusAlleleCount:", caf["locusAlleleCount"])

    # check allele frequency
    expected_allele_freq = 183.0 / 896
    actual_allele_freq = approx(caf["alleleFrequency"], abs=1e-4)
    assert (
        actual_allele_freq == expected_allele_freq
    ), f"incorrect allele frequency, expected {expected_allele_freq} got {actual_allele_freq}"

    # ensure important fields exist
    expected_fields = [
        "focusAlleleCount",
        "locusAlleleCount",
        "alleleFrequency",
        "ancillaryResults",
    ]
    for field in expected_fields:
        assert field in caf, f"expected field {field} in CAF"


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
    print(json.dumps(caf))

    # sanity checks
    assert (
        caf["type"] == "CohortAlleleFrequency"
    ), f"object of type CohortAlleleFrequency not returned, returned {caf['type']} instead"
    assert (
        caf["focusAlleleCount"] <= caf["locusAlleleCount"]
    ), f"Focus allele count ({caf['focusAlleleCount']}) is larger than locus allele count ({caf['locusAlleleCount']})"

    print("focusAlleleCount:", caf["focusAlleleCount"])
    print("locusAlleleCount:", caf["locusAlleleCount"])

    # check allele frequency
    expected_allele_freq = 2.0 / 26
    actual_allele_freq = approx(caf["alleleFrequency"], abs=1e-4)
    assert (
        actual_allele_freq == expected_allele_freq
    ), f"incorrect allele frequency, expected {expected_allele_freq} got {actual_allele_freq}"

    # ensure important fields exist
    expected_fields = [
        "focusAlleleCount",
        "locusAlleleCount",
        "alleleFrequency",
        "ancillaryResults",
    ]
    for field in expected_fields:
        assert field in caf, f"expected field {field} in CAF"


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
    print(f"CAF generated for {caf['focusAllele']}")
    print("focusAlleleCount:", caf["focusAlleleCount"])
    print("locusAlleleCount:", caf["locusAlleleCount"])

    # check allele frequency
    expected_allele_freq = 0.0491
    actual_allele_freq = approx(caf["alleleFrequency"], abs=1e-4)
    assert (
        actual_allele_freq == expected_allele_freq
    ), f"incorrect allele frequency, expected {expected_allele_freq} got {actual_allele_freq}"

    # ensure important fields exist
    expected_fields = [
        "focusAlleleCount",
        "locusAlleleCount",
        "alleleFrequency",
        "ancillaryResults",
    ]
    for field in expected_fields:
        assert field in caf, f"expected field {field} in CAF"


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
    print(json.dumps(caf))

    expected_allele_freq = 0.1034
    actual_allele_freq = approx(caf["alleleFrequency"], abs=1e-4)
    assert (
        actual_allele_freq == expected_allele_freq
    ), f"incorrect allele frequency, expected {expected_allele_freq} got {actual_allele_freq}"

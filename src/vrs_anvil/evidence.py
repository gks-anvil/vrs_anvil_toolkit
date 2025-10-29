import os
from pathlib import Path
import sqlite3

from datetime import datetime
from ga4gh.va_spec.base.core import CohortAlleleFrequencyStudyResult as CAF
from ga4gh.va_spec.base.core import DataSet, StudyGroup
from pysam import VariantFile, VariantRecord
from plugin_system.plugins.base_plugin import BasePlugin

# location to register plugins classes
PLUGIN_MODULE_PATH = "plugin_system.plugins"


def get_cohort_allele_frequency(
    variant_id: str,
    vcf_path: str,
    vcf_index_path: str | None = None,
    participant_list: list[str] | None = None,
    phenotype: str | None = None,
    plugin: BasePlugin | None = None,
) -> CAF:
    """Create a cohort allele frequency for either genotypes or phenotypes

    Args:
        variant_id (str): variant ID (VRS ID)
        vcf_path (str): path to VCF
        vcf_index_path (str): path to VRS to VCF coordinates index (SQLite table)
        phenotype_table (str, optional): where to pull phenotype information from. Defaults to None.
        participant_list (list[str], optional): Subset of participants to use. Defaults to None.
        phenotype (str, optional): Specific phenotype to subset on. Defaults to None.
        plugin (BasePlugin, optional): Plugin object to use for custom processing. Defaults to None, loading in the BasePlugin.

    Returns:
        dict: Cohort Allele Frequency object
    """

    # check variant_id is VRS ID
    assert (
        "ga4gh:VA" in variant_id
    ), "variant ID type not yet supported, use VRS ID instead"

    # use default plugin if none specified
    if plugin is None:
        plugin = BasePlugin()

    # get index of variant to patient
    # in this case, the VCF row of the variant_id
    vcf = VariantFile(vcf_path)
    record = get_vcf_row(variant_id, vcf, vcf_index_path)

    # if multiple alts, get index associated with alt allele
    # if ref has been saved in VRS Alleles IDs, adjust indices to match
    alt_index = record.info["VRS_Allele_IDs"].index(variant_id)
    if "REF" not in vcf.header.info["VRS_Allele_IDs"].description:
        alt_index -= 1

    # create cohort, defaults to all samples listed in VCF
    cohort = (
        set(participant_list) if participant_list is not None else set(record.samples)
    )

    # variables for cohort allele frequency (CAF) object
    focus_allele_count = 0
    locus_allele_count = 0
    cohort_phenotypes = set() if phenotype is None else [phenotype]

    # only relevat if the variant is diploid
    num_homozygotes = 0
    num_hemizygotes = 0

    # aggregate data for CAF so long as...
    for sample_id, genotype in record.samples.items():
        # 1. sample genotype call exist
        alleles = genotype.allele_indices
        if all(a is None for a in alleles):
            continue

        # 2. sample in specified cohort
        if sample_id not in cohort:
            continue

        # 3. matches subcohort criteria
        should_include_sample = plugin.include_sample(sample_id, record, phenotype)
        if phenotype is not None and not should_include_sample:
            continue

        # with these conditions satisfied...
        num_focus_alleles, num_locus_alleles = plugin.process_sample_genotype(
            sample_id, record, alt_index
        )

        # increment allele counts
        focus_allele_count += num_focus_alleles
        locus_allele_count += num_locus_alleles

        # record zygosity
        if num_focus_alleles == 1:
            num_hemizygotes += 1
        elif num_focus_alleles == 2:
            num_homozygotes += 1

        # update phenotypes as necessary
        if phenotype is not None:
            continue
        else:
            # aggregate phenotypes if they exist
            phenotype_index = plugin.get_phenotype_index()
            if phenotype_index is not None and sample_id in phenotype_index:
                cohort_phenotypes.update(phenotype_index[sample_id])

    # format caf fields before populating caf object
    allele_frequency = (
        focus_allele_count * 1.0 / locus_allele_count if locus_allele_count != 0 else 0
    )

    if phenotype is None:
        cohort = StudyGroup(id="ALL", name="Overall")
    else:
        cohort = StudyGroup(id=phenotype, name=phenotype)

    ancillary_results = {
        "homozygotes": num_homozygotes,
        "hemizygotes": num_hemizygotes,
        "phenotypes": list(cohort_phenotypes),
    }

    # populate final caf object according to va-spec-python
    caf = CAF(
        sourceDataSet=DataSet(id=vcf_path, description=f"Created {datetime.now()}"),
        focusAllele=variant_id,
        focusAlleleCount=focus_allele_count,
        focusAlleleFrequency=allele_frequency,
        locusAlleleCount=locus_allele_count,
        cohort=cohort,
        ancillaryResults=ancillary_results,
    )

    return caf


def fetch_by_vrs_ids(
    vrs_ids: list[str], db_location: Path | None = None
) -> list[tuple]:
    """Access index by VRS ID.

    :param vrs_id: VRS allele hash (i.e. everything after ``"ga4gh:VA."``)
    :param db_location: path to sqlite file (assumed to exist)
    :return: location description tuple if available
    """

    trunc_vrs_ids = []
    for vrs_id in vrs_ids:
        trunc_vrs_id = vrs_id[9:] if vrs_id.startswith("ga4gh:VA.") else vrs_id
        trunc_vrs_ids.append(trunc_vrs_id)

    if not db_location.exists():
        raise OSError(f"Index at {db_location} does not exist")

    conn = sqlite3.connect(db_location)

    # have to manually make placeholders for python sqlite API --
    # should still be safe against injection by using parameterized query
    placeholders = ",".join("?" for _ in trunc_vrs_ids)
    result = conn.cursor().execute(
        f"SELECT vrs_id, chr, pos FROM vrs_locations WHERE vrs_id IN ({placeholders})",  # noqa: S608
        trunc_vrs_ids,
    )
    data = result.fetchall()

    if len(data) == 0:
        raise Exception(
            f"No matching rows in the VCF index for the VRS IDs specified \n   - VRS IDs: {vrs_ids} \n   - Index path: {db_location}"
        )

    conn.close()
    return data


def get_vcf_row(
    variant_id: str, vcf: VariantFile, index_path: str = None
) -> VariantRecord:
    """given a variant id and annotated VCF, get the associated VCF row

    Args:
        variant_id (str): VRS ID for the variant of interest
        vcf (VariantFile): Pysam VariantFile
        index_path (str, optional): Index used to speed up search for variant. Defaults to iterating through VCF.

    Raises:
        Exception: outputs if no index is found

    Returns:
        VariantRecord: A Pysam VariantRecord (VCF row)
    """
    if "VRS_Allele_IDs" not in vcf.header.info:
        raise KeyError(
            "no VRS_Allele_IDs key in INFO found, "
            "please ensure that this is an VRS annotated VCF"
        )

    # try to populate from Bash env variable
    if not index_path:
        assert (
            "VRS_VCF_INDEX" in os.environ
        ), "no genotype index specified, no index path was provided nor was a variable name VRS_VCF_INDEX found."
        index_path = Path(os.environ.get("VRS_VCF_INDEX"))

    if index_path:
        # if index provided, use it to get VCF row
        index_path = Path(index_path)

        # find variant of interest
        for _, chr, pos in fetch_by_vrs_ids([variant_id], index_path):
            # TODO [ISSUE-103]: generalize VCF fixtures
            for record in vcf.fetch(chr, pos - 1, pos):
                if variant_id in record.info["VRS_Allele_IDs"]:
                    return record

        raise KeyError(f"no VCF row found matching variant ID {variant_id}")
    else:
        # otherwise, iterate through VCF
        for record in enumerate(vcf.fetch()):
            print(
                "no VCF index specified, iterating through VCF to locate variant of interest"
            )
            if variant_id in record.info["VRS_Allele_IDs"]:
                return record

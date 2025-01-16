import io
import json
import os
import pandas as pd
import pysam

from firecloud import api as fapi
from plugin_system.plugins.base_plugin import BasePlugin
from plugin_system.utils import (
    load_dict,
    csv_to_dataframe,
    terra_data_table_to_dataframe,
)


class GregorPlugin(BasePlugin):
    """
    Plugin for GREGoR U08 release data on Terra
    """

    def __init__(
        self, phenotype_table_path: str | None = None, index_path: str | None = None
    ):
        self.phenotype_index = self.__create_phenotype_index__(
            phenotype_table_path=phenotype_table_path, index_path=index_path
        )

    def __create_phenotype_index__(
        self, phenotype_table_path: str | None = None, index_path: str | None = None
    ) -> dict[str, list[str] | set[str]]:
        """given phenotypical data input specified by the GREGoR Data model (in either tsv/csv/Terra data table),
        return a dictionary mapping from each sample to its list of phenotypes

        Args:
            phenotype_table_path (str, optional): Path to csv/tsv of phenotype data specified by the GREGoR data model.
                    When not specified, defaults to loading from Terra data table in existing workspace titled "phenotypes".
                    For more info on the data model, see https://gregorconsortium.org/data-model
            index_path (str, optional): Path to pre-computed index. Defaults to None.

        Returns:
            dict[str, list[str]]: index of a sample id to sample's phenotypes. For example: {"patient_A": ["lactose intolerance", "anxiety"], "patient_B": ["shortness of breath"]}
        """

        # load index from file if already created
        if index_path is not None:
            return load_dict(index_path)

        # if no path specified, load phenotype table from Terra Data Table by default (must be in Terra workspace)
        if phenotype_table_path is None:
            phenotype_df = terra_data_table_to_dataframe(table_name="phenotype")
        else:  # otherwise load phenotype data table from file
            phenotype_df = csv_to_dataframe(phenotype_table_path)

        # create participant to phenotypes mapping
        phenotype_index = {}
        for participant_id in phenotype_df["participant_id"].unique():
            all_phenotypes = phenotype_df[
                phenotype_df["participant_id"] == participant_id
            ]["term_id"]

            phenotype_index[participant_id] = list(all_phenotypes.unique())

        return phenotype_index

    def include_sample(
        self, sample_id: str, record: pysam.VariantRecord, phenotype: str
    ) -> bool:
        """determine whether to include a sample in the cohort allele frequency based on its variant data and phenotypic traits

        Returns:
            bool: whether to include the sample
        """
        has_specified_phenotype = (
            sample_id in self.phenotype_index
            and phenotype in self.phenotype_index[sample_id]
        )

        # TODO: possibly implement filtering by sex of participant
        return has_specified_phenotype

    def process_sample_genotype(
        self,
        sample_id: str,
        record: pysam.VariantRecord,
        alt_index: int,
    ) -> tuple[int, int]:
        """given a sample's genotype, return focus and locus allele counts

        Args:
            sample_id (str): sample_id used to uniquly identify a sample ID
            record (pysam.VariantRecord): pysam record object representing a VCF row
            sample_phenotype_index (dict[str, list[str]]): mapping from sample IDs to each sample list of phenotypes
            alt_index (int): index matching the variant of interest

        Returns:
            tuple[int, int]: number of focus (specified) alleles, followed by number of locus (total) alleles.
        """

        # FIXME: support for hemizygous regions (chrY / mitochondrial variants)
        # GREGOR uses DRAGEN's continuous allele frequency approach:
        # https://support-docs.illumina.com/SW/DRAGEN_v40/Content/SW/DRAGEN/MitochondrialCalling.htm
        within_hemizygous_region = record.chrom in ["chrM", "chrY"]
        within_x_chr = record.chrom == "chrX"

        # get focus allele count, handling if there are multiple alts
        alleles = record.samples[sample_id].allele_indices
        num_focus_alleles = sum(
            [1 for _, alt_number in enumerate(alleles) if alt_number == alt_index]
        )

        if within_hemizygous_region and num_focus_alleles > 0:
            num_focus_alleles = 1

        # get total allele count
        if within_hemizygous_region:
            num_total_alleles = 1
        elif within_x_chr:
            # FIXME: make use of sex of participant?
            # by default considers all chrX samples as diploid (regardless of sex)
            is_female = True
            num_total_alleles = 2 if is_female else 1
        else:
            num_total_alleles = len(alleles)

        return num_focus_alleles, num_total_alleles

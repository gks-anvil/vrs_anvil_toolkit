import io
import json
import os
import pandas as pd
import pysam

from firecloud import api as fapi
from plugin_system.plugins.base_plugin import BasePlugin
from plugin_system.utils import load_dict, csv_to_dataframe, terra_data_table_to_dataframe


class GregorPlugin(BasePlugin):
    """
    Plugin for GREGoR U08 release data on Terra
    """

    def __init__(self, phenotype_table_path: str | None = None, index_path: str | None = None):
        self.phenotype_index = self.get_patient_phenotype_index(phenotype_table_path=phenotype_table_path, index_path=index_path)

    def get_patient_phenotype_index(
        self, phenotype_table_path: str = None, index_path: str = None
    ) -> dict[str, list | set]:
        """get index of phenotypes associated with a patient""

        Args:
            phenotype_table (str, optional): Path to csv/tsv of phenotype data specified by the GREGoR data model.
                    Defaults to pulling from a Terra data table in existing workspace titled "phenotypes".
                    For more info on the data model, see https://gregorconsortium.org/data-model
            cached_dict (str, optional): Path to cached dictionary to use. Defaults to None.

        Returns:
            dict: phenotypes associated with each participant
        """

        # load from cache if exists
        if index_path is not None:
            return load_dict(index_path)

        # otherwise generate it
        return self.__create_phenotype_index__(phenotype_table_path=phenotype_table_path)

    def __create_phenotype_index__(self, phenotype_table_path: str = None) -> dict[str, list[str] | set[str]]:
        """given phenotypical data input specified by the GREGoR Data model (in either tsv/csv or Terra data table),
        return a dict mapping from each sample id to the sample's list of phenotypes

        Returns:
            dict[str, list[str]]: index of a sample id to sample's phenotypes
        """

        # load table from within workspace using Terra Data Table or from csv/tsv
        if phenotype_table_path is None:
            phenotype_df = terra_data_table_to_dataframe(table_name="phenotype")
        else:
            phenotype_df = csv_to_dataframe(phenotype_table_path)

        # create participant to phenotypes mapping
        phenotype_index = {}
        for participant_id in phenotype_df["participant_id"].unique():
            all_phenotypes = phenotype_df[phenotype_df["participant_id"] == participant_id][
                "term_id"
            ]

            phenotype_index[participant_id] = list(all_phenotypes.unique())

        return phenotype_index

    def include_sample(self, sample_id: str, record: pysam.VariantRecord, phenotype: str) -> bool:
        """given a sample id and its genotype and phenotype of interest, determine whether to include it in the allele counts

        Returns:
            bool: whether to include the sample
        """
        has_specified_phenotype = (
            sample_id in self.phenotype_index and phenotype in self.phenotype_index[sample_id]
        )

        # TODO: possibly implement filtering by sex of participant
        return has_specified_phenotype

    def process_sample_genotype(
        self,
        sample_id: str,
        record: pysam.VariantRecord,
        alt_index: int,
    ) -> tuple[int, int]:
        """given a genotype for a particular sample, count the genotype

        Args:
            sample (pysam.VariantRecordSample): _description_
            patient_phenotype_index (dict[str, list[str]]): _description_

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

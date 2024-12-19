import json
import os
import pandas as pd
import pysam

from plugin_system.base_plugin import BasePlugin


class GregorPlugin(BasePlugin):
    """
    Implementation of methods to aggregrate allele frequency GREGoR data
    """

    def create_sample_phenotype_index(
        self, phenotype_table: str = None, as_set: bool = True, save_path: str = None
    ) -> dict[str, list[str] | set[str]]:
        """given phenotypical data input specified by the GREGoR Data model (in either tsv/csv or Terra data table),
        return a dict mapping from each sample id to the sample's list of phenotypes

        Returns:
            dict[str, list[str]]: index of a sample id to sample's phenotypes
        """

        if phenotype_table is None:
            # if unspecified, ensure valid Terra environment
            for env_key in ["WORKSPACE_NAMESPACE", "WORKSPACE_NAME"]:
                assert (
                    env_key in os.environ
                ), f"ERROR: No {env_key} key found in environmental variables in the Terra workspace. If you are working in a Terra workspace, please ensure both a WORKSPACE_NAMESPACE and a WORKSPACE_NAME are specified."

            # create dataframe from Terra data table
            # https://github.com/broadinstitute/fiss/blob/master/firecloud/api.py
            try:
                response = fapi.get_entities_tsv(
                    os.environ["WORKSPACE_NAMESPACE"],
                    os.environ["WORKSPACE_NAME"],
                    "phenotype",
                    model="flexible",
                )
                response.raise_for_status()
            except Exception as e:
                if response.json() and e in response.json():
                    error_message = response.json()["message"]
                else:
                    error_message = e
                print(
                    f"Error while loading phenotype data table from workspace: \n{error_message}"
                )

            phenotype_tsv = io.StringIO(response.text)
            phenotype_df = pd.read_csv(phenotype_tsv, sep="\t")
        else:
            # table path specified, parse using that table
            with open(phenotype_table, "r") as file:
                if phenotype_table.endswith(".csv"):
                    phenotype_df = pd.read_csv(file)
                elif phenotype_table.endswith(".tsv"):
                    phenotype_df = pd.read_csv(file, sep="\t")
                else:
                    raise Exception(
                        "Only csv and tsv file types implemented for phenotype table"
                    )

        # create patient to phenotype dict
        phenotype_index = {}
        for participant_id in phenotype_df["participant_id"].unique():
            phenotypes = phenotype_df[phenotype_df["participant_id"] == participant_id][
                "term_id"
            ].unique()

            phenotype_index[participant_id] = (
                set(phenotypes) if as_set else list(phenotypes)
            )

        # save to disk if specified
        if save_path:
            if os.path.exists(save_path):
                raise Exception(
                    f"index already exists at path: {save_path}. Please delete it before continuing"
                )
            else:
                with open(save_path, "w") as file:

                    # make it serializable so can be written to disk
                    if as_set:
                        for patient, pheno_set in phenotype_index.items():
                            phenotype_index[patient] = list(pheno_set)

                    json.dump(phenotype_index, file)

        return phenotype_index

    def include_sample(self, sample_id: str) -> bool:
        """given a sample id, determine whether to include it in the allele counts

        Returns:
            bool: whether to include the sample
        """
        raise NotImplementedError(
            "Plugins must implement include include_sample method"
        )

    def process_sample_genotype(
        self,
        record: pysam.VariantRecord,
        sample_id: str,
        sample_phenotype_index: dict[str, list[str]],
        alt_index: int,
    ) -> tuple[int, int]:
        """given a genotype for a particular sample, count the genotype

        Args:
            sample (pysam.VariantRecordSample): _description_
            patient_phenotype_index (dict[str, list[str]]): _description_

        Returns:
            tuple[int, int]: number of focus (specified) alleles, followed by number of locus (total) alleles.
        """

        # FIXME: support for hemizygous regions (chrX / mitochondrial variants)
        # GREGOR makes use of DRAGEN's continuous allele frequency approach:
        # https://support-docs.illumina.com/SW/DRAGEN_v40/Content/SW/DRAGEN/MitochondrialCalling.htm
        within_hemizygous_region = record.chrom in ["chrM", "chrY"]
        within_x_chr = record.chrom == "chrX"

        # increment focus allele count, handling multiple alts edge case
        alleles = record.samples[sample_id].allele_indices
        num_focus_alleles = sum(
            [1 for _, alt_number in enumerate(alleles) if alt_number == alt_index]
        )
        if within_hemizygous_region and num_focus_alleles > 0:
            num_focus_alleles = 1

        # increment total allele count
        if within_hemizygous_region:
            num_total_alleles = 1
        elif within_x_chr:
            # FIXME: make use of sex of participant?
            # by default considers all chrX samples as diploid
            is_female = True
            num_total_alleles = 2 if is_female else 1
        else:
            num_total_alleles = len(alleles)

        return num_focus_alleles, num_total_alleles

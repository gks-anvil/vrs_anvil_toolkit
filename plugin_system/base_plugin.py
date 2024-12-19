import pysam


class BasePlugin:
    """
    Interface for caf generation plugins
    """

    is_plugin = True

    def create_sample_phenotype_index(self) -> dict[str, list[str] | set[str]]:
        """given any arbitrary phenotypical data input, return a dict mapping from each sample id to the sample's list of phenotypes

        Returns:
            dict[str, list[str]]: index of a sample id to sample's phenotypes
        """
        raise NotImplementedError(
            "Plugins must implement sample-to-phenotype index method"
        )

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
        participant_id: str,
        sample_phenotype_index: dict[str, list[str]],
    ) -> tuple[int, int]:
        """given a genotype for a particular sample, count the genotype

        Args:
            sample (pysam.VariantRecordSample): _description_
            patient_phenotype_index (dict[str, list[str]]): _description_

        Returns:
            tuple[int, int]: number of focus (specified) alleles, followed by number of locus (total) alleles.
        """
        raise NotImplementedError(
            "Plugins must implement the process_sample_genotype method"
        )

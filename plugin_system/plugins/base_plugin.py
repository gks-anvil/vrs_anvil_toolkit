import pysam


class BasePlugin:
    """
    Interface for caf generation plugins
    """

    __is_plugin__ = True

    def create_sample_phenotype_index(self) -> dict[str, list[str] | set[str]]:
        """given any arbitrary phenotypical data input, return a dict mapping from each sample id to the sample's list of phenotypes

        Returns:
            dict[str, list[str] | set[str]]: index of a sample id to sample's phenotypes
        """
        raise NotImplementedError(
            "Plugins must implement sample-to-phenotype index method"
        )

    def include_sample(self, sample_id: str, record: pysam.VariantRecord, phenotype: str, sample_phenotype_index: dict[str, list[str]]) -> bool:
        """given a sample id, determine whether to include it in the allele counts

        Args:
            sample_id (str): sample_id used to uniquly identify a sample ID
            record (pysam.VariantRecord): PySam record object representing a VCF row
            phenotype (str): phenotype of interest, matching phenotype codes used in sample_phenotype_)index
            sample_phenotype_index (dict[str, list[str]]): mapping from sample IDs to each sample list of phenotypes

        Returns:
            bool: whether to include the sample
        """
        return True

    def process_sample_genotype(
        self,
        sample_id: str,
        record: pysam.VariantRecord,
        alt_index: int,
    ) -> tuple[int, int]:
        """given a genotype for a particular sample, determine how to sum its alleles

        Args:
            sample_id (str): sample_id used to uniquly identify a sample ID
            record (pysam.VariantRecord): PySam record object representing a VCF row
            sample_phenotype_index (dict[str, list[str]]): mapping from sample IDs to each sample list of phenotypes
            alt_index (int): index matching the variant of interest

        Returns:
            tuple[int, int]: number of focus (specified) alleles, followed by number of locus (total) alleles.
        """

        # increment focus allele count, handling multiple alts edge case
        # for example, if the alt of interest is at index 2, then a genotype of
        # (1,2) would have a 1 focus allele out of 2 total alleles
        alleles = record.samples[sample_id].allele_indices
        num_focus_alleles = sum(
            [1 for _, alt_number in enumerate(alleles) if alt_number == alt_index]
        )
        num_total_alleles = len(alleles)

        return num_focus_alleles, num_total_alleles

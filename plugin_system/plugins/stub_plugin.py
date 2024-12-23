import pysam
from plugin_system.plugins.base_plugin import BasePlugin


class StubPlugin(BasePlugin):
    """
    <description of plugin>
    """

    # FIXME: implement!
    def create_sample_phenotype_index(self) -> dict[str, list[str] | set[str]]:
        """given any arbitrary phenotypical data input, return a dict mapping from each sample id to the sample's list of phenotypes

        Args:
            FIXME: for you to define!

        Returns:
            dict[str, list[str] | set[str]]: index of a sample id to sample's phenotypes
        """
        raise NotImplementedError(
            "Plugins must implement create_sample_phenotype_index method"
        )

    # FIXME: implement!
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
        raise NotImplementedError(
            "Plugins must implement include_sample method"
        )

    # FIXME: implement!
    def process_sample_genotype(
        self,
        sample_id: str,
        record: pysam.VariantRecord,
        sample_phenotype_index: dict[str, list[str]],
        alt_index: int,
    ) -> tuple[int, int]:
        """given a genotype for a particular sample, count the genotype

        Args:
            sample_id (str): sample_id used to uniquly identify a sample ID
            record (pysam.VariantRecord): PySam record object representing a VCF row
            sample_phenotype_index (dict[str, list[str]]): mapping from sample IDs to each sample list of phenotypes
            alt_index (int): index matching the variant of interest

        Returns:
            tuple[int, int]: number of focus (specified) alleles, followed by number of locus (total) alleles.
        """

        raise NotImplementedError(
            "Plugins must implement process_sample_genotype method"
        )

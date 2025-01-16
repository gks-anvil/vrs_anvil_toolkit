import pysam

from abc import ABC


class BasePlugin(ABC):
    """
    Interface for caf generation plugins
    """

    __is_plugin__ = True

    def __init__(self, phenotype_index: dict[str, list[str]] | None = None):
        self.phenotype_index = phenotype_index

    def get_phenotype_index(self) -> dict[str, list[str]]:
        """getter for a dictionary mapping from each sample to the sample's phenotypes

        Returns:
            dict[str, list[str]]: index mapping from sample id to sample's phenotypes. For example: {"patient_A": ["lactose intolerance", "anxiety"], "patient_B": ["shortness of breath"]}
        """
        return self.phenotype_index

    def include_sample(
        self, sample_id: str, record: pysam.VariantRecord, phenotype: str
    ) -> bool:
        """determine whether to include a sample in the cohort allele frequency based on its variant data and phenotypic traits

        Args:
            sample_id (str): sample_id used to uniquly identify a sample ID
            record (pysam.VariantRecord): PySam record object representing a VCF row
            phenotype (str): phenotype of interest, matching phenotype codes used in sample_phenotype_)index

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
        """given a sample's genotype, return focus and locus allele counts

        Args:
            sample_id (str): sample_id used to uniquly identify a sample ID
            record (pysam.VariantRecord): pysam record object representing a VCF row
            sample_phenotype_index (dict[str, list[str]]): mapping from sample IDs to each sample list of phenotypes
            alt_index (int): index matching the variant of interest

        Returns:
            tuple[int, int]: number of focus (specified) alleles, then number of locus (total) alleles.
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

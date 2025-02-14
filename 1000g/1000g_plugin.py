from plugin_system.plugins.base_plugin import BasePlugin
from plugin_system.utils import (
    load_dict,
    csv_to_dataframe,
    terra_data_table_to_dataframe,
)


class ThousandGenomesPlugin(BasePlugin):
    """
    Plugin for AnVIL 1000G PRIMED data release on Terra

    Link: https://anvil.terra.bio/#workspaces/anvil-datastorage/AnVIL_1000G_PRIMED-data-model

    Note that get_phenotype_index is inherited from the parent BasePlugin class.
    """

    def __init__(
        self, phenotype_table_path: str | None = None, index_path: str | None = None
    ):
        """constructor used to set a phenotype index if provided a file path for the index (index_path).
        Otherwise create a phenotype index using a Terra data table (no path specified) or with a csv/tsv filepath.

        Index example: {"sample_A": ["HP:0001263", "HP:0000002"], "sample_B": ["HP:0001263"]}
        Note that we actively do not use super() to invoke the BasePlugin's constructor to create custom functionality.

        Args:
            phenotype_table_path (str, optional): Path to csv/tsv of phenotype data specified by the GREGoR data model.
                When not specified, defaults to loading from Terra data table in existing workspace titled "phenotypes".
                For more info on the data model, see https://gregorconsortium.org/data-model. Defaults to None.
            index_path (str, optional): Path to existing phenotype index. Defaults to None.
        """

        self.phenotype_index = self.__create_phenotype_index(
            phenotype_table_path=phenotype_table_path, index_path=index_path
        )

    def __create_phenotype_index(
        self, phenotype_table_path: str | None = None, index_path: str | None = None
    ) -> dict[str, list[str]]:
        """[private method] given phenotypical data input specified by the GREGoR Data model (in either tsv/csv/Terra data table),
        return a dictionary mapping from each sample to its list of phenotypes

        Args:
            phenotype_table_path (str, optional): Path to csv/tsv of phenotype data specified by the GREGoR data model.
                    When not specified, defaults to loading from Terra data table in existing workspace titled "phenotypes".
                    For more info on the data model, see https://gregorconsortium.org/data-model
            index_path (str, optional): Path to pre-computed index. Defaults to None.

        Returns:
            dict[str, list[str]]: index of a sample id to sample's phenotypes.
        """

        # load index from file if already created
        if index_path is not None:
            return load_dict(index_path)

        # if no path specified, load phenotype table from Terra Data Table by default (must be in Terra workspace)
        if phenotype_table_path is None:
            phenotype_df = terra_data_table_to_dataframe(
                table_name="population_descriptor"
            )
        else:  # otherwise load phenotype data table from file
            phenotype_df = csv_to_dataframe(phenotype_table_path)

        # create participant to phenotypes mapping
        phenotype_index = {}
        for subject_id in phenotype_df["subject_id"].unique():
            all_phenotypes = phenotype_df[phenotype_df["subject_id"] == subject_id][
                "country_of_recruitment"
            ]

            phenotype_index[subject_id] = list(all_phenotypes.unique())

        return phenotype_index

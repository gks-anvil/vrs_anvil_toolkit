import pytest

from plugin_system.utils import terra_data_table_to_dataframe


@pytest.fixture()
def table_name():
    return "phenotype"


def test_terra_loads_gregor_phenotype_table_in_workspace(table_name):
    df = terra_data_table_to_dataframe(table_name)

    for field in ["participant_id", "presence", "term_id"]:
        assert (
            field in df.columns
        ), f"Could not find field named {field} in columns. \nAll columns: {df.columns}"

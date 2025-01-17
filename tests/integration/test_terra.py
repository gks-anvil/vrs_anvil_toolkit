import os
import pytest

from plugin_system.utils import terra_data_table_to_dataframe, WORKSPACE_ENV_KEYS


@pytest.fixture()
def table_name():
    return "phenotype"


@pytest.mark.skipif(
    any([os.getenv(key) is None for key in WORKSPACE_ENV_KEYS]),
    reason="can only be tested on Terra",
)
def test_terra_loads_gregor_phenotype_table_in_workspace(table_name):
    df = terra_data_table_to_dataframe(table_name)

    for field in ["participant_id", "presence", "term_id"]:
        assert (
            field in df.columns
        ), f"Could not find field named {field} in columns. \nAll columns: {df.columns}"

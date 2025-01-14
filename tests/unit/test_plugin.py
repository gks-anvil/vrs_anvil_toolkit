import pytest

from plugin_system.plugin_manager import PluginManager
from vrs_anvil.evidence import PLUGIN_DIR


@pytest.fixture()
def plugin_methods():
    return [
        "get_phenotype_index",
        "include_sample",
        "process_sample_genotype",
    ]


@pytest.fixture()
def plugin_names():
    return ["BasePlugin", "GregorPlugin"]


def test_plugin_manager_can_load_base_and_gregor_plugin(plugin_methods, plugin_names):
    for plugin_name in plugin_names:
        plugin_manager = PluginManager(PLUGIN_DIR)
        plugin = plugin_manager.load_plugin(plugin_name)

        print("plugin:", plugin)

        for method in plugin_methods:
            assert hasattr(
                plugin, method
            ), f"{plugin_name} does not have method '{method}'"

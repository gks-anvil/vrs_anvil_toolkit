# plugin_system/plugin_manager.py
import pkgutil
import importlib
import os

from plugin_system.plugins.base_plugin import BasePlugin
from vrs_anvil.evidence import PLUGIN_MODULE_PATH


class PluginManager:
    def load_plugin(self, plugin_name: str) -> BasePlugin:
        """grab first plugin class matching name `plugin_name` within the `plugin_package` directory

        Args:
            plugin_name (str): name of the plugin class

        Raises:
            OSError: unable to find class named `plugin_name`
            ImportError: unable to import package

        Returns:
            BasePlugin: non-instantiated Plugin class
        """

        # get full path for the default plugin directory
        default_plugin_dir = os.path.dirname(
            importlib.import_module(PLUGIN_MODULE_PATH).__file__
        )

        # look for plugin by name first in default directory then in top-level directory
        for i, iter in enumerate(
            [pkgutil.iter_modules([default_plugin_dir]), pkgutil.iter_modules()]
        ):
            for _, name, _ in iter:
                # only look for specific file names
                if not name.endswith("_plugin"):
                    continue

                # get full module name, use only plugin name if iterating across top-level dirs
                module_name = f"{PLUGIN_MODULE_PATH}.{name}" if i == 0 else name

                try:
                    # dynamically import the module
                    module = importlib.import_module(module_name)

                    # grab first plugin class matching name <plugin_name>
                    for attribute_name in dir(module):
                        if attribute_name == plugin_name:
                            attribute = getattr(module, attribute_name)
                            if isinstance(attribute, type) and hasattr(
                                attribute, "__is_plugin__"
                            ):
                                return attribute
                except ImportError as e:
                    # plugin modules found but unable to load specified plugin
                    print(f"Error loading plugin {name}: {e}")
                    raise

        # if no plugin found, raise error
        raise OSError(
            f"Plugin {plugin_name} not found. Make sure the path is stored in the top-level"
        )

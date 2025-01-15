# plugin_system/plugin_manager.py
import pkgutil
import importlib
import os

from plugin_system.plugins.base_plugin import BasePlugin


class PluginManager:
    def __init__(self, plugin_package: str):
        self.plugin_package = plugin_package

    def load_plugin(self, plugin_name: str) -> BasePlugin:
        # find and import all modules in the specified package
        package_path = os.path.dirname(
            importlib.import_module(self.plugin_package).__file__
        )

        for _, name, _ in pkgutil.iter_modules([package_path]):
            try:
                # dynamically import the module
                full_module_name = f"{self.plugin_package}.{name}"
                module = importlib.import_module(full_module_name)

                # grab first plugin class matching name <plugin_name>
                for attribute_name in dir(module):
                    if attribute_name == plugin_name:
                        attribute = getattr(module, attribute_name)
                        if isinstance(attribute, type) and hasattr(
                            attribute, "__is_plugin__"
                        ):
                            return attribute
            except ImportError as e:
                print(f"Error loading plugin {name}: {e}")
                raise

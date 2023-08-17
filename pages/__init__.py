import os
import importlib

# List all files in the current directory
for module_file in os.listdir(os.path.dirname(__file__)):
    # Check if it's a python file and not this __init__ file
    if module_file.endswith('.py') and module_file != '__init__.py':
        module_name = module_file[:-3]  # remove the .py extension
        module = importlib.import_module('.' + module_name, package=__name__)

        # Add the page function to the current module's namespace
        globals()[module_name] = module.page

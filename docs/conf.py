# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys

from pathlib import Path
from typing import Literal

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))


# To include notebooks from example folder in the documentation
link_from = str(Path(__file__).resolve().parent.parent / "examples")
link_to = str(Path(__file__).resolve().parent / "examples")

os.symlink(link_from, link_to, target_is_directory=True)


project = "StereoMolGraph"
copyright = "2025, Maxim Papusha"
author = "Maxim Papusha"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

nb_execution_mode: Literal["off", "force", "auto", "cache", "inline"] = "off"
nb_execution_show_tb = True  # Show traceback in the output
nb_execution_raise_on_error = True # Critical - makes exceptions fail the build

nb_execution_allow_errors = False  # Don't allow errors in the output
nb_execution_timeout = 300  # Timeout in seconds

extensions = [
    # Core Sphinx extensions first
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    # Type hints should come after autodoc
    "sphinx_autodoc_typehints",
    # MyST extensions
    "myst_nb",
    # Theme and UI extensions
    "sphinx_rtd_theme",
    "sphinx_copybutton",
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
}

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

autodoc_member_order: Literal["alphabetical", "bysource", "groupwise"] = (
    "groupwise"
)
# autodoc_class_signature = "separated"
# Force type hints to be links when possible
typehints_defaults: Literal["comma", "braces", "braces-after"] = "braces-after"

# (Optional) Make the links more readable
typehints_document_rtype_none = False
typehints_use_signature = True
typehints_use_signature_return = True

rst_prolog = """
:github_url: https://github.com/maxim-papusha/StereoMolGraph
"""

# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
from datetime import datetime

from ansys_sphinx_theme import get_version_match

import lamberthub

# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = "lamberthub"
copyright = f"{datetime.now().year} Jorge Martínez Garrido"
author = "Jorge Martínez Garrido"
release = version = lamberthub.__version__
cname = os.getenv("DOCUMENTATION_CNAME", "lamberthub.docs.jorgemartinez.space")

# -- Options for HTML output -------------------------------------------------

html_logo = "_static/logo.png"
html_theme = "ansys_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

html_context = {
    "github_user": "jorgepiloto",
    "github_repo": "lamberthub",
    "github_version": "main",
    "doc_path": "doc/source",
}

html_theme_options = {
    "github_url": "https://github.com/jorgepiloto/lamberthub",
    "use_edit_page_button": True,
    "contact_mail": "contact@jorgemartinez.space",
    "additional_breadcrumbs": [
        ("Jorge Martínez website", "https://jorgemartinez.space"),
    ],
    "switcher": {
        "json_url": f"https://{cname}/release/versions.json",
        "version_match": get_version_match(version),
    },
    "navbar_end": ["version-switcher", "theme-switcher", "navbar-icon-links"],
}

html_short_title = html_title = "lamberthub"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "autoapi.extension",
    "sphinx.ext.autodoc",
    #"sphinx.ext.mathjax",
    "sphinx.ext.autosummary",
    "numpydoc",
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",
    "nbsphinx",
    "sphinx_gallery.load_style",
    "myst_parser",
    "jupyter_sphinx",
    "sphinx_design",
]

# Source files
source_suffix = {
    ".rst": "restructuredtext",
    ".mystnb": "jupyter_notebook",
    ".md": "markdown",
}

# Index document
master_doc = "index"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
if os.getenv("LAMBERTHUB_SKIP_PERFORMANCE") == "true":
    exclude_patterns = ["source/explanations/performance_comparison.mystnb"]
else:
    exclude_patterns = []

# -- Options for Sphinx autoapi ----------------------------------------------

autoapi_type = "python"
autoapi_dirs = ["../../src/lamberthub"]
autodoc_typehints = "none"
autoapi_template_dir = "_autoapi_templates"
exclude_patterns.append("_autoapi_templates/index.rst")
exclude_patterns.append("_autoapi_templates/python/module.rst")

nbsphinx_custom_formats = {
    ".mystnb": ["jupytext.reads", {"fmt": "mystnb"}],
}

# Custom thumbnails for gallery of examples
nbsphinx_thumbnails = {
    "tutorials/gauss1809_solver": "_static/tutorials/gauss_thumbnail.png"
}

# The performance comparison takes a bit long. This avoids the documentation to
# assume something has failed due to long time computations
nb_execution_timeout = 900

# MathJax config
# See https://github.com/spatialaudio/nbsphinx/issues/572#issuecomment-853389268
#mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
mathjax_path="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"


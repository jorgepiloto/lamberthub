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

import lamberthub

# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = lamberthub.__name__
copyright = "2022, Jorge Martínez Garrido"
author = "Jorge Martínez Garrido"

# The full version, including alpha/beta/rc tags
release = lamberthub.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "autoapi.extension",
    "myst_nb",
    "nbsphinx",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_gallery.load_style",
]

# Source files
source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
    ".myst": "myst-nb",
}

# Index document
master_doc = "index"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
if os.environ.get("LAMBERTHUB_SKIP_PERFORMANCE") == "True":
    exclude_patterns = ["source/explanations/performance_comparison.md"]
else:
    exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"
# html_theme_options = {
#        "toc_title": "Sections in this page",
#        "extra_navbar": "",
# }
html_title = "lamberthub"
# html_logo = "_static/lamberts_problem_geometry.png"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["custom.css"]

# -- Options for Sphinx autoapi ----------------------------------------------

autoapi_type = "python"
autoapi_dirs = ["../../src"]
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

# Custom mathjax configuration
myst_update_mathjax = False

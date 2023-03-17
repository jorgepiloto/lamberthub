"""Sphinx configuration file."""

from datetime import datetime
import os

from ansys_sphinx_theme import get_version_match

import lamberthub


# Project information
project = lamberthub.__name__
copyright = f"(c) {datetime.now().year} Jorge Martinez Garrido. All rights reserved"
author = lamberthub.__author__
release = version = lamberthub.__version__
cname = os.getenv("DOCUMENTATION_CNAME", "lamberthub.docs.jorgemartinez.space")


# Extensions and source file configuration
extensions = [
    "autoapi.extension",
    "myst_nb",
    "nbsphinx",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_gallery.load_style",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
    ".myst": "myst-nb",
}

master_doc = "index"
templates_path = ["_templates"]
exclude_patterns = []


# HTML theme configuration
html_theme = "ansys_sphinx_theme"
html_short_title = html_title = project
html_context = {
    "github_user": "lamberthub",
    "github_repo": "lamberthub",
    "github_version": "main",
    "doc_path": "doc/source",
}
html_theme_options = {
    "github_url": "https://github.com/lamberthub/lamberthub",
    "contact_mail": "contact@jorgemartinez.space",
    "additional_breadcrumbs": [
        ("Jorge's Website", "https://jorgemartinez.space"),
    ],
    "switcher": {
        "json_url": f"https://{cname}/versions.json",
        "version_match": get_version_match(version),
    },
    "check_switcher": False,
}
html_static_path = ["_static"]
html_css_files = ["custom.css"]


# Autoapi configuration
autoapi_type = "python"
autoapi_dirs = ["../../src"]
autodoc_typehints = "none"
autoapi_template_dir = "_autoapi_templates"
exclude_patterns.append("_autoapi_templates/index.rst")
exclude_patterns.append("_autoapi_templates/python/module.rst")

# Nbsphinx configuration
nbsphinx_custom_formats = {
    ".mystnb": ["jupytext.reads", {"fmt": "mystnb"}],
}
nbsphinx_thumbnails = {
    "tutorials/gauss1809_solver": "_static/tutorials/gauss_thumbnail.png"
}
nb_execution_timeout = 900
myst_update_mathjax = False

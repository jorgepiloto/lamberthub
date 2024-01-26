"""Sphinx configuration file."""

from datetime import datetime
import os
from pathlib import Path

from ansys_sphinx_theme import get_autoapi_templates_dir_relative_path, get_version_match

import lamberthub

# Project information
project = lamberthub.__name__
copyright = (
    f"(c) {datetime.now().year} Jorge Martinez Garrido. All rights reserved"
)
author = lamberthub.__author__
release = version = lamberthub.__version__
cname = os.getenv("DOCUMENTATION_CNAME", "lamberthub.docs.jorgemartinez.space")


# Extensions and source file configuration
extensions = [
    "autoapi.extension",
    "nbsphinx",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_gallery.load_style",
]

source_suffix = {
    ".rst": "restructuredtext",
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
autoapi_dirs = ["../../src/lamberthub"]
autoapi_root = "api"
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    "special-members",
]
autoapi_template_dir = get_autoapi_templates_dir_relative_path(Path(__file__))
autoapi_keep_files = True
autoapi_render_in_single_page = ["class", "enum", "exception", "function"]

# Nbsphinx configuration
nbsphinx_custom_formats = {
    ".mystnb": ["jupytext.reads", {"fmt": "mystnb"}],
}
nbsphinx_thumbnails = {
    "tutorials/gauss1809_solver": "_static/tutorials/gauss_thumbnail.png"
}
nb_execution_timeout = 900
myst_update_mathjax = False


# -- Declare the Jinja context -----------------------------------------------
exclude_patterns = []
BUILD_API = True if os.environ.get("BUILD_API", "true") == "true" else False
if not BUILD_API:
    exclude_patterns.append("api")

BUILD_EXAMPLES = True if os.environ.get("BUILD_EXAMPLES", "true") == "true" else False
if not BUILD_EXAMPLES:
    exclude_patterns.append("examples/**")
    exclude_patterns.append("examples.rst")

jinja_contexts = {
    "main_toctree": {
        "build_api": BUILD_API,
        "build_examples": BUILD_EXAMPLES,
    },
}

def prepare_jinja_env(jinja_env) -> None:
    """
    Customize the jinja env.

    Notes
    -----
    See https://jinja.palletsprojects.com/en/3.0.x/api/#jinja2.Environment
    """
    jinja_env.globals["project_name"] = project

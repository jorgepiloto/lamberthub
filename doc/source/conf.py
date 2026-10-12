"""Configure the lamberthub documentation build."""

from datetime import UTC, datetime
from importlib.metadata import version as distribution_version
import os

from sphinx_unbloated_theme import get_autoapi_templates_dir

project = "lamberthub"
author = "Jorge Martinez"
copyright = f"{datetime.now(tz=UTC).year}, Jorge Martinez"
release = version = distribution_version("lamberthub")

extensions = [
    "autoapi.extension",
    "numpydoc",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinx_design",
]

intersphinx_mapping = {
    "numpy": ("https://numpy.org/doc/stable/", None),
    "python": ("https://docs.python.org/3", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

autoapi_dirs = ["../../src/lamberthub"]
autoapi_root = "api-reference"
autoapi_template_dir = get_autoapi_templates_dir()
autoapi_add_toctree_entry = False
autoapi_options = [
    "members",
    "show-inheritance",
    "show-module-summary",
    "imported-members",
]
numpydoc_show_class_members = False

html_theme = "sphinx_unbloated_theme"
html_title = "lamberthub"
html_short_title = "lamberthub"
html_baseurl = os.getenv("DOCUMENTATION_BASE_URL", "")
documentation_root_url = os.getenv(
    "DOCUMENTATION_ROOT_URL",
    "https://jorgemartinez.space/projects/lamberthub",
)
html_theme_options = {
    "site_description": "Lambert's problem solvers for orbital mechanics.",
    "footer_text": (
        f"&copy; {datetime.now(tz=UTC).year} Jorge Martinez. "
        "Built with the Sphinx Unbloated Theme."
    ),
    "github_url": "https://github.com/jorgepiloto/lamberthub",
    "versions_url": f"{documentation_root_url}/versions.json",
}
html_additional_pages = {"404": "404.html"}

master_doc = "index"
exclude_patterns = ["_build", "conf.py"]
source_suffix = {
    ".rst": "restructuredtext",
}

autodoc_typehints = "description"
nitpicky = True
nitpick_ignore = [
    ("py:class", "numpy.ndarray"),
]
suppress_warnings = ["ref.citation"]
linkcheck_ignore = [
    r"https://doi\.org/10\.1177/0954410016676847",
]

latex_documents = [
    (
        "index",
        "lamberthub.tex",
        "lamberthub Documentation",
        author,
        "manual",
    ),
]

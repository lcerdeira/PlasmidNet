# Configuration file for Sphinx documentation builder.
project = "PlasmidNet"
copyright = "2026, Louise Cerdeira"
author = "Louise Cerdeira"
release = "1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_parser",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

templates_path = ["_templates"]
exclude_patterns = ["_build"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_logo = "../assets/logo.png"
html_theme_options = {
    "logo_only": False,
    "display_version": True,
    "navigation_depth": 3,
}

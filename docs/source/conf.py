import os
import sys

sys.path.insert(0, os.path.abspath('../..'))

project = 'DenseArrayToolkit'
copyright = '2025, DAT v1.0.0'
author = 'Seismology Group'


extensions = [
    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]

myst_url_schemes = ["http", "https", "mailto"]

html_theme = 'sphinx_rtd_theme'

html_static_path = ['_static']

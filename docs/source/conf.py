import os
import sys

# 添加项目根目录到Python路径
sys.path.insert(0, os.path.abspath('../..'))

# 项目信息
project = 'DenseArrayToolkit'
copyright = '2024, DAT v1.0.0'
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

# 主题设置
html_theme = 'sphinx_rtd_theme'

# 静态文件路径
html_static_path = ['_static']

# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


import os
on_rtd = os.environ.get('READTHEDOCS') == 'True'

if on_rtd:
# import mock
# import sys 
# for mod_name in MOCK_MODULES:
#     sys.modules[mod_name] = mock.Mock()

    MOCK_MODULES = ['numpy', 'scipy', 'matplotlib', 'plotly', 'plotly.subplots', 'plotly.graph_objs', 'seaborn', 'matplotlib.pyplot', 'matplotlib.lines', 'scipy.interpolate', 'scipy.stats', 'streamlit', 'pandas' ]
    autodoc_mock_imports = MOCK_MODULES #[ "plotly", "plotly.graph_objs", "plotly.subplots" ]


# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

else:

    import logging
    import sys
    src = os.path.abspath("..")
    src = os.path.dirname(src)
    logging.critical( f"{src=}" )
    sys.path.insert(0, src)


# -- Project information -----------------------------------------------------

project = 'qpcr'
copyright = '2022, Noah Henrik Kleinschmidt'
author = 'Noah Henrik Kleinschmidt'

# The full version, including alpha/beta/rc tags
release = '4.0.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
                'nbsphinx',
                'sphinx.ext.mathjax',
                'sphinx.ext.napoleon',
                'sphinx.ext.viewcode',
                'sphinx.ext.autodoc',
            ]

napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_admonition_for_notes = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

html_logo = "./qpcr_dark.svg"

html_theme_options = dict( style_nav_header_background = "#00aaffff" ) 

# colorful
pygments_style = "one-dark"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
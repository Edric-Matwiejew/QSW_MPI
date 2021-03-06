# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
#from unittest.mock import MagicMock
import sys
sys.path.insert(0, os.path.abspath("../.."))

autodoc_mock_imports = ["mpi4py", "qsw_mpi.foperators", "qsw_mpi.fMPI", "qsw_mpi.parallel_io", "h5py"]

numfig = True

#class Mock(MagicMock):
#    @classmethod
#    def __getattr__(cls, name):
#        return MagicMock()
#
#MOCK_MODULES = ['mpi4py', 'qsw_mpi.foperators', 'qsw_mpi.fMPI']
#sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)

# -- Project information -----------------------------------------------------

project = 'QSW_MPI'
copyright = '2019, Edric Matwiejew'
author = 'Edric Matwiejew'

# The full version, including alpha/beta/rc tags
release = '1.0.2'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc']

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
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


master_doc = 'index'


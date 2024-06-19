# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html



# -- Project information -----------------------------------------------------

project = 'ngstents'
copyright = '2023, Jay Gopalakrishnan'
author = 'Jay Gopalakrishnan'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["sphinx.ext.autodoc","sphinx.ext.mathjax","sphinx.ext.todo","sphinx.ext.githubpages",
              "IPython.sphinxext.ipython_console_highlighting", "IPython.sphinxext.ipython_directive",
              "jupyter_sphinx",
              "nbsphinx",
              "myst_parser",
              "sphinxemoji.sphinxemoji",
              "myst_nb.sphinx_ext",
              ]

# Run notebook configuration

# The template used when exporting from nbconvert
#   full  - Outputs the full HTML document [Default]
#   basic - Outputs a single div (with no additional resources)
run_notebook_export_template = 'basic'  # Default: 'full'

# Display the source links to the generated evaluated files
run_notebook_display_source_links = False  # Default: True

# Whether or not to evaluate the notebooks prior to embedding them
evaluate_notebooks = True  # Default: True

# START nbsphinx stuff
#increase timeout for cell execution, since some files take long to execute
nbsphinx_timeout = 100000

# If True, the build process is continued even if an exception occurs:
nbsphinx_allow_errors = False

# This is processed by Jinja2 and inserted before each notebook
nbsphinx_prolog = r"""
.. raw:: html

    <style>
        .p-Widget {
            height: 400px;
        }
        .dg.main {
            margin-left: 0px;
        }
        div.p-Widget div div div div.dg ul li {
            list-style: none;
            margin-left: 0px;
        }
        div.p-Widget div div div div.dg ul li div.dg {
            margin-bottom: 0px;
        }

        .MathJax { font-size: 0.9em !important; }
    </style>

.. only:: html
    .. role:: raw-html(raw)
        :format: html
"""

nbsphinx_widgets_path = ""
# END nbsphinx stuff

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'paper', 'env', 'jupyter_execute']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

html_context = {
  'display_github': True,
  'github_user': 'jayggg',
  'github_repo': 'ngstents',
}

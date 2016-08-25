"""Central place for package metadata."""

# NOTE: We use __title__ instead of simply __name__ since the latter would
#       interfere with a global variable __name__ denoting object's name.
__title__ = 'resolwe-bio'
__summary__ = 'Bioinformatics pipelines for the Resolwe platform'
__url__ = 'https://github.com/genialis/resolwe-bio'

# Semantic versioning is used. For more information see:
# https://packaging.python.org/en/latest/distributing/#semantic-versioning-preferred
__version__ = '1.2.1'

__author__ = 'Genialis d.o.o.'
__email__ = 'dev-team@genialis.com'

__license__ = 'Apache License (2.0)'
__copyright__ = '2015-2016, ' + __author__

__all__ = (
    "__title__", "__summary__", "__url__", "__version__", "__author__",
    "__email__", "__license__", "__copyright__",
)

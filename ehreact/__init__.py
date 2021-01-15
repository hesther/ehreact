"""
ehreact
A Python package for extracting and scoring reaction templates based on extended Hasse diagrams.
"""

# Add imports here
import ehreact.arguments
import ehreact.train


# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

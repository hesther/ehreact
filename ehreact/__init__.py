"""
ehreact
A Python package for extracting and scoring reaction templates based on extended Hasse diagrams.
"""

import ehreact.arguments
import ehreact.train
import ehreact.predict
import ehreact.helpers

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

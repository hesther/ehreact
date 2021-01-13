"""
Unit and regression test for the ehreact package.
"""

# Import package, test suite, and other packages as needed
import ehreact
import pytest
import sys

def test_ehreact_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "ehreact" in sys.modules

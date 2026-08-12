import unittest

from isca import GreyCodeBase, IscaCodeBase, DryCodeBase, SocratesCodeBase, SocColumnCodeBase, GFDL_BASE

def test_grey_codebase():
    cb = GreyCodeBase.from_directory(GFDL_BASE)
    cb.compile()

def test_isca_codebase():
    cb = IscaCodeBase.from_directory(GFDL_BASE)
    cb.compile()

def test_dry_codebase():
    cb = DryCodeBase.from_directory(GFDL_BASE)
    cb.compile()

def test_socrates_codebase():
    # Exercises the full pipeline: fetching the vendored/default Socrates
    # version if it isn't already present locally (a sparse partial clone
    # from github.com/MetOffice/socrates), generating its path_names if not
    # already committed, and compiling.
    cb = SocratesCodeBase.from_directory(GFDL_BASE)
    cb.compile()

def test_socrates_column_codebase():
    cb = SocColumnCodeBase.from_directory(GFDL_BASE)
    cb.compile()

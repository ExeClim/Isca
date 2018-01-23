import unittest

from isca import GreyCodeBase, IscaCodeBase, DryCodeBase, GFDL_BASE

def test_grey_codebase():
    cb = GreyCodeBase.from_directory(GFDL_BASE)
    cb.compile()

def test_isca_codebase():
    cb = IscaCodeBase.from_directory(GFDL_BASE)
    cb.compile()

def test_dry_codebase():
    cb = DryCodeBase.from_directory(GFDL_BASE)
    cb.compile()

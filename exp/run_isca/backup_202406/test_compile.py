import unittest

from isca import GreyCodeBase, IscaCodeBase, DryCodeBase, GFDL_BASE, SocratesCodeBase

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
    cb = SocratesCodeBase.from_directory(GFDL_BASE)
    cb.compile(debug=False)

test_socrates_codebase()          

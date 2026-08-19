import unittest

from isca import GreyCodeBase, IscaCodeBase, DryCodeBase, ColumnCodeBase, ShallowCodeBase, BarotropicCodeBase, GFDL_BASE

def test_grey_codebase():
    cb = GreyCodeBase.from_directory(GFDL_BASE)
    cb.compile()

def test_isca_codebase():
    cb = IscaCodeBase.from_directory(GFDL_BASE)
    cb.compile()

def test_dry_codebase():
    cb = DryCodeBase.from_directory(GFDL_BASE)
    cb.compile()

def test_column_codebase():
    cb = ColumnCodeBase.from_directory(GFDL_BASE)
    cb.compile()

def test_shallow_codebase():
    cb = ShallowCodeBase.from_directory(GFDL_BASE)
    cb.compile()

def test_barotropic_codebase():
    cb = BarotropicCodeBase.from_directory(GFDL_BASE)
    cb.compile()

# Not tested here: SocratesCodeBase and SocColumnCodeBase, since both require the
# Socrates source code to be downloaded separately and placed in
# src/atmos_param/socrates/src/trunk - see exp/test_cases/socrates_test/README.md.

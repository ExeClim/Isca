import os

import pytest

from isca.socrates_paths import (
    build_module_index,
    build_procedure_index,
    build_stem_index,
    generate,
    _drop_textually_included_duplicates,
    _prune_leaky_files,
)


def _write(path, content):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        fh.write(content)


# A file that both calls into and textually INCLUDEs another file's full
# source (mirrors the real aux/seaalbedo_driver.f / aux/seaalbedo.f pair in
# Socrates itself -- see _drop_textually_included_duplicates).
FOO_F90 = """\
      SUBROUTINE seaalbedo_like(x)
      REAL x
      x = x + 1.0
      END SUBROUTINE seaalbedo_like
"""

FOO_DRIVER_F90 = """\
      SUBROUTINE foo_driver(x)
      REAL x
      CALL seaalbedo_like(x)
      END SUBROUTINE foo_driver

      INCLUDE 'foo.f90'
"""

# A file whose only path to an out-of-scope dependency is via an INTERFACE
# block restating that dependency's call signature (mirrors the real
# general/solar_intensity_90.f90, whose INTERFACE block for planck_90
# previously made build_procedure_index() treat it as a second definition
# of planck_90 -- see test_build_procedure_index_ignores_interface_block_declarations).
BAR_F90 = """\
      SUBROUTINE bar_impl(x)
      INTERFACE
        SUBROUTINE helper(y)
        REAL y
        END SUBROUTINE helper
      END INTERFACE
      REAL x
      CALL helper(x)
      END SUBROUTINE bar_impl
"""

# The real implementation of `helper`, living outside ALLOWED_SRC_SUBDIRS.
HELPER_F90 = """\
      SUBROUTINE helper(y)
      REAL y
      y = y * 2.0
      END SUBROUTINE helper
"""


@pytest.fixture
def socrates_tree(tmp_path):
    """A miniature Socrates-like source tree (make/Mk_src_aux + src/<dir>/*.f90),
    laid out like a real checkout, used to exercise socrates_paths.py without
    needing one."""
    root = str(tmp_path)
    _write(os.path.join(root, 'make', 'Mk_src_aux'),
           'SRC_AUX = foo.f90 foo_driver.f90 bar.f90\n')
    _write(os.path.join(root, 'src', 'aux', 'foo.f90'), FOO_F90)
    _write(os.path.join(root, 'src', 'aux', 'foo_driver.f90'), FOO_DRIVER_F90)
    _write(os.path.join(root, 'src', 'aux', 'bar.f90'), BAR_F90)
    _write(os.path.join(root, 'src', 'excluded', 'helper.f90'), HELPER_F90)
    return root


def test_build_procedure_index_ignores_interface_block_declarations(socrates_tree):
    src_root = os.path.join(socrates_tree, 'src')
    index = build_procedure_index(src_root)
    # bar.f90's INTERFACE block merely restates helper's call signature --
    # it must not be counted as a second definition of `helper` alongside
    # the real one in excluded/helper.f90. If it were, `helper` would look
    # ambiguous (defined in two places) and get dropped from the index
    # entirely, silently breaking anything that needs to resolve a CALL to
    # it (see build_procedure_index's docstring).
    assert index.get('helper') == os.path.join('excluded', 'helper.f90')
    assert index.get('bar_impl') == os.path.join('aux', 'bar.f90')


def test_build_procedure_index_still_flags_genuine_duplicates(tmp_path):
    # Sanity check that the interface-block fix didn't also paper over a
    # *real* collision: two independent definitions of the same bare
    # procedure name must still be treated as ambiguous and left out.
    root = str(tmp_path)
    _write(os.path.join(root, 'src', 'a', 'one.f90'), 'SUBROUTINE dup(x)\n      END SUBROUTINE dup\n')
    _write(os.path.join(root, 'src', 'b', 'two.f90'), 'SUBROUTINE dup(x)\n      END SUBROUTINE dup\n')
    index = build_procedure_index(os.path.join(root, 'src'))
    assert 'dup' not in index


def test_prune_leaky_files_drops_file_that_calls_out_of_scope_code(socrates_tree):
    src_root = os.path.join(socrates_tree, 'src')
    stem_index = build_stem_index(src_root)
    module_index = build_module_index(src_root)
    procedure_index = build_procedure_index(src_root)
    files = {
        os.path.join('aux', 'foo.f90'),
        os.path.join('aux', 'foo_driver.f90'),
        os.path.join('aux', 'bar.f90'),
    }

    clean, dropped = _prune_leaky_files(socrates_tree, files, stem_index, module_index, procedure_index)

    # bar.f90 calls helper(), which build_procedure_index resolves to
    # excluded/helper.f90 -- outside ALLOWED_SRC_SUBDIRS -- so bar.f90
    # itself must be dropped rather than left in a build that will fail to
    # link against the excluded implementation.
    assert dropped == {os.path.join('aux', 'bar.f90')}
    assert clean == {os.path.join('aux', 'foo.f90'), os.path.join('aux', 'foo_driver.f90')}


def test_drop_textually_included_duplicates(socrates_tree):
    stem_index = build_stem_index(os.path.join(socrates_tree, 'src'))
    files = {os.path.join('aux', 'foo.f90'), os.path.join('aux', 'foo_driver.f90')}

    clean, dropped = _drop_textually_included_duplicates(socrates_tree, files, stem_index)

    # foo_driver.f90 both calls into and textually INCLUDEs the whole of
    # foo.f90; compiling foo.f90 a second time as its own translation unit
    # would duplicate every symbol it defines ("multiple definition of
    # seaalbedo_like_"), so the independently-listed copy must be dropped.
    assert dropped == {os.path.join('aux', 'foo.f90'): os.path.join('aux', 'foo_driver.f90')}
    assert clean == {os.path.join('aux', 'foo_driver.f90')}


def test_drop_textually_included_duplicates_keeps_module_files(tmp_path):
    # A file that defines its own MODULE must never be dropped just
    # because something else happens to INCLUDE it too -- some other file
    # may need to USE it directly, which dropping it would silently break.
    root = str(tmp_path)
    _write(os.path.join(root, 'src', 'aux', 'mod_file.f90'),
           'MODULE real_mod\n      INTEGER :: n\n      END MODULE real_mod\n')
    _write(os.path.join(root, 'src', 'aux', 'includer.f90'),
           "SUBROUTINE includer_sub()\n      END SUBROUTINE includer_sub\n\n      INCLUDE 'mod_file.f90'\n")
    stem_index = build_stem_index(os.path.join(root, 'src'))
    files = {os.path.join('aux', 'mod_file.f90'), os.path.join('aux', 'includer.f90')}

    clean, dropped = _drop_textually_included_duplicates(root, files, stem_index)

    assert dropped == {}
    assert clean == files


def test_generate_end_to_end_drops_both_problem_files(socrates_tree):
    # Combines both fixed bugs in one pass, mirroring the real Socrates
    # checkouts (1703/2207/2411/2026.07.1) this was validated against by
    # actual compilation: an out-of-scope dependency reachable only through
    # an INTERFACE block (bar.f90, standing in for
    # general/solar_intensity_90.f90) and a textual-INCLUDE duplicate
    # (foo.f90/foo_driver.f90, standing in for aux/seaalbedo.f /
    # aux/seaalbedo_driver.f).
    result = generate(socrates_tree)
    assert result == [os.path.join('aux', 'foo_driver.f90')]

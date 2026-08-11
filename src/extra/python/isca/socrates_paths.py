"""Automated generation of per-version Socrates path_names entries for Isca.

Isca's Socrates interface only needs a subset of the upstream Socrates
source tree, and the Met Office restructures directories between releases
often enough that a hand-maintained file list breaks on every new version
(files get moved between directories, occasionally renamed). This module
derives that subset for a given checked-out Socrates source tree by
combining two independent methods and taking the union of what they find:

1. Manifest parsing: the Met Office's own make/Mk_src_*, make/Mk_mod_*
   Makefile fragments list, by bare filename, which files belong to each
   library component (aux, general, core modules, radiance core, nlte).
   We parse those directly rather than re-deriving the same information
   by hand for every release.

2. Dependency closure: starting from Isca's own socrates interface files,
   follow every Fortran ``USE <module>`` and ``INCLUDE 'x.finc'``
   statement recursively through the Socrates tree, resolving module
   names and filenames against an index built by walking the tree. This
   catches files the manifests don't mention because Isca's interface
   uses them even though the standard Socrates build doesn't (e.g.
   qsat_gill.F90), and pulls in include files that are never separate
   compile targets so never appear in a Makefile's source list.

Filenames are resolved to their real location by walking the source tree
at generation time rather than assuming fixed directories, since
directories move between releases too.
"""
import os
import re
import sys

# Met Office Makefile fragments (under <socrates_root>/make/) that between
# them define the "core science" library Isca's Socrates interface calls
# into: general-purpose numerics (gen), core module definitions
# (gencore/mod_gen), the radiative transfer solver itself (radcore), and
# the small NLTE heating utility Isca's interface pulls in. Deliberately
# excludes the correlated-k and scattering *database-generation* tools,
# COSP, flexchem, and the UM/LFRic driver code -- Isca's interface never
# calls into any of those.
CORE_FRAGMENTS = [
    'Mk_src_aux',
    'Mk_src_gen',
    'Mk_mod_gen',
    'Mk_src_gencore',
    'Mk_mod_gencore',
    'Mk_src_radcore',
    'Mk_src_nlte',
]

# Top-level src/ subdirectories the dependency closure is allowed to pull
# files from. Mk_src_aux ropes in a couple of files (Met Office's own
# `socrates_runes`/`socrates_bones` standalone driver-program support code)
# that reach out into correlated-k tools, the COSP satellite simulator,
# illumination/orbit calculations and Met Office's newer "interface_core"
# high-level driver layer -- none of which Isca needs, since Isca's own
# hand-written interface (src/atmos_param/socrates/interface) already does
# the job interface_core does, and Isca supplies its own solar
# geometry via FMS' astronomy_mod rather than Socrates' illumination
# modules. Pulling in an arbitrary subset of e.g. COSP is also a real
# compile-failure risk, since a satellite simulator has its own large
# dependency tree we would only be including a fragment of. Anything only
# reachable through one of these directories is dropped, and files whose
# own references lead there are dropped too (see _prune_leaky_files).
ALLOWED_SRC_SUBDIRS = frozenset({
    'aux', 'general', 'radiance_core', 'modules_core', 'modules_gen', 'nlte',
})

SOURCE_EXTENSIONS = ('.f90', '.finc', '.f')

_MAKEFILE_ASSIGN_RE = re.compile(r'^[A-Za-z0-9_]+\s*[:+]?=')
_MODULE_DECL_RE = re.compile(r'^\s*module\s+(\w+)\s*$', re.IGNORECASE)
_MODULE_SKIP_RE = re.compile(r'^\s*(end\s*module|module\s+procedure)\b', re.IGNORECASE)
_USE_RE = re.compile(r'^\s*use\s*(?:,\s*\w+\s*::)?\s*[:]{0,2}\s*(\w+)', re.IGNORECASE)
_INCLUDE_RE = re.compile(r'''^\s*#?\s*include\s+['"]([^'"]+)['"]''', re.IGNORECASE)
_CALL_RE = re.compile(r'^\s*call\s+(\w+)', re.IGNORECASE)
_PROC_DECL_RE = re.compile(r'\b(?:subroutine|function)\s+(\w+)', re.IGNORECASE)
_PROC_SKIP_RE = re.compile(r'^\s*(end\s|end$|!)', re.IGNORECASE)
_INTERFACE_START_RE = re.compile(r'^\s*(?:abstract\s+)?interface\b', re.IGNORECASE)
_INTERFACE_END_RE = re.compile(r'^\s*end\s*interface\b', re.IGNORECASE)


def _stem(filename):
    """Extension-stripped, lowercased filename, used as the key files are
    matched on. Met Office Makefiles occasionally reference a file with a
    different extension (or none) than the real one on disk, e.g. a
    build-target name like ``foo.f`` for a source file that is really
    ``foo.F90`` after preprocessing -- matching on stem alone is robust to
    that."""
    base = filename.lower()
    for ext in SOURCE_EXTENSIONS:
        if base.endswith(ext):
            return base[:-len(ext)]
    return base


def build_stem_index(src_root):
    """Map extension-stripped, lowercased filename -> list of paths
    (relative to src_root) where a file with that stem is found."""
    index = {}
    for root, _dirs, files in os.walk(src_root):
        for fname in files:
            if not fname.lower().endswith(SOURCE_EXTENSIONS):
                continue
            rel = os.path.relpath(os.path.join(root, fname), src_root)
            index.setdefault(_stem(fname), []).append(rel)
    return index


def build_module_index(src_root):
    """Map lowercased Fortran module name -> path (relative to src_root)
    of the file that declares it."""
    modules = {}
    for root, _dirs, files in os.walk(src_root):
        for fname in files:
            if not fname.lower().endswith(SOURCE_EXTENSIONS):
                continue
            path = os.path.join(root, fname)
            rel = os.path.relpath(path, src_root)
            try:
                with open(path, errors='ignore') as fh:
                    for line in fh:
                        if _MODULE_SKIP_RE.match(line):
                            continue
                        m = _MODULE_DECL_RE.match(line)
                        if m:
                            modules.setdefault(m.group(1).lower(), rel)
            except OSError:
                continue
    return modules


def build_procedure_index(src_root):
    """Map lowercased subroutine/function name -> path (relative to
    src_root) of the file that declares it.

    Used to follow bare ``CALL foo(...)`` statements to external
    (non-module) procedures, which have no ``USE`` statement linking
    caller to callee -- the only static trace of the dependency is the
    call site itself.

    Unlike module names, bare procedure names are not required to be
    unique across the whole Socrates tree (e.g. both the core radiative
    solver and the unrelated, never-linked-with-it COSP simulator define
    a routine called ``two_stream``). A name that resolves to more than
    one file is ambiguous and is deliberately left out rather than
    guessing, since guessing wrong can silently substitute the wrong
    dependency for a real one.

    Lines inside ``INTERFACE``/``END INTERFACE`` blocks are skipped: an
    explicit interface only restates a procedure's call signature (name,
    argument types) for the compiler's benefit, it isn't a second
    definition. Counting it as one would make the name look ambiguous
    (dropping a genuine external dependency out of the index, silently
    defeating both dependency-closure following and _prune_leaky_files
    for that name) or, if the interface is the only "definition" found,
    would wrongly point callers at a file with no actual subroutine body.
    """
    candidates = {}
    for root, _dirs, files in os.walk(src_root):
        for fname in files:
            if not fname.lower().endswith(SOURCE_EXTENSIONS):
                continue
            path = os.path.join(root, fname)
            rel = os.path.relpath(path, src_root)
            try:
                with open(path, errors='ignore') as fh:
                    in_interface = False
                    for line in fh:
                        if _INTERFACE_END_RE.match(line):
                            in_interface = False
                            continue
                        if _INTERFACE_START_RE.match(line):
                            in_interface = True
                            continue
                        if in_interface:
                            continue
                        if _PROC_SKIP_RE.match(line):
                            continue
                        m = _PROC_DECL_RE.search(line)
                        if m:
                            candidates.setdefault(m.group(1).lower(), set()).add(rel)
            except OSError:
                continue
    return {name: next(iter(rels)) for name, rels in candidates.items() if len(rels) == 1}


def parse_makefile_fragment(path):
    """Return the set of source filenames (as written, lowercased)
    referenced anywhere in a Mk_src_*/Mk_mod_* Makefile fragment,
    regardless of which variable(s) it defines."""
    with open(path, errors='ignore') as fh:
        text = fh.read()
    text = text.replace('\\\n', ' ')
    found = set()
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if _MAKEFILE_ASSIGN_RE.match(line):
            line = line.split('=', 1)[1]
        for tok in line.split():
            if tok.lower().endswith(SOURCE_EXTENSIONS):
                found.add(tok.lower())
    return found


def _top_dir(relpath):
    return relpath.split(os.sep, 1)[0]


def _resolve(name, stem_index):
    matches = stem_index.get(_stem(name), [])
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        exact = [m for m in matches if os.path.basename(m).lower() == name.lower()]
        if len(exact) == 1:
            return exact[0]
    return None


def manifest_files(socrates_root, stem_index):
    """Filenames referenced by the Met Office's own core-library Makefile
    fragments, resolved against the tree.

    Returns (resolved relpaths set, unresolved filenames set, out-of-scope
    filename->relpath dict for files that resolved fine but live outside
    ALLOWED_SRC_SUBDIRS).
    """
    make_dir = os.path.join(socrates_root, 'make')
    wanted = set()
    for frag in CORE_FRAGMENTS:
        frag_path = os.path.join(make_dir, frag)
        if os.path.exists(frag_path):
            wanted |= parse_makefile_fragment(frag_path)

    resolved = set()
    unresolved = set()
    out_of_scope = {}
    for fname in wanted:
        target = _resolve(fname, stem_index)
        if not target:
            unresolved.add(fname)
        elif _top_dir(target) in ALLOWED_SRC_SUBDIRS:
            resolved.add(target)
        else:
            out_of_scope[fname] = target
    return resolved, unresolved, out_of_scope


def dependency_closure(socrates_root, stem_index, module_index, procedure_index,
                        seed_relpaths, extra_seed_files=()):
    """Follow USE/INCLUDE/CALL references from a set of starting files
    until no new Socrates files are discovered.

    seed_relpaths: paths already known to be needed, relative to
        <socrates_root>/src.
    extra_seed_files: absolute paths to additional files to scan for
        references *into* the Socrates tree without being part of the
        output themselves (Isca's own interface files).
    """
    found = set(seed_relpaths)
    worklist = list(seed_relpaths)

    def maybe_add(target):
        if target and target not in found and _top_dir(target) in ALLOWED_SRC_SUBDIRS:
            found.add(target)
            worklist.append(target)

    def scan(full_path):
        try:
            with open(full_path, errors='ignore') as fh:
                lines = fh.readlines()
        except OSError:
            return
        for line in lines:
            m = _USE_RE.match(line)
            if m:
                maybe_add(module_index.get(m.group(1).lower()))
                continue
            m = _INCLUDE_RE.match(line)
            if m:
                maybe_add(_resolve(m.group(1), stem_index))
                continue
            m = _CALL_RE.match(line)
            if m:
                maybe_add(procedure_index.get(m.group(1).lower()))

    for extra in extra_seed_files:
        scan(extra)

    seen = set()
    while worklist:
        rel = worklist.pop()
        if rel in seen:
            continue
        seen.add(rel)
        scan(os.path.join(socrates_root, 'src', rel))

    return found


def _prune_leaky_files(socrates_root, files, stem_index, module_index, procedure_index):
    """Drop files that are nominally inside an allowed subdirectory (so
    manifest parsing or closure kept them) but themselves reference a file
    outside ALLOWED_SRC_SUBDIRS -- e.g. a driver-program file that happens
    to live under aux/ but pulls in the COSP simulator. Since that
    out-of-scope target is deliberately excluded, keeping the file that
    needs it would just produce a "module not found" compile error, so we
    drop it too and report it instead.

    Returns (clean_files, dropped_files).
    """
    def references_excluded_dir(relpath):
        full = os.path.join(socrates_root, 'src', relpath)
        try:
            with open(full, errors='ignore') as fh:
                lines = fh.readlines()
        except OSError:
            return False
        for line in lines:
            m = _USE_RE.match(line)
            if m:
                target = module_index.get(m.group(1).lower())
            else:
                m = _INCLUDE_RE.match(line)
                if m:
                    target = _resolve(m.group(1), stem_index)
                else:
                    m = _CALL_RE.match(line)
                    target = procedure_index.get(m.group(1).lower()) if m else None
            if target and _top_dir(target) not in ALLOWED_SRC_SUBDIRS:
                return True
        return False

    dropped = {f for f in files if references_excluded_dir(f)}
    return files - dropped, dropped


# Extensions mkmf gives their own independent compile (and hence link) rule
# to. .finc files never do, so a .finc file being both listed on its own and
# textually INCLUDEd elsewhere is harmless -- only these extensions can
# produce the "multiple definition" collision _drop_textually_included_duplicates
# guards against.
_INCLUDED_COMPILABLE_EXTENSIONS = ('.f90', '.f')


def _drop_textually_included_duplicates(socrates_root, files, stem_index, warn=None):
    """Drop files that are both independently listed as their own compile
    target *and* textually INCLUDEd (full source pasted in at compile time)
    by another file also being compiled.

    Socrates' own Makefile fragments can list a file as a standalone target
    even though some other file in the same fragment both calls into it and
    INCLUDEs its entire source -- observed for aux/seaalbedo.f, which
    aux/seaalbedo_driver.f both calls (see the CALL to ``seaalbedo`` inside
    it) and, in its final line, pulls in wholesale via
    ``INCLUDE 'seaalbedo.f'``. In Socrates' own multi-target build these
    two apparently never end up compiled into the same executable, so nothing
    collides; Isca builds everything into one flat library, so compiling
    aux/seaalbedo.f a second time as its own translation unit would
    duplicate every symbol it defines (``morel``, ``mousse``, ``seaalbedo``,
    ...), which the linker rejects with "multiple definition of ...". Drop
    the independently-listed duplicate and rely on the includer's copy
    instead. A file that defines its own MODULE is never dropped this way,
    since some other file may need to USE it directly.
    """
    if warn is None:
        warn = lambda msg: None

    def defines_module(rel):
        full = os.path.join(socrates_root, 'src', rel)
        try:
            with open(full, errors='ignore') as fh:
                return any(_MODULE_DECL_RE.match(line) for line in fh)
        except OSError:
            return False

    dropped = {}
    for rel in sorted(files):
        full = os.path.join(socrates_root, 'src', rel)
        try:
            with open(full, errors='ignore') as fh:
                lines = fh.readlines()
        except OSError:
            continue
        for line in lines:
            m = _INCLUDE_RE.match(line)
            if not m:
                continue
            target = _resolve(m.group(1), stem_index)
            if (target and target in files and target != rel
                    and target.lower().endswith(_INCLUDED_COMPILABLE_EXTENSIONS)
                    and target not in dropped and not defines_module(target)):
                dropped[target] = rel

    if dropped:
        warn(
            "generate_socrates_path_names: dropped %d file(s) that are "
            "independently listed as their own compile target but are "
            "also textually INCLUDEd (full source pasted in) by another "
            "file also being compiled -- compiling both would duplicate "
            "every symbol the included file defines: %s" %
            (len(dropped), ', '.join('%s (included by %s)' % kv for kv in sorted(dropped.items())))
        )
    return files - set(dropped), dropped


def generate(socrates_root, isca_interface_dir=None, warn=None):
    """Return the sorted list of Socrates source paths (relative to
    <socrates_root>/src) that Isca's Socrates interface needs, for the
    checked-out Socrates source tree at socrates_root.
    """
    if warn is None:
        warn = lambda msg: print(msg, file=sys.stderr)

    src_root = os.path.join(socrates_root, 'src')
    stem_index = build_stem_index(src_root)
    module_index = build_module_index(src_root)
    procedure_index = build_procedure_index(src_root)

    resolved, unresolved, out_of_scope = manifest_files(socrates_root, stem_index)
    if unresolved:
        warn(
            "generate_socrates_path_names: could not resolve %d filename(s) "
            "referenced by Socrates' own Makefile fragments (they may be "
            "offline-tool-only or upstream may have renamed them): %s" %
            (len(unresolved), ', '.join(sorted(unresolved)))
        )
    if out_of_scope:
        warn(
            "generate_socrates_path_names: %d filename(s) from Socrates' "
            "own Makefile fragments resolved to files outside "
            "aux/general/radiance_core/modules_core/modules_gen/nlte, so "
            "were excluded: %s" %
            (len(out_of_scope), ', '.join('%s -> %s' % kv for kv in sorted(out_of_scope.items())))
        )

    extra_seeds = []
    if isca_interface_dir and os.path.isdir(isca_interface_dir):
        for fname in sorted(os.listdir(isca_interface_dir)):
            if fname.lower().endswith(SOURCE_EXTENSIONS):
                extra_seeds.append(os.path.join(isca_interface_dir, fname))

    closure = dependency_closure(socrates_root, stem_index, module_index, procedure_index,
                                  resolved, extra_seeds)

    clean, dropped = _prune_leaky_files(socrates_root, closure, stem_index, module_index, procedure_index)
    if dropped:
        warn(
            "generate_socrates_path_names: excluded %d file(s) that are "
            "nominally in-scope but themselves depend on code outside "
            "aux/general/radiance_core/modules_core/modules_gen/nlte "
            "(correlated-k, scatter, cosp, flexchem, illumination, "
            "interface_core, odepack), which Isca's interface doesn't use: "
            "%s" % (len(dropped), ', '.join(sorted(dropped)))
        )

    clean, _include_dupes = _drop_textually_included_duplicates(socrates_root, clean, stem_index, warn=warn)

    return sorted(clean)


def format_path_names(paths, version):
    return ['atmos_param/socrates/src/%s/src/%s' % (version, p) for p in paths]


def write_path_names(paths, version, output_file):
    lines = format_path_names(paths, version)
    with open(output_file, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    return len(lines)

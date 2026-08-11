#!/usr/bin/env python
"""Generate an Isca socrates_version_paths/<version> file from a checked-out
Socrates source tree.

Usage:

    python generate_socrates_path_names.py /path/to/socrates/checkout 2026.07.1

By default this writes to
src/extra/model/socrates/socrates_version_paths/<version> in this Isca
checkout, and seeds the dependency-closure scan with
src/atmos_param/socrates/interface. See isca.socrates_paths for how the
file list is actually derived.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from isca.socrates_paths import generate, write_path_names  # noqa: E402

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('socrates_root', help='Path to a checked-out Socrates source tree (must contain src/ and make/)')
    parser.add_argument('version', help="Version label used in the generated path_names entries, e.g. '2026.07.1'")
    parser.add_argument('--interface-dir', default=None,
                         help='Path to socrates/interface, used to seed the dependency closure scan. '
                              'Defaults to src/atmos_param/socrates/interface in this Isca checkout.')
    parser.add_argument('--output', default=None,
                         help='Output file. Defaults to '
                              'src/extra/model/socrates/socrates_version_paths/<version> in this Isca checkout.')
    args = parser.parse_args()

    interface_dir = args.interface_dir or os.path.join(
        REPO_ROOT, 'src', 'atmos_param', 'socrates', 'interface')
    output = args.output or os.path.join(
        REPO_ROOT, 'src', 'extra', 'model', 'socrates', 'socrates_version_paths', args.version)

    paths = generate(args.socrates_root, interface_dir)
    os.makedirs(os.path.dirname(output), exist_ok=True)
    n = write_path_names(paths, args.version, output)
    print('Wrote %d entries to %s' % (n, output))


if __name__ == '__main__':
    main()

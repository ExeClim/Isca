This is the location for the Socrates source code, one subdirectory per version, e.g. `2026.07.1/`.

You don't need to populate this by hand: `SocratesCodeBase`/`SocColumnCodeBase` fetch the requested `socrates_version` automatically the first time you compile (or use it as-is if you've already run `git submodule update --init` for the default vendored version -- see `.gitmodules` at the repository root). See `docs/source/modules/socrates.rst` for details.

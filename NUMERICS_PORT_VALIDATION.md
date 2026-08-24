# Numerics/build port validation: `uob_fftw_sit23_numerics_2026`

This documents validation of the FFTW3 spectral-transform core and HPC
build/environment support ported from `uob_fftw_sit23_new_maths2intel` onto
current `master`, done on branch `uob_fftw_sit23_numerics_2026`. Companion to
the physics port, validated separately in `PHYSICS_PORT_VALIDATION.md` on
`uob_fftw_sit23_physics_2026`. Two checks are covered:

1. Does the standard (non-FFTW) build still produce bit-identical results to
   `master`, given this branch touches several shared-physics files outside
   any `#ifdef FFTW3` guard?
2. Does the new FFTW3 spectral core actually work — does it compile, link,
   run, and produce physically correct output?

## What was ported

Branch `uob_fftw_sit23_numerics_2026` (commit `85bd52f1`) carries 9 commits
from `uob_fftw_sit23_new_maths2intel` onto current `master`, with original
authorship preserved (cherry-pick / diff+3-way-apply as appropriate):

- **FFTW3 spectral transform core** (`grid_fourier.F90`,
  `spherical_fourier.F90`, `spherical_fourier_fftw.F90` (new), `fft.F90`,
  `fftw.F90` (new), plus `spectral_dynamics.F90`'s dispatch to it) —
  George Lancaster, 2019.
- **`cb.enable_fftw3()` codebase API**, `DryCodeBaseFFTW` export, T213
  resolution, slurm multi-node support — sit23, 2020-2021.
- **`dry_fftw`/`grey_fftw` model directories** and the `held_suarez_fftw`
  test case — qv18258, 2019.
- **Isambard/BluePebble/Cray-GNU HPC support**: env files, mkmf templates,
  job scripts (kept per Stephen's explicit instruction — "just in case
  they're useful for others" — even though these specific HPC systems may
  be retired) and a set of compiler-compatibility fixes to shared physics
  code needed to build on that toolchain — George Lancaster / Lancasterg,
  2019.
- **Postprocessing script improvements** — `run_plevel.py`'s `all_vars`
  variable-subsetting option (reconciled against master's own since-modernised
  scratch config), a `plevel_fn.py` fix for variables with a `scalar_axis`
  dimension, and small fixes to `mppnccombine_run.sh` /
  `compile_plev_interpolation.sh` / `modified_time_script.py` — sit23,
  2020-2021.

**Deliberately excluded**, per Stephen's explicit decision:
- The branch's own incomplete single-precision attempt — deferred to the
  separate, more complete `single_prec` branch, to be merged later as its
  own project.
- `maths2` (env file + mkmf template) — turns out to **already be on
  `master`**, added at commit `9e78f165` (this branch's own merge-base), so
  it predates the branch entirely and was never part of the port. This
  resolves the open question of whether `maths2` is the same lineage as the
  currently-used `maths5`/`maths-gpu` nodes: it's moot, since master already
  carries whatever `maths2` support it has independently of this branch.

## The concern: two files change outside any `#ifdef FFTW3` guard

Unlike the physics port (fully additive, all gated behind new namelist
options), this branch modifies shared code that runs on *every* build,
FFTW3 or not:

- **`transforms.F90`**: `trans_spherical_to_grid_2d`/`_3d` and
  `reverse_transpose_fourier` gain an optional `block` argument (defaults
  preserve old behaviour); the `#ifdef INTERNAL_FILE_NML` namelist-reading
  branch was rewritten to open `input.nml` directly via
  `open_namelist_file()`/`close_file()` instead of reading from the
  preloaded `input_nml_file` internal-file array.
- **`spectral_dynamics.F90`**: the same `INTERNAL_FILE_NML` namelist-reading
  rewrite.
- The same `open_namelist_file()`/`close_file()` rewrite, applied
  consistently, also appears in `diffusivity.F90`, `hs_forcing.F90`,
  `lscale_cond.F90`, `monin_obukhov_kernel.F90`, `vert_turb_driver.F90`,
  `vert_advection.F90`, `atmos_model.F90`, `atmosphere.F90`,
  `polvani_2007.F90`, `spectral_initialize_fields.F90`,
  `vert_coordinate.F90`, and `surface_flux.F90` (part of the Cray/Isambard
  compiler-compatibility commit) — evidently a deliberate, systematic fix
  for that toolchain, not an accidental one-off.

None of these files add or rename any namelist variables (checked directly:
no diff in either module's `namelist /..._nml/` variable list), so unlike
the physics port's `update_land_mask_from_ice` case, `trip_test`'s
shared-namelist limitation ([[infra-trip-test-shared-namelist-limitation]])
doesn't apply here — a direct `trip_test` comparison is valid.

## Bit-reproducibility of the standard build

Ran `trip_test_command_line` comparing `master` (`d1321ac3`) against this
branch's tip (`85bd52f1`) on `held_suarez` and `frierson` — chosen because
between them they exercise the dynamical core (`transforms.F90`,
`spectral_dynamics.F90`), the boundary-layer/turbulence physics
(`monin_obukhov_kernel.F90`, `vert_turb_driver.F90`, `diffusivity.F90`),
large-scale condensation (`lscale_cond.F90`), and the driver/coupler layer
(`atmos_model.F90`, `atmosphere.F90`, `surface_flux.F90`) — i.e. essentially
every non-FFTW-gated file this branch touches.

**Result: pass on both.**

```
Test passed for held_suarez. Commit 85bd52f1 gives the same answer as commit d1321ac3
Test passed for frierson. Commit 85bd52f1 gives the same answer as commit d1321ac3
held_suarez : pass
frierson : pass
Congratulations, all tests have passed
```

Bit-identical output confirms the `open_namelist_file()` rewrite and the new
optional `block` argument are behaviourally neutral for the standard build —
consistent with what they look like on inspection (reading the same
`input.nml` content by a different mechanism; an optional argument that
defaults to the old behaviour everywhere it isn't explicitly passed).

## Does the FFTW3 spectral core actually work?

Compiling `held_suarez_fftw` with `-DFFTW3` was already confirmed to compile
cleanly during the port itself (both without and with FFTW3 enabled, the
latter needing the `dry_fftw`/`grey_fftw` `path_names` fix documented in the
commit history — stale manifests missing `cloud_simple`/
`frierson_monin_obukhov.F90`). That only proves it compiles, not that it
runs correctly, so a further check was run here: **an actual 5-day
`held_suarez` integration**, comparing the FFTW3 spectral core against the
default (Temperton) FFT, same resolution (T42L25), same namelist, same
initial conditions.

(This machine's `isca_env_26` conda environment doesn't ship FFTW3 dev
headers/libraries, unlike the target Isambard/BluePebble/Cray systems where
this would come from an environment module. A local-only, uncommitted
`-I`/`-L`/`-Wl,-rpath`/`-lfftw3` addition to `mkmf.template.ubuntu_conda`,
pointing at an existing FFTW3 install in the `isca_analysis` conda
environment, was used to make this testable here; reverted immediately
after, confirmed via `git status`/`git diff` afterward.)

**Result: agrees with the reference FFT to floating-point roundoff, not
"close" — essentially exact.**

| Field | max&#124;diff&#124; (FFTW3 &minus; default FFT), day 5 |
|---|---|
| ps | 0.0 |
| temp | 0.0 |
| ucomp | 5.96&times;10<sup>-8</sup> m/s |
| vcomp | 2.98&times;10<sup>-8</sup> m/s |
| vor | 1.14&times;10<sup>-13</sup> s<sup>-1</sup> |

These differences are all at or below double-precision rounding noise
(compare to typical field magnitudes: `ucomp`/`vcomp` range &plusmn;2&ndash;4
m/s, `vor` range &plusmn;3&times;10<sup>-6</sup> s<sup>-1</sup>) — several
orders of magnitude below the day-1 chaotic-noise-floor differences seen
between two independently-compiled binaries in the physics port's validation
(order 10<sup>-2</sup> in physical units). Two different FFT algorithms
computing the same spherical harmonic transform, to this level of agreement
after 5 days of nonlinear integration, is strong evidence the FFTW3 dispatch
in `fft.F90`/`fftw.F90` is correctly wired up — not just "doesn't crash",
but numerically equivalent to the reference implementation it replaces.

## Bottom line

- Every file this branch touches outside an `#ifdef FFTW3` guard is
  confirmed bit-reproducibility-neutral via `trip_test` (`held_suarez`,
  `frierson`, both pass against `master`).
- The FFTW3 spectral core, once actually run (not just compiled), agrees
  with the reference FFT implementation to floating-point roundoff over a
  5-day integration — the transform is correctly implemented, not just
  syntactically valid.
- Both the standard and FFTW3-enabled builds compile cleanly, including
  after fixing a real bug surfaced during the port itself (stale
  `dry_fftw`/`grey_fftw` `path_names` manifests missing modules current
  master's `dry`/`grey` models depend on).
- `maths2` and single-precision support are out of scope for this branch
  (see above) — nothing further needed here for either.

## Appendix: how this was tested

`trip_test_command_line -e held_suarez frierson -n 8 -r
<numerics-2026-worktree-path> d1321ac3 85bd52f1`, run locally against this
worktree as the `-r` repo argument (a plain local git clone source — no
GitHub round-trip needed, since all worktrees on this machine share one
`.git` object store).

The FFTW3-vs-default-FFT comparison used a standalone script (not
`trip_test`, since it needs to compile two different `CodeBase` subclasses —
`DryCodeBaseFFTW` vs `DryCodeBase` — against the *same* commit, which
`trip_test` isn't set up to do): both compiled from this branch's worktree,
run for 5 days at T42L25 with identical namelists/diagnostics/initial
conditions, output compared field-by-field with `xarray`.

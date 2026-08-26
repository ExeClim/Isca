# Single-precision port validation: `single_prec_2026`

This documents validation of the single-precision compilation support ported
from `single_prec` onto current `master`, done on branch `single_prec_2026`.
Companion to the physics and numerics/FFTW ports, validated separately on
`uob_fftw_sit23_physics_2026`/`uob_fftw_sit23_numerics_2026` — an unrelated
branch, not part of that merge. Two checks are covered:

1. Does the standard (double-precision) build still produce bit-identical
   results to `master`, given this branch touches shared code
   (`lcl.F90`, `cloud_simple/*.F90`, `codebase.py`, `compile.sh`) outside
   any precision-specific guard?
2. Does `cb.use_single_precision()` actually work — does it compile, link,
   run, and produce physically sane output?

## What was ported, and what wasn't finished on the original branch

Branch `single_prec_2026` (commit `faaa6d03`) carries the 5 commits from
`single_prec` onto current `master` (all sit23, April 2024), plus one new
commit completing what the original branch left unfinished:

- `cb.use_single_precision()`: a new `CodeBase` method that switches the
  precision-related compile macros from `-DOVERLOAD_C8` (double, the
  default) to `-DOVERLOAD_C4 -DOVERLOAD_R4` (single), renames the
  executable/build dir, and switches to a `<env>_single` environment file
  so a separate, precision-specific mkmf template can be selected (needed
  because e.g. Intel's `-r8`/gfortran's `-fdefault-real-8` flags must be
  present for double precision and *absent* for single — there's no way to
  toggle that from compile flags alone).
- `lcl.F90` (the simple-cloud scheme's lifting-condensation-level /
  Lambert-W solver): every `double precision` declaration changed to
  plain `real`, so the type responds to whichever precision the rest of
  the build is using, instead of being hardcoded to 8 bytes regardless.
- `cloud_simple.F90`/`cloud_cover_diags.F90`/`large_scale_cloud.F90`/
  `marine_strat_cloud.F90`: add `file_exist` to each file's `fms_mod`
  import list.
- `held_suarez_test_case_single.py`: a single-precision Held-Suarez test
  case, identical to the standard one plus `cb.use_single_precision()`.
- `mkmf.template.maths2_single`: a copy of `mkmf.template.maths2` with
  `-r8` removed.

**One bug fixed during the port, in the very first commit carried over**
(`b09d638f`): removing the hardcoded `-DOVERLOAD_C8` from `compile.sh`'s
`cppDefs` (correct — it's now supplied dynamically via
`precision_compile_flags`, which defaults to `['-DOVERLOAD_C8']`) also
dropped a letter from `-DINTERNAL_FILE_NML`, mangling it to
`-DINTERNAL_FILE_NM` — evidently a typo, not intentional, since nothing in
any of the 5 commits' messages mentions `INTERNAL_FILE_NML` at all. Left
as-is, this would have silently forced *every* build (single or double
precision) onto the non-`INTERNAL_FILE_NML` codepath everywhere in the
model. Fixed to read `-DINTERNAL_FILE_NML` again, while keeping the
`-DOVERLOAD_C8` removal and reconciling with master's own independent
addition of a `${CDEFS}` variable to this line since 2024.

**One thing genuinely never finished, completed in a new commit here**
(`faaa6d03`): `use_single_precision()` appends `_single` to whatever
`GFDL_ENV` file is loaded (so e.g. `maths2` &rarr; `maths2_single`), but no
`src/extra/env/maths2_single` file was ever created in the original branch
— its own last commit message ("Needs new env file to point to new mkmf
file") describes exactly this gap but the commit actually only adds the
mkmf *template*, not the env file that would select it. The feature was
therefore never actually runnable, even on the branch's own target machine.
Added the missing `src/extra/env/maths2_single` here, plus an equivalent
`ubuntu_conda_single` env file and `mkmf.template.ubuntu_conda_single` (this
machine doesn't use `maths2`/Intel; `ubuntu_conda`/gfortran needs the same
`-fdefault-real-8`/`-fdefault-double-8` removal as `maths2_single` removes
`-r8`) — needed to test any of this at all, and generically useful since
`ubuntu_conda` is the common development environment, not specific to this
machine.

**Also still true, independent of this branch:** the `file_exist` import
fix (`27342490`) fixes a real, still-present bug on current `master` —
all 4 `cloud_simple` files call `file_exist('input.nml')` in their
`#else` (non-`INTERNAL_FILE_NML`) branch without importing it. Currently
invisible because `master`'s `compile.sh` always hardcodes
`-DINTERNAL_FILE_NML`, so that branch is dead code under every normal
build — would only surface compiling without that flag, exactly the
situation this branch's own typo (above) would have created if left
unfixed.

## Bit-reproducibility of the standard (double-precision) build

Ran `trip_test_command_line` comparing `master` (`d1321ac3`) against this
branch's tip (`faaa6d03`) on two cases: `held_suarez` (general dynamical
core / build regression check) and `socrates_aquaplanet_cloud` (chosen
specifically because it's the one trip_test case that actually turns on
the simple-cloud scheme, exercising `cloud_simple.F90`/`lcl.F90` at
runtime, not just at compile time).

**Result: pass on both.**

```
Test passed for held_suarez. Commit faaa6d03... gives the same answer as commit d1321ac3
Test passed for socrates_aquaplanet_cloud. Commit faaa6d03... gives the same answer as commit d1321ac3
Congratulations, all tests have passed
```

Bit-identical output confirms the `lcl.F90` type change (`double
precision` &rarr; `real`) and the `file_exist` import additions are
behaviourally neutral for the default (double-precision) build — expected,
since the whole point of the `real` change is that it resolves to 8 bytes
under this build's `-fdefault-real-8`/`-fdefault-double-8` flags, i.e. the
same representation `double precision` already forced.

## Does `cb.use_single_precision()` actually work?

A standalone script ran `held_suarez` for 5 days at T42L25, once with
`cb.use_single_precision()` and once without (both from this branch's
worktree, identical namelist/diagnostics/initial conditions), and compared
output.

**Result: both runs complete successfully, no NaNs, no crash.** Differences
between single and double precision after 5 days:

| Field | max&#124;diff&#124; | relative |
|---|---|---|
| ps | 87.6 Pa | 0.09% |
| temp | 0.58 K | 0.2% |
| ucomp | 0.48 m/s | 13.1% |
| vcomp | 0.066 m/s | 2.9% |
| vor | 2.3&times;10<sup>-6</sup> s<sup>-1</sup> | 75%* |

(*`vor` is near a zero-crossing at day 5 in this configuration — small
absolute values make the relative figure large and not very meaningful on
its own; the absolute difference is the same order as the field itself,
consistent with the other fields' behaviour, not an outlier.)

This is a materially different comparison from the FFTW/numerics port's
validation (two *double-precision* compiles, agreeing to floating-point
roundoff, ~10<sup>-8</sup> relative) — here, one side is genuinely
single-precision (~7 significant decimal digits vs ~15), so a much larger
initial perturbation is expected from step one, then amplified by 5 days of
ordinary chaotic error growth in the same way as any other Isca comparison
in this project's validation work. The magnitudes above are consistent with
that: single precision is doing what it's supposed to (introducing
roundoff at the ~10<sup>-7</sup> relative level per operation), not
silently corrupting the run.

## Bottom line

- `lcl.F90`/`cloud_simple`'s changes are confirmed bit-reproducibility-
  neutral for the default double-precision build via `trip_test`
  (`held_suarez`, and `socrates_aquaplanet_cloud` which actually exercises
  the changed code at runtime).
- `cb.use_single_precision()` compiles, links, and runs a real experiment
  to completion, producing physically sane output whose divergence from
  the double-precision reference is exactly what single-precision roundoff
  amplified by 5 days of chaotic dynamics should look like.
- Two real gaps in the original branch were fixed as part of this port:
  the `compile.sh` `INTERNAL_FILE_NML` typo, and the missing
  `maths2_single`/`ubuntu_conda_single` env files that meant the feature
  was never actually runnable on any machine before now.

## Appendix: how this was tested

`trip_test_command_line -e held_suarez socrates_aquaplanet_cloud -n 8 -r
<single-prec-2026-worktree-path> d1321ac3 <branch-tip-hash>`, run locally
against this worktree as the `-r` repo argument.

The single-vs-double comparison used a standalone script (not `trip_test`,
since it needs to compile the *same* commit twice with different
`CodeBase` configuration, which `trip_test` isn't set up to do): both
compiled from this branch's worktree, run for 5 days at T42L25 with
identical namelists/diagnostics/initial conditions, output compared
field-by-field with `xarray`. Needed `GFDL_BASE` explicitly overridden to
this worktree's path for the duration of the test script — the installed
`isca` Python package resolves environment files via a module-level
`GFDL_BASE` constant read from `os.environ['GFDL_BASE']` at import time,
which on this machine normally points at the canonical `/home/links/sit204/Isca`
checkout, not whichever worktree `CodeBase.from_directory()` was pointed
at — without the override, `get_env_file()` looked for
`ubuntu_conda_single` in the wrong location, silently failed to source it
(`compile.sh` doesn't `set -e`), and fell through to the unrelated `ia64`
Intel-flavoured default template, which doesn't work with gfortran.

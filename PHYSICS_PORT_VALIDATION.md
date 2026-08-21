# Physics port validation: `uob_fftw_sit23_physics_2026`

This documents validation of the sea-ice/land-SST physics port from
`uob_fftw_sit23_new_maths2intel` onto current `master`, done on branch
`uob_fftw_sit23_physics_2026`. Two checks are covered:

1. Does an existing production experiment that depends on this physics
   (the Maunder Minimum q-flux run) still compile and run on the new branch?
2. Does it produce the same output as the original, un-ported branch?

## What was ported

Branch `uob_fftw_sit23_physics_2026` (commit `75044a0a`) carries 6 commits
cherry-picked from `uob_fftw_sit23_new_maths2intel` (commit `a9dc5390`) onto
current `master`, with original authorship preserved:

- 4 commits adding sea-ice/land-SST patching options to
  `src/atmos_spectral/driver/solo/mixed_layer.F90` (`specify_sst_over_sea_ice`,
  `linearly_interpolate_sea_ice_temp_and_sst`, `ice_sst_file`,
  `specify_sst_over_land_from_separate_file_to_ocean_sst`, `land_sst_file`,
  `allow_qflux_over_land`, `update_land_mask_from_ice`), hand-reconciled
  against master's independent ice-albedo work added since the branch forked.
- 1 compatibility commit, preserving existing bit-reproducibility for
  `realistic_continents`/`socrates_simcloud` (see below).
- 1 commit porting `calculate_qflux.py`'s Socrates support and
  `ignore_ice_for_calculation` option from the original branch.

One deliberate change during the port: the original branch's `binary_ice_albedo`
switch was dropped, because master independently added an equivalent (and more
general) `ice_albedo_method` option (`'step_function'` / `'ramp_function'`) in
the years since the branch forked. The two are algebraically identical for the
non-default case (`binary_ice_albedo=False` &equiv; `ice_albedo_method='ramp_function'`),
so this is a rename, not a behaviour change.

## Bit-reproducibility of existing test cases

Isca's `trip_test` suite couldn't be used directly to validate the
`update_land_mask_from_ice` compatibility fix, because it feeds one shared
namelist file to two different compiled binaries, and the fix specifically
requires a namelist key the old binary doesn't know about (see the
"how this was tested" note in the appendix). A standalone comparison giving
each commit its own correct namelist confirmed **bit-identical** output
(all fields, 0.0 max abs diff) against `master` for `realistic_continents_fixed_sst`
and `realistic_continents_variable_qflux` with `update_land_mask_from_ice`
set explicitly, as well as `frierson`, `held_suarez`, and `bucket_model`
passing `trip_test` unmodified.

## Does the Maunder Minimum experiment run?

**Yes.** `exp/socrates_fast_land_ice_cmpi5_ozone/maunder_min/socrates_fast_land_fft_new_ice_sbm_do_simple_false_bucket_qflux_mm.py`
compiles (`SocratesCodeBase`) and runs to completion on
`uob_fftw_sit23_physics_2026`, continuing from the experiment's existing
restart file (`input/res0240.tar.gz`, run 240 &rarr; 241). One namelist
adaptation was needed, matching the `binary_ice_albedo` retirement above:

```diff
- 'binary_ice_albedo': False,
+ 'ice_albedo_method': 'ramp_function',
```

`cb.enable_fftw3()` remains commented out in the script, as it is on the
original branch checkout used here — the FFTW/numerics port is a separate,
not-yet-started piece of work (`uob_fftw_sit23_numerics_2026`), so this
validation is physics-only and doesn't exercise FFTW either way.

## Does it produce the same output as the original branch?

A direct comparison was run: both `uob_fftw_sit23_new_maths2intel` (the
original, un-ported branch) and `uob_fftw_sit23_physics_2026` compiled and ran
**the same experiment, from the same restart file, with the same namelist**
(modulo the `binary_ice_albedo`&rarr;`ice_albedo_method` substitution above)
for one month (run 241, 30 days).

**Short answer: not bit-identical, and that's expected, not a bug.** Isca is a
chaotic dynamical system — even a rounding-level difference anywhere in the
computation (here, unavoidable given the two branches were compiled against
today's toolchain/library versions rather than whatever was current when the
original branch was last touched in 2023) amplifies over a month of
integration into large-looking differences in the instantaneous weather
state, with no error in the physics involved. Two pieces of evidence support
that this is exactly what's happening here, not a real regression:

**1. Fields that don't evolve chaotically match exactly.** `albedo`,
`ice_conc`, `flux_oceanq`, `bk`, `pk`, and `zsurf` — the fields most directly
exercised by the ported physics (ice-concentration interpolation, the q-flux
input file) or otherwise fixed/diagnostic — have **0.0 max absolute
difference** after 30 days. If the port had gotten any of this physics wrong,
these are exactly the fields that would show it, and they don't.

**2. The dynamical fields diverge with a textbook chaotic-error-growth
signature, not a step change.** A separate 5-day run with daily output shows
the difference between the two branches starting essentially at the noise
floor and growing roughly exponentially:

| Day | max&#124;&Delta;ps&#124; (Pa) | max&#124;&Delta;T&#124; (K) | max&#124;&Delta;u&#124; (m/s) | max&#124;&Delta;v&#124; (m/s) |
|---|---|---|---|---|
| 1 | 0.61  | 0.050 | 0.061 | 0.056 |
| 2 | 7.6   | 0.454 | 0.663 | 0.403 |
| 3 | 3.5   | 0.443 | 1.007 | 1.021 |
| 4 | 6.3   | 0.370 | 1.513 | 1.212 |
| 5 | 12.0  | 0.465 | 3.129 | 2.885 |

Day 1's differences (order 10<sup>-2</sup>&ndash;10<sup>-1</sup> in physical
units, i.e. parts-per-million relative to typical field magnitudes) are
consistent with ordinary floating-point rounding differences between two
different compiles, not a physics discrepancy. By day 30 this has amplified,
via normal atmospheric chaos, into the full-field differences summarised
below — large in places, but the same kind of divergence you'd see comparing
any two bit-non-identical runs of the same physical configuration.

### Full 30-day field comparison (`new_2026` &minus; `old_branch`)

| Field | max&#124;diff&#124; | relative | | Field | max&#124;diff&#124; | relative |
|---|---|---|---|---|---|---|
| ps | 916.9 Pa | 0.9% | | flux_lhe | 52.9 W/m&sup2; | 7.4% |
| temp | 5.03 K | 1.6% | | flux_t | 45.6 W/m&sup2; | 15.8% |
| ucomp | 14.26 m/s | 18.9% | | soc_surf_flux_lw | 13.80 W/m&sup2; | 8.0% |
| vcomp | 11.54 m/s | 24.7% | | soc_surf_flux_sw | 2.94 W/m&sup2; | 0.9% |
| t_surf | 1.74 K | 0.6% | | soc_olr | 4.88 W/m&sup2; | 1.5% |
| temp_2m | 1.85 K | 0.6% | | soc_toa_sw | 2.49 W/m&sup2; | 0.6% |
| rh | 17.2 % | 15.8% | | soc_ozone | 6.3&times;10<sup>-8</sup> | 0.4% |
| precipitation | 4.8&times;10<sup>-5</sup> | 19.7% | | vor | 2.1&times;10<sup>-5</sup> | 33.0% |
| **albedo** | **0.0** | **0.0%** | | div | 1.5&times;10<sup>-5</sup> | 49.8% |
| **ice_conc** | **0.0** | **0.0%** | | **bk / pk / zsurf** | **0.0** | **0.0%** |
| **flux_oceanq** | **0.0** | **0.0%** | | | | |

(Bolded rows are the ones most directly tied to the ported physics; all
match exactly.)

## Where the divergence actually comes from

The chaotic-growth explanation above is correct as far as it goes, but it
doesn't answer the sharper question: *is the divergence a generic property of
comparing any two Isca compiles, or is it caused by something specific?* This
was checked directly, using a plain `socrates_aquaplanet` control case (no
ice/land-SST physics involved at all) run the same way (1 month +
a 5-day daily-output growth check) across several pairings:

| Comparison | Result |
|---|---|
| `master` vs. `uob_fftw_sit23_physics_2026` | **Bit-identical.** 0.0 difference in every field, every day, for the full month. |
| Merge-base (2023) vs. current `master` | **Bit-identical**, bar one ~10<sup>-7</sup>-relative rounding blip in a diagnostic-only flux field that never reaches any dynamical field. |
| Original `uob_fftw_sit23_new_maths2intel` vs. `physics_2026` | **Diverges** — same signature as the Maunder Minimum comparison above. |

So `master` has been essentially perfectly reproducible for three years and
372 commits, and the physics port doesn't disturb that at all. The divergence
is real, and it's specific to the original branch.

**Root cause, isolated by bisection:** every one of the ~35 source files the
original branch modifies (relative to its merge-base) was reverted to the
merge-base version in a scratch copy of the branch, one at a time / in
batches, re-testing against `physics_2026` after each step, until the
divergence disappeared and reappeared cleanly. Two environment-only patches
(the `affinity.c` glibc fix and the `spectral_dynamics.F90` format-string fix,
both already on `master`) were re-applied throughout purely so the branch
would compile on this machine's current toolchain - they're unrelated to the
result. This narrowed the entire divergence to a single file:
**`src/atmos_param/socrates/interface/socrates_interface.F90`** — the
branch's own (never-merged, unrelated to this port's `mixed_layer.F90`
changes) modifications to the Socrates radiation coupling code. Reverting
just that one file, with everything else at the branch's tip, restores
bit-identical output; restoring just that one file, with everything else at
merge-base, reproduces the full divergence on its own.

The individual changes within that file (syntax cleanups from `.eqv. .true.`
to bare booleans, an explicit `real(kind(r_def))` cast on `coszen`/`rrsun`
that's a no-op given `r_def` resolves to the same double-precision kind as
this build's default real, six new gated-off diagnostic registrations, and a
reordered spectral-filename validity check) don't individually look like they
should change any computed value under this build's configuration — yet the
whole file, as a unit, does. That residual is a genuine loose end in the
original branch's own Socrates interface work, not yet pinned to one exact
line. It doesn't affect this PR: **`uob_fftw_sit23_physics_2026` never
touches `socrates_interface.F90` at all** — confirmed byte-identical to
`master`'s current version — so this file's issue, whatever it precisely is,
isn't carried into the physics port and isn't part of what's being merged
here.

## Bottom line

- The physics port compiles and runs the real Maunder Minimum experiment
  successfully on top of current `master`.
- The specific physics carried over (sea-ice/land SST patching, ice albedo,
  q-flux masking) is reproduced correctly — every field tied to it matches
  the original branch exactly.
- `master` and `physics_2026` are bit-identical to each other on cases the
  port doesn't touch, and master has been bit-identical to its own 2023
  merge-base for three years — so the port introduces no detectable
  regression of its own.
- The Maunder Minimum comparison's divergence from the *original* branch is
  real, but traced to that branch's own unported `socrates_interface.F90`
  changes — a file this port never touches - not to chaos alone and not to
  anything in the physics being merged here.

## Appendix: how this was tested

Both branches were compiled and run directly (not via `trip_test`, whose
shared-namelist-across-commits design can't represent a fix that adds a new
namelist key only one side understands) via a small standalone script that:
imports the experiment script's namelist/diagnostics/input-file configuration
once, builds a `SocratesCodeBase` against each branch's worktree, applies the
`ice_albedo_method` substitution only for the new branch, compiles, and runs
`exp.run(241, use_restart=True, ...)` from the existing `res0240.tar.gz`
restart file, then diffs the resulting `atmos_monthly.nc`/`atmos_daily.nc`
with `xarray`.

Getting `uob_fftw_sit23_new_maths2intel` to compile at all on this machine's
current toolchain (gfortran 12.4.0, three years newer than the branch's last
commit) required two environment-compatibility patches, applied locally to
the test worktree only, not committed to either branch, and unrelated to
physics:

- `src/shared/mpp/affinity.c`: the branch's own `static gettid()` declaration
  conflicts with glibc &ge;2.30's own `gettid()` — master already carries a
  version-guarded fix for this.
- `src/atmos_spectral/model/spectral_dynamics.F90`: two `FORMAT` statements
  (json-logging labels 300/400) are missing a comma between edit descriptors,
  tolerated by older compilers but a hard runtime error under this gfortran —
  master already carries the fix (commit `36e7398b`).

Both fixes are already present on `master`/`uob_fftw_sit23_physics_2026`
through independent, unrelated work; they only needed re-applying to the old
branch to get a fair comparison running on this specific machine.

# AGENTS.md

Guidance for AI coding agents (and humans in a hurry) working on
[ExeClim/Isca](https://github.com/ExeClim/Isca). It summarises how the project
is developed and points to the full docs; where this file and `docs/` disagree,
`docs/` and the maintainers win. `CLAUDE.md` just imports this file.

## What Isca is

An idealised global atmospheric model framework. The model is Fortran (an
outgrowth of GFDL's FMS-based models); a Python package, `isca`, compiles,
configures and runs it. Scientists rely on Isca giving the **same answers**
for the same configuration across versions, so reproducibility is the
project's first concern when changing anything.

## Isca's philosophy

- Isca is built to support a **hierarchy of models**. Many of its
  parameterisations and settings can be changed and configured, so the same
  framework covers set-ups from very simple to more complex.
- It is therefore best suited to mechanism-style experiments, testing
  hypotheses, and similar idealised work.
- It should **not** be used as an operational weather-forecasting tool, or
  as an operational climate model for policy-relevant predictions.

For agents: keep this purpose in mind when proposing changes or experiments.
For example, prefer configurable, switchable options over hard-wiring extra
"realism" into the model.

## Repository layout

| Path | Contents |
|---|---|
| `src/atmos_*`, `src/coupler`, `src/shared`, ... | Fortran model source |
| `src/extra/python/isca/` | Python front end (`CodeBase`, `Experiment`, `Namelist`, `DiagTable`) |
| `src/extra/model/<codebase>/path_names` | List of Fortran files compiled for each codebase (`isca`, `dry`, `grey`, `column`, `socrates`, `socrates_column`, `shallow`, `barotropic`) |
| `src/extra/env/<name>` | Site/compiler environment files, selected by `GFDL_ENV` |
| `src/extra/python/isca/templates/` | mkmf compiler templates |
| `exp/test_cases/` | Reference experiments; also the inputs to the trip tests |
| `exp/test_cases/trip_test/` | Regression test harness (`trip_test_command_line`, `trip_test_functions.py`) |
| `exp/site_specific/` | Notes for running on particular HPC systems |
| `docs/source/` | Sphinx/reST docs (published at execlim.github.io/Isca) |
| `ci/*.yml` | Conda environment files (`environment-py3.12_frozen.yml` is the tested one) |
| `test/` | pytest compile tests run by GitHub Actions |
| `postprocessing/` | Post-processing tools (e.g. pressure-level interpolation) |

The default branch is **`master`**.

## Setting up and building

Full instructions: `ReadMe.md` and `docs/source/install.rst`.

```bash
conda env create -f ci/environment-py3.12_frozen.yml   # gfortran + OpenMPI + netCDF
conda activate isca_env
pip install -e src/extra/python                         # editable install of `isca`

export GFDL_BASE=/path/to/this/checkout   # source tree that gets compiled
export GFDL_ENV=ubuntu_conda              # file in src/extra/env/
export GFDL_WORK=/path/to/isca_work       # build + run scratch
export GFDL_DATA=/path/to/isca_data       # model output
```

- Experiments compile on first use: `cd exp/test_cases/held_suarez && python held_suarez_test_case.py`.
- `GFDL_BASE` decides which source tree is built. If you work in a git
  worktree or a second clone, point `GFDL_BASE` at it, or you will silently
  build and test a different checkout. Env files are looked up relative to
  `GFDL_BASE`, so copy any untracked env file you rely on into the worktree.
- **Intel compilers / oneAPI.** Add or reuse a `GFDL_ENV` file that sets
  `F90`/`CC` (for example `mpiifort`/`mpiicc`) and `GFDL_MKMF_TEMPLATE`. See
  `isca_intel_2021`, `maths2` and `emps-gv` in `src/extra/env/` for examples.
  Then `source <oneapi>/setvars.sh` before building. Things that commonly go
  wrong in **non-interactive** scripts (batch jobs, agent shells):
  - `conda`/`mamba activate` can silently no-op unless the shell hook is
    sourced first (`etc/profile.d/conda.sh`, plus `mamba.sh` for mamba).
  - Link errors like `cannot find -lsz` / `-lcurl` mean the link-time path is
    missing: `export LIBRARY_PATH=$CONDA_PREFIX/lib:$LIBRARY_PATH`.
  - Undefined `GLIBCXX` symbols mean a stale library from another conda env
    is being picked up: put `$CONDA_PREFIX/lib` *first* in `LD_LIBRARY_PATH`,
    after sourcing `setvars.sh`.
  - Intel MPI failing in `PMPI_Init` on a single node (fabric/RDMA errors):
    `export I_MPI_FABRICS=shm`.
- **Socrates** (Met Office radiation scheme) is not shipped with Isca. On
  `master`, `SocratesCodeBase` symlinks `$GFDL_SOC` to
  `src/atmos_param/socrates/src/trunk`; see
  `exp/test_cases/socrates_test/README.md`. If `GFDL_SOC` is unset or points
  at a missing directory, any Socrates codebase raises `OSError`.
- Building writes `Makefile`, `*.o` and `*.x` into `$GFDL_WORK`, not the
  repo. Full compiles take minutes, so don't recompile in a loop.

## Testing

There are two layers. **A green GitHub Actions run does not validate a
Fortran change.**

1. **`python -m pytest`** (what CI runs, see `.github/workflows/tests.yml`):
   checks that each non-Socrates codebase compiles. Nothing is run.
2. **Trip tests** (`docs/source/testing/trip_test.rst`): build two commits,
   run each test case briefly with both, and check the output is
   bit-for-bit identical. These are **required for any change to model
   (Fortran) code**. Python-only front-end changes are exempt. CI does not
   run them; the contributor does.

```bash
cd $GFDL_BASE/exp/test_cases/trip_test
./trip_test_command_line <master_HEAD> <branch_HEAD> -n 8 -r <repo>
```

Points that trip people (and agents) up:

- **`-r` defaults to an SSH GitHub URL.** In shells without GitHub SSH keys,
  pass a local path (`-r $GFDL_BASE`, which works if both commits exist
  locally) or an HTTPS URL. A permission error here is easy to misread as
  "commit not found".
- **`-e` selects which cases to run** (default `all`); it is an
  include list, not an exclude list. If Socrates is not set up, list the
  cases explicitly, omitting `socrates_aquaplanet`,
  `socrates_aquaplanet_cloud` and `ape_aquaplanet`. Otherwise the uncaught
  `OSError` stops the whole run, not just those cases. Say in the PR which
  cases you skipped and why. The full list is in
  `list_all_test_cases_implemented_in_trip_test()` in `trip_test_functions.py`.
- **Read the printed summary, not the exit status.** The script exits 0
  even when cases fail. A case that crashes on *both* commits is a
  pre-existing or environment problem, not your regression, but report it.
- **Stale output is reused silently.** Run directories are named
  `<case>_trip_test_21_<commit>` under `$GFDL_DATA`, and `Experiment.run`
  will not overwrite existing data. If you rerun the same commit pair after
  changing anything else (`GFDL_BASE`, env, inputs), delete those
  directories first, or you will be comparing old output.
- **Namelists are shared between the two commits.** The harness imports each
  test case's namelist once, from the current `GFDL_BASE`, and feeds it to
  both binaries. So a *new* namelist variable can only be tested at its
  default. Setting it explicitly makes the old binary fail with
  `Unknown namelist ... IOSTAT = 5010`. To show that a new option reproduces
  old behaviour when set explicitly, write a small script that builds each
  commit with `CodeBase.from_repo(repo=..., commit=...)`, gives each commit
  its own namelist, runs both with `Experiment`, and diffs the output with
  xarray.
- Trip tests compare two commits **on the same machine and compiler**.
  There is no stored reference output, so don't add one.

## Coding conventions

- **New features are off by default.** Existing configurations must give
  bit-identical results after your change unless the change is a deliberate,
  documented bug fix.
- **Namelist options:** add the variable to the module's `namelist /..._nml/`
  block with a default that preserves current behaviour, and read it with the
  existing `#ifdef INTERNAL_FILE_NML` / `check_nml_error` pattern used across
  `src/atmos_param/`. Document new options in the module docs
  (`docs/source/modules/`).
- **New Fortran files** must be added to the relevant
  `src/extra/model/<codebase>/path_names` file(s), or they will not be
  compiled.
- **New modules or major options** need a test case in `exp/test_cases/`.
  Add it to the trip tests too, via `get_nml_diag()` and
  `list_all_test_cases_implemented_in_trip_test()`.
- **Do not change existing test cases** (scripts, namelists, inputs). They
  are both user documentation and the trip-test reference. Fixing a test case
  that is broken is fine, as a separate, clearly described change.
- **Python:** black is the stated code style. Keep the public `isca` API
  (`CodeBase`, `Experiment`, `Namelist`, `DiagTable`) backward compatible,
  because users' experiment scripts depend on it.
- **Docs:** reST under `docs/source/`. Start new pages from `_template.rst`,
  link them from the local `index.rst`, and add citations to
  `references.rst` (see `docs/source/contributing.rst`).
- **Third-party code** (FMS/GFDL, RRTM, Socrates interface code under Crown
  copyright): keep licence headers, and keep changes to vendored code minimal
  and clearly marked.
- New site/compiler setups go in `src/extra/env/<site>`, with no personal
  paths, usernames or credentials.

## Branches, forks and pull requests

- Contributions arrive as pull requests to `ExeClim/Isca:master`, usually
  from a personal fork (`<user>/Isca`) on a topic branch. Never commit
  directly to `master`.
- One purpose per PR. Check that the diff against `master` contains only
  what you intend to merge (no stray experiment scripts, data or env files).
- For model-code PRs, paste the trip-test summary (the two commit IDs and the
  pass/fail per case) into the PR. Then label it `trip tests passing`, or
  `changes previous results` with an explanation of each expected change.
- Unexpected trip-test failures are investigated before review, not
  explained away.
- Reviewers check: defaults preserve answers; trip-test evidence is present
  and matches the claimed commits; new options are documented and have a test
  case; no generated files or data are included.
- Commits and PRs carry the human contributor as the author. Do not add AI
  attribution (`Co-Authored-By: <AI tool>`, "Generated with ..." footers or
  similar); the person submitting the change is accountable for it.
- Maintainers are research scientists, so reviews can take a while.
  Smaller, well-described PRs get reviewed faster.

## Things an agent must NOT do

- Push to `master`, force-push shared branches, rewrite published history,
  or merge PRs. Leave pushing and PR creation to the human you work for
  unless they have explicitly asked you to do it.
- Add AI attribution lines (`Co-Authored-By: <AI tool>`, "Generated with
  ..." footers) to commits or PR descriptions.
- Change model answers silently. That includes flipping a default, changing
  a constant, or "tidying" numerics in a way that alters results without
  saying so in the PR.
- Edit test cases, trip-test case lists or comparison tolerances to make
  trip tests pass.
- Claim trip tests pass without having run them, or report a partial run
  (some cases skipped, stale output) as a full pass.
- Commit build or run products (`*.x`, `*.o`, `Makefile`s, `exec.*`,
  `compilation*/`, `$GFDL_WORK`/`$GFDL_DATA` contents), large NetCDF files,
  Socrates source, or anything containing personal paths, credentials or
  tokens. New small input `.nc` files for a test case need a clear reason in
  the PR.
- Modify someone else's working tree or local branches. Use a separate clone
  or `git worktree` with its own `GFDL_BASE`.
- Run large or many parallel model jobs on shared machines without being
  asked. Use a few cores (`-n`) and short runs.

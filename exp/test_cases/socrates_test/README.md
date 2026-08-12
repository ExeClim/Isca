# Using Isca with the Socrates radiation scheme

The Socrates radiation scheme is a highly-flexible scheme written and maintained by the UK Met Office. It has many significant advantages over RRTM, notably the ease with the properties of the radiation scheme can be changed using input-files, e.g. different numbers of spectral bands, different atmospheric compositions, etc.

Socrates is open source (BSD-3-Clause) and hosted at [github.com/MetOffice/socrates](https://github.com/MetOffice/socrates).

## Getting started

### 1. Get the Socrates source code

All the code needed to interface Isca and the Socrates radiation scheme is provided in the folder

`$GFDL_BASE/src/atmos_param/socrates/interface`,

and this is packed with Isca itself. The Socrates source code itself is **not** packaged with this repository, but you don't need to download it by hand: `SocratesCodeBase`/`SocColumnCodeBase` fetch it automatically the first time you compile, e.g.

```python
from isca import SocratesCodeBase
cb = SocratesCodeBase.from_directory(GFDL_BASE)  # uses the default, tested Socrates version
cb.compile()
```

The first time this runs it fetches a sparse, partial clone of the matching Socrates source (a few MB, not the whole ~300MB repository) from GitHub into a cache directory, and symlinks it into `$GFDL_BASE/src/atmos_param/socrates/src/<version>`. If you'd rather have it as a proper git checkout, you can instead run:

```
git submodule update --init
```

from the top of this repository, which populates the default version's directory the standard git way (see `.gitmodules`). Either way, once the source is present at `src/atmos_param/socrates/src/<version>`, Isca will use it as-is and won't try to fetch it again.

If you need a specific Socrates release, pass its GitHub tag as `socrates_version`, e.g. `SocratesCodeBase.from_directory(GFDL_BASE, socrates_version='2025.12.1')`. If your machine doesn't have outbound network access (e.g. an HPC compute node), fetch the version you need in advance somewhere that does, and either set `GFDL_SOC_DIR` to point at a directory containing a `<version>/` subfolder for it, or set `GFDL_SOC` to point directly at a single checkout. See `docs/source/modules/socrates.rst` for the full details of how version selection and fetching work.

### 2. Edit the number of angles in the phase function

From this point, the Isca test case for Socrates should run without issue. However, it is strongly advised to first make the following changes to the Socrates source code.

Open the file

`$GFDL_BASE/src/atmos_param/socrates/src/<version>/src/modules_core/dimensions_spec_ucf.F90`

and make the following changes:

```
npd_k_term=14
npd_scale_variable = 4
npd_continuum = 2
npd_drop_type = 5
npd_cloud_parameter = 30
npd_phase_term = 1
```

The final one of these is the most important, as a large value for this term significantly increases Socrates' memory usage, and will make Isca much slower.

### 3. Run the Socrates test case.

Navigate to

`$GFDL_BASE/exp/test_cases/socrates_test/`

and run the test-case `socrates_aquaplanet.py`. This will compile and run Isca with a default Earth-like aquaplanet.


### 4. If Socrates seems slow, here are some hints

* Isca is set up to pass socrates a certain number of vertical profiles for each time the Socrates is called. This number is set as `chunk_size` in the `socrates_rad_nml`. A value of 16 was found to be optimal on large linux-server-type machines in Exeter, but it is worth playing with this number to find the optimal number on your machine.

* Socrates reads external input files that tell it the number of spectral bands to use, with one file setting the short-wave options, and another file setting the long-wave options. Some spectral files have lots of bands, which will make the model run slowly. The default files used in the Met Office's Unified Model-GA7, and also in Isca, can be found here (`<version>` is whatever `cb.socrates_version` was used to compile):
	* `$GFDL_BASE/src/atmos_param/socrates/src/<version>/data/spectra/ga7/sp_lw_ga7` for the long-wave
	* `$GFDL_BASE/src/atmos_param/socrates/src/<version>/data/spectra/ga7/sp_sw_ga7` for the short-wave
* Other spectral files are available directly from the [Socrates GitHub repository](https://github.com/MetOffice/socrates/tree/main/data/spectra) (Isca's automatic fetch only pulls the `ga7` and `ga3_1` sets its own test cases use -- widen your local sparse-checkout, or fetch a specific file by hand, if you need a different one), and a useful set of other spectral files are provided via this [webpage](https://simplex.giss.nasa.gov/gcm/ROCKE-3D/).

### 5. If you find Socrates doesn't compile

* Isca derives the list of Socrates files to compile for a given version automatically (see `docs/source/modules/socrates.rst` and `isca.socrates_paths`), from Socrates' own build manifests plus a scan of what Isca's interface code actually uses. A handful of tested versions have this list already committed at `$GFDL_BASE/src/extra/model/socrates/socrates_version_paths/<version>`; for any other version it's generated on the fly the first time you compile, and logged as it does so.

* If compilation still fails for a version, it's most likely because Socrates changed something the automatic generation doesn't know how to follow statically (e.g. a genuinely new dependency reached only through a code path the generator's static scan doesn't cover), rather than a missing/renamed file, which the generator already handles by re-deriving the list from the checked-out source. The compiler error will name the missing symbol or module; search the Socrates source tree for it, and either add the containing file directly to `cb.extra_path_names` before calling `cb.compile()`, or fix it in the generator itself (`src/extra/python/isca/socrates_paths.py`) and regenerate with `src/extra/python/scripts/generate_socrates_path_names.py` so future users of that version don't hit the same issue.

* You may also find errors of the type `multiple definition of MAIN`. This happens because various pieces of code within Socrates are designed to be run offline, meaning they have MAIN sections. These MAIN sections cause conflicts with Isca's MAIN section, and so the relevant Socrates files cannot be part of Isca's compiled version. The automatic generator already excludes Socrates' own standalone driver programs and unrelated components (correlated-k/scattering database tools, COSP, flexchem, illumination, its newer high-level driver layer) for exactly this kind of reason -- if you still hit this, remove the offending file from `cb.extra_path_names` before compiling.

* If you still can't get Isca to compile with Socrates, then check your compiler options against ours in `$GFDL_BASE/src/extra/python/isca/templates/mkmf.template.ia64`.

# Topography Source Files
## ERA-5

Isca provides a number of preconstructed topography files for use with common resolutions, which can be found under `Isca/input/era5_smoothed_topography_land_masks`. These were constructed with the script `Isca/src/extra/python/scripts/create_era5_topography.py` and make use of the ERA-5 invariant fields as an input. The ERA-5 invariant fields are not included as part of Isca, however they may be obtained by using either directly from the [CEDA](https://catalogue.ceda.ac.uk/uuid/2c8f38fac04945b89cf12d6e9c928c6f/) archive (recommended) or from the ERA-5 [CDS](https://cds.climate.copernicus.eu/datasets). Note that a login is required for both websites.

If downloading from CEDA the following files are required:
- Topography:    `ecmwf-era5_oper_an_sfc_200001010000.z.inv`
- Land-sea mask: `ecmwf-era5_oper_an_sfc_200001010000.lsm.inv`

## Other Sources/Planets
In theory the script `create_era5_topography.py` may be used to create smoothed/filtered topography input files from other data sources. You may need to change input variable names or disable the land-sea mask, but in principle the topography of other planets (e.g. Mars) can be used.


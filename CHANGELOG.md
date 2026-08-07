
# Version X (Date)


## Major changes:

### General changes
- Reworked the scaling workflow, which should now be broadly automated rather than relying on two pipeline runs
- Updated syntax to support the latest Nextflow (26.04), which enforces strict syntax by default; required a lot of refactoring
- Added a validation routine for input files
- Moved processes downstream of finalize from storeDir to publishDir. Folders that are populated this way (i.e. their contents are symlinked back to the cached files in the workdir, rather than being independent copies) are prefixed with `rr__`. This means these will now be re-run if their input changes. Old caching style is maintained for cellpose, decon, images, flatfields and registration, so these will still need to manually be removed to force a rerun. 
- Removed the subcell workflow, it was out of date


### Input & parameter changes
- Dropped `-entry` (no longer supported by NF > 26.04), added `--workflow` to control which workflow gets executed
- Changed how PSFs are specified in the manifest, now in `<channel>=<path>` format; dropped `dc_channels` as it is no longer needed
- Dropped `channels` column in manifest, it wasn't used
- Several configurations have been renamed, please see the parameters header below.


## Bugfixes
- Fixed index_images & index_cellcrops re-triggering unnecessarily due to `.last()`; now uses a fingerprint-based approach that checks modification times without needing to stage every input file

## Minor changes:
- Plates are now automatically renamed during staging to the plate name provided in the manifest, instead of the name from the PE index file
- Updated flatfield fitting to drop the requirement for a basicpy install
- Fixed `-stub-run` to fire properly after recent updates
- Added a script to perform a stub run, useful for testing that the wiring is working properly
- Updated the cellcrops process to produce a single parquet file per field instead of one CSV per field
- Cellcrop indexing now aggregates into a single parquet file instead of a CSV
- Added new configuration for the Sanger CUB cluster
- Updated the GPU queue for the Sanger farm22 cluster
- Renamed the manifest parameter `bp_channels` to `ff_channels`
- General code cleanup, removed deperacted scripts and comments

## Parameters:
- Added `--workflow`
- Removed subcell parameters
- `sc_` prefix now used for scaling parameters
- `bp_` prefix (basicpy) renamed to `ff_` (flatfield)
- `tg_` prefix (tglow) renamed to `rn_` (run); applies to `tg_conda_env` and `tg_container`
- Removed `cp_cell_power`, `cp_nucl_power`, `rn_dummy_mode`, `rn_threshold` and `cp_dont_postprocess` options
- Removed `executor.poolSize` from the LSF definition
- Removed `dc_channels` from the manifest

## Documentation:
- Fixed small documentation issues that didn't match actual behaviour


# Version pre X
Was not keeping a changelog at this time

# Version X (Date)


## Major changes:
- Removed subcell workflow, it was out of date
- Re-worked scaling workflow which should now be broadly automated rather then relying on two pipeline runs
- Updated syntax to allow use of latest nextflow 26.04 (which default forces strict syntax, lots of refactors)
- Dropped -entry, added --workflow to control entrypoint
- Added validation routine for input files
- Changed handling of PSFs in manifest, now in <channel>=<path> format, dropped dc_channels
- Moved processes downstream of finalize to publishDir instead of storeDir

## Minor changes:
- Plates are now automatically renamed during staging to the plate name provided in the manifest instead of the name provided in the PE index file.
- Updated flatfield fitting to drop requerement for basicpy install. 
- Updated -stub-run to fire properly after updates
- Added script to perform stub run, usefull for testing the wiring is working properly
- Updated cellcrop process to produce a single parquet file instead of one csv per field
- Cellcrop indexing now produces a parquet instead of csv file
- Added new configuration for Sanger Cub cluster
- Updated GPU queue for Sanger farm22 cluster
- Updated parameter name in manifest from bp_channels to ff_channels

## Parameters:
- Added --workflow
- Removed subcell parameters
- sc_ prefix now used for scaling parameters
- bp_ prefix (basicpy) to ff_ (flatfield)
- tg_ prefix (tglow) now rn_ (run). Applies to tg_conda_env and tg_container
- removed cp_cell_power, cp_nucl_power, rn_dummy_mode, rn_threshold and cp_dont_postprocess options
- removed 'executor.poolSize' from lsf definition
- removed dc_channels from manifest

## Documentation:
- Updates to fix small issues in documentation not matching real behaviour



# Version  pre X
Was not keeping a changelog at this time
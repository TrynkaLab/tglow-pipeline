
# Version X (Date)


## Major changes:
- Removed subcell workflow, it was out of date
- Re-worked scaling workflow which should now be broadly automated rather then relying on two pipeline runs
- Updated syntax to allow use of latest nextflow 26.04 which default forces strict syntax
- Dropped -entry, added --workflow to control entrypoint

## Minor changes:
- Plates are now automatically renamed during staging to the plate name provided in the manifest instead of the name provided in the PE index file.
- Updated flatfield fitting to drop requerement for basicpy install. 
- Updated -stub-run to fire properly after updates
- Added script to perform stub run, usefull for testing the wiring is working properly
- Updated cellcrop process to produce a single parquet file instead of one csv per field
- Cellcrop indexing now produces a parquet instead of csv file
- Added new configuration for Sanger Cub cluster
- Updated GPU queue for Sanger farm22 cluster

## Parameters:
- Added --workflow
- Removed subcell parameters
- sc_ prefix now used for scaling parameters
- bp_ prefix (basicpy) to ff_ (flatfield)

## Documentation:
- Updates to fix small issues in documentation not matching real behaviour



# Version  pre X
Was not keeping a changelog at this time
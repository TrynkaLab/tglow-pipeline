
# Version X (Date)


## Major changes:
- Removed subcell workflow, it was out of date
- Re-worked scaling workflow which should now be broadly automated
- Updated syntax to allow use of latest nextflow 26.04 which default forces strict syntax
- Dropped -entry, added --workflow to control entrypoint

## Minor changes:
- Plates are now automatically renamed during staging to the plate name provided in the manifest instead of the name provided in the PE index file.
- Added new configuration for Sanger Cub cluster
- Updated GPU queue for Sanger farm22 cluster
- Updated flatfield fitting to drop requerement for basicpy install. 


## Parameters:
- Added --workflow
- Removed subcell parameters



Documentation:



# Version  pre X
Was not keeping a changelog at this time
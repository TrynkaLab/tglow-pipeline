# Tglow: Nextflow pipeline for analyzing HCI tglow data


# Installation & dependencies
The pipeline currently uses conda as package manager.

This repo relies on the core component of tglow which should be installed into python
and the runners should be availabkle.
https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-core

# Running the pipeline

The pipeline runs in two stages. One where data is staged, this can be done manually, or through the workflow 'stage' if the raw data has a perkin elmer index.xml or index.idx.xml file

Then the pipeline can be configured and run using workflow 'run_pipeline'

## Options
See nextflow.config for available options and their descriptions

## Workflows

### 1. stage:

This workflow takes a perkin elmer (currently only works for phenix) index xml and stages the files into a plate/row/col/field.ome.tiff file structure. IO is handled through AICSImageio, channel names and pixel sizes are extracted from the index.xml and set as metadata items in the ome tiff. Some additional files are staged for reproducabillity and ease of reading channel orders etc. This output is also intended to be backed up to iRODS

Proccess:
1. Create a manifest with the wells to run to re-use later as a nextflow channel
2. Read the files from /nfs and save directly into the above format

### 2. run_pipeline:

This takes as input the images produced during staging and then runs the following processes.

Processes:
1. Basicpy [optional]: parallelized on the plate + channel level, flatfields are saved, no images are stored
2. Registration (pystackreg) [optional]: parallelized on the well level (all fields). Only registration matrices are saved, no images are stored
3. Cellpose: parallelized on the well level (all fields), runs on GPU. If registering, currently only runs on the reference plate. Must run a cell segmentation, can optionally provide nucleus channel as well.
4. Deconvolition [optional]: Will spin out another copy of the data, runs on GPU
5. feature extraction / cellprofiler: 
    1. Stage the files into a cellprofiler compatable format, apply flatfields, bring together imaging cylces and apply registration if applicable
    2. run feature extraction (cellprofiler) or other custom script
<br>

## Quick explanation

Install as indicated above

Adjust the required configurations either in a config file or on the nextflow command (see nettflow docs for more detaills). See nextflow/nextflow.config for available parameters.

prepare a manifest.tsv with one line per plate:
```
plate	index_xml	channels	bp_channels	cp_nucl_channel cp_cell_channel
ref_plate   index.xml   1,2,4,5 3,5 5   2  
qry_plate1  index.xml  2,3,4 none   none    none
qry_plate2  index.xml  3,4 3,4 3   4  

```
- plate: The exact plate name in the export
- index_xml: path the perkin elmer (or Martin's version) index xml file
- channels: list of available channels
- bp_channels: channels to run basicpy for
- cp_nucl_channel:optional channel for nuclei segmentation
- cp_cell_channel: cellpose channel for cell segmentation


First stage the data so its accesible on the lustre. IMPORTANT, as it needs access to /nfs to stage the index xml, run this through a job in the imaging queue, otherwise it will crash.

```
OUTPUT=./output

nextflow run tglow.nf \
-profile lsf \
-w ${OUTPUT}/workdir \
-resume \
-entry stage \
--rn_manifest manifest.tsv \
--rn_publish_dir ${OUTPUT}/results \
--rn_image_dir ${OUTPUT}/results/images
```

After the data is staged into rn_image_dir, build the cellprofiler pipeline. To find the channel order check the images/plate/aquistion_info.txt file which contains a channel list and voxel resolutions. 

For extracting the nucleus masks configure pattern as such:
`^(?P<field>\d+)_(?P<mask_type>.*_mask)_d\d+_ch\d+_cp_masks.tif`

Then for the channels configure an additonal pattern only matching to ome.tiff files not containing the _masks keyword
`^(?P<field>\d+)_(?P<plate>.*)_(?P<well>\w\d\d)_ch(?P<channel>\d+)`

Channel order of registation channels is explained below in the registration section.

Then run the pipeline
```
OUTPUT=./output

nextflow run tglow.nf \
-profile lsf \
-w ${OUTPUT}/workdir \
-resume \
-entry run_pipeline \
--rn_manifest manifest.tsv \
--rn_publish_dir ${OUTPUT}/results \
--rn_image_dir ${OUTPUT}/results/images
```


## Blacklisting troublesome wells

Sometimes wells might have issues or have different aquisition parameters (e.q. number of channels or planes). As currently the pipeline assumes all image stacks have the same shape, this is not supported. Wells which not match can be excluded in two ways. This is also usefull for excluding single stain controls which don't make sense to run through the pipeline.

1. Edit the manifest.tsv in the `rn_image_dir` of the plate to remove the well
2. Specify `--rn_blacklist` which is  a tsv file with header `plate<tab>well` and then specifying a plate well combination, one per line, to exclude from the run. 

I reccomend the blacklist prior to staging, especially if many wells need to be blacklisted

## Registering multiple cycles of imaging

To register imaging cycles additionally provide a registration manifest which links plates together. This should have the form:

```
reference_plate	reference_channel	query_plates	query_channels
ref_plate   5   qry_plate1,qry_plate2   3,4
```

- reference_plate: Plate name of reference plate
- reference_channel: 1 based index of channel to use in registering (nucleus)
- query_plates: comma seperated list of plate names to register against reference
- query_channels: comma seperated list of 1 based channel indices of channel to register (nucleus)


IMPORTANT: this assumes that all wells in the ref_plate manifest are available in the query plate. Currently this is NOT automatched. To fix this you can:

1. Go into the image folder and override the manifest.tsv to only include the wells that overlap between the plates. It's suggested to make a copy of the original to avoid having to re-run the stage workflow.
2. Supply the wells that DO NOT match to `--rn_blacklist`. See detaills above

If this is provided, in downstream tasks (cellprofiler/feature extraction) data is treated as one plate and the query plates are treated as extra channels, with their indices increasing seqeuntially in the order specified in the manifest. 

Query plates must be provided in the manifest.tsv!
Currently cellpose is only run on reference plates, even if the channels are provided in the registration manifest. If needed, will update cellpose to run on the registered data so other cycle channels can be used in segmenting.

## Extra detaills
The nextflow pipline runs in two stages. Some of it is not done through fully "proper" nextflow, as nextflow is very storage heavy and can easily store redundant copies of data which is not great for imaging. 

Hence most processes rely on the storeDir option which serves as a permanent cache between runs. Keep in mind that as long as nextflow finds the files in this storeDir, processes are NOT re-run, so you will have to remove them manually if you want to re-do some steps! This is the case for the stageing of the raw images, the deconvolution, registration and basicpy. This is not "proper" nextflow in the sense everything is localized to the workdir, but it made the most sense in this case as:
1. We wanted to save on as much redundant stroing of images as possible (and automated cleanup of workdirs is not currently available in nextflow)
2. Many of the processes rely on eachother yet don't have a natural  a > b > c structure but represent a complicated tree. 

Parallelization is done on the per well level to not overload the system with 20 second jobs wasting a lot of time on overheads (loading conda, initializing python etc). This makes some of the processes implicitly assume all the fields are supposed to be run. Given the same well always has the same number of fields this should work fine. But nextflow itself is not aware of the fields! This pipeline doesn't currently support seperate fields between cycles, but this should ideally never happen anyway. If it does happen, first stage the data, then remove or rename the fields so they match manually in the storeDir, then it should work fine.

 
Note to self: If the pipeline does need to be 'field aware' most of the runners implement the --fields argument, as does the io reader, so it should be trivial to implement should the need arise by adding a field to the plate level manifest configuring which fields to run.


The first workflow involves staging the data from an Harmony export on NFS. This is done with the -entry stage

```
nextflow <path/to/tglow.nf> -entry stage <params>
```

There are a bunch of configurations to be set, have a look at the nextflow.config for detaills. In general, you make a manifest in the following form

```
<plate> <index.xml> <channels> <basicpy_channels> <cellpose_nucleus> <cellpose_cell>
```
that tells the pipeline where the data lives, and which channels are which and how to apply them. The file can run multiple plates at the same time by adding more lines.




# Known issues

## Basicpy not producing proper flatfields
Basicpy is quite sensitive to cell density, and I found it needs a proper coverage of foreground signal to work. It might take a little tweaking to get it to work well (see the options available). Alternatively, if the dataset is very sparse (10-20 cells per field) it might be better to only fit one model for all plates in an experiment run manually, rather then fitting one model per plate. This can easily be done, by just putting the basicpy models in the correct format in the results folder. See example bash scripts here (tbd) to manually fit a basicpy model per plate. On the TODO list is to generalize the correction options to use the PE flatfields if available (but these also stuggle with low density from quick inspection), or use some other strategy (averaging etc).

## Stage Cellprofiler not finding decon results when registering
If you get crashes of some cellprofiler processes when deconveluting with multiple cycles, this is due to the fact the decon output is only syncronized to the reference plate, but currently ignores the 2nd cycle. This is a bug which needs fixing. The workarround is to first run with --cpr_run false to make sure all the decons are done, then run with --cpr_run true to run cellprofiler after that instance has completed.

## 'file not found' exceptions on the Index.xml files.
Sometimes upon running the pipeline for the first time, the run crashes, with 'file not found' exceptions on the Index.xml files. I have no clue why this happens, maybe its something with the specific node in the imaging queue not having the mount, but usually when I re-run without changing, it seems to run fine the 2nd time. I have had wierd issues like this with other NF pipelines in the past, so think its more a NF issue then anything else. Update: I have not had this recently so perhaps it was a FARM related issue



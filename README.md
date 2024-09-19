# Tglow: Nextflow pipeline for analyzing HCI data


This repo contains the nextflow pipeline and binaries and scripts to run a tglow-pipeline instance for the analysis of high content imaging data.
A detailled walkthrough of the steps, installation and configuration is given on the [wiki](https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-pipeline/-/wikis/home)

There are three components to the overall workflow
1. [tglow-core](https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-core) A python library with IO, parsing and convience functions based on AICSImageIO.
2. [tglow-pipeline](https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-pipeline) Nextflow files and binary python scripts for running pipeline processes
3. [tglow-r-core](https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-r-core) A Seurat like R package for analyzing HCI features at scale


# Installation & dependencies

See [here](https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-pipeline/-/wikis/home/1_installation) for install instructions

The pipeline currently uses conda as package manager through Nextflow.

This repo relies on the core tglow python library of tglow which should be installed into python [tglow-core](https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-core)

# Quick overview of pipeline flow

The pipeline runs in two stages. One where data is staged, this can be done manually, or through the workflow 'stage' if the raw data has a perkin elmer index.xml or index.idx.xml file. Then the pipeline can be configured and run using workflow 'run_pipeline'

## Options
See nextflow.config or the [wiki](https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-pipeline/-/wikis/home/options) for available options and their descriptions.

## Workflows
See [wiki](https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-pipeline/-/wikis/home/2_setup_and_stage) for more detaills.

### 1. stage:

This workflow takes a perkin elmer (currently only works for phenix) index xml and stages the files into a plate_name/row/col/field.ome.tiff file structure. IO is handled through AICSImageio, channel names and pixel sizes are extracted from the index.xml and set as metadata items in the ome tiff. Some additional files are staged for reproducabillity and ease of reading channel orders etc. This output is also intended to be backed up to iRODS.

Nextflow proccess:
1. prepare_manifest: Create a manifest with the wells to run to re-use later as a nextflow channel
2. fetch_raw: Read the files from /nfs and save directly into the above format

### 2. run_pipeline:

This takes as input the images produced during staging and then runs the following processes.

Nextflow processes (in order):
1. basicpy [optional]: parallelized on the plate + channel level, flatfields are saved, no images are stored
2. register (pystackreg) [optional]: parallelized on the well level (all fields). Only registration matrices are saved, no images are stored
3. cellpose: parallelized on the well level (all fields), runs on GPU. If registering, currently only runs on the reference plate. Must run a cell segmentation, can optionally provide nucleus channel as well.
4. deconvolute [optional]: Will spin out another copy of the data, runs on GPU
5. cellprofiler (feature extraction): 
    1. Stage the files into a cellprofiler compatable format, apply flatfields, bring together imaging cylces and apply registration if applicable
    2. run feature extraction (cellprofiler) or other custom script


## Known issues
See the [wiki](https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-pipeline/-/wikis/home) for known issues. If you find an issue please raise it on the git.
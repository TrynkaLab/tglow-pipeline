# Tglow: Nextflow pipeline for analyzing HCI data

This repo contains the nextflow pipeline and binaries and scripts to run a tglow-pipeline instance for the analysis of high content imaging data.
A detailed walkthrough of the steps, installation and configuration is given on the [wiki - tbd]()

There are three components to the overall workflow
1. [tglow-pipeline] - this repo - Nextflow files and Python scripts for running pipeline processes
2. [tglow-core](https://github.com/TrynkaLab/tglow-core) A python library with IO, parsing and convenience functions based on AICSImageIO.
3. [tglow-r](https://github.com/TrynkaLab/tglow-r) A Seurat-like R package for analyzing HCI features at scale


# Installation & dependencies

The pipeline currently uses conda as package manager and relies on pre-creating the conda environments needed. In future we will streamline the install process.

See [here - tbd]() for full install instructions of all pipeline components.

# Pipeline overview

The following readme gives a high level overview, for more detailed guide please see the wiki. The pipeline consists of two main stages:
- stage: prepare and standardize raw images into a well/field-organized OME-TIFF layout with metadata.
- run_pipeline: perform image processing and feature extraction on the staged images.

Both stages are implemented as Nextflow workflows and can be run independently using `-entry stage|run_pipeline`. 

> Some steps in the pipeline require GPU's to be available. These are semgmentation and deconvolution. Deconvolution will not run without GPU. Segmentation (CellPose) will run, but we only reccomend this in cases where you are generating masks in 2D. In 3d the computational burden for large datasets will be too much for CPU. 

> The pipeline is intended to run on high performance compute (HPC) clusters, and bundled resource profiles should work for most HPC, but some tweaks to queue names and GPU settings may be required as flags differ between vendors and HPC configurations. Go to conf/processes.config and search for `queue` and `clusterOptions` to update. Furthermore each HPC is different, with different machines and resource limits. You may need to add a profile for your HPC enviroment in the conf folder. The nf-core config directory may be of help for your HPC: https://nf-co.re/configs/. If something is unclear, feel free to raise an issue on github.  

> If you dont want to run the pipeline on HPC but run it locally, supply `-profile local`.

## 1) stage
Purpose: Stage Revity/PerkinElmer (currently Phenix or Operetta) acquisitions into a reproducible plate/row/col/field.ome.tiff structure and capture metadata (channel names, pixel sizes, channel order, original index files).

> If you don't have a Phenix or Operetta export, you can skip this step, but will need to organize the images using your own script. See more details here: [manual staging - tbd]()

Input:
- PerkinElmer index.xml / index.idx.xml and raw instrument files (or manually organized raw files).

Output:
- plate_name/row/col/field.ome.tiff (with metadata)
- manifest listing wells to process (used as a Nextflow channel)
- auxiliary files to capture provenance (index.xml, channel maps, etc.) (optional)

Nextflow processes:
1. prepare_manifest — create a manifest with wells/fields to run (re-usable Nextflow channel)
2. fetch_raw — read raw files and write standardized OME-TIFFs and metadata

## 2) run_pipeline
Purpose: Run the core image-processing and feature-extraction steps on the staged images. The workflow is modular — many steps are optional or configurable.

Input:
- Staged OME-TIFFs (from `stage`) and optional per-plate/field metadata (flatfields, registration references, etc.)

Output:
- Segmentation outputs, registration matrices, flatfields, extracted feature tables, and logs/artifacts needed for downstream analysis.

Main processing steps (in typical execution order — each step can be enabled/disabled via config):
1. estimate flatfield (Polynomial / BaSiCPY) (optional)
   - Parallelization: per-plate + channel or single flatfield for all plates + channels.
   - Output: flatfield images only (no transformed images saved).
2. register (cross correlation / pystackreg) (optional)
   - Parallelization: per-well
   - Output: registration matrices (no transformed images saved).
3. cellpose segmentation
   - Parallelization: per-well, GPU-enabled
   - Notes: If registration is used, segmentation currently runs on the reference plate. Nucleus channel optional but segmentation is required.
   - Output: 2D or 3D cell & nucleus masks as tiffs
4. deconvolute with CLIJ2-fft (optional)
   - Parallelization: per-well, GPU-enabled
   - Output: deconvolved images (creates a data copy)
5. finalizing images
   - Parallelization: per-well
   - Applies all the registration, flatfields, scaling, max projection to the (deconvolved) images and collects the masks
   - Output: Analysis reade OME-TIFFs
6. feature extraction with CellProfiler
   - Parallelization: per-well
   - Stage images into a CellProfiler-compatible layout, apply flatfields and registration (if enabled), and run feature extraction.
   - Outputs: CellProfiler artifacts as a zip archive per well
7. cellcrops (optional)
   - Parallelization: per-well
   - Produces a HDF5 file for each field where each h5 group is a cell
   - Outputs: h5 file with fully processed cellcrops   

# Options
See nextflow.config or the [wiki - tbd]() for available options and their descriptions.

# Quick Usage

Prerequisites:
- Nextflow and conda
- Create required conda environments (see install instructions)
- Ensure GPU drivers are available in the conda env when running GPU steps (cellpose, deconvolution)

I strongly reccomend to configure through a configuration file. Altough parameters can be overridden on the commandline. I would reccomend a project structure as follows:

- my_project
  - results: By default this is where the pipeline stores outputs
  - scripts
    - logs
    - my_config.config
    - run_pipeline.sh
  - workdir: By default this is the Nextflow workdir

Quick examples:

Stage PerkinElmer data from a raw export:
```
nextflow \
-log logs/stage.nextflow.log \
run </path/to/main.nf> \
-profile <your profile> \
-w ../workdir \
-resume \
-entry stage \
-with-report logs/stage.nextflow.html \
-with-trace logs/stage.nextflow.trace \
-c my_config.config"
```

Run the main pipeline on staged images:
```
nextflow \
-log logs/run_pipeline.nextflow.log \
run </path/to/main.nf> \
-profile <your profile> \
-w ../workdir \
-resume \
-entry run_pipeline \
-with-report logs/run_pipeline.nextflow.html \
-with-trace logs/run_pipeline.nextflow.trace \
-c my_config.config"
```

# Known issues
We will put known issues here, or in the issue tracker. If you find an issue please raise it on the git or contact us directly.

# Authors:
- Olivier Bakker
- Francesco Cisterno

# References
- https://github.com/clij/clij2-fft
- https://cellprofiler.org/
- https://scikit-image.org/
- https://github.com/MouseLand/cellpose
- https://github.com/glichtner/pystackreg/tree/master
- https://basicpy.readthedocs.io/en/latest/
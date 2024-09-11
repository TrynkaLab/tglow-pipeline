# Tglow: Nextflow pipeline for analyzing HCI tglow data


# Installation & dependencies
The pipeline currently uses conda as package manager.


This repo relies on the core component of tglow which should be installed into python
and the runners should be availabkle.
https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-core

# Running the pipeline

The pipeline runs in two stages. One where data is staged, this can be done manually, or through the workflow 'stage' if the raw data has a perkin elmer index.xml or index.idx.xml file

Then the pipeline can be configured and run using workflow 'run_pipeline'

# Setting up the pipeline



# Options
See nextflow.config for available options and their descriptions

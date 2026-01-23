#  pipeline parameters



## Pipeline



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `rn_manifest` | Core manifest set bp_channels to none to not run basicpy set cp_cell_channel to none to not run cellpose set cp_nucl_channel to none to run cellpose without a nucleus channel If a registration manifest is provided, only reference plates will be run through cellpose <plate> <index_xml> <channels> <bp_channels> <cp_nucl_channel> <cp_cell_channel> <dc_channels> <dc_psfs> <mask_channels> <scale_factors> | `string` |  | True |  |
| `rn_blacklist` | Well level manifest (blacklist), supply path to file. Specifying which wells to ignore. This is useful to deal with edgecases in a plate or single stain controls etc. you don't want to analyze. Leave at null to ignore. <plate> <well> plate1 A09 plate2 A09 plate2 D12 | `string` |  |  |  |
| `rn_control_list` | Control file to use to calculate plate offsets based on control cells for each channel <plate> <well> <channels> <name> plate1 A10 1,2,3 controlA plate1 B10 1,2,3 controlA plate1 C10 1,2,3 controlA plate1 D11 5,6,7 controlB plate1 D12 5,6,7 controlB | `string` |  |  |  |
| `rn_manifest_well` | Well level manifest (whitelist) This is a file with <well> <plate> <pe xml> structure indicating which wells on which plates to run Set to null if auto generated from perkinelmer XML If not null, supply a list /path/to/plate1_manifest.csv,/path/to/plate2_manifest.csv,... | `string` |  |  |  |
| `rn_manifest_registration` | Registration manifest This dictates which plates will be merged and registered. If left to null no merging or registration is performed | `string` |  |  |  |
| `rn_publish_dir` | Directory to store and cache the output. Generally name results | `string` | ../results |  |  |
| `rn_image_dir` | Permanent cache for images [optional] Defaults to: ${rn_publish_dir}/images | `string` | ../results/images |  |  |
| `rn_decon_dir` | Permanent cache for deconvolutions [optional] Defaults to: ${rn_publish_dir}/decon | `string` | ../results/decon |  |  |
| `rn_max_project` | Max project prior to running segmentation and cellprofiler Decon is still done in 3d, but results are saved as max projections. | `boolean` |  |  |  |
| `rn_hybrid` | Run in hybrid 2d/3d mode. Masks and decon run and saved in 3d but only cellprofiler is run using max projections. In true ignores rn_max_project. Must be true to enable demultiplexing of nuclear and non-nuclear signals. | `boolean` |  |  |  |
| `rn_wells` | Select only these wells, useful for testing [optional] comma separated string of well ids: A06,B19,C22 | `string` |  |  |  |
| `rn_scratch` | Use scratch space for most workdir operations, which saves IO load on networked filesystems. But this makes debugging harder as tmp results are not available | `boolean` | True |  |  |
| `rn_cache_images` | Cache the final output images (flatfield, decon, demultiplexed, registered, scaled, max_projected) They are saved as plate/row/col/field.ome.tiff in CZYX with additional cycle channels sequentially added. If storage is a concern, disable this and enable rn_scratch to not save large intermediates. NOTE: When mode is hybrid or max project, the storage this generates will not be that bad, so it’s enabled by default. NOTE: This must be true when running subcell, but as subcell only works with max projected images, the storage overhead should be ok NOTE: These are not directly compatible with cellprofiler when working in 3d as they are saved in .ome.tiff for pipeline compatibility, which CPR does not handle. Prior to running cellprofiler they are split up into <field><plate><well>_ch<channel>.tiff which is compatible with both 2d and 3d formats. Use this pattern to set up your pipeline. | `boolean` | True |  |  |
| `rn_make_cellcrops` | Produce h5 files with cellcrops | `boolean` | True |  |  |
| `rn_max_per_field` | When making cellcrops, if there are more then this many cells per field, skip the field | `integer` | 1000 |  |  |
| `tg_conda_env` | Tglow conda env path | `string` | /software/teamtrynka/installs/tglow |  |  |
| `tg_container` | Container for the tglow environment | `string` |  |  |  |

## Staging



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `st_label` | Resource label for staging | `string` | small_img |  |  |

## Flatfield



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `bp_run` | Run flatfield estimation or not. Overrides the manifest | `boolean` | True |  |  |
| `bp_global_flatfield` | Fit one flatfield per plate-channel, or one global flatfield per channel shared among the plates of a cycle | `boolean` |  |  |  |
| `bp_label` | Resource label for flatfield estimation | `string` | himem |  |  |
| `bp_mode` | Mode, one of BASICPY, POLY or PE | `string` | POLY |  |  |
| `bp_threshold` | Should images be multi-otsu thresholded prior to fitting. Two top tiers treated as foreground. Helpful in sparse images but can throw off valuation if strong flatfield not handled well or outliers exist. In BASICPY mode mask is supplied as fitting_weights to fit only actual signal avoiding background fitting. In POLY the polynomial is fit on foreground pixels only. | `boolean` |  |  |  |
| `bp_degree` | Degree for polynomial (2-4 recommended). If <=0 special polynomial model used same as PE Harmony fits, else numpy.polynomial.polynomial.polyvander2d used to generate design matrix with all combinations of x^i*y^i. Higher degrees more likely to overfit. | `integer` | 0 |  |  |
| `bp_channels` | Channels to fit basicpy models on, leave null to run all channels specified in manifest (recommended). To not run basicpy for a plate, set channels to "none". Otherwise specify [[<plate>,<channel>,<index_xml>],...] to override manifest. | `string` |  |  |  |
| `bp_nimg` | Number of random images to read into memory. Sampled with replacement. If no other option specified, this is number of images flatfield is trained on. | `integer` | 200 |  |  |
| `bp_nimg_test` | Number of images for independent sampling used for testing flatfield. Set 0 to skip flatfield evaluation. | `integer` | 100 |  |  |
| `bp_merge_n` | Number of images to max project into a compound — nimg times. If >1 basicpy run on nimg images each compound of merge_n images. Set null to run vanilla basicpy with no merging. Useful if low density images and flatfields tend to background signal rather than foreground. Recommended bp_nimg=100 and bp_merge_n=50 starting point but mileage varies. WARNING: samples same images into different compounds so some overlapping. | `integer` |  |  |  |
| `bp_pseudoreplicates` | Pseudoreplicate in memory related to merge_n but instead of disk I/O only nimg images read then pseudoreplicate compound images of size merge_n made. Sampling with replacement and overlaps possible. Goal similar to merge_n but avoids major IO load. Recommended to set nimg high to reduce overlap. | `integer` |  |  |  |
| `bp_pseudoreplicates_test` | Same as bp_pseudoreplicates but for flatfield evaluation | `integer` |  |  |  |
| `bp_use_ridge` | Use ridge regression instead of OLS to fit polynomial. Uses RidgeCV and 10 fold CV to find optimal alpha | `boolean` |  |  |  |
| `bp_all_planes` | Instead of randomly picking one plane for a stack, use all planes. Can cause issues with basicpy as it assumes random variation between images. When rn_max_project true all planes are read and max projected so not an issue. Not recommended. | `boolean` |  |  |  |
| `bp_autosegment` | Apply basicpy autosegment option, opposite of threshold, applies mask erosion. Not recommended, basicpy only. | `boolean` |  |  |  |
| `bp_no_tune` | Do not tune basicpy model. Not recommended, basicpy only. | `boolean` |  |  |  |

## Registration



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `rg_mode` | Mode use skimage phase cross correlation or pystackreg with translation. CROSS or STACKREG. Currently STACKREG is disabled. | `string` | CROSS |  |  |
| `rg_label` | Resource label for registration | `string` | small |  |  |
| `rg_offset_x` | Offset in X. Positive shifts down (scipy.ndimage.shift convention) | `string` |  |  |  |
| `rg_offset_y` | Offset in Y. Positive shifts down (scipy.ndimage.shift convention) | `string` |  |  |  |
| `rg_eval` | Evaluate registration results with pearson correlation between registration channels. Useful if many signals as images mostly noise. TODO: Threshold image first then correlate | `boolean` | True |  |  |

## Deconvolution



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `dc_label` | Resource label for deconvolution | `string` | gpu_normal |  |  |
| `dc_run` | Run deconvolution or not | `boolean` |  |  |  |
| `dc_psfs` | Point spread functions used for deconvolution. Leave null to not run deconvolution, better specified in manifest | `string` |  |  |  |
| `dc_niter` | Number of iterations to run Richardson Lucy deconvolution for | `integer` | 100 |  |  |
| `dc_psf_crop_z` | Number of planes around PSF center to use for decon. Defaults to all in PSF | `integer` |  |  |  |
| `dc_psf_subsample_z` | Select every x planes starting from PSF center. Useful if PSF image has higher z resolution than actual image. For example PSF 100nm spacing and image 500nm spacing set this to 5. dc_psf_crop_z is applied after this. | `integer` |  |  |  |
| `dc_mode` | Implementation of Richardson Lucy to use. Options: clij2 - clij2-fft richardson_lucy clij2_nc - clij2-fft non circulant richardson_lucy rlf_cpu - RedLionFish CPU mode rlf_gpu - RedLionFish GPU mode RedLionFish not recommended when strong edges in data | `string` | clij2_nc |  |  |
| `dc_regularization` | Regularization parameter for Richardson Lucy. Only used for clij2 and clij2_nc modes. Default recommended. See https://forum.image.sc/t/richardson-lucy-deconvolution-large-images/85325/7 and https://pubmed.ncbi.nlm.nih.gov/24436314/ | `number` | 0.0002 |  |  |
| `dc_clip_max` | Pixel value clip after deconvolution in 32->16 bit conversion. new_intensity=round((intensity/dc_clip_max)65535) Default 565535=327675. Values lower preserved but may lose precision. Value above clip to 65535 in 16bit output. Set 65535 to clip everything outside 16bit range but preserve dynamic range. Set > 65535 to scale values down but keep dynamic range upper end. Keep consistent between runs to interpret intensities correctly. | `integer` | 327675 |  |  |

## Segmentation



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `cp_run` | Run cellpose or not (will stall the pipeline at downstream steps that need masks) | `boolean` | True |  |  |
| `cp_label` | Resource label for cellpose | `string` | gpu_normal_plus |  |  |
| `cp_cell_size` | Estimated cell size for cellpose in pixels. Default is a reasonable estimates for T cells at 0.149um pixel size | `integer` | 75 |  |  |
| `cp_nucl_size` | Estimated nucleus size for cellpose in pixels. Default is a reasonable estimates for T cells at 0.149um pixel size. | `integer` | 60 |  |  |
| `cp_min_cell_area` | Cellpose min cell area in pixels. Minimal area of ROIs. If null defaults to π (1/6th cp_cell_size)^2 for 2d or 4/3π (1/6th cp_cell_size)^3 for 3d. | `string` |  |  |  |
| `cp_min_nucl_area` | Cellpose min nucleus area in pixels. Minimal area of ROIs. If null defaults to π (1/6th cp_nucl_size)^2 for 2d or 4/3π (1/6th cp_nucl_size)^3 for 3d. | `string` |  |  |  |
| `cp_model` | Cellpose model. util built on cyto2 model and if possible will use nucleus channel. Other models should work in principle. | `string` | cyto2 |  |  |
| `cp_dont_use_nucl_for_declump` | Fit nucleus mask but do not use nuclei to declump objects in mask creation. | `boolean` |  |  |  |
| `cp_cell_power` | Raise cell images to power then scale. Deprecated/experimental. Soft thresholding if cellpose fits masks too strongly to noise. Option if tweaking cellprob threshold or post-processing insufficient. | `string` |  |  | True |
| `cp_nucl_power` | Raise nucleus images to power then scale. Deprecated/experimental. Soft thresholding if cellpose fits masks too strongly to noise. Option if tweaking cellprob threshold or post-processing insufficient. | `string` |  |  | True |
| `cp_downsample` | Downsample images in YX prior to running cellpose. Scales diameter, anisotropy, min cell area, min nucl area to match. Improves speed at cost of mask resolution. Masks scaled up by nearest neighbour interpolation. Recommended integer values producing whole number in YX, 2 is good. | `number` |  |  |  |
| `cp_dont_post_process` | To not post process masks in 3d mode set true. In 3d some post-processing on nuclei applied: local otsu thresholding, hole closing, mask multiplied with cellpose masks to contain region with signal. Prevents incorrect masks due to z-bleedover from PSF residual signal. Deprecated/experimental. | `boolean` | True |  | True |
| `cp_cell_flow_thresh` | Cellpose flow threshold for cells. See cellpose docs. | `number` | 0.4 |  |  |
| `cp_nucl_flow_thresh` | Cellpose flow threshold for nuclei. See cellpose docs. | `number` | 0.4 |  |  |
| `cp_cell_prob_threshold` | Cellpose cellprob threshold for cells. Between -6 and 6. Higher is tighter masks, lower looser masks. See cellpose docs for details. | `integer` | 0 |  |  |
| `cp_nucl_prob_threshold` | Cellpose cellprob threshold for nuclei. Between -6 and 6. Higher is tighter masks, lower looser masks. See cellpose docs for details. | `integer` | 0 |  |  |
| `rg_plot` | Run only if registration manifest provided. Make before/after images of registration results | `boolean` | True |  |  |

## Scaling



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `rn_manualscale` | Path to scaling_factors.txt with line formatted: <plate>_ch<channel>=<scale> <plate>_ch<channel>=<scale> <plate>_ch<channel>=<scale> Space separated Channels zero indexed <plate> is the reference plate <channel> after merging cycles <scale> factor by which channel in plate is divided | `string` |  |  |  |
| `rn_scale_slope` | Path to scaling_<slope/bias>.txt with one line formatted: <plate>_ch<channel>=<slope/bias> <plate>_ch<channel>=<slope/bias> <plate>_ch<channel>=<slope/bias> Sets shape of sigmoid curve used to weigh scaling factors differently in intensity ranges. Background pixels remain unscaled, smooth transition scales foreground pixels. Pipeline supports only pre-calculated slope and bias, cannot auto-estimate. Must supply both slope and bias if used. If null, scaling_factors applied uniformly. Works with --rn_manualscale and --rn_autoscale. | `string` |  |  |  |
| `rn_scale_bias` | See rn_scale_slope (must supply with rn_scale_slope if supplied) | `string` |  |  |  |
| `rn_autoscale` | Automatically determine scaling factors based on all images in manifest (excluding blacklist) Overrides rn_manualscale Runs after dc_clip_max applied Waits till deconvolution queue empty No feature extraction jobs submitted until autoscale completes Interacts with rn_controllist so dynamic range is optimal after plate offset factors. If rn_controllist is provided, rn_autoscale is on by default. | `boolean` |  |  |  |
| `rn_autoscale_q1` | Controls value to scale to within an image (quantiles of pixels). Valid options: 'q0', 'q0.1', 'q1', 'q5', 'q25', 'q50', 'q75', 'q90', 'q99', 'q99.9', 'q99.99', 'q99.999', 'q99.9999', 'q99.99999', 'q100', 'mean' Not all options practically sensible, recommend not below q99. | `string` | q99.9999 |  |  |
| `rn_autoscale_q2` | Controls which quantile is taken across all images for chosen q1 Valid options: 0 - 1 | `integer` | 95 |  |  |
| `rn_dummy_mode` | Generate plate offsets in dummy mode (returns all 1’s for equal scaling, just for testing) | `boolean` |  |  | True |
| `rn_threshold` | When calculating plate offsets, threshold channel images prior to calculating mean object intensities. Background regions ignored. Otsu threshold on whole image used. | `boolean` |  |  | True |

## Cellprofiler



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `cpr_label` |  | `string` | normal |  |  |
| `cpr_conda_env` |  | `string` | /software/teamtrynka/installs/cellprofiler |  |  |
| `cpr_plugins` |  | `string` | $projectDir/bin/cellprofiler/plugins |  |  |
| `cpr_run` |  | `boolean` | True |  |  |
| `cpr_pipeline_2d` |  | `string` |  |  |  |
| `cpr_pipeline_3d` |  | `string` |  |  |  |
| `cpr_no_zip` |  | `boolean` |  |  |  |
| `cpr_container` | Container for the cellprofiler enviroment | `string` |  |  |  |

## SubCell



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `sc_label` | Resource label for subcell | `string` | gpu_short |  | True |
| `sc_gpu` | Should GPU be used, make sure to set sc_label to gpu_<x> | `boolean` | True |  | True |
| `sc_dont_mask` | Should the cell crops not be masked by the cell object prior to embedding | `boolean` |  |  | True |
| `sc_conda_env` | Subcell conda env path | `string` |  |  | True |
| `sc_dl_conda_env` | Subcell model download conda env path | `string` |  |  | True |
| `sc_model_ref_channels` | Subcell model reference channels string like "rybg" | `string` | rybg |  | True |
| `sc_model` | Subcell model string like "mae_contrast_supcon_model" | `string` | mae_contrast_supcon_model |  | True |
| `sc_channels` | Channels to find localization for should be string '<name>=<channel> <name>=<channel>' | `string` |  |  | True |
| `sc_nucl` | Nucleus channel | `integer` |  |  | True |
| `sc_tub` | Tubulin channel | `integer` |  |  | True |
| `sc_er` | Endoplasmic reticulum channel | `integer` |  |  | True |
| `sc_scale` | Scale factor to get pixel size to 80nm. In Phenix 1px=149nm at 40x, so sc_scale=149/80. TODO set properly using physical pixel size attribute if available | `number` | 1.8625 |  | True |

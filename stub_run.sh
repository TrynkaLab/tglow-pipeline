#!/usr/bin/env bash
# Runs the pipeline with `-stub-run` against fabricated, minimal inputs so the
# workflow wiring can be exercised without real images, conda envs or GPUs.
# Every process below already ships a `stub:` block (see processes/*.nf) which
# just touches placeholder outputs instead of running the real script.
#
# Usage: ./stub_run.sh [stage|run_pipeline]   (default: run_pipeline)

prev_wd="$(pwd)"

set -euo pipefail

WORKFLOW="${1:-run_pipeline}"
case "$WORKFLOW" in
    stage|run_pipeline) ;;
    *) echo "Usage: $0 [stage|run_pipeline]"; exit 1 ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUB_DIR="$SCRIPT_DIR/.stub_run"

rm -rf "$STUB_DIR"
mkdir -p "$STUB_DIR/results/images" "$STUB_DIR/workdir"

# --- Fabricated inputs -----------------------------------------------------

# Dummy PSF for deconvolution (dc_psfs is indexed by 0-based channel index, so a
# single dc_channels=1 entry only ever needs psf index 0 - one file is enough).
touch "$STUB_DIR/dummy_psf.tif"

# Core manifest: one plate, two channels, nucleus=1 & cell=2 (so cellpose has
# something to key off), everything else disabled ("none").
cat > "$STUB_DIR/manifest.tsv" <<EOF
plate	index_xml	channels	ff_channels	cp_nucl_channel	cp_cell_channel	dc_channels	dc_psfs	mask_channels
PLATE1	dummy_index.xml	1,2	1	1	2	1	$STUB_DIR/dummy_psf.tif	none
PLATE2	dummy_index.xml	1,2	1	1	2	1	$STUB_DIR/dummy_psf.tif	none
EOF
touch "$STUB_DIR/dummy_index.xml"

# sc_autoscale needs a channel map and a control list to derive scaling factors
# from - contents don't matter under -stub-run (the real script never executes),
# they just need to exist so Nextflow can stage them.
touch "$STUB_DIR/dummy_channel_map.tsv"
printf 'PLATE1\tA1\t1,2\tcontrolA\n' > "$STUB_DIR/dummy_control_list.tsv"

# Dummy image + per-plate well manifest, so `setup`'s well_channel is
# actually populated (index_imagedir's stub block just touches an empty
# manifest.tsv regardless of what's on disk — pre-existing manifest.tsv is
# needed so its storeDir cache-hit is used instead).
mkdir -p "$STUB_DIR/results/images/PLATE1/A/1"
touch "$STUB_DIR/results/images/PLATE1/A/1/1.ome.tiff"
printf 'A1\tA\t1\tPLATE1\tnone\n' > "$STUB_DIR/results/images/PLATE1/manifest.tsv"

# PLATE2 gets the same well as PLATE1 - registration merges PLATE2 (query) onto
# PLATE1 (reference) below, and with dc_run on, the ref/query grouping in
# run_pipeline.nf waits on a decon output from both plates for the same well,
# so PLATE2 needs a matching A1 well to avoid hanging.
mkdir -p "$STUB_DIR/results/images/PLATE2/A/1"
touch "$STUB_DIR/results/images/PLATE2/A/1/1.ome.tiff"
printf 'A1\tA\t1\tPLATE2\tnone\n' > "$STUB_DIR/results/images/PLATE2/manifest.tsv"

# Registration manifest: merge PLATE2 (query) onto PLATE1 (reference), both on channel 1.
printf 'reference_plate\treference_channel\tquery_plates\tquery_channels\nPLATE1\t1\tPLATE2\t1\n' > "$STUB_DIR/dummy_registration_manifest.tsv"

# CellProfiler pipeline: only needs to exist, cellprofiler runs its stub block.
touch "$STUB_DIR/dummy.cppipe"

# Local override config: the repo config enables conda globally and points at
# environments under /software/teamtrynka which won't exist outside that HPC.
cat > "$STUB_DIR/stub.config" <<'EOF'
conda.enabled = false
docker.enabled = false
singularity.enabled = false

// conf/processes.config sizes labels (e.g. gpu_normal_plus = 24GB) for the
// real HPC cluster. A stub run only touches placeholder files, so cap every
// label down to something any laptop can satisfy.
process {
    withLabel: '.*' {
        cpus = 1
        memory = 1.GB
        time = 10.m
    }
}
EOF

# --- Run --------------------------------------------------------------------
cd $STUB_DIR

nextflow run "$SCRIPT_DIR/main.nf" \
    -profile local \
    -w "$STUB_DIR/workdir" \
    -c "$STUB_DIR/stub.config" \
    -stub-run \
    --workflow "$WORKFLOW" \
    --rn_skip_checks true \
    --rn_manifest "$STUB_DIR/manifest.tsv" \
    --rn_publish_dir "$STUB_DIR/results" \
    --rn_image_dir "$STUB_DIR/results/images" \
    --rn_decon_dir "$STUB_DIR/results/decon" \
    --cpr_pipeline_3d "$STUB_DIR/dummy.cppipe" \
    --dc_run true \
    --sc_autoscale true \
    --sc_channel_map "$STUB_DIR/dummy_channel_map.tsv" \
    --sc_control_list "$STUB_DIR/dummy_control_list.tsv" \
    --rn_manifest_registration "$STUB_DIR/dummy_registration_manifest.tsv" \
    -with-report "$STUB_DIR/report.html" \
    -with-trace "$STUB_DIR/trace.txt"

echo "Stub run complete. Outputs under: $STUB_DIR"

cd ${prev_wd}
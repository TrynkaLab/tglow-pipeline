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

# Core manifest: one plate, two channels, nucleus=1 & cell=2 (so cellpose has
# something to key off), everything else disabled ("none").
cat > "$STUB_DIR/manifest.tsv" <<'EOF'
plate	index_xml	channels	bp_channels	cp_nucl_channel	cp_cell_channel	dc_channels	dc_psfs	mask_channels
PLATE1	dummy_index.xml	1,2	none	1	2	none	none	none
EOF
touch "$STUB_DIR/dummy_index.xml"

# CellProfiler pipeline: only needs to exist, cellprofiler runs its stub block.
touch "$STUB_DIR/dummy.cppipe"

# Local override config: the repo config enables conda globally and points at
# environments under /software/teamtrynka which won't exist outside that HPC.
cat > "$STUB_DIR/stub.config" <<'EOF'
conda.enabled = false
docker.enabled = false
singularity.enabled = false
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
    --cpr_pipeline_3d "$STUB_DIR/dummy.cppipe" \
    -with-report "$STUB_DIR/report.html" \
    -with-trace "$STUB_DIR/trace.txt"

echo "Stub run complete. Outputs under: $STUB_DIR"

cd ${prev_wd}
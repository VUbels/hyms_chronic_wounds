#!/usr/bin/env bash
#
# run_proseg.sh — Batch proseg pipeline for Xenium spatial transcriptomics
#
# Recursively finds all transcripts.parquet files under a target directory,
# runs proseg on each, and produces outputs compatible with:
#   - Python (spatialdata zarr)
#   - R / Seurat v5 (mtx + csv.gz metadata)
#   - Xenium Explorer (proseg-to-baysor conversion)
#
# Usage:
#   ./run_proseg.sh [OPTIONS]
#
# Options:
#   -i, --input DIR        Root directory containing Xenium region folders
#                          (default: current working directory)
#   -o, --output DIR       Output root directory
#                          (default: current working directory)
#   -t, --threads N        Number of threads per proseg run (default: all available)
#   -p, --parallel N       Run N regions in parallel (default: 1 = sequential)
#   -n, --dry-run          Print what would be run without executing
#   --proseg PATH          Path to proseg binary (default: searches PATH,
#                          then ./proseg/target/release/proseg)
#   --no-diffusion         Disable transcript diffusion modeling
#   --voxel-size S         Voxel size in microns (proseg default: 0.5)
#   --ncomponents N        Number of expression mixture components (proseg default: 10)
#   -h, --help             Show this help message
#
# Requirements:
#   - proseg (built via cargo install proseg or from source)
#   - proseg-to-baysor (installed alongside proseg)
#
# Output structure per region:
#   <output_dir>/<region_name>/
#     ├── proseg-output.zarr              # spatialdata format
#     ├── counts.mtx.gz                   # sparse count matrix
#     ├── cell-metadata.csv.gz            # cell centroids, volume, etc.
#     ├── gene-metadata.csv.gz            # per-gene summary stats
#     ├── transcript-metadata.csv.gz
#     ├── cell-polygons.geojson.gz
#     ├── proseg.log                      # stdout + stderr log
#     └── xenium_explorer/
#         ├── transcript-metadata.csv     # baysor-compatible transcript table
#         ├── cell-polygons.geojson       # baysor-compatible polygons
#         └── import_to_xenium_ranger.sh
#

set -euo pipefail

# ============================================================================
# Defaults
# ============================================================================
INPUT_DIR="$(pwd)"
OUTPUT_DIR="$(pwd)"
THREADS=""
PARALLEL_JOBS=1
DRY_RUN=false
PROSEG_BIN=""
EXTRA_PROSEG_ARGS=()

# ============================================================================
# Color output helpers
# ============================================================================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info()  { echo -e "${BLUE}[INFO]${NC}  $(date '+%Y-%m-%d %H:%M:%S') $*"; }
log_ok()    { echo -e "${GREEN}[OK]${NC}    $(date '+%Y-%m-%d %H:%M:%S') $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC}  $(date '+%Y-%m-%d %H:%M:%S') $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') $*"; }

# ============================================================================
# Argument parsing
# ============================================================================
show_help() {
    sed -n '/^# Usage:/,/^# Requirements:/{ /^# Requirements:/d; s/^# \?//p }' "$0"
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)       INPUT_DIR="$2";                   shift 2 ;;
        -o|--output)      OUTPUT_DIR="$2";                  shift 2 ;;
        -t|--threads)     THREADS="$2";                     shift 2 ;;
        -p|--parallel)    PARALLEL_JOBS="$2";               shift 2 ;;
        -n|--dry-run)     DRY_RUN=true;                     shift   ;;
        --proseg)         PROSEG_BIN="$2";                  shift 2 ;;
        --no-diffusion)   EXTRA_PROSEG_ARGS+=("--no-diffusion"); shift ;;
        --voxel-size)     EXTRA_PROSEG_ARGS+=("--voxel-size" "$2"); shift 2 ;;
        --ncomponents)    EXTRA_PROSEG_ARGS+=("--ncomponents" "$2"); shift 2 ;;
        -h|--help)        show_help ;;
        *)                log_error "Unknown argument: $1"; exit 1 ;;
    esac
done

# ============================================================================
# Resolve proseg binary
# ============================================================================
resolve_proseg() {
    if [[ -n "$PROSEG_BIN" ]]; then
        if [[ -x "$PROSEG_BIN" ]]; then
            echo "$PROSEG_BIN"
            return
        else
            log_error "Specified proseg binary not found or not executable: $PROSEG_BIN"
            exit 1
        fi
    fi

    if command -v proseg &>/dev/null; then
        command -v proseg
        return
    fi

    local local_bin="./proseg/target/release/proseg"
    if [[ -x "$local_bin" ]]; then
        echo "$local_bin"
        return
    fi

    log_error "proseg binary not found. Install with 'cargo install proseg' or build from source."
    exit 1
}

PROSEG_BIN="$(resolve_proseg)"
log_info "Using proseg binary: $PROSEG_BIN"

PROSEG_TO_BAYSOR=""
if command -v proseg-to-baysor &>/dev/null; then
    PROSEG_TO_BAYSOR="$(command -v proseg-to-baysor)"
elif [[ -x "./proseg/target/release/proseg-to-baysor" ]]; then
    PROSEG_TO_BAYSOR="./proseg/target/release/proseg-to-baysor"
else
    log_warn "proseg-to-baysor not found. Xenium Explorer output will be skipped."
fi

# ============================================================================
# Validate input directory
# ============================================================================
if [[ ! -d "$INPUT_DIR" ]]; then
    log_error "Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# ============================================================================
# Discover all transcripts.parquet files
# ============================================================================
log_info "Scanning for transcripts.parquet files in: $INPUT_DIR"

mapfile -t PARQUET_FILES < <(find "$INPUT_DIR" -name "transcripts.parquet" -type f | sort)

if [[ ${#PARQUET_FILES[@]} -eq 0 ]]; then
    log_error "No transcripts.parquet files found under $INPUT_DIR"
    exit 1
fi

log_info "Found ${#PARQUET_FILES[@]} region(s):"
for f in "${PARQUET_FILES[@]}"; do
    region_dir="$(dirname "$f")"
    region_name="$(basename "$region_dir")"
    echo "         - $region_name"
done

# ============================================================================
# Process a single region
# ============================================================================
process_region() {
    local parquet_path="$1"
    local region_dir
    local region_name
    local out_base

    region_dir="$(dirname "$parquet_path")"
    region_name="$(basename "$region_dir")"

    local safe_name
    safe_name="$(echo "$region_name" | tr ' ' '_')"

    out_base="${OUTPUT_DIR}/${safe_name}"

    log_info "========================================"
    log_info "Processing: $region_name"
    log_info "  Input:  $parquet_path"
    log_info "  Output: $out_base"
    log_info "========================================"

    if $DRY_RUN; then
        log_info "[DRY RUN] Would process $region_name"
        return 0
    fi

    # Only xenium_explorer gets its own subdirectory
    mkdir -p "$out_base"
    mkdir -p "$out_base/xenium_explorer"

    # ------------------------------------------------------------------
    # Step 1: Run proseg — all outputs flat in out_base
    # ------------------------------------------------------------------
    local proseg_cmd=(
        "$PROSEG_BIN"
        --xenium
        --overwrite
        --output-spatialdata "$out_base/proseg-output.zarr"
        --output-counts "$out_base/counts.mtx.gz"
        --output-cell-metadata "$out_base/cell-metadata.csv.gz"
        --output-gene-metadata "$out_base/gene-metadata.csv.gz"
        --output-transcript-metadata "$out_base/transcript-metadata.csv.gz"
        --output-cell-polygons "$out_base/cell-polygons.geojson.gz"
    )

    if [[ -n "$THREADS" ]]; then
        proseg_cmd+=(--nthreads "$THREADS")
    fi

    if [[ ${#EXTRA_PROSEG_ARGS[@]} -gt 0 ]]; then
        proseg_cmd+=("${EXTRA_PROSEG_ARGS[@]}")
    fi

    proseg_cmd+=("$parquet_path")

    log_info "Running proseg..."
    log_info "  Command: ${proseg_cmd[*]}"

    local start_time
    start_time=$(date +%s)

    if "${proseg_cmd[@]}" 2>&1 | tee "$out_base/proseg.log"; then
        local end_time
        end_time=$(date +%s)
        local elapsed=$(( end_time - start_time ))
        log_ok "proseg completed for $region_name in ${elapsed}s"
    else
        local end_time
        end_time=$(date +%s)
        local elapsed=$(( end_time - start_time ))
        log_error "proseg FAILED for $region_name after ${elapsed}s"
        log_error "  Check log: $out_base/proseg.log"
        return 1
    fi

    # ------------------------------------------------------------------
    # Step 2: Convert for Xenium Explorer (proseg-to-baysor)
    # ------------------------------------------------------------------
    if [[ -n "$PROSEG_TO_BAYSOR" ]]; then
        log_info "Converting to Xenium Explorer format..."

        local baysor_cmd=(
            "$PROSEG_TO_BAYSOR"
            "$out_base/proseg-output.zarr"
            --output-transcript-metadata "$out_base/xenium_explorer/transcript-metadata.csv"
            --output-cell-polygons "$out_base/xenium_explorer/cell-polygons.geojson"
        )

        if "${baysor_cmd[@]}" 2>&1 | tee -a "$out_base/proseg.log"; then
            log_ok "Xenium Explorer conversion complete for $region_name"

            cat > "$out_base/xenium_explorer/import_to_xenium_ranger.sh" <<XENIUM_EOF
#!/usr/bin/env bash
# Run this to import proseg segmentation into Xenium Explorer.
# Requires xeniumranger to be installed.
# Edit the --xenium-bundle path to point to the ORIGINAL Xenium output bundle.

xeniumranger import-segmentation \\
    --id="${safe_name}_proseg" \\
    --xenium-bundle="${region_dir}" \\
    --transcript-assignment="$out_base/xenium_explorer/transcript-metadata.csv" \\
    --viz-polygons="$out_base/xenium_explorer/cell-polygons.geojson" \\
    --units=microns \\
    --localcores=32 \\
    --localmem=128
XENIUM_EOF
            chmod +x "$out_base/xenium_explorer/import_to_xenium_ranger.sh"
            log_info "  Xenium Ranger import script written to:"
            log_info "  $out_base/xenium_explorer/import_to_xenium_ranger.sh"
        else
            log_warn "proseg-to-baysor conversion failed for $region_name"
        fi
    fi

    log_ok "Region complete: $region_name"
    return 0
}

# ============================================================================
# Execute across all regions
# ============================================================================
TOTAL=${#PARQUET_FILES[@]}
SUCCESS=0
FAILED=0
FAILED_REGIONS=()

log_info "Starting proseg pipeline"
log_info "  Mode: $([ "$PARALLEL_JOBS" -gt 1 ] && echo "parallel ($PARALLEL_JOBS jobs)" || echo "sequential")"
log_info "  Regions: $TOTAL"
log_info "  Output root: $OUTPUT_DIR/"
echo ""

if [[ "$PARALLEL_JOBS" -le 1 ]]; then
    for i in "${!PARQUET_FILES[@]}"; do
        log_info "Region $(( i + 1 )) / $TOTAL"
        if process_region "${PARQUET_FILES[$i]}"; then
            (( SUCCESS++ )) || true
        else
            (( FAILED++ )) || true
            FAILED_REGIONS+=("$(basename "$(dirname "${PARQUET_FILES[$i]}")")")
        fi
        echo ""
    done
else
    export INPUT_DIR OUTPUT_DIR THREADS DRY_RUN PROSEG_BIN PROSEG_TO_BAYSOR
    export -a EXTRA_PROSEG_ARGS
    export -f process_region log_info log_ok log_warn log_error

    if command -v parallel &>/dev/null; then
        log_info "Using GNU parallel with $PARALLEL_JOBS jobs"
        printf '%s\n' "${PARQUET_FILES[@]}" | \
            parallel -j "$PARALLEL_JOBS" --halt soon,fail=1 \
            process_region {}
    else
        log_warn "GNU parallel not found. Falling back to background jobs."
        log_warn "  Install with: sudo apt install parallel"

        local_running=0
        for parquet_path in "${PARQUET_FILES[@]}"; do
            process_region "$parquet_path" &
            (( local_running++ )) || true
            if (( local_running >= PARALLEL_JOBS )); then
                wait -n
                (( local_running-- )) || true
            fi
        done
        wait
    fi

    for f in "${PARQUET_FILES[@]}"; do
        region_dir="$(dirname "$f")"
        region_name="$(basename "$region_dir")"
        safe_name="$(echo "$region_name" | tr ' ' '_')"
        if [[ -d "${OUTPUT_DIR}/${safe_name}/proseg-output.zarr" ]]; then
            (( SUCCESS++ )) || true
        else
            (( FAILED++ )) || true
            FAILED_REGIONS+=("$region_name")
        fi
    done
fi

# ============================================================================
# Summary
# ============================================================================
echo ""
log_info "========================================"
log_info "Pipeline complete"
log_info "  Total:     $TOTAL"
log_ok   "  Succeeded: $SUCCESS"
if [[ $FAILED -gt 0 ]]; then
    log_error "  Failed:    $FAILED"
    for r in "${FAILED_REGIONS[@]}"; do
        log_error "    - $r"
    done
fi
log_info "  Output at: $OUTPUT_DIR/"
log_info "========================================"
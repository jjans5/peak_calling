#!/bin/bash

# ============================================
# CONFIGURATION - MODIFY THESE PATHS AS NEEDED
# ============================================

# Base directory
BASE_DIR="/cluster/project/treutlein/USERS/jjans/analysis/adult_intestine/dl_models/enterocytes_v0"

# Input: subsampled fragments directory
INPUT_BASE="${BASE_DIR}/subsampled_fragments"

# Output: bigwig directory
OUTPUT_BASE="${BASE_DIR}/subsampled_bigwigs"

# BigWig type: "coverage" or "cut-sites"
BIGWIG_TYPE="cut-sites"  # Change to "cut-sites" for cut-site bigwigs

# Number of parallel jobs (reduce if OOM errors occur)
PARALLEL_JOBS=4

# Skip existing valid bigWigs (set to false to overwrite all)
SKIP_EXISTING=true

# Chromosome sizes per species
declare -A CHROM_SIZES
CHROM_SIZES[human]="/cluster/home/jjanssens/jjans/analysis/cerebellum/genomes_new/homo_sapiens/hg38.chrom.sizes"
CHROM_SIZES[gorilla]="/cluster/work/treutlein/jjans/data/intestine/nhp_atlas/genomes/reference_/gorGor4/fasta/gorGor4.chrom.sizes"
CHROM_SIZES[chimpanzee]="/cluster/work/treutlein/jjans/data/intestine/nhp_atlas/genomes/reference_/panTro3/fasta/panTro3.chrom.sizes"
CHROM_SIZES[bonobo]="/cluster/work/treutlein/jjans/data/intestine/nhp_atlas/genomes/reference_/panPan1/fasta/panPan1.chrom.sizes"
CHROM_SIZES[macaque]="/cluster/work/treutlein/jjans/data/intestine/nhp_atlas/genomes/reference_/Mmul10/fasta/Mmul10.chrom.sizes"
CHROM_SIZES[marmoset]="/cluster/work/treutlein/jjans/data/intestine/nhp_atlas/genomes/reference_/calJac1_mito/fasta/calJac1_mito.chrom.sizes"

# ============================================
# END CONFIGURATION
# ============================================

# Activate environment
mamba activate scatac_fragment_tools

# Determine bigwig suffix and flags
if [[ "$BIGWIG_TYPE" == "cut-sites" ]]; then
    BW_SUFFIX=".cutsite.bw"
    BW_FLAG="--cut-sites"
    echo "🔪 Mode: Cut-sites bigWig"
else
    BW_SUFFIX=".cov.bw"
    BW_FLAG="-n"
    echo "📊 Mode: Coverage bigWig (normalized)"
fi

echo "🔍 Scanning for subsampled fragments..."
echo "   Input: ${INPUT_BASE}"
echo "   Output: ${OUTPUT_BASE}"
echo ""

# Find all fragment files
fragment_files=$(find "${INPUT_BASE}" -name "*.fragments_sub.tsv.gz" -type f)

if [[ -z "$fragment_files" ]]; then
    echo "❌ No fragment files found in ${INPUT_BASE}"
    exit 1
fi

file_count=$(echo "$fragment_files" | wc -l)
echo "📝 Found ${file_count} fragment files to convert"

# Validate chromosome sizes files
echo "🔍 Validating chromosome sizes files..."
missing_chrom=0
for species in "${!CHROM_SIZES[@]}"; do
    if [[ ! -f "${CHROM_SIZES[$species]}" ]]; then
        echo "   ❌ Missing: $species -> ${CHROM_SIZES[$species]}"
        missing_chrom=1
    fi
done

if [[ $missing_chrom -eq 1 ]]; then
    echo "❌ Some chromosome sizes files are missing"
    exit 1
fi
echo "   ✅ All chromosome sizes files found"
echo ""

# Create commands for parallel execution
commands_file=$(mktemp)
trap "rm -f $commands_file" EXIT

while IFS= read -r input_frag; do
    # Extract relative path and filename
    relative_path="${input_frag#${INPUT_BASE}/}"
    dir_name=$(dirname "$relative_path")
    file_name=$(basename "$relative_path" .fragments_sub.tsv.gz)
    
    # Extract species from filename (e.g., bonobo_enterocytes -> bonobo)
    species=$(echo "$file_name" | cut -d'_' -f1)
    
    # Get chromosome sizes for this species
    chrom_file="${CHROM_SIZES[$species]}"
    
    if [[ -z "$chrom_file" || ! -f "$chrom_file" ]]; then
        echo "⚠️  Unknown species or missing chrom sizes for: $species (file: $file_name) — skipping"
        continue
    fi
    
    # Create output path
    output_dir="${OUTPUT_BASE}/${dir_name}"
    output_bw="${output_dir}/${file_name}${BW_SUFFIX}"
    log_file="${output_dir}/${file_name}.bigwig.log"
    
    # Create output directory
    mkdir -p "$output_dir"
    
    # Skip if valid bigwig already exists (only if SKIP_EXISTING=true)
    if [[ "$SKIP_EXISTING" == true && -s "$output_bw" ]]; then
        if command -v bigWigInfo &>/dev/null; then
            if bigWigInfo "$output_bw" >/dev/null 2>&1; then
                continue
            fi
        else
            continue
        fi
    fi
    
    # Create command with logging
    cmd="scatac_fragment_tools bigwig -i \"$input_frag\" -c \"$chrom_file\" -o \"$output_bw\" $BW_FLAG > \"$log_file\" 2>&1 && echo \"✅ $file_name ($species)\" || echo \"❌ $file_name ($species)\""
    
    echo "$cmd" >> "$commands_file"
done <<< "$fragment_files"

cmd_count=$(wc -l < "$commands_file")

if [[ $cmd_count -eq 0 ]]; then
    if [[ "$SKIP_EXISTING" == true ]]; then
        echo "✅ All bigWig files already exist and are valid"
    else
        echo "✅ No files to process"
    fi
    exit 0
fi

echo "🚀 Converting ${cmd_count} fragments to bigWig (${PARALLEL_JOBS} parallel jobs)..."
echo "   Logs will be saved alongside bigWig files"
echo ""

# Run conversions in parallel
parallel --no-notice -j "$PARALLEL_JOBS" < "$commands_file"

echo ""
echo "🎉 BigWig conversion complete!"
echo "   Output: ${OUTPUT_BASE}"
echo "   Type: ${BIGWIG_TYPE}"
echo "   Check individual .bigwig.log files for details"
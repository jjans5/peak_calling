#!/bin/bash

# batch_subsample_fragments.sh
# Subsample all fragment files while preserving directory structure
#
# Usage:
#   batch_subsample_fragments.sh --input-dir <dir> --output-dir <dir> [OPTIONS]
#
# Options:
#   --input-dir, -i     Input directory (e.g., fragment_files/)
#   --output-dir, -o    Output directory (e.g., subsampled_fragment_files/)
#   --ncells, -c        Number of cells to sample (optional)
#   --nfrags, -f        Total number of fragments to sample (optional)
#   --frags-per-cell, -p  Number of fragments per cell (optional)
#   --seed, -s          Random seed for reproducibility (default: 42)
#   --test              Test mode: create 1M read test file first
#   --help, -h          Show this help message

set -euo pipefail

# === Default values ===
input_dir=""
output_dir=""
ncells=""
nfrags=""
frags_per_cell=""
seed=42
test_mode=false

# === Parse command line arguments ===
while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input-dir)
      input_dir="$2"
      shift 2
      ;;
    -o|--output-dir)
      output_dir="$2"
      shift 2
      ;;
    -c|--ncells)
      ncells="$2"
      shift 2
      ;;
    -f|--nfrags)
      nfrags="$2"
      shift 2
      ;;
    -p|--frags-per-cell)
      frags_per_cell="$2"
      shift 2
      ;;
    -s|--seed)
      seed="$2"
      shift 2
      ;;
    --test)
      test_mode=true
      shift
      ;;
    -h|--help)
      grep '^#' "$0" | grep -v '#!/bin/bash' | sed 's/^# //'
      exit 0
      ;;
    *)
      echo "❌ ERROR: Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

# === Validate required arguments ===
if [[ -z "$input_dir" ]]; then
  echo "❌ ERROR: --input-dir is required" >&2
  exit 1
fi

if [[ -z "$output_dir" ]]; then
  echo "❌ ERROR: --output-dir is required" >&2
  exit 1
fi

if [[ ! -d "$input_dir" ]]; then
  echo "❌ ERROR: Input directory not found: $input_dir" >&2
  exit 1
fi

# Check that at least one subsampling option is specified
if [[ -z "$ncells" && -z "$nfrags" && -z "$frags_per_cell" ]] && [[ "$test_mode" == false ]]; then
  echo "❌ ERROR: At least one of --ncells, --nfrags, or --frags-per-cell must be specified" >&2
  exit 1
fi

# Get absolute path to subsample_fragments.sh
script_dir="$(cd "$(dirname "$0")" && pwd)"
subsample_script="$script_dir/subsample_fragments.sh"

if [[ ! -f "$subsample_script" ]]; then
  echo "❌ ERROR: subsample_fragments.sh not found in $script_dir" >&2
  exit 1
fi

# === TEST MODE: Create small test file ===
if [[ "$test_mode" == true ]]; then
  echo "🧪 TEST MODE: Creating 1M read test file..."
  
  # Find first fragment file
  test_input=$(find "$input_dir" -name "*.fragments.tsv.gz" -type f | head -1)
  
  if [[ -z "$test_input" ]]; then
    echo "❌ ERROR: No fragment files found in $input_dir" >&2
    exit 1
  fi
  
  echo "   Using: $test_input"
  
  # Create test file with 1M reads
  test_output="test_fragments_1M.tsv.gz"
  echo "   Creating $test_output..."
  
  gunzip -c "$test_input" | head -1000000 | gzip > "$test_output"
  
  test_size=$(gunzip -c "$test_output" | wc -l)
  echo "✅ Test file created: $test_output ($test_size fragments)"
  echo ""
  echo "   Test it with:"
  echo "   bash subsample_fragments.sh -i $test_output -o test_subsampled.tsv.gz --nfrags 100000"
  echo ""
  exit 0
fi

# === Build subsample command arguments ===
subsample_args=""
[[ -n "$ncells" ]] && subsample_args="$subsample_args --ncells $ncells"
[[ -n "$nfrags" ]] && subsample_args="$subsample_args --nfrags $nfrags"
[[ -n "$frags_per_cell" ]] && subsample_args="$subsample_args --frags-per-cell $frags_per_cell"
subsample_args="$subsample_args --seed $seed"

echo "🔬 Batch subsampling fragment files..."
echo "   Input directory: $input_dir"
echo "   Output directory: $output_dir"
[[ -n "$ncells" ]] && echo "   Cells to sample: $ncells"
[[ -n "$nfrags" ]] && echo "   Total fragments to sample: $nfrags"
[[ -n "$frags_per_cell" ]] && echo "   Fragments per cell: $frags_per_cell"
echo "   Random seed: $seed"
echo ""

# === Create output directory ===
mkdir -p "$output_dir"

# === Find all fragment files ===
fragment_files=($(find "$input_dir" -name "*.fragments.tsv.gz" -type f))
total_files=${#fragment_files[@]}

if [[ $total_files -eq 0 ]]; then
  echo "❌ ERROR: No fragment files found in $input_dir" >&2
  exit 1
fi

echo "📁 Found $total_files fragment files"
echo ""

# === Process each file ===
processed=0
failed=0

for input_file in "${fragment_files[@]}"; do
  ((processed++))
  
  # Get relative path from input_dir
  rel_path="${input_file#$input_dir/}"
  
  # Create output path preserving directory structure
  output_file="$output_dir/$rel_path"
  output_subdir="$(dirname "$output_file")"
  
  echo "[$processed/$total_files] Processing: $rel_path"
  
  # Create subdirectory if needed
  mkdir -p "$output_subdir"
  
  # Run subsampling
  if bash "$subsample_script" --input "$input_file" --output "$output_file" $subsample_args; then
    echo "   ✅ Success"
  else
    echo "   ❌ Failed"
    ((failed++))
  fi
  
  echo ""
done

# === Summary ===
echo "================================"
echo "✅ Batch processing complete!"
echo "   Total files: $total_files"
echo "   Successful: $((total_files - failed))"
echo "   Failed: $failed"
echo "   Output directory: $output_dir"

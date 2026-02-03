#!/bin/bash
# =============================================================================
# Subsample Fragments
# =============================================================================
# Subsample a fragments file by cells and/or fragments for testing or analysis.
#
# USAGE:
#   bash subsample_fragments.sh --input <file> --output <file> [OPTIONS]
#
# OPTIONS:
#   --input, -i         Input fragments file (.tsv.gz or .tsv)
#   --output, -o        Output fragments file (.tsv.gz or .tsv)
#   --ncells, -c        Number of cells to sample
#   --nfrags, -f        Total number of fragments to sample
#   --frags-per-cell, -p  Number of fragments per cell
#   --min-frags, -m     Minimum fragments per cell to include
#   --filter-chroms     Keep only main chromosomes (default)
#   --no-filter-chroms  Keep all chromosomes including scaffolds
#   --seed, -s          Random seed (default: 42)
#   --help, -h          Show this help message
#
# EXAMPLES:
#   # Subsample to 1000 cells
#   bash subsample_fragments.sh -i input.tsv.gz -o output.tsv.gz -c 1000
#
#   # Subsample to 1M fragments total
#   bash subsample_fragments.sh -i input.tsv.gz -o output.tsv.gz -f 1000000
#
#   # Subsample 500 cells with at least 1000 fragments each
#   bash subsample_fragments.sh -i input.tsv.gz -o output.tsv.gz -c 500 -m 1000
#
# =============================================================================

set -euo pipefail

# === Default values ===
input_file=""
output_file=""
ncells=""
nfrags=""
min_frags=""
frags_per_cell=""
seed=42
FILTER_CHROMS=true

# === Parse arguments ===
while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input)
      input_file="$2"
      shift 2
      ;;
    -o|--output)
      output_file="$2"
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
    -m|--min-frags)
      min_frags="$2"
      shift 2
      ;;
    --filter-chroms)
      FILTER_CHROMS=true
      shift
      ;;
    --no-filter-chroms)
      FILTER_CHROMS=false
      shift
      ;;
    -s|--seed)
      seed="$2"
      shift 2
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

# === Validate arguments ===
if [[ -z "$input_file" ]]; then
  echo "❌ ERROR: --input is required" >&2
  exit 1
fi

if [[ -z "$output_file" ]]; then
  echo "❌ ERROR: --output is required" >&2
  exit 1
fi

if [[ ! -f "$input_file" ]]; then
  echo "❌ ERROR: Input file not found: $input_file" >&2
  exit 1
fi

if [[ -z "$ncells" && -z "$nfrags" && -z "$frags_per_cell" ]]; then
  echo "❌ ERROR: At least one of --ncells, --nfrags, or --frags-per-cell required" >&2
  exit 1
fi

# === Setup ===
echo "📂 Input: $input_file"
echo "📤 Output: $output_file"
echo "🎲 Seed: $seed"

# Create temp directory
tmpdir=$(mktemp -d -t subsample_XXXXXXXX)
trap 'rm -rf "$tmpdir"' EXIT

# === Decompress input ===
if [[ "$input_file" == *.gz ]]; then
  zcat "$input_file" > "$tmpdir/input.tsv"
else
  cp "$input_file" "$tmpdir/input.tsv"
fi

# === Filter chromosomes ===
if [[ "$FILTER_CHROMS" == "true" ]]; then
  echo "🧬 Filtering to main chromosomes..."
  grep -E '^(chr)?([0-9]+|[XY]|2[AB])\t' "$tmpdir/input.tsv" > "$tmpdir/filtered.tsv" || true
  if [[ -s "$tmpdir/filtered.tsv" ]]; then
    mv "$tmpdir/filtered.tsv" "$tmpdir/input.tsv"
  fi
fi

total_frags=$(wc -l < "$tmpdir/input.tsv")
echo "📊 Total fragments: $total_frags"

# === Get unique cells and counts ===
echo "📊 Counting cells..."
cut -f4 "$tmpdir/input.tsv" | sort | uniq -c | sort -rn > "$tmpdir/cell_counts.txt"
total_cells=$(wc -l < "$tmpdir/cell_counts.txt")
echo "   Found $total_cells unique cells"

# === Apply minimum fragment filter ===
if [[ -n "$min_frags" ]]; then
  echo "🔍 Filtering cells with at least $min_frags fragments..."
  awk -v min="$min_frags" '$1 >= min {print $2}' "$tmpdir/cell_counts.txt" > "$tmpdir/valid_cells.txt"
  valid_count=$(wc -l < "$tmpdir/valid_cells.txt")
  echo "   $valid_count cells pass minimum threshold"
else
  cut -d' ' -f2 "$tmpdir/cell_counts.txt" | sed 's/^ *//' > "$tmpdir/valid_cells.txt"
fi

# === Sample cells ===
if [[ -n "$ncells" ]]; then
  echo "🎲 Sampling $ncells cells..."
  shuf --random-source=<(yes "$seed") "$tmpdir/valid_cells.txt" | head -n "$ncells" > "$tmpdir/sampled_cells.txt"
else
  cp "$tmpdir/valid_cells.txt" "$tmpdir/sampled_cells.txt"
fi

sampled_cell_count=$(wc -l < "$tmpdir/sampled_cells.txt")
echo "   Selected $sampled_cell_count cells"

# === Filter fragments to sampled cells ===
echo "🔀 Extracting fragments for sampled cells..."
# Create lookup file for grep
sort "$tmpdir/sampled_cells.txt" > "$tmpdir/sampled_cells_sorted.txt"

# Use awk for faster filtering
awk 'NR==FNR {cells[$1]=1; next} $4 in cells' \
    "$tmpdir/sampled_cells_sorted.txt" "$tmpdir/input.tsv" > "$tmpdir/cell_filtered.tsv"

filtered_frags=$(wc -l < "$tmpdir/cell_filtered.tsv")
echo "   Fragments after cell filtering: $filtered_frags"

# === Sample fragments ===
if [[ -n "$nfrags" && "$nfrags" -lt "$filtered_frags" ]]; then
  echo "🎲 Sampling $nfrags fragments..."
  shuf --random-source=<(yes "$seed") "$tmpdir/cell_filtered.tsv" | head -n "$nfrags" > "$tmpdir/output.tsv"
elif [[ -n "$frags_per_cell" ]]; then
  echo "🎲 Sampling $frags_per_cell fragments per cell..."
  # Group by cell and sample
  awk -v n="$frags_per_cell" -v seed="$seed" '
    BEGIN { srand(seed) }
    {
      cell = $4
      if (!(cell in count)) count[cell] = 0
      if (count[cell] < n) {
        print
        count[cell]++
      }
    }
  ' "$tmpdir/cell_filtered.tsv" > "$tmpdir/output.tsv"
else
  mv "$tmpdir/cell_filtered.tsv" "$tmpdir/output.tsv"
fi

# === Sort output ===
echo "📝 Sorting output..."
sort -k1,1 -k2,2n "$tmpdir/output.tsv" > "$tmpdir/sorted.tsv"

# === Write output ===
mkdir -p "$(dirname "$output_file")"

if [[ "$output_file" == *.gz ]]; then
  gzip -c "$tmpdir/sorted.tsv" > "$output_file"
else
  mv "$tmpdir/sorted.tsv" "$output_file"
fi

final_frags=$(wc -l < "$tmpdir/sorted.tsv")
final_cells=$(cut -f4 "$tmpdir/sorted.tsv" | sort -u | wc -l)

# === Summary ===
echo ""
echo "================================"
echo "✅ SUBSAMPLING COMPLETED"
echo "================================"
echo "Input:  $total_frags fragments, $total_cells cells"
echo "Output: $final_frags fragments, $final_cells cells"
echo "File:   $output_file"
echo "================================"

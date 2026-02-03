#!/bin/bash

# subsample_fragments.sh
# Subsample a fragments file by cells and/or fragments
#
# Usage:
#   subsample_fragments.sh --input <file> --output <file> [OPTIONS]
#
# Options:
#   --input, -i         Input fragments file (.tsv.gz or .tsv)
#   --output, -o        Output fragments file (.tsv.gz or .tsv)
#   --ncells, -c        Number of cells to sample (optional)
#   --nfrags, -f        Total number of fragments to sample (optional)
#   --frags-per-cell, -p  Number of fragments per cell (optional)
#   --min-frags, -m     Minimum fragments per cell to include (optional)
#   --seed, -s          Random seed for reproducibility (default: 42)
#   --help, -h          Show this help message

set -euo pipefail

# === Default values ===
input_file=""
output_file=""
ncells=""
nfrags=""
min_frags=""
frags_per_cell=""
seed=42

# === Parse command line arguments ===
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
      echo "Use --help for usage information" >&2
      exit 1
      ;;
  esac
done

# === Validate required arguments ===
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

# Check that at least one subsampling option is specified
if [[ -z "$ncells" && -z "$nfrags" && -z "$frags_per_cell" ]]; then
  echo "❌ ERROR: At least one of --ncells, --nfrags, or --frags-per-cell must be specified" >&2
  exit 1
fi

# === Calculate and validate parameter relationships ===
# If all three are given, check consistency: c × p should equal f
if [[ -n "$ncells" && -n "$frags_per_cell" && -n "$nfrags" ]]; then
  expected_frags=$((ncells * frags_per_cell))
  if [[ $expected_frags -ne $nfrags ]]; then
    echo "⚠️  WARNING: --ncells ($ncells) × --frags-per-cell ($frags_per_cell) = $expected_frags, but --nfrags is $nfrags" >&2
    echo "   Using --ncells and --frags-per-cell, ignoring --nfrags" >&2
    nfrags=""  # Ignore nfrags when c and p are provided
  fi
fi

# If two parameters given, calculate the third
if [[ -n "$ncells" && -n "$frags_per_cell" && -z "$nfrags" ]]; then
  nfrags=$((ncells * frags_per_cell))
  echo "ℹ️  Calculated --nfrags = $nfrags (from $ncells cells × $frags_per_cell frags/cell)"
elif [[ -n "$ncells" && -n "$nfrags" && -z "$frags_per_cell" ]]; then
  frags_per_cell=$((nfrags / ncells))
  echo "ℹ️  Calculated --frags-per-cell = $frags_per_cell (from $nfrags frags ÷ $ncells cells)"
elif [[ -n "$frags_per_cell" && -n "$nfrags" && -z "$ncells" ]]; then
  ncells=$((nfrags / frags_per_cell))
  echo "ℹ️  Calculated --ncells = $ncells (from $nfrags frags ÷ $frags_per_cell frags/cell)"
fi

# === Validate command dependencies ===
for cmd in awk sort gzip gunzip shuf; do
  if ! command -v $cmd &> /dev/null; then
    echo "❌ ERROR: Required command not found: $cmd" >&2
    exit 1
  fi
done

# === Setup temporary directory ===
tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT

echo "🔬 Subsampling fragments file..."
echo "   Input: $input_file"
echo "   Output: $output_file"
if [[ -n "$ncells" && -n "$frags_per_cell" ]]; then
  echo "   Target: $ncells cells × $frags_per_cell frags/cell = $((ncells * frags_per_cell)) total frags"
elif [[ -n "$ncells" ]]; then
[[ -n "$min_frags" ]] && echo "   Minimum fragments per cell: $min_frags"
  echo "   Cells to sample: $ncells"
elif [[ -n "$frags_per_cell" ]]; then
  echo "   Fragments per cell: $frags_per_cell"
fi
[[ -n "$nfrags" ]] && echo "   Total fragments limit: $nfrags"
[[ -n "$min_frags" ]] && echo "   Minimum fragments per cell: $min_frags"
echo "   Random seed: $seed"

# === Determine if input is compressed ===
if [[ "$input_file" == *.gz ]]; then
  cat_cmd="gunzip -c"
else
  cat_cmd="cat"
fi

# === Count total fragments and cells in input ===
echo "📊 Analyzing input file..."

# Count fragments per cell and optionally filter by minimum
if [[ -n "$min_frags" ]]; then
  # First count all cells
  all_cells_count=$($cat_cmd "$input_file" | awk '{print $4}' | sort -u | wc -l)
  
  # Then filter by minimum fragments
  $cat_cmd "$input_file" | awk -v min="$min_frags" '
    { count[$4]++ }
    END {
      for (cell in count) {
        if (count[cell] >= min) {
          print cell
        }
      }
    }' | sort -u > "$tmpdir/all_cells.txt"
  
  filtered_cells=$(wc -l < "$tmpdir/all_cells.txt")
  total_frags=$($cat_cmd "$input_file" | wc -l)
  
  echo "   Total cells: $all_cells_count (before filtering)"
  echo "   Cells with ≥$min_frags fragments: $filtered_cells"
  echo "   Total fragments: $total_frags"
  
  if [[ $filtered_cells -eq 0 ]]; then
    echo "❌ ERROR: No cells have ≥$min_frags fragments" >&2
    exit 1
  fi
  
  total_cells=$filtered_cells
else
  $cat_cmd "$input_file" | awk '{print $4}' | sort -u > "$tmpdir/all_cells.txt"
  total_cells=$(wc -l < "$tmpdir/all_cells.txt")
  total_frags=$($cat_cmd "$input_file" | wc -l)
  echo "   Total cells: $total_cells"
  echo "   Total fragments: $total_frags"
fi

# === Step 1: Subsample cells if requested ===
if [[ -n "$ncells" ]]; then
  if [[ $ncells -ge $total_cells ]]; then
    echo "⚠️  WARNING: Requested $ncells cells, but only $total_cells available. Using all cells."
    ncells=$total_cells
  fi
  
  echo "🎲 Sampling $ncells cells..."
  shuf --random-source=<(yes $seed) -n "$ncells" "$tmpdir/all_cells.txt" > "$tmpdir/selected_cells.txt"
  
  # Filter fragments to only selected cells
  $cat_cmd "$input_file" | awk 'NR==FNR{cells[$1]=1; next} $4 in cells' \
    "$tmpdir/selected_cells.txt" - > "$tmpdir/fragments_by_cells.tsv"
  
  input_for_next_step="$tmpdir/fragments_by_cells.tsv"
  cat_next="cat"
else
  input_for_next_step="$input_file"
  cat_next="$cat_cmd"
fi

# === Step 2: Subsample by fragments per cell if requested ===
if [[ -n "$frags_per_cell" ]]; then
  echo "🎲 Sampling $frags_per_cell fragments per cell..."
  
  $cat_next "$input_for_next_step" | \
    awk -v seed="$seed" -v n="$frags_per_cell" '
    BEGIN { srand(seed) }
    {
      cell = $4
      # Store fragment with random number for shuffling
      frags[cell, ++count[cell]] = $0 "\t" rand()
    }
    END {
      for (cell in count) {
        # Collect all fragments for this cell
        n_cell = count[cell]
        delete arr
        for (i = 1; i <= n_cell; i++) {
          arr[i] = frags[cell, i]
        }
        
        # Sort by random number (last field)
        asort(arr)
        
        # Output up to n fragments
        limit = (n < n_cell) ? n : n_cell
        for (i = 1; i <= limit; i++) {
          # Remove random number before printing
          sub(/\t[0-9.]+$/, "", arr[i])
          print arr[i]
        }
      }
    }' > "$tmpdir/fragments_per_cell.tsv"
  
  input_for_next_step="$tmpdir/fragments_per_cell.tsv"
  cat_next="cat"
fi

# === Step 3: Subsample by total fragments if requested (and not already satisfied) ===
# Skip if we already have exactly the right number from ncells × frags_per_cell
if [[ -n "$nfrags" ]]; then
  current_frags=$($cat_next "$input_for_next_step" | wc -l)
  
  # Only subsample if current doesn't match target
  if [[ $current_frags -ne $nfrags ]]; then
    if [[ $nfrags -ge $current_frags ]]; then
      echo "⚠️  WARNING: Requested $nfrags fragments, but only $current_frags available. Using all fragments."
    else
      echo "🎲 Sampling $nfrags total fragments (from $current_frags)..."
      $cat_next "$input_for_next_step" | \
        shuf --random-source=<(yes $seed) -n "$nfrags" > "$tmpdir/fragments_total.tsv"
      
      input_for_next_step="$tmpdir/fragments_total.tsv"
      cat_next="cat"
    fi
  fi
fi

# === Sort output by genomic coordinates ===
echo "📑 Sorting output by coordinates..."
$cat_next "$input_for_next_step" | \
  sort -k1,1 -k2,2n > "$tmpdir/sorted_fragments.tsv"

# === Write output (compress if needed) ===
if [[ "$output_file" == *.gz ]]; then
  echo "📦 Compressing output..."
  gzip -c "$tmpdir/sorted_fragments.tsv" > "$output_file"
else
  cp "$tmpdir/sorted_fragments.tsv" "$output_file"
fi

# === Report final statistics ===
final_cells=$(awk '{print $4}' "$tmpdir/sorted_fragments.tsv" | sort -u | wc -l)
final_frags=$(wc -l < "$tmpdir/sorted_fragments.tsv")
avg_frags_per_cell=$(awk "BEGIN {printf \"%.1f\", $final_frags / $final_cells}")

echo ""
echo "✅ Subsampling complete!"
echo "   Final cells: $final_cells"
echo "   Final fragments: $final_frags"
echo "   Avg fragments/cell: $avg_frags_per_cell"
echo "   Output written to: $output_file"

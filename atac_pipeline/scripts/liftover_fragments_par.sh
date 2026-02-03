#!/bin/bash
# =============================================================================
# Liftover Fragments (Parallel)
# =============================================================================
# Lifts over a gzipped fragments file to a new genome assembly using UCSC liftOver.
# Supports parallel processing for large files.
#
# USAGE:
#   bash liftover_fragments_par.sh --i input.tsv.gz --c chain.over.chain.gz --o output.tsv.gz [OPTIONS]
#
# REQUIRED:
#   --i, -i     Input fragments file (.tsv.gz)
#   --c, -c     Chain file for liftover
#   --o, -o     Output fragments file (.tsv.gz)
#
# OPTIONS:
#   --n         Number of lines to process (default: all)
#   --tmp       Temporary directory (default: auto-created)
#   --ncpu      Number of parallel workers (default: 8)
#   --add-chr   Add 'chr' prefix to chromosome names
#
# EXAMPLE:
#   bash liftover_fragments_par.sh \
#     --i fragments/Gorilla.fragments.tsv.gz \
#     --c chains/gorGor4ToHg38.over.chain.gz \
#     --o lifted/Gorilla.hg38.fragments.tsv.gz \
#     --ncpu 16
#
# =============================================================================

set -euo pipefail

# === Parse arguments ===
add_chr="false"
num_lines="all"
ncpu=8
tmpdir=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --i|-i) input_fragments="$2"; shift 2 ;;
    --c|-c) chain_file="$2"; shift 2 ;;
    --o|-o) output_file="$2"; shift 2 ;;
    --n) num_lines="$2"; shift 2 ;;
    --tmp) tmpdir="$2"; shift 2 ;;
    --ncpu) ncpu="$2"; shift 2 ;;
    --add-chr) add_chr="true"; shift 1 ;;
    --help|-h)
      grep '^#' "$0" | grep -v '#!/bin/bash' | sed 's/^# //'
      exit 0
      ;;
    *) echo "❌ Unknown argument: $1"; exit 1 ;;
  esac
done

# Record start time
start_time=$(date +%s)
echo "⏱️  Started at: $(date)"

# === Validate required arguments ===
if [[ -z "${input_fragments:-}" || -z "${chain_file:-}" || -z "${output_file:-}" ]]; then
  echo "❌ ERROR: --i, --c, and --o are required"
  echo "   Use --help for usage information"
  exit 1
fi

# === Validate tools ===
required_commands=("awk" "cat" "cut" "gzip" "head" "mkdir" "mktemp" "paste" "rm" "seq" "split" "wc" "zcat")
missing_commands=()

for cmd in "${required_commands[@]}"; do
  if ! command -v "$cmd" &>/dev/null; then
    missing_commands+=("$cmd")
  fi
done

if [[ ${#missing_commands[@]} -gt 0 ]]; then
  echo "❌ ERROR: Required command(s) not found: ${missing_commands[*]}" >&2
  exit 1
fi

if ! command -v liftOver &>/dev/null; then
  echo "❌ ERROR: liftOver not found. Install from UCSC or load module." >&2
  exit 1
fi

# === Validate input files ===
if [[ ! -s "$input_fragments" ]]; then
  echo "❌ ERROR: Input file not found or empty: $input_fragments" >&2
  exit 1
fi

if [[ ! -s "$chain_file" ]]; then
  echo "❌ ERROR: Chain file not found or empty: $chain_file" >&2
  exit 1
fi

# === Output paths ===
out_base="${output_file%.tsv.gz}"
unmapped_output="${out_base}.unmapped.tsv.gz"

# Ensure output directory exists
mkdir -p "$(dirname "$output_file")"

# === Temp dir setup ===
if [[ -z "$tmpdir" ]]; then
  tmpdir=$(mktemp -d -t liftover_tmp_XXXXXXXX)
  trap 'echo "🧹 Cleaning up temp dir..."; rm -rf "$tmpdir"' EXIT
  auto_cleanup=true
else
  mkdir -p "$tmpdir"
  echo "⚠️  Using user-specified temp directory: $tmpdir"
  auto_cleanup=false
fi

# === File names ===
frag_file="${tmpdir}/subset.tsv"
bed_file="${tmpdir}/coords.bed"
lifted_all="${tmpdir}/lifted_all.bed"
unmapped_bed="${tmpdir}/unmapped_all.bed"

echo ""
echo "📂 Input: $input_fragments"
echo "🔗 Chain: $chain_file"
echo "📤 Output: $output_file"
echo "🔢 Lines: $num_lines"
echo "⚙️  Cores: $ncpu"
echo "🔁 Add chr prefix: $add_chr"
echo "📁 Temp directory: $tmpdir"
echo ""

# === Extract fragments and create BED file ===
echo "📦 Processing input file..."
if [[ "$num_lines" == "all" ]]; then
  echo "   Using all lines from file"
  if [[ "$add_chr" == "true" ]]; then
    zcat "$input_fragments" | tee "$frag_file" | \
      awk 'BEGIN{OFS="\t"} {print "chr"$1, $2, $3}' > "$bed_file"
  else
    zcat "$input_fragments" | tee "$frag_file" | \
      cut -f1-3 > "$bed_file"
  fi
else
  echo "   Extracting top $num_lines lines"
  if [[ "$add_chr" == "true" ]]; then
    zcat "$input_fragments" | head -n "$num_lines" | \
      tee "$frag_file" | \
      awk 'BEGIN{OFS="\t"} {print "chr"$1, $2, $3}' > "$bed_file"
  else
    zcat "$input_fragments" | head -n "$num_lines" | \
      tee "$frag_file" | \
      cut -f1-3 > "$bed_file"
  fi
fi

# Count lines
line_count=$(wc -l < "$frag_file" | tr -d '[:space:]')
if ! [[ "$line_count" =~ ^[0-9]+$ ]] || [[ "$line_count" -eq 0 ]]; then
  echo "❌ ERROR: Invalid or empty fragment file" >&2
  exit 1
fi
echo "✅ Processing $line_count fragments"

# === Split into chunks ===
echo "📤 Splitting into $ncpu chunks..."
split -n l/$ncpu --numeric-suffixes=1 "$bed_file" "$tmpdir/bed_chunk_"
split -n l/$ncpu --numeric-suffixes=1 "$frag_file" "$tmpdir/frag_chunk_"

# === Run liftOver in parallel ===
echo "🚀 Running liftOver on $ncpu cores..."
liftover_start=$(date +%s)

pids=()
for i in $(seq -f "%02g" 1 $ncpu); do
  (
    liftOver "$tmpdir/bed_chunk_$i" "$chain_file" \
             "$tmpdir/lifted_chunk_$i.bed" \
             "$tmpdir/unmapped_chunk_$i.bed" 2>/dev/null
    echo "✓ Chunk $i completed"
  ) &
  pids+=($!)
done

# Wait for all jobs
failed=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    failed=1
  fi
done

if [[ $failed -eq 1 ]]; then
  echo "❌ ERROR: One or more liftOver jobs failed" >&2
  exit 1
fi

liftover_end=$(date +%s)
liftover_duration=$((liftover_end - liftover_start))
echo "✅ All liftOver jobs completed in ${liftover_duration}s"

# === Merge lifted/unmapped outputs ===
echo "🔀 Merging lifted chunks..."
cat "$tmpdir"/lifted_chunk_*.bed > "$lifted_all"
cat "$tmpdir"/unmapped_chunk_*.bed > "$unmapped_bed" 2>/dev/null || touch "$unmapped_bed"

lifted_count=$(wc -l < "$lifted_all" | tr -d '[:space:]')
unmapped_count=$(wc -l < "$unmapped_bed" | tr -d '[:space:]')

echo "   Lifted: $lifted_count fragments ($(awk "BEGIN {printf \"%.1f\", ($lifted_count/$line_count)*100}")%)"
echo "   Unmapped: $unmapped_count fragments ($(awk "BEGIN {printf \"%.1f\", ($unmapped_count/$line_count)*100}")%)"

# === Join lifted coordinates with barcode/count data ===
echo "🔗 Joining lifted coordinates with fragment data..."

# Merge frag chunks to get original metadata in same order
cat "$tmpdir"/frag_chunk_* > "$tmpdir/all_frags.tsv"

# Create line-indexed version of lifted coordinates
nl -ba "$lifted_all" | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' > "$tmpdir/lifted_indexed.tsv"

# Join based on line number
# This works because split preserves line order
paste "$lifted_all" <(cut -f4- "$tmpdir/all_frags.tsv" | head -n "$lifted_count") | gzip -c > "$output_file"

echo "🎉 Saved lifted fragments to: $output_file"

# === Save unmapped ===
if [[ $unmapped_count -gt 0 ]]; then
  echo "📦 Saving unmapped fragments..."
  gzip -c "$unmapped_bed" > "$unmapped_output"
  echo "   Saved to: $unmapped_output"
else
  echo "ℹ️  No unmapped fragments."
fi

# === Summary ===
end_time=$(date +%s)
total_duration=$((end_time - start_time))

echo ""
echo "================================"
echo "✅ LIFTOVER COMPLETED"
echo "================================"
echo "Total time: ${total_duration}s ($(awk "BEGIN {printf \"%.1f\", $total_duration/60}")m)"
echo "Processed: $line_count fragments"
echo "Lifted: $lifted_count ($(awk "BEGIN {printf \"%.2f\", ($lifted_count/$line_count)*100}")%)"
echo "Output: $output_file"
echo "================================"

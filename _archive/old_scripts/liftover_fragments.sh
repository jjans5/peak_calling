#!/bin/bash

# USAGE:
#   bash liftover_fragments.sh --i input.tsv.gz --c chain.over.chain.gz --o output.tsv.gz [--n 100000] [--tmp tmp_dir]
# Flags:
#   --i     input fragment file (TSV.GZ, 5-column)
#   --c     chain file (e.g., to hg38)
#   --o     output lifted fragment file (TSV.GZ)
#   --n     number of lines (default: all)
#   --tmp   optional temp directory to use (default: mktemp)

# === Parse arguments ===
while [[ $# -gt 0 ]]; do
  case "$1" in
    --i) input_fragments="$2"; shift 2 ;;
    --c) chain_file="$2"; shift 2 ;;
    --o) output_file="$2"; shift 2 ;;
    --n) num_lines="$2"; shift 2 ;;
    --tmp) tmpdir="$2"; shift 2 ;;
    *) echo "❌ Unknown argument: $1"; exit 1 ;;
  esac
done

# === Set defaults ===
num_lines="${num_lines:-all}"

# === Validate inputs ===
if [[ -z "${input_fragments:-}" || -z "${chain_file:-}" || -z "${output_file:-}" ]]; then
  echo "❌ ERROR: --i, --c, and --o are required"
  exit 1
fi

if [[ ! -s "$input_fragments" ]]; then
  echo "❌ ERROR: Input fragments file not found or empty: $input_fragments" >&2
  exit 1
fi

if [[ ! -s "$chain_file" ]]; then
  echo "❌ ERROR: Chain file not found or empty: $chain_file" >&2
  exit 1
fi

if ! command -v liftOver &>/dev/null; then
  echo "❌ ERROR: liftOver not found. Load with: module load ucsc" >&2
  exit 1
fi

# === Output files ===
out_base="${output_file%.tsv.gz}"
unmapped_output="${out_base}.unmapped.tsv.gz"

# === Temp dir ===
if [[ -z "${tmpdir:-}" ]]; then
  tmpdir=$(mktemp -d -t liftover_tmp_XXXXXXXX)
  trap 'rm -rf "$tmpdir"' EXIT
else
  mkdir -p "$tmpdir"
  echo "⚠️  Using user-specified temp directory: $tmpdir"
fi

# === Temp file paths ===
frag_file="${tmpdir}/subset.tsv"
bed_file="${tmpdir}/coords.bed"
lifted_file="${tmpdir}/lifted.bed"
unmapped_bed="${tmpdir}/unmapped.bed"
filtered_file="${tmpdir}/filtered.tsv"

# === Logging ===
echo "📂 Input: $input_fragments"
echo "🔗 Chain: $chain_file"
echo "📤 Output: $output_file"
echo "🔢 Lines: $num_lines"

# === Extract lines ===
if [[ "$num_lines" == "all" ]]; then
  echo "📦 Using all lines..."
  zcat "$input_fragments" > "$frag_file"
else
  echo "📦 Extracting top $num_lines lines..."
  zcat "$input_fragments" | head -n "$num_lines" > "$frag_file"
  line_count=$(wc -l < "$frag_file" | tr -d '[:space:]')
  if ! [[ "$line_count" =~ ^[0-9]+$ ]]; then
    echo "❌ ERROR: Invalid line count: '$line_count'" >&2
    exit 1
  elif [[ "$line_count" -eq 0 ]]; then
    echo "❌ ERROR: No lines extracted from $input_fragments" >&2
    exit 1
  fi
  echo "✅ Extracted $line_count lines"
fi

# === Create BED file with chr prefix ===
awk 'BEGIN{OFS="\t"} {print "chr"$1, $2, $3}' "$frag_file" > "$bed_file"

# === Run liftOver ===
echo "🚀 Running liftOver..."
liftOver "$bed_file" "$chain_file" "$lifted_file" "$unmapped_bed"

if [[ ! -s "$lifted_file" ]]; then
  echo "❌ ERROR: liftOver failed — no lifted output" >&2
  exit 1
fi

lifted_count=$(wc -l < "$lifted_file" | tr -d '[:space:]')
echo "✅ LiftOver complete: $lifted_count fragments lifted"

# === Join lifted coordinates with barcode + count ===
awk 'NR==FNR {a[NR]=$0; next} FNR in a {print a[FNR]}' "$frag_file" "$lifted_file" > "$filtered_file"
paste "$lifted_file" <(cut -f4-5 "$filtered_file") | gzip > "$output_file"
echo "🎉 Saved lifted fragments to: $output_file"

# === Save unmapped ===
if [[ -s "$unmapped_bed" ]]; then
  echo "🧩 Creating unmapped fragment file..."
  awk 'NR==FNR {a[NR]=$0; next} FNR in a {print a[FNR]}' "$frag_file" "$unmapped_bed" > "${tmpdir}/unmapped.tsv"
  gzip -c "${tmpdir}/unmapped.tsv" > "$unmapped_output"
  echo "📦 Saved unmapped fragments to: $unmapped_output"
else
  echo "ℹ️ No unmapped fragments."
fi
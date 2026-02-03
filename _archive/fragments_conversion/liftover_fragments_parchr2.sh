#!/bin/bash

# USAGE:
# bash liftover_fragments_par.sh --i input.tsv.gz --c chain.over.chain.gz --o output.tsv.gz [--n 100000] [--tmp tmp_dir] [--ncpu 8] [--add-chr] [--sort-chr]

# === Parse arguments ===
add_chr="false"
sort_chr="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --i) input_fragments="$2"; shift 2 ;;
    --c) chain_file="$2"; shift 2 ;;
    --o) output_file="$2"; shift 2 ;;
    --n) num_lines="$2"; shift 2 ;;
    --tmp) tmpdir="$2"; shift 2 ;;
    --ncpu) ncpu="$2"; shift 2 ;;
    --add-chr) add_chr="true"; shift 1 ;;
    --sort-chr) sort_chr="true"; shift 1 ;;
    *) echo "❌ Unknown argument: $1"; exit 1 ;;
  esac
done

# === Set defaults ===
num_lines="${num_lines:-all}"
ncpu="${ncpu:-1}"

# === Validate required arguments ===
if [[ -z "${input_fragments:-}" || -z "${chain_file:-}" || -z "${output_file:-}" ]]; then
  echo "❌ ERROR: --i, --c, and --o are required"
  exit 1
fi

# === Validate tools and files ===
if ! command -v liftOver &>/dev/null; then
  echo "❌ ERROR: liftOver not found. Load with: module load ucsc" >&2
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

# === Output paths ===
out_base="${output_file%.tsv.gz}"
unmapped_output="${out_base}.unmapped.tsv.gz"

# === Temp dir setup ===
if [[ -z "${tmpdir:-}" ]]; then
  tmpdir=$(mktemp -d -t liftover_tmp_XXXXXXXX)
  trap 'rm -rf "$tmpdir"' EXIT
else
  mkdir -p "$tmpdir"
  echo "⚠️  Using user-specified temp directory: $tmpdir"
fi

# === File names ===
frag_file="${tmpdir}/subset.tsv"
bed_file="${tmpdir}/coords.bed"
lifted_all="${tmpdir}/lifted_all.bed"
unmapped_bed="${tmpdir}/unmapped_all.bed"
filtered_file="${tmpdir}/filtered.tsv"

echo "📂 Input: $input_fragments"
echo "🔗 Chain: $chain_file"
echo "📤 Output: $output_file"
echo "🔢 Lines: $num_lines"
echo "⚙️  Cores: $ncpu"
echo "🔁 Add chr prefix: $add_chr"
echo "🧪 Sort by chromosome: $sort_chr"

# === Extract fragments ===
if [[ "$num_lines" == "all" ]]; then
  echo "📦 Using all lines..."
  zcat "$input_fragments" > "$frag_file"
else
  echo "📦 Extracting top $num_lines lines..."
  zcat "$input_fragments" | head -n "$num_lines" > "$frag_file"
  line_count=$(wc -l < "$frag_file" | tr -d '[:space:]')
  if ! [[ "$line_count" =~ ^[0-9]+$ ]] || [[ "$line_count" -eq 0 ]]; then
    echo "❌ ERROR: Invalid or empty fragment subset" >&2
    exit 1
  fi
  echo "✅ Extracted $line_count lines"
fi

# === Create BED file ===
echo "🧬 Creating BED file ${add_chr:+(with chr prefix)}..."
if [[ "$add_chr" == "true" ]]; then
  awk 'BEGIN{OFS="\t"} {print "chr"$1, $2, $3}' "$frag_file" > "$bed_file"
else
  cut -f1-3 "$frag_file" > "$bed_file"
fi

# === Optional chromosome sort ===
if [[ "$sort_chr" == "true" ]]; then
  echo "🧮 Sorting by chromosome..."
  paste "$bed_file" "$frag_file" | LC_ALL=C sort -k1,1 -k2,2n --parallel=4 -S 2G > "${tmpdir}/sorted.tsv"
  #paste "$bed_file" "$frag_file" | sort -k1,1 -k2,2n > "${tmpdir}/sorted.tsv"
  cut -f1-3 "${tmpdir}/sorted.tsv" > "$bed_file"
  cut -f4- > "${tmpdir}/sorted.tsv" > "$frag_file"
fi

# === Split into chunks ===
echo "📤 Splitting into $ncpu chunks..."
split -n l/$ncpu --numeric-suffixes=1 "$bed_file" "$tmpdir/bed_chunk_"
split -n l/$ncpu --numeric-suffixes=1 "$frag_file" "$tmpdir/frag_chunk_"

# === Run liftOver in parallel with timing ===
echo "🚀 Running liftOver in parallel with timing..."
for i in $(seq -f "%02g" 1 $ncpu); do
  (
    echo "🕒 Chunk $i started at $(date)"
    start_time=$(date +%s)

    liftOver "$tmpdir/bed_chunk_$i" "$chain_file" \
             "$tmpdir/lifted_chunk_$i.bed" \
             "$tmpdir/unmapped_chunk_$i.bed"

    end_time=$(date +%s)
    duration=$((end_time - start_time))
    echo "✅ Chunk $i finished in ${duration}s at $(date)"
  ) &
done
wait
echo "✅ All liftOver jobs completed"

# === Merge lifted/unmapped outputs ===
cat "$tmpdir"/lifted_chunk_*.bed > "$lifted_all"
cat "$tmpdir"/unmapped_chunk_*.bed > "$unmapped_bed"

# === Join lifted BED with barcode/count ===
awk 'NR==FNR {a[NR]=$0; next} FNR in a {print a[FNR]}' "$frag_file" "$lifted_all" > "$filtered_file"
paste "$lifted_all" <(cut -f4-5 "$filtered_file") | gzip > "$output_file"
echo "🎉 Saved lifted fragments to: $output_file"

# === Save unmapped ===
if [[ -s "$unmapped_bed" ]]; then
  awk 'NR==FNR {a[NR]=$0; next} FNR in a {print a[FNR]}' "$frag_file" "$unmapped_bed" > "$tmpdir/unmapped.tsv"
  gzip -c "$tmpdir/unmapped.tsv" > "$unmapped_output"
  echo "📦 Saved unmapped fragments to: $unmapped_output"
else
  echo "ℹ️ No unmapped fragments."
fi
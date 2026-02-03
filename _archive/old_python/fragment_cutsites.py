#!/usr/bin/env python3
"""
Convert ATAC-seq fragment files to cut-site BED files.

Fragment files contain paired-end read positions (chr, start, end).
This script extracts the Tn5 cut sites at both ends of each fragment,
outputting a BED6 file with + strand for 5' cuts and - strand for 3' cuts.
"""

import os
import re
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# --- CONFIGURATION ---
SPECIES = "Human"  # Change this to switch species

# Input: directory containing fragment files (*.fragments.tsv.gz)
INPUT_DIR = f"../atac/consensus_peak_calling_{SPECIES}_filter/pseudobulk_bed_files/"

# Output: directory for cut-site BED files (saved in current working directory)
OUTPUT_DIR = os.path.join(os.getcwd(), f"fragment_files/{SPECIES}")

# Number of parallel workers (I/O bound, so can be higher than CPU cores)
MAX_WORKERS = 8


def convert_fragments_to_cutsites(input_fragments: str, output_bed: str) -> dict:
    """
    Convert a paired-end fragments file into a BED file of Tn5 cut sites.
    
    For ATAC-seq, each fragment has two Tn5 insertion sites:
    - 5' end (start position) → + strand
    - 3' end (end position - 1) → - strand
    
    Input format:  chr  start  end  [barcode]  [count]
    Output format: BED6 (chr, start, end, name, score, strand)
    
    Args:
        input_fragments: Path to gzipped fragment file
        output_bed: Path for output gzipped BED file
    
    Returns:
        dict with status info
    """
    sample_name = Path(input_fragments).name.split('.')[0]
    
    # Optimized awk command:
    # - Uses OFS for cleaner output
    # - Processes both cut sites per line in single pass
    # - BED is 0-based half-open: start is inclusive, end is exclusive
    awk_cmd = r"""awk -v OFS='\t' '{
        print $1, $2, $2+1, ".", ".", "+";
        print $1, $3-1, $3, ".", ".", "-"
    }'"""
    
    # Pipeline: zcat → awk → gzip
    # Using pigz for parallel compression if available, fallback to gzip
    cmd = f"zcat {input_fragments} | {awk_cmd} | gzip > {output_bed}"
    
    try:
        subprocess.run(
            cmd, 
            shell=True, 
            check=True, 
            executable='/bin/bash',
            capture_output=True,
            text=True
        )
        
        # Get output file size for reporting
        output_size = os.path.getsize(output_bed) / (1024 * 1024)  # MB
        
        return {
            "sample": sample_name,
            "status": "success",
            "output_size_mb": round(output_size, 2),
            "message": f"✅ {sample_name}: {output_size:.1f} MB"
        }
    except subprocess.CalledProcessError as e:
        return {
            "sample": sample_name,
            "status": "error",
            "output_size_mb": 0,
            "message": f"❌ {sample_name}: {e.stderr}"
        }


def process_all_fragments(
    input_dir: str,
    output_dir: str,
    max_workers: int = MAX_WORKERS,
    pattern: str = r"\.fragments\.tsv\.gz$"
) -> list:
    """
    Process all fragment files in a directory in parallel.
    
    Args:
        input_dir: Directory containing fragment files
        output_dir: Directory for output cut-site files
        max_workers: Number of parallel workers
        pattern: Regex pattern to match fragment files
    
    Returns:
        List of result dicts
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    
    # Create output directory
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find all fragment files
    fragment_files = []
    for entry in os.scandir(input_path):
        if entry.is_file() and re.search(pattern, entry.name):
            fragment_files.append(entry)
    
    if not fragment_files:
        print(f"⚠️ No fragment files found in {input_dir}")
        return []
    
    print(f"📂 Found {len(fragment_files)} fragment files")
    print(f"📁 Output directory: {output_dir}")
    print(f"👷 Workers: {max_workers}")
    print("-" * 60)
    
    # Build job list: (input_path, output_path) tuples
    jobs = []
    for entry in fragment_files:
        input_file = entry.path
        output_file = str(output_path / entry.name)
        jobs.append((input_file, output_file))
    
    # Process in parallel using ThreadPoolExecutor (I/O bound task)
    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(convert_fragments_to_cutsites, inp, out): inp 
            for inp, out in jobs
        }
        
        for future in as_completed(futures):
            result = future.result()
            print(result["message"])
            results.append(result)
    
    return results


# --- EXECUTION ---
if __name__ == "__main__":
    results = process_all_fragments(
        input_dir=INPUT_DIR,
        output_dir=OUTPUT_DIR,
        max_workers=MAX_WORKERS
    )
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    successful = [r for r in results if r["status"] == "success"]
    failed = [r for r in results if r["status"] == "error"]
    total_size = sum(r["output_size_mb"] for r in successful)
    
    print(f"Total files processed: {len(results)}")
    print(f"  ✅ Successful: {len(successful)}")
    print(f"  ❌ Failed: {len(failed)}")
    print(f"  💾 Total output size: {total_size:.1f} MB")
    
    if failed:
        print("\nFailed files:")
        for r in failed:
            print(f"  - {r['sample']}")

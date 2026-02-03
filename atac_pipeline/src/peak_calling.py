"""
Peak Calling Functions
======================

Functions for MACS3 peak calling on ATAC-seq data.
Includes fragment-to-cutsite conversion and parallel execution.
"""

import os
import re
import subprocess
import json
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from functools import partial
from typing import Optional, Dict, List, Any

# =============================================================================
# CONSTANTS
# =============================================================================

EFFECTIVE_GENOME_SIZES = {
    'Bonobo': 2595269547,
    'Macaque': 2653677440,
    'Chimpanzee': 2792339170,
    'Gorilla': 2661668758,
    'Marmoset': 2597026658,
    'Human': 2913022398,  # From MACS3 documentation (deeptools)
}

DEFAULT_MACS3_PARAMS = {
    "format": "BED",
    "qvalue": 0.01,
    "shift": -73,       # Negative shift for ATAC-seq to center on Tn5 cut site
    "extsize": 146,     # Extension size for reads
    "keep_dup": "all",  # Keep all duplicates for ATAC-seq
    "min_length": 200,  # Minimum peak length
    "nomodel": True,    # Use shift/extsize directly, skip model building
    "call_summits": True,  # Call peak summits (required for narrowPeak)
    "nolambda": True,   # Use fixed background lambda
}

# =============================================================================
# FRAGMENT CONVERSION
# =============================================================================

def convert_fragments_to_cutsites(input_fragments: str, output_bed: str) -> Dict[str, Any]:
    """
    Convert a paired-end fragments file into a BED file of Tn5 cut sites.
    
    For ATAC-seq, each fragment has two Tn5 insertion sites:
    - 5' end (start position) → + strand
    - 3' end (end position - 1) → - strand
    
    Parameters
    ----------
    input_fragments : str
        Path to gzipped fragment file (.fragments.tsv.gz)
    output_bed : str
        Path for output gzipped BED file
    
    Returns
    -------
    dict
        Result containing sample name, status, output size, and message
    """
    sample_name = Path(input_fragments).name.split('.')[0]
    
    # awk command to extract both cut sites per fragment
    # BED is 0-based half-open: start is inclusive, end is exclusive
    awk_cmd = r"""awk -v OFS='\t' '{
        print $1, $2, $2+1, ".", ".", "+";
        print $1, $3-1, $3, ".", ".", "-"
    }'"""
    
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
    max_workers: int = 8,
    pattern: str = r"\.fragments\.tsv\.gz$"
) -> List[Dict[str, Any]]:
    """
    Process all fragment files in a directory in parallel.
    
    Parameters
    ----------
    input_dir : str
        Directory containing fragment files
    output_dir : str
        Directory for output cut-site files
    max_workers : int
        Number of parallel workers (I/O bound, so can be higher than CPU cores)
    pattern : str
        Regex pattern to match fragment files
    
    Returns
    -------
    list
        List of result dictionaries
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    
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
    
    # Build job list
    jobs = [(entry.path, str(output_path / entry.name)) for entry in fragment_files]
    
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


# =============================================================================
# MACS3 PEAK CALLING
# =============================================================================

def build_macs3_command(
    sample_name: str,
    fragment_path: str,
    species: str,
    out_dir: str,
    macs3_path: str,
    params: Dict[str, Any]
) -> List[str]:
    """
    Build MACS3 command with configurable parameters.
    
    Parameters
    ----------
    sample_name : str
        Name for output files
    fragment_path : str
        Path to input fragment/cut-site file
    species : str
        Species name (must be in EFFECTIVE_GENOME_SIZES)
    out_dir : str
        Output directory
    macs3_path : str
        Path to MACS3 executable
    params : dict
        MACS3 parameters
    
    Returns
    -------
    list
        Command as list of strings
    """
    gsize = EFFECTIVE_GENOME_SIZES.get(species)
    if gsize is None:
        raise ValueError(f"Unknown species: {species}. Available: {list(EFFECTIVE_GENOME_SIZES.keys())}")
    
    cmd = [
        macs3_path, "callpeak",
        "--treatment", fragment_path,
        "--name", sample_name,
        "--outdir", out_dir,
        "--format", params["format"],
        "--gsize", str(gsize),
        "--qvalue", str(params["qvalue"]),
        "--shift", str(params["shift"]),
        "--extsize", str(params["extsize"]),
        "--keep-dup", str(params["keep_dup"]),
        "--min-length", str(params["min_length"]),
    ]
    
    # Add boolean flags
    if params.get("nomodel"):
        cmd.append("--nomodel")
    if params.get("call_summits"):
        cmd.append("--call-summits")
    if params.get("nolambda"):
        cmd.append("--nolambda")
    
    return cmd


def run_macs3_worker(
    job: tuple,
    species: str,
    out_dir: str,
    macs3_path: str,
    params: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Worker function for parallel MACS3 execution.
    
    Parameters
    ----------
    job : tuple
        (sample_name, fragment_path)
    species : str
        Species name
    out_dir : str
        Output directory
    macs3_path : str
        Path to MACS3 executable
    params : dict
        MACS3 parameters
    
    Returns
    -------
    dict
        Result with sample_name, status, peak_count, and message
    """
    sample_name, fragment_path = job
    
    cmd = build_macs3_command(sample_name, fragment_path, species, out_dir, macs3_path, params)
    
    print(f"🚀 Starting: {sample_name}")
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Count peaks from the narrowPeak file
        narrowpeak_file = os.path.join(out_dir, f"{sample_name}_peaks.narrowPeak")
        peak_count = 0
        if os.path.exists(narrowpeak_file):
            with open(narrowpeak_file, 'r') as f:
                peak_count = sum(1 for _ in f)
        
        return {
            "sample_name": sample_name,
            "status": "success",
            "peak_count": peak_count,
            "message": f"✅ Finished: {sample_name} ({peak_count:,} peaks)"
        }
    except subprocess.CalledProcessError as e:
        return {
            "sample_name": sample_name,
            "status": "error",
            "peak_count": 0,
            "message": f"❌ Error in {sample_name}: {e.stderr}"
        }


def run_peak_calling(
    species: str,
    frag_dir: str,
    out_dir: str,
    macs3_path: str = "macs3",
    max_workers: int = 15,
    params: Optional[Dict[str, Any]] = None,
    **param_overrides
) -> List[Dict[str, Any]]:
    """
    Run MACS3 peak calling in parallel for all fragment files.
    
    Parameters
    ----------
    species : str
        Species name (must be in EFFECTIVE_GENOME_SIZES)
    frag_dir : str
        Directory with cut-site fragment files
    out_dir : str
        Output directory for peaks
    macs3_path : str
        Path to MACS3 executable
    max_workers : int
        Number of parallel workers/cores (default: 15)
    params : dict, optional
        Full parameter dict (if None, uses DEFAULT_MACS3_PARAMS)
    **param_overrides
        Individual parameters to override (e.g., qvalue=0.05)
    
    Returns
    -------
    list
        List of result dicts containing sample info and peak counts
    
    Notes
    -----
    Output files created by MACS3 per sample:
        - {sample}_peaks.narrowPeak: BED6+4 format with peak information
        - {sample}_peaks.xls: Spreadsheet with peak info
        - {sample}_summits.bed: Peak summit positions
    """
    # Build final parameters
    final_params = (params if params is not None else DEFAULT_MACS3_PARAMS).copy()
    final_params.update(param_overrides)
    
    # Ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)
    
    # Find all fragment files
    fragment_files = [f for f in os.listdir(frag_dir) if f.endswith(".fragments.tsv.gz")]
    
    if not fragment_files:
        print(f"⚠️ No fragment files found in {frag_dir}")
        return []
    
    print(f"📂 Found {len(fragment_files)} fragment files for {species}")
    print(f"📁 Output directory: {out_dir}")
    print(f"🧬 Genome size: {EFFECTIVE_GENOME_SIZES[species]:,}")
    print(f"⚙️ Parameters: qvalue={final_params['qvalue']}, shift={final_params['shift']}, "
          f"extsize={final_params['extsize']}, min_length={final_params['min_length']}")
    print(f"👷 Workers: {max_workers}")
    print("-" * 60)
    
    # Save parameters to file
    params_file = os.path.join(out_dir, "macs3_parameters.json")
    params_to_save = {
        "species": species,
        "genome_size": EFFECTIVE_GENOME_SIZES[species],
        "macs3_path": macs3_path,
        "max_workers": max_workers,
        "frag_dir": frag_dir,
        "out_dir": out_dir,
        "run_date": datetime.now().isoformat(),
        "macs3_params": final_params
    }
    with open(params_file, 'w') as f:
        json.dump(params_to_save, f, indent=2)
    print(f"💾 Parameters saved to: {params_file}")
    print("-" * 60)
    
    # Create jobs list
    jobs = [(f.split('.')[0], os.path.join(frag_dir, f)) for f in fragment_files]
    
    # Create worker with fixed arguments
    worker = partial(
        run_macs3_worker,
        species=species,
        out_dir=out_dir,
        macs3_path=macs3_path,
        params=final_params
    )
    
    # Run in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(worker, jobs))
    
    # Generate peak count report
    report_file = os.path.join(out_dir, "peak_counts_report.tsv")
    with open(report_file, 'w') as f:
        f.write("cell_type\tpeak_count\tstatus\n")
        for result in results:
            f.write(f"{result['sample_name']}\t{result['peak_count']}\t{result['status']}\n")
    print(f"\n📊 Peak count report saved to: {report_file}")
    
    return results

import os
import subprocess
import json
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
from functools import partial

# --- SPECIES AND GENOME SIZES ---
EFFECTIVE_GENOME_SIZES = {
    'Bonobo': 2595269547,
    'Macaque': 2653677440,
    'Chimpanzee': 2792339170,
    'Gorilla': 2661668758,
    'Marmoset': 2597026658,
    'Human': 2913022398  # value from macs3 site (deeptools)
}

# --- CONFIGURATION ---
SPECIES = "Human"  # Change this to switch species

# Input/Output directories
BASE_DIR = "/cluster/home/jjanssens/jjans/analysis/adult_intestine/peaks"
FRAG_DIR = os.path.join(BASE_DIR, f"fragment_files/{SPECIES}")
# Output will be saved in current working directory under consensus_peak_calling_{SPECIES}/
OUT_DIR_BASE = os.path.join(os.getcwd(), f"consensus_peak_calling_{SPECIES}")
MACS3_PATH = "/cluster/project/treutlein/jjans/software/miniforge3/envs/scenicplus/bin/macs3"

# --- MACS3 PARAMETERS (easy to modify) ---
# Full documentation: https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md
MACS3_PARAMS = {
    # Input format: BED, BAM, SAM, BEDPE, etc.
    "format": "BED",
    # q-value (minimum FDR) cutoff for peak detection
    "qvalue": 0.01,
    # Shift reads by this amount (negative for ATAC-seq)
    "shift": -73,
    # Extend reads to this fragment size
    "extsize": 146,
    # How to handle duplicate reads: "auto", "all", or integer
    "keep_dup": "all",
    # Minimum length of peak region
    "min_length": 200,
    # Flags (set to True to enable)
    "nomodel": True,       # Skip model building, use shift/extsize directly
    "call_summits": True,  # Call peak summits (required for narrowPeak output)
    "nolambda": True,      # Use fixed background lambda
}

# Number of parallel workers (cores)
MAX_WORKERS = 15


def build_macs3_command(sample_name, fragment_path, species, out_dir, macs3_path, params):
    """Build MACS3 command with configurable parameters."""
    
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


def run_macs3_worker(job, species, out_dir, macs3_path, params):
    """Worker function for parallel execution.
    
    Returns:
        dict with sample_name, status, peak_count, and error message if any
    """
    sample_name, fragment_path = job
    
    cmd = build_macs3_command(sample_name, fragment_path, species, out_dir, macs3_path, params)
    
    print(f"🚀 Starting: {sample_name}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
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
    species,
    frag_dir=None,
    out_dir=None,
    macs3_path=MACS3_PATH,
    max_workers=MAX_WORKERS,
    params=None,
    **param_overrides
):
    """
    Run MACS3 peak calling in parallel for all fragment files.
    
    Args:
        species: Species name (must be in EFFECTIVE_GENOME_SIZES)
        frag_dir: Directory with fragment files (default: BASE_DIR/fragment_files/{species})
        out_dir: Output directory (default: ./consensus_peak_calling_{species}/)
        macs3_path: Path to macs3 executable
        max_workers: Number of parallel workers/cores (default: 15)
        params: Full parameter dict (if None, uses MACS3_PARAMS)
        **param_overrides: Individual parameters to override (e.g., qvalue=0.05)
    
    Returns:
        List of result dicts containing sample info and peak counts
    
    Output files created by MACS3:
        - {sample}_peaks.narrowPeak: BED6+4 format with peak information
        - {sample}_peaks.xls: Spreadsheet with peak info
        - {sample}_summits.bed: Peak summit positions
        - {sample}_treat_pileup.bdg: Treatment pileup (if --bdg flag)
        - {sample}_control_lambda.bdg: Local lambda (if --bdg flag)
    """
    # Set defaults based on species
    if frag_dir is None:
        frag_dir = os.path.join(BASE_DIR, f"fragment_files/{species}")
    if out_dir is None:
        # Save in current working directory
        out_dir = os.path.join(os.getcwd(), f"consensus_peak_calling_{species}")
    
    # Build final parameters
    final_params = (params if params is not None else MACS3_PARAMS).copy()
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
    print(f"⚙️ Parameters: qvalue={final_params['qvalue']}, shift={final_params['shift']}, extsize={final_params['extsize']}, min_length={final_params['min_length']}")
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


# --- EXECUTION ---
if __name__ == "__main__":
    results = run_peak_calling(
        species=SPECIES,
        max_workers=MAX_WORKERS,
        # Easy parameter overrides - uncomment/modify as needed:
        # qvalue=0.05,
        # min_length=150,
        # shift=-100,
        # extsize=200,
    )
    
    print("\n" + "=" * 60)
    print("SUMMARY - Peak counts per cell type")
    print("=" * 60)
    
    total_peaks = 0
    successful = 0
    failed = 0
    
    for result in results:
        print(result["message"])
        if result["status"] == "success":
            successful += 1
            total_peaks += result["peak_count"]
        else:
            failed += 1
    
    print("\n" + "-" * 60)
    print(f"Total samples processed: {len(results)}")
    print(f"  ✅ Successful: {successful}")
    print(f"  ❌ Failed: {failed}")
    print(f"  📊 Total peaks called: {total_peaks:,}")
    print(f"\nOutput saved to: {os.path.join(os.getcwd(), f'consensus_peak_calling_{SPECIES}')}")
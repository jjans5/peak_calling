"""
Liftover Functions
==================

Functions for lifting over genomic coordinates between assemblies.
Supports both peak files and fragment files.
"""

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Dict, Any, List
import pandas as pd

# =============================================================================
# CHAIN FILES
# =============================================================================

# Default chain file directory (Treutlein lab shared location)
DEFAULT_CHAIN_DIR = "/cluster/work/treutlein/jjans/data/intestine/nhp_atlas/genomes/chain_files"

# Chain file mapping (species -> hg38)
CHAIN_FILES = {
    "Bonobo": "panPan2ToHg38.over.chain",
    "Chimpanzee": "panTro5ToHg38.over.chain",
    "Gorilla": "gorGor4ToHg38.over.chain",
    "Macaque": "rheMac10ToHg38.over.chain",
    # Marmoset requires two-step liftover (calJac1 -> calJac4 -> hg38)
    "Marmoset_step1": "calJac1ToCalJac4.over.chain",
    "Marmoset_step2": "calJac4ToHg38.over.chain",
}


def print_chain_info(chain_dir: str = None):
    """Print available chain files and their locations."""
    chain_dir = chain_dir or DEFAULT_CHAIN_DIR
    print("=" * 70)
    print("CHAIN FILES FOR LIFTOVER")
    print("=" * 70)
    print(f"📁 Chain file directory: {chain_dir}")
    print()
    print("Available species -> hg38 chain files:")
    for species, chain_file in CHAIN_FILES.items():
        full_path = os.path.join(chain_dir, chain_file)
        exists = "✅" if os.path.exists(full_path) else "❌"
        print(f"  {exists} {species}: {chain_file}")
    print("=" * 70)


def get_chain_file(species: str, chain_dir: str = None) -> str:
    """
    Get the chain file path for a species.
    
    Parameters
    ----------
    species : str
        Species name (Bonobo, Chimpanzee, Gorilla, Macaque, Marmoset)
    chain_dir : str, optional
        Directory containing chain files. Defaults to DEFAULT_CHAIN_DIR.
        
    Returns
    -------
    str
        Full path to chain file
    """
    chain_dir = chain_dir or DEFAULT_CHAIN_DIR
    
    if species not in CHAIN_FILES:
        raise ValueError(f"Unknown species: {species}. Available: {list(CHAIN_FILES.keys())}")
    
    chain_file = os.path.join(chain_dir, CHAIN_FILES[species])
    
    print(f"🔗 Using chain file: {chain_file}")
    
    if not os.path.exists(chain_file):
        print(f"⚠️  Warning: Chain file not found at {chain_file}")
    
    return chain_file


# =============================================================================
# PEAK LIFTOVER
# =============================================================================

def liftover_peaks(
    input_bed: str,
    output_bed: str,
    chain_file: str,
    liftover_path: str = "liftOver",
    min_match: float = 0.95,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Liftover a BED file to a new assembly using UCSC liftOver.
    
    Parameters
    ----------
    input_bed : str
        Path to input BED file
    output_bed : str
        Path for output BED file
    chain_file : str
        Path to chain file for liftover
    liftover_path : str
        Path to liftOver executable
    min_match : float
        Minimum ratio of bases that must remap (default 0.95)
    verbose : bool
        Print chain file info
    
    Returns
    -------
    dict
        Result with counts of lifted and unmapped peaks
    """
    if verbose:
        print(f"🔗 Chain file: {chain_file}")
    
    # Create unmapped output path
    unmapped_bed = output_bed.replace(".bed", ".unmapped.bed")
    
    cmd = [
        liftover_path,
        "-minMatch=" + str(min_match),
        input_bed,
        chain_file,
        output_bed,
        unmapped_bed
    ]
    
    try:
        # Check if input file exists
        if not os.path.exists(input_bed):
            print(f"❌ ERROR: Input file not found: {input_bed}")
            return {
                "status": "error",
                "lifted": 0,
                "unmapped": 0,
                "message": f"❌ Input file not found: {input_bed}"
            }
        
        # Check if chain file exists
        if not os.path.exists(chain_file):
            print(f"❌ ERROR: Chain file not found: {chain_file}")
            return {
                "status": "error",
                "lifted": 0,
                "unmapped": 0,
                "message": f"❌ Chain file not found: {chain_file}"
            }
        
        # Check input file content
        input_count = sum(1 for _ in open(input_bed))
        if verbose:
            print(f"📄 Input file has {input_count:,} lines")
            # Show first few chromosome names
            with open(input_bed) as f:
                first_chroms = [line.split('\t')[0] for line in [next(f, '') for _ in range(3)] if line]
            print(f"📄 First chromosomes in input: {first_chroms}")
        
        # Check chain file chromosome naming
        if verbose:
            with open(chain_file) as f:
                for line in f:
                    if line.startswith('chain'):
                        parts = line.split()
                        if len(parts) > 2:
                            print(f"📄 Chain file expects chromosomes like: {parts[2]}")
                        break
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_bed) or '.', exist_ok=True)
        
        if verbose:
            print(f"🔧 Running: {' '.join(cmd)}")
        
        # Run liftOver
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        
        # Check for errors
        if result.returncode != 0:
            error_msg = result.stderr if result.stderr else result.stdout
            print(f"❌ liftOver failed with return code {result.returncode}")
            print(f"   stderr: {result.stderr}")
            print(f"   stdout: {result.stdout}")
            return {
                "status": "error",
                "lifted": 0,
                "unmapped": 0,
                "message": f"❌ liftOver failed: {error_msg}",
                "command": " ".join(cmd)
            }
        
        # Count results
        lifted_count = sum(1 for _ in open(output_bed)) if os.path.exists(output_bed) else 0
        unmapped_count = sum(1 for _ in open(unmapped_bed)) if os.path.exists(unmapped_bed) else 0
        
        # If both are 0, something went wrong
        if lifted_count == 0 and unmapped_count == 0:
            print(f"⚠️  WARNING: No peaks were processed!")
            print(f"   Input had {input_count} lines")
            print(f"   This usually means chromosome names don't match between BED and chain file")
            print(f"   Check if your BED has 'chr' prefix and chain file expects it (or vice versa)")
            
            # Show what chromosomes are in the input
            with open(input_bed) as f:
                chroms = set()
                for i, line in enumerate(f):
                    if i >= 100:
                        break
                    chroms.add(line.split('\t')[0])
            print(f"   Chromosomes in input (first 100 lines): {sorted(chroms)[:10]}")
            
            return {
                "status": "warning",
                "lifted": 0,
                "unmapped": 0,
                "output_file": output_bed,
                "unmapped_file": unmapped_bed,
                "message": f"⚠️  No peaks processed. Input had {input_count} lines. Chromosome naming mismatch?",
                "command": " ".join(cmd)
            }
        
        return {
            "status": "success",
            "lifted": lifted_count,
            "unmapped": unmapped_count,
            "output_file": output_bed,
            "unmapped_file": unmapped_bed,
            "message": f"✅ Lifted {lifted_count:,} peaks, {unmapped_count:,} unmapped"
        }
    except subprocess.CalledProcessError as e:
        print(f"❌ Subprocess error: {e.stderr}")
        return {
            "status": "error",
            "lifted": 0,
            "unmapped": 0,
            "message": f"❌ Error: {e.stderr}",
            "command": " ".join(cmd)
        }
    except Exception as e:
        print(f"❌ Unexpected error: {str(e)}")
        import traceback
        traceback.print_exc()
        return {
            "status": "error",
            "lifted": 0,
            "unmapped": 0,
            "message": f"❌ Unexpected error: {str(e)}",
            "command": " ".join(cmd)
        }


def liftover_narrowpeak(
    input_narrowpeak: str,
    output_narrowpeak: str,
    chain_file: str,
    liftover_path: str = "liftOver",
) -> Dict[str, Any]:
    """
    Liftover a narrowPeak file, preserving all columns.
    
    narrowPeak format: chr, start, end, name, score, strand, signalValue, pValue, qValue, summit
    
    Parameters
    ----------
    input_narrowpeak : str
        Path to input narrowPeak file
    output_narrowpeak : str
        Path for output narrowPeak file
    chain_file : str
        Path to chain file
    liftover_path : str
        Path to liftOver executable
    
    Returns
    -------
    dict
        Result dictionary
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Extract BED3 for liftover
        bed3_file = os.path.join(tmpdir, "coords.bed")
        extra_file = os.path.join(tmpdir, "extra.tsv")
        lifted_bed3 = os.path.join(tmpdir, "lifted.bed")
        unmapped_file = os.path.join(tmpdir, "unmapped.bed")
        
        # Read narrowPeak
        df = pd.read_csv(input_narrowpeak, sep="\t", header=None)
        
        # Save BED3 (chr, start, end) with line numbers for rejoining
        df[[0, 1, 2]].to_csv(bed3_file, sep="\t", header=False, index=True)
        
        # Save extra columns
        df.iloc[:, 3:].to_csv(extra_file, sep="\t", header=False, index=True)
        
        # Run liftOver
        cmd = [liftover_path, bed3_file, chain_file, lifted_bed3, unmapped_file]
        subprocess.run(cmd, capture_output=True, check=True)
        
        # Read lifted coordinates
        lifted_df = pd.read_csv(lifted_bed3, sep="\t", header=None,
                                names=["idx", "Chromosome", "Start", "End"])
        
        # Read extra columns
        extra_df = pd.read_csv(extra_file, sep="\t", header=None, index_col=0)
        
        # Join on index
        result_df = lifted_df.set_index("idx").join(extra_df)
        result_df = result_df.reset_index(drop=True)
        
        # Recalculate summit relative to new start
        # Summit is the 10th column (index 9), relative to Start
        # This is an approximation - summit might not lift correctly
        
        # Save output
        result_df.to_csv(output_narrowpeak, sep="\t", header=False, index=False)
        
        return {
            "status": "success",
            "lifted": len(result_df),
            "original": len(df),
            "output_file": output_narrowpeak,
        }


# =============================================================================
# FRAGMENT LIFTOVER
# =============================================================================

def liftover_fragments_parallel(
    input_fragments: str,
    output_fragments: str,
    chain_file: str,
    ncpu: int = 8,
    add_chr: bool = False,
    liftover_script: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Liftover a fragments file using parallel processing.
    
    This wraps the liftover_fragments_par.sh script for efficient
    parallel liftover of large fragment files.
    
    Parameters
    ----------
    input_fragments : str
        Path to input fragments file (.tsv.gz)
    output_fragments : str
        Path for output fragments file (.tsv.gz)
    chain_file : str
        Path to chain file
    ncpu : int
        Number of parallel workers
    add_chr : bool
        Whether to add 'chr' prefix to chromosome names
    liftover_script : str, optional
        Path to liftover_fragments_par.sh script
    
    Returns
    -------
    dict
        Result dictionary
    """
    if liftover_script is None:
        # Try to find script in same directory as this module
        module_dir = Path(__file__).parent.parent
        liftover_script = str(module_dir / "scripts" / "liftover_fragments_par.sh")
    
    if not os.path.exists(liftover_script):
        raise FileNotFoundError(f"Liftover script not found: {liftover_script}")
    
    cmd = [
        "bash", liftover_script,
        "--i", input_fragments,
        "--c", chain_file,
        "--o", output_fragments,
        "--ncpu", str(ncpu),
    ]
    
    if add_chr:
        cmd.append("--add-chr")
    
    print(f"🚀 Running liftover for {Path(input_fragments).name}...")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Parse output for statistics
        output_lines = result.stdout.split("\n")
        lifted_count = 0
        unmapped_count = 0
        
        for line in output_lines:
            if "Lifted:" in line:
                lifted_count = int(line.split()[1])
            elif "Unmapped:" in line:
                unmapped_count = int(line.split()[1])
        
        return {
            "status": "success",
            "lifted": lifted_count,
            "unmapped": unmapped_count,
            "output_file": output_fragments,
            "message": f"✅ Lifted {lifted_count:,} fragments"
        }
    except subprocess.CalledProcessError as e:
        return {
            "status": "error",
            "lifted": 0,
            "unmapped": 0,
            "message": f"❌ Error: {e.stderr}"
        }


def run_liftover_all_species(
    fragments_dir: str,
    output_dir: str,
    chain_dir: str,
    ncpu: int = 30,
    species_list: Optional[List[str]] = None,
) -> Dict[str, Dict[str, Any]]:
    """
    Run liftover for all species to human (hg38).
    
    Parameters
    ----------
    fragments_dir : str
        Directory containing species fragment files
    output_dir : str
        Output directory for lifted fragments
    chain_dir : str
        Directory containing chain files
    ncpu : int
        Number of parallel workers
    species_list : list, optional
        List of species to process (default: all non-Human species)
    
    Returns
    -------
    dict
        Results for each species
    """
    if species_list is None:
        species_list = ["Bonobo", "Chimpanzee", "Gorilla", "Macaque"]
    
    os.makedirs(output_dir, exist_ok=True)
    results = {}
    
    for species in species_list:
        print(f"\n{'='*60}")
        print(f"Processing {species}")
        print('='*60)
        
        # Find input file
        input_pattern = f"{species}*.fragments.tsv.gz"
        input_files = list(Path(fragments_dir).glob(input_pattern))
        
        if not input_files:
            print(f"⚠️  No fragment files found for {species}")
            continue
        
        input_file = str(input_files[0])
        chain_file = os.path.join(chain_dir, CHAIN_FILES.get(species, ""))
        
        if not os.path.exists(chain_file):
            print(f"⚠️  Chain file not found for {species}: {chain_file}")
            continue
        
        output_file = os.path.join(output_dir, f"{species}.hg38.fragments.tsv.gz")
        
        result = liftover_fragments_parallel(
            input_fragments=input_file,
            output_fragments=output_file,
            chain_file=chain_file,
            ncpu=ncpu,
        )
        
        results[species] = result
        print(result["message"])
    
    return results

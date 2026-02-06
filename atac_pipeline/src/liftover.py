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

def _liftover_chunk(args):
    """Worker function to liftover a single chunk."""
    chunk_bed, chain_file, liftover_path, min_match, min_blocks, multiple, output_bed, unmapped_bed = args
    
    cmd = [
        liftover_path,
        f"-minMatch={min_match}",
    ]
    if min_blocks is not None:
        cmd.append(f"-minBlocks={min_blocks}")
    if multiple:
        cmd.append("-multiple")
    cmd.extend([
        chunk_bed,
        chain_file,
        output_bed,
        unmapped_bed
    ])
    
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    return result.returncode == 0


def liftover_peaks(
    input_bed: str,
    output_bed: str,
    chain_file: str,
    liftover_path: str = "liftOver",
    min_match: float = 0.95,
    min_blocks: float = None,
    multiple: bool = False,
    verbose: bool = True,
    auto_chr: bool = True,
    ncpu: int = 1,
    unmapped_bed: str = None,
) -> Dict[str, Any]:
    """
    Liftover a BED file to a new assembly using UCSC liftOver.
    
    Handles BED files with extra columns (beyond chr, start, end) by
    temporarily extracting coordinates, running liftOver, and rejoining.
    
    Automatically handles chromosome naming mismatches (chr1 vs 1).
    
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
    min_blocks : float, optional
        Minimum ratio of alignment blocks that must map (default: not used)
    multiple : bool
        Allow multiple output regions for a single input (default False)
    verbose : bool
        Print chain file info
    auto_chr : bool
        Automatically fix chr prefix mismatches (default True)
    ncpu : int
        Number of parallel workers (default 1 = sequential)
    unmapped_bed : str, optional
        Path for unmapped peaks (default: output_bed with .unmapped.bed suffix)
    
    Returns
    -------
    dict
        Result with counts of lifted and unmapped peaks
    """
    if verbose:
        print(f"🔗 Chain file: {chain_file}")
        if ncpu > 1:
            print(f"⚙️  Using {ncpu} parallel workers")
    
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
    
    # Get first chromosome from input
    with open(input_bed) as f:
        first_line = f.readline().strip()
        input_chrom = first_line.split('\t')[0] if first_line else ""
    input_has_chr = input_chrom.startswith('chr')
    
    # Get expected chromosome format from chain file
    chain_has_chr = False
    with open(chain_file) as f:
        for line in f:
            if line.startswith('chain'):
                parts = line.split()
                if len(parts) > 2:
                    chain_has_chr = parts[2].startswith('chr')
                break
    
    if verbose:
        print(f"📄 Input file has {input_count:,} lines")
        print(f"📄 Input chromosomes: {'chr1, chr2...' if input_has_chr else '1, 2...'}")
        print(f"📄 Chain expects: {'chr1, chr2...' if chain_has_chr else '1, 2...'}")
    
    # Determine if we need to fix chromosome naming
    need_add_chr = auto_chr and chain_has_chr and not input_has_chr
    need_remove_chr = auto_chr and not chain_has_chr and input_has_chr
    
    if need_add_chr and verbose:
        print(f"🔧 Auto-adding 'chr' prefix to match chain file")
    elif need_remove_chr and verbose:
        print(f"🔧 Auto-removing 'chr' prefix to match chain file")
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_bed) or '.', exist_ok=True)
    
    # Create unmapped output path
    if unmapped_bed is None:
        unmapped_bed = output_bed.replace(".bed", ".unmapped.bed")
    
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Read input BED and check how many columns it has
            with open(input_bed) as f:
                first_line = f.readline().strip()
                num_cols = len(first_line.split('\t'))
            
            if verbose:
                print(f"📄 Input has {num_cols} columns")
            
            # Prepare files for liftover
            bed3_file = os.path.join(tmpdir, "coords.bed")
            extra_file = os.path.join(tmpdir, "extra.tsv")
            lifted_bed3 = os.path.join(tmpdir, "lifted.bed")
            unmapped_tmp = os.path.join(tmpdir, "unmapped.bed")
            
            if verbose:
                print(f"📄 Splitting BED file for liftOver (preserving extra columns)...")
            
            # Split input: BED3 with index, and extra columns with index
            # Also handle chr prefix if needed
            with open(input_bed) as fin, \
                 open(bed3_file, 'w') as fbed, \
                 open(extra_file, 'w') as fextra:
                for idx, line in enumerate(fin):
                    parts = line.rstrip('\n').split('\t')
                    chrom = parts[0]
                    
                    # Fix chromosome naming if needed
                    if need_add_chr and not chrom.startswith('chr'):
                        chrom = 'chr' + chrom
                    elif need_remove_chr and chrom.startswith('chr'):
                        chrom = chrom[3:]
                    
                    # Write chr, start, end, and use index as name (4th column)
                    fbed.write(f"{chrom}\t{parts[1]}\t{parts[2]}\t{idx}\n")
                    # Write extra columns (if any beyond the first 3)
                    if len(parts) > 3:
                        fextra.write(f"{idx}\t" + "\t".join(parts[3:]) + "\n")
                    else:
                        fextra.write(f"{idx}\n")
            
            # Run liftOver (parallel or sequential)
            if ncpu > 1 and input_count > 1000:
                # Parallel processing: split into chunks
                if verbose:
                    print(f"🚀 Running liftOver in parallel ({ncpu} workers)...")
                
                from multiprocessing import Pool
                import math
                
                # Read all lines from bed3 file
                with open(bed3_file) as f:
                    all_lines = f.readlines()
                
                # Split into chunks
                chunk_size = math.ceil(len(all_lines) / ncpu)
                chunks = [all_lines[i:i + chunk_size] for i in range(0, len(all_lines), chunk_size)]
                
                # Create chunk files and args
                chunk_args = []
                for i, chunk in enumerate(chunks):
                    chunk_bed = os.path.join(tmpdir, f"chunk_{i}.bed")
                    chunk_lifted = os.path.join(tmpdir, f"lifted_{i}.bed")
                    chunk_unmapped = os.path.join(tmpdir, f"unmapped_{i}.bed")
                    
                    with open(chunk_bed, 'w') as f:
                        f.writelines(chunk)
                    
                    chunk_args.append((chunk_bed, chain_file, liftover_path, min_match, min_blocks, multiple, chunk_lifted, chunk_unmapped))
                
                # Run in parallel
                with Pool(ncpu) as pool:
                    results = pool.map(_liftover_chunk, chunk_args)
                
                if not all(results):
                    return {
                        "status": "error",
                        "lifted": 0,
                        "unmapped": 0,
                        "message": "❌ One or more liftOver chunks failed",
                    }
                
                # Merge lifted chunks
                with open(lifted_bed3, 'w') as fout:
                    for i in range(len(chunks)):
                        chunk_lifted = os.path.join(tmpdir, f"lifted_{i}.bed")
                        if os.path.exists(chunk_lifted):
                            with open(chunk_lifted) as fin:
                                fout.write(fin.read())
                
                # Merge unmapped chunks
                with open(unmapped_tmp, 'w') as fout:
                    for i in range(len(chunks)):
                        chunk_unmapped = os.path.join(tmpdir, f"unmapped_{i}.bed")
                        if os.path.exists(chunk_unmapped):
                            with open(chunk_unmapped) as fin:
                                for line in fin:
                                    if not line.startswith('#'):
                                        fout.write(line)
            else:
                # Sequential processing
                cmd = [
                    liftover_path,
                    f"-minMatch={min_match}",
                ]
                if min_blocks is not None:
                    cmd.append(f"-minBlocks={min_blocks}")
                if multiple:
                    cmd.append("-multiple")
                cmd.extend([
                    bed3_file,
                    chain_file,
                    lifted_bed3,
                    unmapped_tmp
                ])
                
                if verbose:
                    print(f"🔧 Running liftOver...")
                
                result = subprocess.run(cmd, capture_output=True, text=True, check=False)
                
                if result.returncode != 0:
                    print(f"❌ liftOver failed with return code {result.returncode}")
                    print(f"   stderr: {result.stderr}")
                    print(f"   stdout: {result.stdout}")
                    return {
                        "status": "error",
                        "lifted": 0,
                        "unmapped": 0,
                        "message": f"❌ liftOver failed: {result.stderr}",
                        "command": " ".join(cmd)
                    }
            
            # Read extra columns into dict (idx -> extra columns)
            extra_cols = {}
            with open(extra_file) as f:
                for line in f:
                    parts = line.rstrip('\n').split('\t')
                    idx = int(parts[0])
                    extra_cols[idx] = parts[1:] if len(parts) > 1 else []
            
            # Read lifted coordinates and rejoin with extra columns
            # Output will have the TARGET assembly chromosome names (from chain file)
            lifted_count = 0
            with open(lifted_bed3) as fin, open(output_bed, 'w') as fout:
                for line in fin:
                    parts = line.rstrip('\n').split('\t')
                    chrom, start, end = parts[0], parts[1], parts[2]
                    idx = int(parts[3])  # The index we stored as name
                    
                    # Get extra columns for this index
                    extra = extra_cols.get(idx, [])
                    
                    # Write output: chr, start, end, extra columns...
                    if extra:
                        fout.write(f"{chrom}\t{start}\t{end}\t" + "\t".join(extra) + "\n")
                    else:
                        fout.write(f"{chrom}\t{start}\t{end}\n")
                    lifted_count += 1
            
            # Process unmapped - rejoin with extra columns
            unmapped_count = 0
            if os.path.exists(unmapped_tmp):
                with open(unmapped_tmp) as fin, open(unmapped_bed, 'w') as fout:
                    for line in fin:
                        if line.startswith('#'):
                            fout.write(line)
                            continue
                        parts = line.rstrip('\n').split('\t')
                        if len(parts) >= 4:
                            chrom, start, end = parts[0], parts[1], parts[2]
                            idx = int(parts[3])  # The index we stored as name
                            
                            # Get extra columns for this index
                            extra = extra_cols.get(idx, [])
                            
                            # Write output with extra columns
                            if extra:
                                fout.write(f"{chrom}\t{start}\t{end}\t" + "\t".join(extra) + "\n")
                            else:
                                fout.write(f"{chrom}\t{start}\t{end}\n")
                            unmapped_count += 1
        
        # Check results
        if lifted_count == 0 and unmapped_count == 0:
            print(f"⚠️  WARNING: No peaks were processed!")
            print(f"   Input had {input_count} lines")
            print(f"   This is unexpected - check your input file format")
            
            return {
                "status": "warning",
                "lifted": 0,
                "unmapped": 0,
                "output_file": output_bed,
                "unmapped_file": unmapped_bed,
                "message": f"⚠️  No peaks processed.",
            }
        
        if verbose:
            pct = lifted_count / (lifted_count + unmapped_count) * 100
            print(f"✅ Successfully lifted {lifted_count:,} peaks ({pct:.1f}%), {unmapped_count:,} unmapped")
        
        return {
            "status": "success",
            "lifted": lifted_count,
            "unmapped": unmapped_count,
            "output_file": output_bed,
            "unmapped_file": unmapped_bed,
            "message": f"✅ Lifted {lifted_count:,} peaks, {unmapped_count:,} unmapped"
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
        }


def liftover_two_step(
    input_bed: str,
    output_bed: str,
    chain_file_1: str,
    chain_file_2: str,
    unmapped_bed: str = None,
    liftover_path: str = "liftOver",
    min_match: float = 0.95,
    min_blocks: float = None,
    multiple: bool = False,
    auto_chr: bool = True,
    verbose: bool = True,
    ncpu: int = 1,
) -> Dict[str, Any]:
    """
    Perform two-step liftover (e.g., calJac1 -> calJac4 -> hg38).
    
    Uses a temporary intermediate file between the two liftover steps.
    
    Args:
        input_bed: Path to input BED file
        output_bed: Path to output BED file (final result)
        chain_file_1: Path to first chain file (e.g., calJac1ToCalJac4)
        chain_file_2: Path to second chain file (e.g., calJac4ToHg38)
        unmapped_bed: Path for unmapped regions (combined from both steps)
        liftover_path: Path to liftOver executable
        min_match: Minimum match ratio (0.0-1.0) for both steps
        min_blocks: Minimum ratio of alignment blocks that must map (optional)
        multiple: Allow multiple output regions for a single input (default False)
        auto_chr: Automatically detect and fix chr prefix mismatch
        verbose: Print detailed progress
        ncpu: Number of parallel workers (default 1)
    
    Returns:
        Dictionary with status, counts for each step, and file paths
    """
    import tempfile
    
    if verbose:
        print(f"🔄 Two-step liftover:")
        print(f"   Input: {os.path.basename(input_bed)}")
        print(f"   Step 1: {os.path.basename(chain_file_1)}")
        print(f"   Step 2: {os.path.basename(chain_file_2)}")
        print(f"   Output: {os.path.basename(output_bed)}")
        if ncpu > 1:
            print(f"   Workers: {ncpu}")
    
    # Set up unmapped file path
    if unmapped_bed is None:
        unmapped_bed = output_bed.replace('.bed', '.unmapped.bed')
    
    # Create temporary files for intermediate step
    with tempfile.TemporaryDirectory() as tmpdir:
        intermediate_bed = os.path.join(tmpdir, "intermediate.bed")
        unmapped_step1 = os.path.join(tmpdir, "unmapped_step1.bed")
        unmapped_step2 = os.path.join(tmpdir, "unmapped_step2.bed")
        
        # Step 1: First liftover
        if verbose:
            print(f"\n📍 Step 1: {os.path.basename(chain_file_1)}")
        
        result1 = liftover_peaks(
            input_bed=input_bed,
            output_bed=intermediate_bed,
            chain_file=chain_file_1,
            unmapped_bed=unmapped_step1,
            liftover_path=liftover_path,
            min_match=min_match,
            min_blocks=min_blocks,
            multiple=multiple,
            auto_chr=auto_chr,
            verbose=verbose,
            ncpu=ncpu,
        )
        
        if result1["status"] == "error":
            return {
                "status": "error",
                "step1": result1,
                "step2": None,
                "message": f"❌ Step 1 failed: {result1['message']}"
            }
        
        # Step 2: Second liftover
        if verbose:
            print(f"\n📍 Step 2: {os.path.basename(chain_file_2)}")
        
        result2 = liftover_peaks(
            input_bed=intermediate_bed,
            output_bed=output_bed,
            chain_file=chain_file_2,
            unmapped_bed=unmapped_step2,
            liftover_path=liftover_path,
            min_match=min_match,
            min_blocks=min_blocks,
            multiple=multiple,
            auto_chr=auto_chr,
            verbose=verbose,
            ncpu=ncpu,
        )
        
        if result2["status"] == "error":
            return {
                "status": "error",
                "step1": result1,
                "step2": result2,
                "message": f"❌ Step 2 failed: {result2['message']}"
            }
        
        # Combine unmapped files
        with open(unmapped_bed, 'w') as fout:
            fout.write(f"# Unmapped from step 1 ({os.path.basename(chain_file_1)})\n")
            if os.path.exists(unmapped_step1):
                with open(unmapped_step1) as f:
                    for line in f:
                        if not line.startswith('#'):
                            fout.write(line)
            
            fout.write(f"\n# Unmapped from step 2 ({os.path.basename(chain_file_2)})\n")
            if os.path.exists(unmapped_step2):
                with open(unmapped_step2) as f:
                    for line in f:
                        if not line.startswith('#'):
                            fout.write(line)
        
        # Calculate overall statistics
        original_count = result1["lifted"] + result1["unmapped"]
        final_lifted = result2["lifted"]
        total_unmapped = result1["unmapped"] + result2["unmapped"]
        
        if original_count > 0:
            overall_pct = final_lifted / original_count * 100
        else:
            overall_pct = 0
        
        if verbose:
            print(f"\n✅ Two-step liftover complete:")
            print(f"   Original peaks: {original_count:,}")
            print(f"   After step 1:   {result1['lifted']:,} ({result1['lifted']/original_count*100:.1f}%)" if original_count > 0 else f"   After step 1:   {result1['lifted']:,}")
            print(f"   After step 2:   {final_lifted:,} ({overall_pct:.1f}%)")
            print(f"   Total unmapped: {total_unmapped:,}")
        
        return {
            "status": "success",
            "original": original_count,
            "lifted": final_lifted,
            "unmapped_step1": result1["unmapped"],
            "unmapped_step2": result2["unmapped"],
            "total_unmapped": total_unmapped,
            "step1": result1,
            "step2": result2,
            "output_file": output_bed,
            "unmapped_file": unmapped_bed,
            "message": f"✅ Lifted {final_lifted:,}/{original_count:,} peaks ({overall_pct:.1f}%) through two steps"
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

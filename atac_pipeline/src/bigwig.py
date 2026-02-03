"""
BigWig generation from fragment files.

Inspired by scatac_fragment_tools (aertslab/scatac_fragment_tools).
Creates genome coverage bigWig files from ATAC-seq fragment files.
"""

from __future__ import annotations

import gzip
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Generator, Literal

import numpy as np

try:
    import numba
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

try:
    import polars as pl
    HAS_POLARS = True
except ImportError:
    HAS_POLARS = False

try:
    import pyBigWig
    HAS_PYBIGWIG = True
except ImportError:
    HAS_PYBIGWIG = False

try:
    import pybigtools
    HAS_PYBIGTOOLS = True
except ImportError:
    HAS_PYBIGTOOLS = False


def get_chromosome_sizes(chrom_sizes_file: str) -> dict[str, int]:
    """
    Read chromosome sizes from a file.
    
    Parameters
    ----------
    chrom_sizes_file : str
        Path to chromosome sizes file (*.chrom.sizes or *.fa.fai).
        
    Returns
    -------
    dict[str, int]
        Dictionary mapping chromosome names to sizes.
    """
    chrom_sizes = {}
    
    with open(chrom_sizes_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    
    return chrom_sizes


def read_fragments_to_df(
    fragments_file: str,
    min_columns: int = 4,
) -> "pl.DataFrame":
    """
    Read fragments file to a Polars DataFrame.
    
    Parameters
    ----------
    fragments_file : str
        Path to fragments file (TSV, optionally gzipped).
    min_columns : int
        Minimum number of columns required.
        
    Returns
    -------
    pl.DataFrame
        DataFrame with columns: Chromosome, Start, End, Name (barcode), Score (optional).
    """
    if not HAS_POLARS:
        raise ImportError("polars is required for reading fragments. Install with: pip install polars")
    
    # Detect if gzipped
    open_fn = gzip.open if fragments_file.endswith('.gz') else open
    
    # Get column count from first data line
    column_count = 0
    skip_rows = 0
    with open_fn(fragments_file, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                skip_rows += 1
                continue
            column_count = len(line.split('\t'))
            break
    
    if column_count < min_columns:
        raise ValueError(f"Fragments file needs at least {min_columns} columns, got {column_count}")
    
    column_names = ['Chromosome', 'Start', 'End', 'Name', 'Score'][:column_count]
    
    df = pl.read_csv(
        fragments_file,
        separator='\t',
        has_header=False,
        skip_rows=skip_rows,
        new_columns=column_names,
        schema_overrides={
            'Chromosome': pl.Categorical,
            'Start': pl.Int32,
            'End': pl.Int32,
            'Name': pl.Categorical,
        },
    )
    
    # If no score column or score is '.', compute from duplicates
    if 'Score' not in df.columns or df.schema.get('Score') == pl.Utf8:
        df = df.group_by(['Chromosome', 'Start', 'End', 'Name']).agg(
            pl.len().cast(pl.Int32).alias('Score')
        )
    else:
        df = df.with_columns(pl.col('Score').cast(pl.Int32))
    
    return df


if HAS_NUMBA:
    @numba.njit
    def _calculate_depth(chrom_size: int, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
        """Calculate depth per basepair using numba for speed."""
        depth = np.zeros(chrom_size, dtype=np.uint32)
        for i in range(starts.shape[0]):
            depth[starts[i]:ends[i]] += np.uint32(1)
        return depth
    
    @numba.njit
    def _collapse_consecutive(X: np.ndarray):
        """Collapse consecutive identical values into ranges."""
        n = X.shape[0]
        idx = np.empty(n + 1, dtype=np.uint32)
        idx[0] = 0
        j = 1
        
        for i in range(1, n):
            if X[i - 1] != X[i]:
                idx[j] = i
                j += 1
        
        idx[j] = n
        values = X[idx[:j]].astype(np.float32)
        lengths = idx[1:j + 1] - idx[:j]
        
        return idx[:j].copy(), values, lengths
else:
    def _calculate_depth(chrom_size: int, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
        """Calculate depth per basepair (pure numpy fallback)."""
        depth = np.zeros(chrom_size, dtype=np.uint32)
        for start, end in zip(starts, ends):
            depth[start:end] += 1
        return depth
    
    def _collapse_consecutive(X: np.ndarray):
        """Collapse consecutive identical values (pure numpy fallback)."""
        if len(X) == 0:
            return np.array([], dtype=np.uint32), np.array([], dtype=np.float32), np.array([], dtype=np.uint32)
        
        # Find where consecutive values differ
        diff_mask = np.concatenate([[True], X[1:] != X[:-1]])
        idx = np.where(diff_mask)[0].astype(np.uint32)
        values = X[idx].astype(np.float32)
        
        # Calculate lengths
        idx_with_end = np.concatenate([idx, [len(X)]])
        lengths = (idx_with_end[1:] - idx_with_end[:-1]).astype(np.uint32)
        
        return idx, values, lengths


def fragments_to_coverage(
    fragments_df: "pl.DataFrame",
    chrom_sizes: dict[str, int],
    normalize: bool = True,
    scaling_factor: float = 1.0,
    cut_sites: bool = False,
    verbose: bool = False,
) -> Generator[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray], None, None]:
    """
    Calculate genome coverage from fragments.
    
    Parameters
    ----------
    fragments_df : pl.DataFrame
        Polars DataFrame with fragment data.
    chrom_sizes : dict[str, int]
        Chromosome sizes.
    normalize : bool
        Normalize by RPM (reads per million).
    scaling_factor : float
        Additional scaling factor.
    cut_sites : bool
        Use Tn5 cut sites instead of full fragments.
    verbose : bool
        Print progress.
        
    Yields
    ------
    tuple
        (chroms, starts, ends, values) arrays for each chromosome.
    """
    # Initialize depth arrays
    chrom_arrays = {chrom: np.zeros(size, dtype=np.uint32) for chrom, size in chrom_sizes.items()}
    n_fragments = 0
    
    # Partition by chromosome
    per_chrom = {
        str(chrom): df 
        for (chrom,), df in fragments_df.partition_by(['Chromosome'], as_dict=True).items()
    }
    
    if verbose:
        print(f"Processing {len(per_chrom)} chromosomes...")
    
    # Calculate depth per chromosome
    for chrom, chrom_df in per_chrom.items():
        if chrom not in chrom_sizes:
            if verbose:
                print(f"  Skipping {chrom} (not in chrom sizes)")
            continue
        
        starts, ends = chrom_df.select(['Start', 'End']).to_numpy().T
        
        if cut_sites:
            # Convert to Tn5 cut sites (1bp at start and end of each fragment)
            starts, ends = (
                np.hstack((starts, ends - 1)),
                np.hstack((starts + 1, ends)),
            )
        
        chrom_arrays[chrom] = _calculate_depth(chrom_sizes[chrom], starts, ends)
        n_fragments += chrom_df.height
    
    if n_fragments == 0:
        return
    
    rpm_factor = n_fragments / 1_000_000.0
    
    if verbose:
        print(f"Total fragments: {n_fragments:,}")
        print("Collapsing and writing coverage...")
    
    # Yield coverage per chromosome
    for chrom in chrom_sizes:
        if chrom not in chrom_arrays:
            continue
        
        idx, values, lengths = _collapse_consecutive(chrom_arrays[chrom])
        non_zero = np.flatnonzero(values)
        
        if non_zero.shape[0] == 0:
            continue
        
        chroms = np.repeat(chrom, len(non_zero))
        starts_out = idx[non_zero]
        ends_out = idx[non_zero] + lengths[non_zero]
        values_out = values[non_zero]
        
        if normalize:
            values_out = values_out / rpm_factor * scaling_factor
        elif scaling_factor != 1.0:
            values_out *= scaling_factor
        
        yield chroms, starts_out, ends_out, values_out


def fragments_to_bigwig(
    fragments_file: str,
    chrom_sizes_file: str,
    output_file: str,
    normalize: bool = True,
    scaling_factor: float = 1.0,
    cut_sites: bool = True,
    chrom_prefix: str | None = None,
    writer: Literal['pybigwig', 'pybigtools'] = 'pybigwig',
    verbose: bool = False,
) -> dict:
    """
    Convert fragment file to bigWig.
    
    Parameters
    ----------
    fragments_file : str
        Path to fragments file (TSV, optionally gzipped).
    chrom_sizes_file : str
        Path to chromosome sizes file.
    output_file : str
        Output bigWig file path.
    normalize : bool
        Normalize coverage by RPM.
    scaling_factor : float
        Additional scaling factor.
    cut_sites : bool
        Use Tn5 cut sites (recommended for ATAC-seq).
    chrom_prefix : str, optional
        Prefix to add to chromosome names.
    writer : str
        BigWig writer: 'pybigwig' or 'pybigtools'.
    verbose : bool
        Print progress.
        
    Returns
    -------
    dict
        Result with status and metadata.
    """
    if not HAS_POLARS:
        raise ImportError("polars is required. Install with: pip install polars")
    
    if writer == 'pybigwig' and not HAS_PYBIGWIG:
        raise ImportError("pyBigWig is required. Install with: pip install pyBigWig")
    elif writer == 'pybigtools' and not HAS_PYBIGTOOLS:
        raise ImportError("pybigtools is required. Install with: pip install pybigtools")
    
    result = {
        'input': fragments_file,
        'output': output_file,
        'status': 'success',
        'n_fragments': 0,
    }
    
    try:
        # Read data
        if verbose:
            print(f"Reading fragments: {fragments_file}")
        
        chrom_sizes = get_chromosome_sizes(chrom_sizes_file)
        fragments_df = read_fragments_to_df(fragments_file)
        
        # Add chromosome prefix if needed
        if chrom_prefix:
            fragments_df = fragments_df.with_columns(
                (pl.lit(chrom_prefix) + pl.col('Chromosome').cast(pl.Utf8))
                .cast(pl.Categorical)
                .alias('Chromosome')
            )
        
        result['n_fragments'] = fragments_df.height
        
        if verbose:
            print(f"Loaded {fragments_df.height:,} fragments")
        
        # Generate coverage
        coverage_iter = fragments_to_coverage(
            fragments_df=fragments_df,
            chrom_sizes=chrom_sizes,
            normalize=normalize,
            scaling_factor=scaling_factor,
            cut_sites=cut_sites,
            verbose=verbose,
        )
        
        # Write bigWig
        os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
        
        if writer == 'pybigwig':
            with pyBigWig.open(output_file, 'wb') as bw:
                bw.addHeader(list(chrom_sizes.items()))
                
                for chroms, starts, ends, values in coverage_iter:
                    bw.addEntries(
                        chroms=chroms.tolist(),
                        starts=starts.astype(int).tolist(),
                        ends=ends.astype(int).tolist(),
                        values=values.tolist(),
                    )
        else:  # pybigtools
            def entry_iter():
                for chroms, starts, ends, values in coverage_iter:
                    chrom = chroms[0]
                    for s, e, v in zip(starts.tolist(), ends.tolist(), values.tolist()):
                        yield chrom, s, e, v
            
            bw = pybigtools.open(output_file, 'w')
            bw.write(chrom_sizes, entry_iter())
        
        if verbose:
            print(f"✅ Written: {output_file}")
        
    except Exception as e:
        result['status'] = 'error'
        result['error'] = str(e)
        if verbose:
            print(f"❌ Error: {e}")
    
    return result


def fragments_to_bigwig_cli(
    fragments_file: str,
    chrom_sizes_file: str,
    output_file: str,
    cut_sites: bool = True,
    normalize: bool = True,
) -> dict:
    """
    Create bigWig using scatac_fragment_tools CLI (if available).
    
    Falls back to Python implementation if CLI not available.
    
    Parameters
    ----------
    fragments_file : str
        Input fragments file.
    chrom_sizes_file : str
        Chromosome sizes file.
    output_file : str
        Output bigWig file.
    cut_sites : bool
        Use cut sites mode (-x flag).
    normalize : bool
        Normalize coverage.
        
    Returns
    -------
    dict
        Result dictionary.
    """
    # Check if scatac_fragment_tools is available
    try:
        result = subprocess.run(
            ['scatac_fragment_tools', '--version'],
            capture_output=True,
            text=True,
        )
        has_cli = result.returncode == 0
    except FileNotFoundError:
        has_cli = False
    
    if has_cli:
        # Use CLI
        cmd = [
            'scatac_fragment_tools', 'bigwig',
            '-i', fragments_file,
            '-c', chrom_sizes_file,
            '-o', output_file,
        ]
        if cut_sites:
            cmd.append('-x')
        if normalize:
            cmd.append('-n')
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        return {
            'input': fragments_file,
            'output': output_file,
            'status': 'success' if result.returncode == 0 else 'error',
            'method': 'cli',
            'error': result.stderr if result.returncode != 0 else None,
        }
    else:
        # Fall back to Python implementation
        return fragments_to_bigwig(
            fragments_file=fragments_file,
            chrom_sizes_file=chrom_sizes_file,
            output_file=output_file,
            normalize=normalize,
            cut_sites=cut_sites,
        )


def process_all_fragments_to_bigwig(
    input_dir: str,
    output_dir: str,
    chrom_sizes_file: str,
    pattern: str = '*.tsv.gz',
    cut_sites: bool = True,
    normalize: bool = True,
    max_workers: int = 4,
    verbose: bool = False,
) -> list[dict]:
    """
    Convert all fragment files in a directory to bigWig.
    
    Parameters
    ----------
    input_dir : str
        Directory with fragment files.
    output_dir : str
        Output directory for bigWig files.
    chrom_sizes_file : str
        Chromosome sizes file.
    pattern : str
        Glob pattern for fragment files.
    cut_sites : bool
        Use Tn5 cut sites.
    normalize : bool
        Normalize coverage.
    max_workers : int
        Number of parallel workers.
    verbose : bool
        Print progress.
        
    Returns
    -------
    list[dict]
        Results for each file.
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    fragment_files = list(input_path.glob(pattern))
    
    if not fragment_files:
        print(f"No fragment files matching '{pattern}' in {input_dir}")
        return []
    
    print(f"Found {len(fragment_files)} fragment files")
    
    results = []
    
    # Process sequentially (bigWig writing is I/O bound and memory intensive)
    for frag_file in fragment_files:
        output_file = output_path / frag_file.name.replace('.tsv.gz', '.bw').replace('.tsv', '.bw')
        
        if verbose:
            print(f"\nProcessing: {frag_file.name}")
        
        result = fragments_to_bigwig(
            fragments_file=str(frag_file),
            chrom_sizes_file=chrom_sizes_file,
            output_file=str(output_file),
            normalize=normalize,
            cut_sites=cut_sites,
            verbose=verbose,
        )
        results.append(result)
        
        status = "✅" if result['status'] == 'success' else "❌"
        print(f"{status} {frag_file.name}")
    
    # Summary
    successful = sum(1 for r in results if r['status'] == 'success')
    print(f"\nCompleted: {successful}/{len(results)} files")
    
    return results


# Convenience wrapper for the most common use case
def create_bigwig(
    fragments: str,
    chromsizes: str,
    output: str,
    cut_sites: bool = True,
    normalize: bool = True,
    verbose: bool = False,
) -> dict:
    """
    Simple wrapper to create a bigWig from fragments.
    
    Parameters
    ----------
    fragments : str
        Fragment file path.
    chromsizes : str
        Chromosome sizes file.
    output : str
        Output bigWig path.
    cut_sites : bool
        Use Tn5 cut sites (recommended for ATAC-seq).
    normalize : bool
        Normalize by RPM.
    verbose : bool
        Print progress.
        
    Returns
    -------
    dict
        Result dictionary.
        
    Example
    -------
    >>> from src.bigwig import create_bigwig
    >>> result = create_bigwig(
    ...     fragments='sample.fragments.tsv.gz',
    ...     chromsizes='hg38.chrom.sizes',
    ...     output='sample.bw',
    ... )
    """
    return fragments_to_bigwig_cli(
        fragments_file=fragments,
        chrom_sizes_file=chromsizes,
        output_file=output,
        cut_sites=cut_sites,
        normalize=normalize,
    )

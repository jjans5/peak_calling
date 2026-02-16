"""
Quantification module for ATAC-seq peak analysis.

Quantify fragment files, bigWigs, or Tn5 insertion BED files over a set of peaks.
Supports parallel processing and memory-efficient chunked operations.

Public API
----------
quantify : Quantify a single file over peaks.
quantify_matrix : Quantify multiple files into a peaks × samples matrix.
save_matrix : Save a quantification matrix to disk.
load_matrix : Load a quantification matrix from disk.
"""

import os
import tempfile
import subprocess
from pathlib import Path
from typing import List, Dict, Optional, Union, Literal
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

from .utils import load_peaks, clean_sample_name

# Optional imports
try:
    import pyBigWig
    _HAS_PYBIGWIG = True
except ImportError:
    _HAS_PYBIGWIG = False

try:
    import pysam
    _HAS_PYSAM = True
except ImportError:
    _HAS_PYSAM = False

try:
    import polars as pl
    _HAS_POLARS = True
except ImportError:
    _HAS_POLARS = False


# =============================================================================
# INTERNAL: CHROMOSOME NAME HARMONIZATION
# =============================================================================

def _detect_input_chroms(input_file: str, input_type: str) -> set:
    """Read chromosome names from the first lines of an input file or bigWig header."""
    import gzip

    if input_type == "bigwig":
        if _HAS_PYBIGWIG:
            bw = pyBigWig.open(input_file)
            chroms = set(bw.chroms().keys())
            bw.close()
            return chroms
        return set()

    # fragments or tn5 – read first N lines to collect chrom names
    open_fn = gzip.open if input_file.endswith('.gz') else open
    chroms = set()
    try:
        with open_fn(input_file, 'rt') as fh:
            for i, line in enumerate(fh):
                if i >= 1000:
                    break
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                chroms.add(line.split('\t')[0])
    except Exception:
        pass
    return chroms


def _harmonize_chr_prefix(peaks_df: pd.DataFrame, input_file: str, input_type: str) -> pd.DataFrame:
    """
    Ensure peak chromosome names match the input file's naming convention.

    If peaks use 'chr1' but the input uses '1' (or vice versa), the DataFrame
    is adjusted in-place so that bedtools / pyBigWig lookups succeed.

    Returns the (potentially modified) peaks DataFrame.
    """
    input_chroms = _detect_input_chroms(input_file, input_type)
    if not input_chroms:
        return peaks_df

    peak_chroms = set(peaks_df['Chromosome'].unique())
    input_has_chr = any(c.startswith('chr') for c in input_chroms if not c.startswith('chrUn'))
    peaks_has_chr = any(c.startswith('chr') for c in peak_chroms if not c.startswith('chrUn'))

    if input_has_chr == peaks_has_chr:
        return peaks_df  # already matching

    df = peaks_df.copy()
    if input_has_chr and not peaks_has_chr:
        # peaks lack 'chr' – add it
        df['Chromosome'] = 'chr' + df['Chromosome'].astype(str)
    elif not input_has_chr and peaks_has_chr:
        # peaks have 'chr' – remove it
        df['Chromosome'] = df['Chromosome'].astype(str).str.replace(r'^chr', '', regex=True)

    return df


# =============================================================================
# INTERNAL: LOW-LEVEL QUANTIFICATION
# =============================================================================

def _build_peak_index(peaks_df: pd.DataFrame) -> dict:
    """Build a per-chromosome sorted interval index for fast overlap counting.

    Returns a dict mapping each chromosome to (starts, ends, row_indices),
    where starts/ends are sorted numpy arrays and row_indices map back to
    the original peaks_df row positions.

    RAM: ~20 MB for 900K peaks (3 arrays of int64).
    """
    index = {}
    chroms = peaks_df['Chromosome'].values.astype(str)
    starts = peaks_df['Start'].values.astype(np.int64)
    ends = peaks_df['End'].values.astype(np.int64)

    for chrom in np.unique(chroms):
        mask = chroms == chrom
        rows = np.where(mask)[0]
        s = starts[mask]
        e = ends[mask]
        # Sort by start position for searchsorted
        order = np.argsort(s)
        index[chrom] = (s[order], e[order], rows[order])
    return index


def _count_1bp_hits(positions: np.ndarray, peak_starts: np.ndarray,
                    peak_ends: np.ndarray, peak_rows: np.ndarray,
                    counts: np.ndarray) -> None:
    """Vectorized counting of 1bp positions falling within sorted peak intervals.

    For each position p, find peaks where start <= p < end.
    Uses searchsorted to find the nearest candidate peak, then checks the end
    boundary. For non-overlapping peaks (the common case) this is a single pass.
    A depth loop handles overlapping peaks when present.

    This is the inner kernel -- called once per chromosome per chunk.
    """
    if len(positions) == 0 or len(peak_starts) == 0:
        return

    # idx[i] = number of peaks whose start <= positions[i]
    idx = np.searchsorted(peak_starts, positions, side='right')

    # Depth 1: the most common case (handles all non-overlapping peaks)
    candidate_idx = idx - 1
    valid = candidate_idx >= 0
    if not valid.any():
        return
    ci = candidate_idx[valid]
    vp = positions[valid]
    hits = peak_ends[ci] > vp  # start <= pos is guaranteed by searchsorted
    if hits.any():
        np.add.at(counts, peak_rows[ci[hits]], 1)

    # Depth > 1: only needed when peaks overlap. Check a few more candidates
    # but bail out aggressively since overlapping peaks are rare.
    max_depth = min(20, len(peak_starts))
    for depth in range(2, max_depth + 1):
        candidate_idx = idx - depth
        still_valid = candidate_idx >= 0
        if not still_valid.any():
            break
        ci = candidate_idx[still_valid]
        vp = positions[still_valid]
        in_range = peak_starts[ci] <= vp
        if not in_range.any():
            break
        hits = in_range & (peak_ends[ci] > vp)
        if hits.any():
            np.add.at(counts, peak_rows[ci[hits]], 1)


def _count_fragments_coverage(fragments_file: str, peaks_df: Optional[pd.DataFrame] = None,
                               *, peak_index: Optional[dict] = None,
                               n_peaks: Optional[int] = None) -> np.ndarray:
    """Count fragment/cutsite coverage over peaks.

    Uses polars (preferred, ~2x faster) or pandas to parse the fragment file,
    then groups by chromosome and uses vectorized searchsorted for overlap
    counting. No temp files, near-zero RAM overhead.

    RAM: peak index (~20 MB) + one file's worth of chrom + start columns.
    """
    if peak_index is None:
        peak_index = _build_peak_index(peaks_df)
    if n_peaks is None:
        n_peaks = len(peaks_df) if peaks_df is not None else sum(
            len(v[0]) for v in peak_index.values())
    counts = np.zeros(n_peaks, dtype=np.int64)

    if _HAS_POLARS:
        df = pl.read_csv(
            fragments_file, separator='\t', has_header=False,
            columns=[0, 1], new_columns=['c', 's'],
            schema_overrides={'c': pl.Utf8, 's': pl.Int64},
            comment_prefix='#',
        )
        for (ch_name,), group in df.group_by('c'):
            entry = peak_index.get(ch_name)
            if entry is None:
                continue
            positions = group['s'].to_numpy()
            _count_1bp_hits(positions, entry[0], entry[1], entry[2], counts)
    else:
        CHUNK = 500_000
        reader = pd.read_csv(
            fragments_file, sep='\t', header=None,
            usecols=[0, 1], names=['c', 's'],
            dtype={'c': str, 's': np.int64},
            comment='#', engine='c',
            chunksize=CHUNK,
        )
        for chunk in reader:
            chroms = chunk['c'].values
            starts = chunk['s'].values
            for ch in np.unique(chroms):
                entry = peak_index.get(ch)
                if entry is None:
                    continue
                positions = starts[chroms == ch]
                _count_1bp_hits(positions, entry[0], entry[1], entry[2], counts)

    return counts


def _count_fragments_cutsites(fragments_file: str, peaks_df: Optional[pd.DataFrame] = None,
                               *, peak_index: Optional[dict] = None,
                               n_peaks: Optional[int] = None) -> np.ndarray:
    """Count Tn5 cut sites (both fragment ends) within peaks.

    Extracts both start and end-1 positions from each fragment,
    then uses the same polars/pandas + searchsorted approach.
    """
    if peak_index is None:
        peak_index = _build_peak_index(peaks_df)
    if n_peaks is None:
        n_peaks = len(peaks_df) if peaks_df is not None else sum(
            len(v[0]) for v in peak_index.values())
    counts = np.zeros(n_peaks, dtype=np.int64)

    if _HAS_POLARS:
        df = pl.read_csv(
            fragments_file, separator='\t', has_header=False,
            columns=[0, 1, 2], new_columns=['c', 's', 'e'],
            schema_overrides={'c': pl.Utf8, 's': pl.Int64, 'e': pl.Int64},
            comment_prefix='#',
        )
        for (ch_name,), group in df.group_by('c'):
            entry = peak_index.get(ch_name)
            if entry is None:
                continue
            starts = group['s'].to_numpy()
            ends = group['e'].to_numpy() - 1
            positions = np.concatenate([starts, ends])
            _count_1bp_hits(positions, entry[0], entry[1], entry[2], counts)
    else:
        CHUNK = 500_000
        reader = pd.read_csv(
            fragments_file, sep='\t', header=None,
            usecols=[0, 1, 2], names=['c', 's', 'e'],
            dtype={'c': str, 's': np.int64, 'e': np.int64},
            comment='#', engine='c',
            chunksize=CHUNK,
        )
        for chunk in reader:
            chroms = chunk['c'].values
            starts = chunk['s'].values
            ends = chunk['e'].values - 1
            all_chroms = np.concatenate([chroms, chroms])
            all_positions = np.concatenate([starts, ends])
            for ch in np.unique(all_chroms):
                entry = peak_index.get(ch)
                if entry is None:
                    continue
                positions = all_positions[all_chroms == ch]
                _count_1bp_hits(positions, entry[0], entry[1], entry[2], counts)

    return counts


def _are_fragments_already_cutsites(fragments_file: str, n_check: int = 500) -> bool:
    """Check if a fragment file contains pre-computed cut-sites (all entries 1bp wide)."""
    import gzip
    open_fn = gzip.open if fragments_file.endswith('.gz') else open
    try:
        with open_fn(fragments_file, 'rt') as fh:
            checked = 0
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 3:
                    length = int(parts[2]) - int(parts[1])
                    if length != 1:
                        return False
                    checked += 1
                    if checked >= n_check:
                        break
    except Exception:
        return False
    return checked > 0


def _quantify_fragments(fragments_file, peaks_df, method="coverage",
                         peak_index=None, n_peaks=None):
    """Quantify a fragment file over peaks."""
    if method == "coverage":
        return _count_fragments_coverage(fragments_file, peaks_df,
                                         peak_index=peak_index, n_peaks=n_peaks)
    elif method == "cutsites":
        return _count_fragments_cutsites(fragments_file, peaks_df,
                                         peak_index=peak_index, n_peaks=n_peaks)
    else:
        raise ValueError(f"Unknown method '{method}'. Use 'coverage' or 'cutsites'.")


def _quantify_tn5_bed(tn5_bed, peaks_df, peak_index=None, n_peaks=None):
    """Count Tn5 insertions (single-bp BED) within peaks.
    Uses the same streaming coverage approach (1bp entries)."""
    return _count_fragments_coverage(tn5_bed, peaks_df, peak_index=peak_index,
                                     n_peaks=n_peaks)


def _quantify_bigwig(bigwig_file, peaks_df, stat="mean"):
    """Quantify bigWig signal over peaks."""
    if not _HAS_PYBIGWIG:
        raise ImportError(
            "pyBigWig is required for bigWig quantification. "
            "Install with: pip install pyBigWig"
        )
    bw = pyBigWig.open(bigwig_file)
    values = []
    for _, row in peaks_df.iterrows():
        chrom, start, end = str(row['Chromosome']), int(row['Start']), int(row['End'])
        try:
            val = bw.stats(chrom, start, end, type=stat)[0]
            values.append(val if val is not None else 0.0)
        except RuntimeError:
            values.append(0.0)
    bw.close()
    return np.array(values)


def _dispatch(input_file, peaks_df, input_type, method, stat,
              peak_index=None, n_peaks=None):
    """Route to the correct backend based on input_type."""
    if input_type == "fragments":
        return _quantify_fragments(input_file, peaks_df, method=method,
                                    peak_index=peak_index, n_peaks=n_peaks)
    elif input_type == "tn5":
        return _quantify_tn5_bed(input_file, peaks_df, peak_index=peak_index,
                                 n_peaks=n_peaks)
    elif input_type == "bigwig":
        return _quantify_bigwig(input_file, peaks_df, stat=stat)
    else:
        raise ValueError(
            f"Unknown input_type '{input_type}'. "
            "Use 'fragments', 'tn5', or 'bigwig'."
        )


# =============================================================================
# INTERNAL: PARALLEL WORKERS
# =============================================================================

def _worker_chunk(args):
    """Worker: quantify one chunk of peaks for a single file."""
    chunk_df, input_file, input_type, method, stat = args
    return _dispatch(input_file, chunk_df, input_type, method, stat, peak_index=None)


def _worker_file(args):
    """Worker: quantify one file over all peaks. Returns (name, counts)."""
    (input_file, peaks_df, input_type, method, stat,
     name_pattern, name_replacement, use_stem, peak_index, n_peaks) = args
    counts = _dispatch(input_file, peaks_df, input_type, method, stat,
                       peak_index=peak_index, n_peaks=n_peaks)
    name = clean_sample_name(
        input_file, pattern=name_pattern,
        replacement=name_replacement, use_stem=use_stem,
    )
    return name, counts


# =============================================================================
# PUBLIC: quantify  (single file)
# =============================================================================

def quantify(
    input_file: str,
    peak_file: str,
    input_type: Literal["fragments", "tn5", "bigwig"] = "fragments",
    method: Literal["coverage", "cutsites"] = "coverage",
    stat: Literal["mean", "sum", "max", "min"] = "mean",
    n_chunks: int = 1,
    n_workers: int = 1,
    verbose: bool = False,
) -> pd.Series:
    """
    Quantify a single input file over a set of peaks.

    Parameters
    ----------
    input_file : str
        Path to fragment file, Tn5 insertion BED, or bigWig.
    peak_file : str
        Path to BED file with peaks (≥3 columns; col4 used as peak name).
    input_type : {'fragments', 'tn5', 'bigwig'}
        Type of the input file.
    method : {'coverage', 'cutsites'}
        Only for ``input_type='fragments'``.
        *coverage* – count overlapping fragments.
        *cutsites* – count Tn5 cut sites (both fragment ends).
    stat : {'mean', 'sum', 'max', 'min'}
        Only for ``input_type='bigwig'``.
    n_chunks : int
        Split peaks into this many chunks for parallel processing.
        Useful to speed up quantification of a single large file.
    n_workers : int
        Number of parallel workers (only when ``n_chunks > 1``).
    verbose : bool
        Print progress.

    Returns
    -------
    pd.Series
        Values indexed by peak name.
    """
    peaks_df = load_peaks(peak_file)

    # Auto-detect cutsites fragments (check once, not per-chunk)
    if input_type == "fragments" and method == "cutsites":
        if _are_fragments_already_cutsites(input_file):
            if verbose:
                print("  Auto-detected 1bp cut-site fragments, "
                      "using 'coverage' instead of 'cutsites'")
            method = "coverage"

    # Harmonize chr prefix once on the peaks
    peaks_df = _harmonize_chr_prefix(peaks_df, input_file, input_type)

    # Build peak index once (avoids rebuilding in each worker/call)
    peak_index = _build_peak_index(peaks_df)

    if verbose:
        print(f"Quantifying {Path(input_file).name} ({input_type}) "
              f"over {len(peaks_df):,} peaks...")

    if n_chunks <= 1 or n_workers <= 1:
        counts = _dispatch(input_file, peaks_df, input_type, method, stat,
                           peak_index=peak_index, n_peaks=len(peaks_df))
    else:
        chunks = np.array_split(peaks_df, n_chunks)
        args_list = [(c, input_file, input_type, method, stat) for c in chunks]
        ordered = [None] * len(chunks)
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = {pool.submit(_worker_chunk, a): i
                       for i, a in enumerate(args_list)}
            for fut in as_completed(futures):
                ordered[futures[fut]] = fut.result()
        counts = np.concatenate(ordered)

    return pd.Series(counts, index=peaks_df['Name'].values,
                     name=Path(input_file).stem)


# =============================================================================
# PUBLIC: quantify_matrix  (multiple files → matrix)
# =============================================================================

def quantify_matrix(
    input_files: List[str],
    peak_file: str,
    input_type: Literal["fragments", "tn5", "bigwig"] = "fragments",
    method: Literal["coverage", "cutsites"] = "coverage",
    stat: Literal["mean", "sum", "max", "min"] = "mean",
    n_workers: int = 4,
    name_pattern: Optional[str] = None,
    name_replacement: str = "",
    use_stem: bool = True,
    sample_names: Optional[List[str]] = None,
    output_file: Optional[str] = None,
    output_format: Literal["tsv", "feather", "parquet"] = "feather",
    chunk_size: Optional[int] = None,
    verbose: bool = True,
) -> Optional[pd.DataFrame]:
    """
    Quantify multiple files into a peaks × samples matrix.

    Parameters
    ----------
    input_files : list of str
        Paths to input files (all must be the same *input_type*).
    peak_file : str
        Path to BED file with peaks.
    input_type : {'fragments', 'tn5', 'bigwig'}
        Type of the input files.
    method : {'coverage', 'cutsites'}
        For ``input_type='fragments'`` only.
    stat : {'mean', 'sum', 'max', 'min'}
        For ``input_type='bigwig'`` only.
    n_workers : int
        Number of parallel workers.
    name_pattern : str, optional
        Regex applied via ``re.sub`` to clean sample names from filenames.
    name_replacement : str
        Replacement string for *name_pattern*.
    use_stem : bool
        Use filename stem (without extension) as base sample name.
    sample_names : list of str, optional
        Explicit sample names. Overrides automatic name cleaning.
        Must have the same length as *input_files*.
    output_file : str, optional
        Save matrix to this path. If ``None``, the DataFrame is returned.
    output_format : {'tsv', 'feather', 'parquet'}
        Disk format (only when *output_file* is set).
    chunk_size : int, optional
        Process files in batches of this size and write intermediate results
        to avoid holding the full matrix in memory. Requires *output_file*.
    verbose : bool
        Print progress.

    Returns
    -------
    pd.DataFrame or None
        Peaks × samples matrix, or ``None`` when written to *output_file*.
    """
    peaks_df = load_peaks(peak_file)
    n_peaks, n_files = len(peaks_df), len(input_files)

    if sample_names is not None and len(sample_names) != n_files:
        raise ValueError(
            f"sample_names length ({len(sample_names)}) != "
            f"input_files length ({n_files})"
        )

    # Auto-detect cutsites fragments (check first file once)
    if input_type == "fragments" and method == "cutsites" and n_files > 0:
        if _are_fragments_already_cutsites(input_files[0]):
            method = "coverage"
            if verbose:
                print("  Auto-detected 1bp cut-site fragments, "
                      "using 'coverage' instead of 'cutsites'")

    # Harmonize chr prefix once on the peaks (sample first file)
    if n_files > 0:
        peaks_df = _harmonize_chr_prefix(peaks_df, input_files[0], input_type)

    if verbose:
        print(f"Quantifying {n_files} files over {n_peaks:,} peaks")
        extra = (f" | method: {method}" if input_type == "fragments"
                 else f" | stat: {stat}" if input_type == "bigwig" else "")
        print(f"   Input type: {input_type}{extra}")
        print(f"   Workers: {n_workers}")

    # Build peak index once — passed to all workers to avoid rebuilding
    peak_index = _build_peak_index(peaks_df)

    # --- chunked mode (memory-efficient) ---
    if chunk_size is not None and output_file is not None:
        _quantify_matrix_chunked(
            input_files, peaks_df, input_type, method, stat,
            n_workers, name_pattern, name_replacement, use_stem,
            sample_names, output_file, output_format, chunk_size, verbose,
        )
        return None

    # --- standard mode ---
    # For fragments/tn5, workers only need peak_index + n_peaks (skip pickling peaks_df).
    # For bigwig, workers need the actual peaks_df.
    worker_df = None if input_type in ("fragments", "tn5") else peaks_df
    args_list = [
        (f, worker_df, input_type, method, stat,
         name_pattern, name_replacement, use_stem, peak_index, n_peaks)
        for f in input_files
    ]

    results = {}
    completed = 0

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(_worker_file, a): i
                   for i, a in enumerate(args_list)}
        for fut in as_completed(futures):
            idx = futures[fut]
            auto_name, counts = fut.result()
            col_name = sample_names[idx] if sample_names else auto_name
            results[col_name] = counts
            completed += 1
            if verbose:
                print(f"   [{completed}/{n_files}] {col_name}")

    df = pd.DataFrame(results, index=peaks_df['Name'].values)
    df = df.reindex(sorted(df.columns), axis=1)

    if verbose:
        print(f"\nMatrix shape: {df.shape}")

    if output_file is not None:
        save_matrix(df, output_file, output_format, verbose=verbose)
        return None
    return df


# =============================================================================
# INTERNAL: CHUNKED PROCESSING
# =============================================================================

def _quantify_matrix_chunked(
    input_files, peaks_df, input_type, method, stat,
    n_workers, name_pattern, name_replacement, use_stem,
    sample_names, output_file, output_format, chunk_size, verbose,
):
    """Process files in batches; write intermediate chunks to limit memory."""
    import shutil

    n_files = len(input_files)
    n_chunks = (n_files + chunk_size - 1) // chunk_size

    n_peaks = len(peaks_df)

    # Build peak index once for all batches
    peak_index = _build_peak_index(peaks_df)
    worker_df = None if input_type in ("fragments", "tn5") else peaks_df

    if verbose:
        print(f"   Processing in {n_chunks} chunk(s) of {chunk_size} files max")

    temp_dir = tempfile.mkdtemp(prefix="quant_chunks_")
    temp_paths = []

    try:
        for batch_start in range(0, n_files, chunk_size):
            batch_files = input_files[batch_start:batch_start + chunk_size]
            batch_idx = batch_start // chunk_size + 1

            if verbose:
                print(f"\n   Chunk {batch_idx}/{n_chunks} "
                      f"({len(batch_files)} files)...")

            args_list = [
                (f, worker_df, input_type, method, stat,
                 name_pattern, name_replacement, use_stem, peak_index, n_peaks)
                for f in batch_files
            ]

            chunk_results = {}
            with ProcessPoolExecutor(max_workers=n_workers) as pool:
                futures = {pool.submit(_worker_file, a): batch_start + i
                           for i, a in enumerate(args_list)}
                for fut in as_completed(futures):
                    idx = futures[fut]
                    auto_name, counts = fut.result()
                    col_name = sample_names[idx] if sample_names else auto_name
                    chunk_results[col_name] = counts

            chunk_df = pd.DataFrame(chunk_results, index=peaks_df['Name'].values)
            chunk_path = os.path.join(temp_dir, f"chunk_{batch_idx:04d}.feather")
            chunk_df.reset_index().to_feather(chunk_path)
            temp_paths.append(chunk_path)

            if verbose:
                print(f"      Chunk {batch_idx} saved "
                      f"({chunk_df.shape[1]} samples)")

        # merge
        if verbose:
            print(f"\n   Merging {len(temp_paths)} chunk(s)...")

        merged = None
        for tp in temp_paths:
            part = pd.read_feather(tp).set_index('index')
            merged = part if merged is None else pd.concat([merged, part], axis=1)

        merged = merged.reindex(sorted(merged.columns), axis=1)
        save_matrix(merged, output_file, output_format, verbose=verbose)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


# =============================================================================
# PUBLIC: MATRIX I/O
# =============================================================================

def save_matrix(
    df: pd.DataFrame,
    output_file: str,
    output_format: Literal["tsv", "feather", "parquet"] = "feather",
    verbose: bool = True,
) -> str:
    """
    Save a quantification matrix to disk.

    Parameters
    ----------
    df : pd.DataFrame
        Matrix (peaks × samples).
    output_file : str
        Output path (extension adjusted to match *output_format*).
    output_format : {'tsv', 'feather', 'parquet'}
        - **feather** – fast columnar I/O (recommended).
        - **parquet** – compressed columnar (best for very large files).
        - **tsv** – human-readable tab-separated text.
    verbose : bool
        Print the saved path and file size.

    Returns
    -------
    str
        Actual path written (with correct extension).
    """
    base = os.path.splitext(output_file)[0]

    if output_format == "tsv":
        out = base + ".tsv"
        df.to_csv(out, sep='\t')
    elif output_format == "feather":
        out = base + ".feather"
        df.reset_index().to_feather(out)
    elif output_format == "parquet":
        out = base + ".parquet"
        df.reset_index().to_parquet(out, engine='pyarrow')
    else:
        raise ValueError(
            f"Unknown format '{output_format}'. "
            "Use 'tsv', 'feather', or 'parquet'."
        )

    if verbose:
        size_mb = os.path.getsize(out) / (1024 * 1024)
        print(f"   Saved to {out} ({size_mb:.1f} MB)")
    return out


def load_matrix(input_file: str) -> pd.DataFrame:
    """
    Load a quantification matrix from disk.

    Automatically detects format from the file extension
    (.tsv / .txt, .feather, .parquet).

    Parameters
    ----------
    input_file : str
        Path to the matrix file.

    Returns
    -------
    pd.DataFrame
        Peaks × samples matrix.
    """
    ext = os.path.splitext(input_file)[1].lower()

    if ext in ('.tsv', '.txt'):
        df = pd.read_csv(input_file, sep='\t', index_col=0)
    elif ext == '.feather':
        df = pd.read_feather(input_file)
        if 'index' in df.columns:
            df = df.set_index('index')
    elif ext == '.parquet':
        df = pd.read_parquet(input_file)
        if 'index' in df.columns:
            df = df.set_index('index')
    else:
        raise ValueError(f"Unknown file extension '{ext}'.")
    return df

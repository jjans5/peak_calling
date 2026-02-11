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

def _write_temp_bed(peaks_df: pd.DataFrame, suffix: str = ".bed") -> str:
    """Write peaks to temporary BED file."""
    fd, path = tempfile.mkstemp(suffix=suffix)
    peaks_df[['Chromosome', 'Start', 'End', 'Name']].to_csv(
        path, sep='\t', header=False, index=False
    )
    os.close(fd)
    return path


def _count_fragments_coverage(fragments_file: str, peaks_df: pd.DataFrame) -> np.ndarray:
    """Count fragment coverage over peaks (bedtools intersect -c)."""
    temp_peaks = _write_temp_bed(peaks_df)
    try:
        cmd = f"bedtools intersect -a {temp_peaks} -b {fragments_file} -c"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"bedtools error: {result.stderr}")
        counts = [int(line.split('\t')[-1]) for line in result.stdout.strip().split('\n') if line]
        return np.array(counts)
    finally:
        os.unlink(temp_peaks)


def _count_fragments_cutsites(fragments_file: str, peaks_df: pd.DataFrame) -> np.ndarray:
    """Count Tn5 cut sites (fragment ends) within peaks."""
    temp_peaks = _write_temp_bed(peaks_df)
    fd, temp_cuts = tempfile.mkstemp(suffix=".bed")
    os.close(fd)
    try:
        # Extract both ends of each fragment as cut sites
        awk_cmd = (
            f"""awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$2+1; print $1,$3-1,$3}}' """
            f"""{fragments_file} > {temp_cuts}"""
        )
        subprocess.run(awk_cmd, shell=True, check=True)

        cmd = f"bedtools intersect -a {temp_peaks} -b {temp_cuts} -c"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"bedtools error: {result.stderr}")
        counts = [int(line.split('\t')[-1]) for line in result.stdout.strip().split('\n') if line]
        return np.array(counts)
    finally:
        os.unlink(temp_peaks)
        if os.path.exists(temp_cuts):
            os.unlink(temp_cuts)


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


def _quantify_fragments(fragments_file, peaks_df, method="coverage"):
    """Quantify a fragment file over peaks.

    If ``method='cutsites'`` but the fragments are already 1bp cut-sites,
    automatically falls back to ``coverage`` to avoid double-processing.
    """
    if method == "cutsites" and _are_fragments_already_cutsites(fragments_file):
        print(f"  Auto-detected 1bp cut-site fragments, using 'coverage' instead of 'cutsites'")
        method = "coverage"

    if method == "coverage":
        return _count_fragments_coverage(fragments_file, peaks_df)
    elif method == "cutsites":
        return _count_fragments_cutsites(fragments_file, peaks_df)
    else:
        raise ValueError(f"Unknown method '{method}'. Use 'coverage' or 'cutsites'.")


def _quantify_tn5_bed(tn5_bed, peaks_df):
    """Count Tn5 insertions (single-bp BED) within peaks."""
    temp_peaks = _write_temp_bed(peaks_df)
    try:
        cmd = f"bedtools intersect -a {temp_peaks} -b {tn5_bed} -c"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"bedtools error: {result.stderr}")
        counts = [int(line.split('\t')[-1]) for line in result.stdout.strip().split('\n') if line]
        return np.array(counts)
    finally:
        os.unlink(temp_peaks)


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


def _dispatch(input_file, peaks_df, input_type, method, stat):
    """Route to the correct backend based on input_type."""
    # Ensure peak chromosome names match the input file
    peaks_df = _harmonize_chr_prefix(peaks_df, input_file, input_type)

    if input_type == "fragments":
        return _quantify_fragments(input_file, peaks_df, method=method)
    elif input_type == "tn5":
        return _quantify_tn5_bed(input_file, peaks_df)
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
    return _dispatch(input_file, chunk_df, input_type, method, stat)


def _worker_file(args):
    """Worker: quantify one file over all peaks. Returns (name, counts)."""
    (input_file, peaks_df, input_type, method, stat,
     name_pattern, name_replacement, use_stem) = args
    counts = _dispatch(input_file, peaks_df, input_type, method, stat)
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

    if verbose:
        print(f"Quantifying {Path(input_file).name} ({input_type}) "
              f"over {len(peaks_df):,} peaks...")

    if n_chunks <= 1 or n_workers <= 1:
        counts = _dispatch(input_file, peaks_df, input_type, method, stat)
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

    if verbose:
        print(f"Quantifying {n_files} files over {n_peaks:,} peaks")
        extra = (f" | method: {method}" if input_type == "fragments"
                 else f" | stat: {stat}" if input_type == "bigwig" else "")
        print(f"   Input type: {input_type}{extra}")
        print(f"   Workers: {n_workers}")

    # --- chunked mode (memory-efficient) ---
    if chunk_size is not None and output_file is not None:
        _quantify_matrix_chunked(
            input_files, peaks_df, input_type, method, stat,
            n_workers, name_pattern, name_replacement, use_stem,
            sample_names, output_file, output_format, chunk_size, verbose,
        )
        return None

    # --- standard mode ---
    args_list = [
        (f, peaks_df, input_type, method, stat,
         name_pattern, name_replacement, use_stem)
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
        print(f"\n📦 Matrix shape: {df.shape}")

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

    if verbose:
        print(f"   📦 Processing in {n_chunks} chunk(s) of ≤{chunk_size} files")

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
                (f, peaks_df, input_type, method, stat,
                 name_pattern, name_replacement, use_stem)
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
            print(f"\n   🔗 Merging {len(temp_paths)} chunk(s)...")

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

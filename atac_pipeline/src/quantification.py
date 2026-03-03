"""
Quantification module for ATAC-seq peak analysis.

Two quantification backends:

**BigWig-based (recommended)**
    Convert fragments → bigWig once, then quantify using pyBigWig.stats().
    Inspired by CREsted (aertslab/CREsted).  Near-instant per file.

**Fragment-based (legacy)**
    Parse raw fragment / Tn5 files directly with polars/pandas + searchsorted.
    Works without pre-computing bigWigs but is slower on large files.

Public API
----------
quantify_bigwig        : Quantify a single bigWig over peaks.
quantify_bigwig_matrix : Quantify multiple bigWigs into a peaks × samples matrix.
fragments_to_bigwigs   : Convert fragment files → bigWigs in parallel.
quantify               : Quantify a single file (fragments/tn5/bigwig) over peaks.
quantify_matrix        : Quantify multiple files into a peaks × samples matrix.
save_matrix            : Save a quantification matrix to disk.
load_matrix            : Load a quantification matrix from disk.
"""

import os
import re
import gzip
import tempfile
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
    import polars as pl
    _HAS_POLARS = True
except ImportError:
    _HAS_POLARS = False


# =============================================================================
# INTERNAL HELPERS
# =============================================================================

def _clean_sample_name(filepath, pattern=None, replacement="", use_stem=True):
    """Clean up sample name from filepath."""
    if use_stem:
        name = Path(filepath).stem
        for ext in ['.fragments', '.tn5', '.insertions', '.sorted', '.filtered']:
            if name.endswith(ext):
                name = name[:-len(ext)]
    else:
        name = Path(filepath).name
    if pattern:
        name = re.sub(pattern, replacement, name)
    return name


def _detect_chr_convention(bw_chroms: set, peak_chroms: set):
    """Detect whether chr prefix needs adding/removing."""
    bw_has_chr = any(c.startswith('chr') for c in bw_chroms
                     if not c.startswith('chrUn'))
    pk_has_chr = any(c.startswith('chr') for c in peak_chroms
                     if not c.startswith('chrUn'))
    return bw_has_chr, pk_has_chr


def _harmonize_peak_chroms(peaks_df, bw_chroms):
    """Ensure peak chromosome names match the bigWig naming convention."""
    peak_chroms = set(peaks_df['Chromosome'].unique())
    bw_has_chr, pk_has_chr = _detect_chr_convention(bw_chroms, peak_chroms)

    if bw_has_chr == pk_has_chr:
        return peaks_df

    df = peaks_df.copy()
    if bw_has_chr and not pk_has_chr:
        df['Chromosome'] = 'chr' + df['Chromosome'].astype(str)
    elif not bw_has_chr and pk_has_chr:
        df['Chromosome'] = (df['Chromosome'].astype(str)
                            .str.replace(r'^chr', '', regex=True))
    return df


# =============================================================================
# CORE: QUANTIFY BIGWIG OVER PEAKS
# =============================================================================

def _quantify_single_bigwig(bigwig_file, chroms, starts, ends, n_peaks,
                             stat="mean"):
    """
    Quantify one bigWig over peak regions using pyBigWig.stats().

    Groups peaks by chromosome and calls bw.stats() once per chromosome
    which is much faster than row-by-row Python iteration.

    Parameters
    ----------
    bigwig_file : str
        Path to bigWig file.
    chroms, starts, ends : np.ndarray
        Peak coordinates (string, int64, int64).
    n_peaks : int
        Total number of peaks.
    stat : str
        Statistic: 'mean', 'sum', 'max', 'min', 'coverage'.

    Returns
    -------
    np.ndarray
        Array of length n_peaks with quantified values.
    """
    values = np.zeros(n_peaks, dtype=np.float64)
    bw = pyBigWig.open(bigwig_file)
    bw_chroms_dict = bw.chroms()
    bw_chrom_set = set(bw_chroms_dict.keys())

    unique_chroms = np.unique(chroms)
    for chrom in unique_chroms:
        chrom_str = str(chrom)
        if chrom_str not in bw_chrom_set:
            continue

        mask = chroms == chrom
        idxs = np.where(mask)[0]
        ch_starts = starts[mask]
        ch_ends = ends[mask]

        # Clamp ends to chromosome size to avoid pyBigWig errors
        chrom_size = bw_chroms_dict[chrom_str]
        ch_ends = np.minimum(ch_ends, chrom_size)
        # Skip regions where start >= end after clamping
        valid = ch_starts < ch_ends
        if not valid.any():
            continue
        if not valid.all():
            idxs = idxs[valid]
            ch_starts = ch_starts[valid]
            ch_ends = ch_ends[valid]

        # Query one peak at a time (pyBigWig.stats batch mode with lists
        # is only available in some builds; single-call is still fast at C level)
        for i in range(len(idxs)):
            try:
                val = bw.stats(chrom_str, int(ch_starts[i]), int(ch_ends[i]),
                               type=stat, exact=True)
                values[idxs[i]] = val[0] if val[0] is not None else 0.0
            except RuntimeError:
                values[idxs[i]] = 0.0

    bw.close()
    return values


def _worker_quantify_bigwig(args):
    """Worker function for parallel bigWig quantification."""
    (bigwig_file, chroms, starts, ends, n_peaks, stat,
     name_pattern, name_replacement, use_stem) = args
    values = _quantify_single_bigwig(bigwig_file, chroms, starts, ends,
                                      n_peaks, stat)
    name = _clean_sample_name(bigwig_file, pattern=name_pattern,
                               replacement=name_replacement, use_stem=use_stem)
    return name, values


# =============================================================================
# PUBLIC: QUANTIFY SINGLE BIGWIG
# =============================================================================

def quantify_bigwig(
    bigwig_file: str,
    peak_file: str,
    stat: Literal["mean", "sum", "max", "min", "coverage"] = "sum",
    verbose: bool = False,
) -> pd.Series:
    """
    Quantify a single bigWig file over a set of peaks.

    Parameters
    ----------
    bigwig_file : str
        Path to bigWig file.
    peak_file : str
        Path to BED file with peaks (col4 = peak name).
    stat : {'mean', 'sum', 'max', 'min', 'coverage'}
        Summary statistic per peak region.
    verbose : bool
        Print progress.

    Returns
    -------
    pd.Series
        Values indexed by peak name.
    """
    if not _HAS_PYBIGWIG:
        raise ImportError(
            "pyBigWig is required for bigWig quantification. "
            "Install with: pip install pyBigWig"
        )

    peaks_df = load_peaks(peak_file)

    # Harmonize chr naming
    bw = pyBigWig.open(bigwig_file)
    bw_chroms = set(bw.chroms().keys())
    bw.close()
    peaks_df = _harmonize_peak_chroms(peaks_df, bw_chroms)

    chroms = peaks_df['Chromosome'].values.astype(str)
    starts = peaks_df['Start'].values.astype(np.int64)
    ends = peaks_df['End'].values.astype(np.int64)

    if verbose:
        print(f"Quantifying {Path(bigwig_file).name} over "
              f"{len(peaks_df):,} peaks (stat={stat})...")

    values = _quantify_single_bigwig(bigwig_file, chroms, starts, ends,
                                      len(peaks_df), stat)

    if verbose:
        nz = (values > 0).sum()
        print(f"   Non-zero peaks: {nz:,} / {len(peaks_df):,}")

    return pd.Series(values, index=peaks_df['Name'].values,
                     name=Path(bigwig_file).stem)


# =============================================================================
# PUBLIC: QUANTIFY BIGWIG MATRIX
# =============================================================================

def quantify_bigwig_matrix(
    bigwig_files: List[str],
    peak_file: str,
    stat: Literal["mean", "sum", "max", "min", "coverage"] = "sum",
    n_workers: int = 4,
    name_pattern: Optional[str] = None,
    name_replacement: str = "",
    use_stem: bool = True,
    sample_names: Optional[List[str]] = None,
    output_file: Optional[str] = None,
    output_format: Literal["tsv", "feather", "parquet"] = "feather",
    verbose: bool = True,
) -> Optional[pd.DataFrame]:
    """
    Quantify multiple bigWig files into a peaks × samples matrix.

    Uses ProcessPoolExecutor for parallel quantification — each bigWig
    is processed independently.  This is very fast since pyBigWig.stats()
    works at the C level.

    Parameters
    ----------
    bigwig_files : list of str
        Paths to bigWig files.
    peak_file : str
        Path to BED file with peaks.
    stat : {'mean', 'sum', 'max', 'min', 'coverage'}
        Summary statistic per peak region.
    n_workers : int
        Number of parallel workers.
    name_pattern : str, optional
        Regex to clean sample names from filenames.
    name_replacement : str
        Replacement string for name_pattern.
    use_stem : bool
        Use filename stem as base sample name.
    sample_names : list of str, optional
        Explicit sample names (overrides auto-detection).
    output_file : str, optional
        Save matrix to this path.
    output_format : {'tsv', 'feather', 'parquet'}
        Disk format.
    verbose : bool
        Print progress.

    Returns
    -------
    pd.DataFrame or None
        Peaks × samples matrix, or None when written to output_file.
    """
    if not _HAS_PYBIGWIG:
        raise ImportError(
            "pyBigWig is required. Install with: pip install pyBigWig"
        )

    peaks_df = load_peaks(peak_file)
    n_peaks = len(peaks_df)
    n_files = len(bigwig_files)

    if sample_names is not None and len(sample_names) != n_files:
        raise ValueError(
            f"sample_names length ({len(sample_names)}) != "
            f"bigwig_files length ({n_files})"
        )

    # Harmonize chr naming with first bigWig
    if n_files > 0:
        bw = pyBigWig.open(bigwig_files[0])
        bw_chroms = set(bw.chroms().keys())
        bw.close()
        peaks_df = _harmonize_peak_chroms(peaks_df, bw_chroms)

    chroms = peaks_df['Chromosome'].values.astype(str)
    starts = peaks_df['Start'].values.astype(np.int64)
    ends = peaks_df['End'].values.astype(np.int64)

    if verbose:
        print(f"Quantifying {n_files} bigWig files over {n_peaks:,} peaks")
        print(f"   Statistic: {stat}")
        print(f"   Workers: {n_workers}")

    args_list = [
        (f, chroms, starts, ends, n_peaks, stat,
         name_pattern, name_replacement, use_stem)
        for f in bigwig_files
    ]

    results = {}
    completed = 0

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(_worker_quantify_bigwig, a): i
                   for i, a in enumerate(args_list)}
        for fut in as_completed(futures):
            idx = futures[fut]
            auto_name, values = fut.result()
            col_name = sample_names[idx] if sample_names else auto_name
            results[col_name] = values
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
# HELPERS: CUT-SITE DETECTION & SINGLE-FILE CONVERTER
# =============================================================================

def _are_fragments_already_cutsites(fragments_file, n_check=500):
    """Check if a fragment file contains pre-computed cut-sites (all 1bp wide)."""
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
                    if int(parts[2]) - int(parts[1]) != 1:
                        return False
                    checked += 1
                    if checked >= n_check:
                        break
    except Exception:
        return False
    return checked > 0


def _convert_one_fragment(args):
    """Top-level worker for ProcessPoolExecutor (must be picklable)."""
    from .bigwig import fragments_to_bigwig as _frag2bw

    frag_file, chrom_sizes_file, output_dir, cut_sites, normalize = args

    stem = Path(frag_file).name
    for ext in ['.fragments.tsv.gz', '.tsv.gz', '.tsv']:
        if stem.endswith(ext):
            stem = stem[:-len(ext)]
            break
    out_bw = os.path.join(output_dir, stem + '.bw')

    if os.path.exists(out_bw):
        return frag_file, out_bw, 'skipped'

    result = _frag2bw(
        fragments_file=frag_file,
        chrom_sizes_file=chrom_sizes_file,
        output_file=out_bw,
        normalize=normalize,
        cut_sites=cut_sites,
        verbose=False,
    )
    return frag_file, out_bw, result.get('status', 'error')


# =============================================================================
# PUBLIC: CONVERT FRAGMENTS → BIGWIGS
# =============================================================================

def fragments_to_bigwigs(
    fragment_files: List[str],
    chrom_sizes_file: str,
    output_dir: str,
    cut_sites: bool = True,
    normalize: bool = False,
    n_workers: int = 4,
    verbose: bool = True,
) -> Dict[str, str]:
    """
    Convert fragment files to bigWig files in parallel.

    If the input fragments are already 1-bp cut-sites (auto-detected from
    the first file), ``cut_sites`` is automatically set to ``False`` to
    avoid double-expanding.

    Parameters
    ----------
    fragment_files : list of str
        Paths to fragment .tsv.gz files.
    chrom_sizes_file : str
        Path to chrom.sizes file for the genome assembly.
    output_dir : str
        Output directory for bigWig files.
    cut_sites : bool
        Use Tn5 cut-site mode (recommended for ATAC-seq).
        Automatically set to False when input is already cut-sites.
    normalize : bool
        Normalize by RPM.  Set False for raw counts (recommended
        for DESeq2 downstream).
    n_workers : int
        Number of parallel workers.
    verbose : bool
        Print progress.

    Returns
    -------
    dict
        Mapping of input fragment file → output bigWig file path.
    """
    os.makedirs(output_dir, exist_ok=True)

    # --- Auto-detect already-cutsite input ---
    if cut_sites and fragment_files:
        if _are_fragments_already_cutsites(fragment_files[0]):
            cut_sites = False
            if verbose:
                print("   ⚠️  Fragments are already 1-bp cut-sites — "
                      "setting cut_sites=False to avoid double-expansion")

    if verbose:
        print(f"Converting {len(fragment_files)} fragment files → bigWig")
        print(f"   Output:    {output_dir}")
        print(f"   cut_sites: {cut_sites}")
        print(f"   Workers:   {n_workers}")

    args_list = [
        (f, chrom_sizes_file, output_dir, cut_sites, normalize)
        for f in fragment_files
    ]

    mapping = {}
    completed = 0
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(_convert_one_fragment, a): a[0]
                   for a in args_list}
        for fut in as_completed(futures):
            frag, bw, status = fut.result()
            mapping[frag] = bw
            completed += 1
            icon = "✅" if status in ('success', 'skipped') else "❌"
            skip_note = " (cached)" if status == 'skipped' else ""
            if verbose:
                print(f"   [{completed}/{len(fragment_files)}] "
                      f"{icon} {Path(frag).name}{skip_note}")

    if verbose:
        ok = sum(1 for f, b in mapping.items() if os.path.exists(b))
        print(f"\n   {ok}/{len(fragment_files)} bigWig files ready")

    return mapping


# =============================================================================
# FRAGMENT-BASED QUANTIFICATION (legacy – works without bigWigs)
# =============================================================================

def _detect_input_chroms(input_file: str, input_type: str) -> set:
    """Read chromosome names from the first lines of an input file or bigWig header."""
    if input_type == "bigwig":
        if _HAS_PYBIGWIG:
            bw = pyBigWig.open(input_file)
            chroms = set(bw.chroms().keys())
            bw.close()
            return chroms
        return set()

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


def _harmonize_chr_prefix(peaks_df, input_file, input_type):
    """Ensure peak chromosome names match the input file's naming convention."""
    input_chroms = _detect_input_chroms(input_file, input_type)
    if not input_chroms:
        return peaks_df

    peak_chroms = set(peaks_df['Chromosome'].unique())
    input_has_chr = any(c.startswith('chr') for c in input_chroms
                        if not c.startswith('chrUn'))
    peaks_has_chr = any(c.startswith('chr') for c in peak_chroms
                        if not c.startswith('chrUn'))

    if input_has_chr == peaks_has_chr:
        return peaks_df

    df = peaks_df.copy()
    if input_has_chr and not peaks_has_chr:
        df['Chromosome'] = 'chr' + df['Chromosome'].astype(str)
    elif not input_has_chr and peaks_has_chr:
        df['Chromosome'] = (df['Chromosome'].astype(str)
                            .str.replace(r'^chr', '', regex=True))
    return df


def _build_peak_index(peaks_df):
    """Build a per-chromosome sorted interval index for fast overlap counting."""
    index = {}
    chroms = peaks_df['Chromosome'].values.astype(str)
    starts = peaks_df['Start'].values.astype(np.int64)
    ends = peaks_df['End'].values.astype(np.int64)

    for chrom in np.unique(chroms):
        mask = chroms == chrom
        rows = np.where(mask)[0]
        s = starts[mask]
        e = ends[mask]
        order = np.argsort(s)
        index[chrom] = (s[order], e[order], rows[order])
    return index


def _count_1bp_hits(positions, peak_starts, peak_ends, peak_rows, counts):
    """Vectorized counting of 1bp positions falling within sorted peak intervals."""
    if len(positions) == 0 or len(peak_starts) == 0:
        return

    idx = np.searchsorted(peak_starts, positions, side='right')
    candidate_idx = idx - 1
    valid = candidate_idx >= 0
    if not valid.any():
        return
    ci = candidate_idx[valid]
    vp = positions[valid]
    hits = peak_ends[ci] > vp
    if hits.any():
        np.add.at(counts, peak_rows[ci[hits]], 1)

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


def _count_fragments_coverage(fragments_file, peaks_df=None, *,
                               peak_index=None, n_peaks=None):
    """Count fragment coverage over peaks using polars or pandas."""
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
            dtypes={'c': pl.Utf8, 's': pl.Int64},
            comment_prefix='#',
        )
        for ch_name, group in df.group_by('c'):
            if isinstance(ch_name, tuple):
                ch_name = ch_name[0]
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
            comment='#', engine='c', chunksize=CHUNK,
        )
        for chunk in reader:
            chroms_arr = chunk['c'].values
            starts_arr = chunk['s'].values
            for ch in np.unique(chroms_arr):
                entry = peak_index.get(ch)
                if entry is None:
                    continue
                positions = starts_arr[chroms_arr == ch]
                _count_1bp_hits(positions, entry[0], entry[1], entry[2], counts)

    return counts


def _count_fragments_cutsites(fragments_file, peaks_df=None, *,
                               peak_index=None, n_peaks=None):
    """Count Tn5 cut sites (both fragment ends) within peaks."""
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
            dtypes={'c': pl.Utf8, 's': pl.Int64, 'e': pl.Int64},
            comment_prefix='#',
        )
        for ch_name, group in df.group_by('c'):
            if isinstance(ch_name, tuple):
                ch_name = ch_name[0]
            entry = peak_index.get(ch_name)
            if entry is None:
                continue
            starts_arr = group['s'].to_numpy()
            ends_arr = group['e'].to_numpy() - 1
            positions = np.concatenate([starts_arr, ends_arr])
            _count_1bp_hits(positions, entry[0], entry[1], entry[2], counts)
    else:
        CHUNK = 500_000
        reader = pd.read_csv(
            fragments_file, sep='\t', header=None,
            usecols=[0, 1, 2], names=['c', 's', 'e'],
            dtype={'c': str, 's': np.int64, 'e': np.int64},
            comment='#', engine='c', chunksize=CHUNK,
        )
        for chunk in reader:
            chroms_arr = chunk['c'].values
            starts_arr = chunk['s'].values
            ends_arr = chunk['e'].values - 1
            all_chroms = np.concatenate([chroms_arr, chroms_arr])
            all_positions = np.concatenate([starts_arr, ends_arr])
            for ch in np.unique(all_chroms):
                entry = peak_index.get(ch)
                if entry is None:
                    continue
                positions = all_positions[all_chroms == ch]
                _count_1bp_hits(positions, entry[0], entry[1], entry[2], counts)

    return counts


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
    """Count Tn5 insertions (single-bp BED) within peaks."""
    return _count_fragments_coverage(tn5_bed, peaks_df, peak_index=peak_index,
                                     n_peaks=n_peaks)


def _quantify_bigwig_legacy(bigwig_file, peaks_df, stat="mean"):
    """Quantify bigWig signal over peaks (row-by-row, legacy fallback)."""
    if not _HAS_PYBIGWIG:
        raise ImportError("pyBigWig is required. Install with: pip install pyBigWig")
    bw = pyBigWig.open(bigwig_file)
    values = []
    for _, row in peaks_df.iterrows():
        chrom = str(row['Chromosome'])
        start, end = int(row['Start']), int(row['End'])
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
        return _quantify_bigwig_legacy(input_file, peaks_df, stat=stat)
    else:
        raise ValueError(
            f"Unknown input_type '{input_type}'. "
            "Use 'fragments', 'tn5', or 'bigwig'."
        )


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
# PUBLIC: quantify  (single file – fragments / tn5 / bigwig)
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
        Split peaks into chunks for parallel processing.
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

    if input_type == "fragments" and method == "cutsites":
        if _are_fragments_already_cutsites(input_file):
            if verbose:
                print("  Auto-detected 1bp cut-site fragments, "
                      "using 'coverage' instead of 'cutsites'")
            method = "coverage"

    peaks_df = _harmonize_chr_prefix(peaks_df, input_file, input_type)
    peak_index = _build_peak_index(peaks_df)

    if verbose:
        print(f"Quantifying {Path(input_file).name} ({input_type}) "
              f"over {len(peaks_df):,} peaks...")

    counts = _dispatch(input_file, peaks_df, input_type, method, stat,
                       peak_index=peak_index, n_peaks=len(peaks_df))

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
    verbose: bool = True,
) -> Optional[pd.DataFrame]:
    """
    Quantify multiple files into a peaks × samples matrix.

    Parameters
    ----------
    input_files : list of str
        Paths to input files (fragments, tn5, or bigwig).
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
    name_pattern, name_replacement, use_stem
        Controls for automatic sample name cleaning.
    sample_names : list of str, optional
        Explicit sample names.
    output_file : str, optional
        Save matrix to this path.
    output_format : {'tsv', 'feather', 'parquet'}
        Disk format.
    verbose : bool
        Print progress.

    Returns
    -------
    pd.DataFrame or None
    """
    peaks_df = load_peaks(peak_file)
    n_peaks, n_files = len(peaks_df), len(input_files)

    if sample_names is not None and len(sample_names) != n_files:
        raise ValueError(
            f"sample_names length ({len(sample_names)}) != "
            f"input_files length ({n_files})"
        )

    if input_type == "fragments" and method == "cutsites" and n_files > 0:
        if _are_fragments_already_cutsites(input_files[0]):
            method = "coverage"
            if verbose:
                print("  Auto-detected 1bp cut-site fragments, "
                      "using 'coverage' instead of 'cutsites'")

    if n_files > 0:
        peaks_df = _harmonize_chr_prefix(peaks_df, input_files[0], input_type)

    if verbose:
        print(f"Quantifying {n_files} files over {n_peaks:,} peaks")
        extra = (f" | method: {method}" if input_type == "fragments"
                 else f" | stat: {stat}" if input_type == "bigwig" else "")
        print(f"   Input type: {input_type}{extra}")
        print(f"   Workers: {n_workers}")

    peak_index = _build_peak_index(peaks_df)
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
        Output path (extension adjusted to match output_format).
    output_format : {'tsv', 'feather', 'parquet'}
        Format.
    verbose : bool
        Print the saved path and file size.

    Returns
    -------
    str
        Actual path written.
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
        raise ValueError(f"Unknown format '{output_format}'.")

    if verbose:
        size_mb = os.path.getsize(out) / (1024 * 1024)
        print(f"   Saved to {out} ({size_mb:.1f} MB)")
    return out


def load_matrix(input_file: str) -> pd.DataFrame:
    """
    Load a quantification matrix from disk.

    Automatically detects format from the file extension.

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

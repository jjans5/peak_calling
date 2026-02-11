"""
Utility Functions
=================

Helper functions for configuration, file I/O, and common operations.
"""

import os
import json
import yaml
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any, Union, List

try:
    import pyranges as pr
    HAS_PYRANGES = True
except ImportError:
    HAS_PYRANGES = False


# =============================================================================
# CONFIGURATION
# =============================================================================

def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load configuration from a YAML file.
    
    Parameters
    ----------
    config_path : str
        Path to config.yaml file
    
    Returns
    -------
    dict
        Configuration dictionary
    """
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def save_parameters(
    params: Dict[str, Any],
    output_path: str,
    include_timestamp: bool = True,
) -> str:
    """
    Save parameters to a JSON file.
    
    Parameters
    ----------
    params : dict
        Parameters to save
    output_path : str
        Output file path
    include_timestamp : bool
        Whether to add timestamp to parameters
    
    Returns
    -------
    str
        Path to saved file
    """
    if include_timestamp:
        params['saved_at'] = datetime.now().isoformat()
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(params, f, indent=2)
    
    return output_path


# =============================================================================
# CHROMOSOME SIZES
# =============================================================================

# Default paths for chromosome sizes files
CHROMSIZES_PATHS = {
    "Human": "hg38.chrom.sizes",
    "Gorilla": "Gorilla_gorilla.gorGor4.dna.toplevel.chrom.sizes",
    "Chimpanzee": "Pan_troglodytes.Pan_tro_3.0.dna.toplevel.chrom.sizes",
    "Bonobo": "Pan_paniscus.panpan1.1.dna.toplevel.chrom.sizes",
    "Macaque": "Macaca_mulatta.Mmul_10.dna.toplevel.chrom.sizes",
    "Marmoset": "calJac1.chrom.sizes",
}


def get_chromsizes(
    species: str,
    chromsizes_dir: Optional[str] = None,
    chromsizes_path: Optional[str] = None,
    as_pyranges: bool = True,
    add_chr_prefix: bool = False,
) -> Union[pd.DataFrame, 'pr.PyRanges']:
    """
    Load chromosome sizes for a species.
    
    Parameters
    ----------
    species : str
        Species name
    chromsizes_dir : str, optional
        Directory containing chromsizes files
    chromsizes_path : str, optional
        Direct path to chromsizes file (overrides species lookup)
    as_pyranges : bool
        Whether to return as PyRanges (default) or DataFrame
    add_chr_prefix : bool
        Whether to add 'chr' prefix to chromosome names
    
    Returns
    -------
    PyRanges or DataFrame
        Chromosome sizes with columns [Chromosome, Start, End]
    """
    if chromsizes_path is None:
        if chromsizes_dir is None:
            raise ValueError("Either chromsizes_dir or chromsizes_path must be provided")
        chromsizes_path = os.path.join(chromsizes_dir, CHROMSIZES_PATHS.get(species, ""))
    
    if not os.path.exists(chromsizes_path):
        raise FileNotFoundError(f"Chromsizes file not found: {chromsizes_path}")
    
    # Load file
    df = pd.read_csv(chromsizes_path, sep="\t", header=None)
    
    # Handle different formats
    if df.shape[1] == 2:
        df.columns = ['Chromosome', 'End']
        df['Start'] = 0
    else:
        df.columns = ['Chromosome', 'End', 'Start']
    
    df = df[['Chromosome', 'Start', 'End']]
    
    # Add chr prefix if needed
    if add_chr_prefix and not str(df['Chromosome'].iloc[0]).startswith('chr'):
        df['Chromosome'] = 'chr' + df['Chromosome'].astype(str)
    
    if as_pyranges and HAS_PYRANGES:
        return pr.PyRanges(df)
    
    return df


def filter_main_chromosomes(
    chromsizes: Union[pd.DataFrame, 'pr.PyRanges'],
    species: str = "Human",
) -> Union[pd.DataFrame, 'pr.PyRanges']:
    """
    Filter to keep only main chromosomes (1-22, X, Y for human).
    
    Parameters
    ----------
    chromsizes : DataFrame or PyRanges
        Input chromsizes
    species : str
        Species (affects which chromosomes to keep)
    
    Returns
    -------
    DataFrame or PyRanges
        Filtered chromsizes
    """
    is_pyranges = HAS_PYRANGES and isinstance(chromsizes, pr.PyRanges)
    
    if is_pyranges:
        df = chromsizes.df
    else:
        df = chromsizes
    
    # Define main chromosomes per species
    if species == "Human":
        main_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    else:
        # For other primates, keep numbered chromosomes
        main_chroms = [f"chr{i}" for i in range(1, 30)] + ["chrX", "chrY"]
        # Also try without chr prefix
        main_chroms += [str(i) for i in range(1, 30)] + ["X", "Y"]
    
    df_filtered = df[df['Chromosome'].isin(main_chroms)]
    
    if is_pyranges:
        return pr.PyRanges(df_filtered)
    
    return df_filtered


# =============================================================================
# FILE HANDLING
# =============================================================================

def ensure_dir(path: str) -> str:
    """Create directory if it doesn't exist."""
    os.makedirs(path, exist_ok=True)
    return path


def find_files(
    directory: str,
    pattern: str = "*.narrowPeak",
    recursive: bool = False,
) -> list:
    """
    Find files matching a pattern in a directory.
    
    Parameters
    ----------
    directory : str
        Directory to search
    pattern : str
        Glob pattern to match
    recursive : bool
        Whether to search recursively
    
    Returns
    -------
    list
        List of matching file paths
    """
    p = Path(directory)
    if recursive:
        return sorted([str(f) for f in p.rglob(pattern)])
    return sorted([str(f) for f in p.glob(pattern)])


def count_lines(filepath: str) -> int:
    """Count lines in a file (handles gzipped files)."""
    import gzip
    
    if filepath.endswith('.gz'):
        with gzip.open(filepath, 'rt') as f:
            return sum(1 for _ in f)
    else:
        with open(filepath, 'r') as f:
            return sum(1 for _ in f)


# =============================================================================
# REPORTING
# =============================================================================

def generate_run_report(
    results: list,
    output_path: str,
    title: str = "Pipeline Run Report",
) -> str:
    """
    Generate a markdown report from pipeline results.
    
    Parameters
    ----------
    results : list
        List of result dictionaries
    output_path : str
        Path to save markdown report
    title : str
        Report title
    
    Returns
    -------
    str
        Path to saved report
    """
    lines = [
        f"# {title}",
        f"\n**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n",
        "## Summary\n",
    ]
    
    # Calculate statistics
    total = len(results)
    successful = sum(1 for r in results if r.get('status') == 'success')
    failed = total - successful
    
    lines.append(f"- **Total samples:** {total}")
    lines.append(f"- **Successful:** {successful}")
    lines.append(f"- **Failed:** {failed}")
    
    if 'peak_count' in results[0]:
        total_peaks = sum(r.get('peak_count', 0) for r in results)
        lines.append(f"- **Total peaks:** {total_peaks:,}")
    
    lines.append("\n## Details\n")
    lines.append("| Sample | Status | Details |")
    lines.append("|--------|--------|---------|")
    
    for r in results:
        name = r.get('sample_name', r.get('sample', 'Unknown'))
        status = '✅' if r.get('status') == 'success' else '❌'
        details = r.get('message', str(r.get('peak_count', '')))
        lines.append(f"| {name} | {status} | {details} |")
    
    content = "\n".join(lines)
    
    with open(output_path, 'w') as f:
        f.write(content)
    
    return output_path


# =============================================================================
# CHROMOSOME NAME UTILITIES
# =============================================================================

def modify_chr_prefix(
    input_bed: str,
    output_bed: str,
    add_chr: bool = True,
) -> Dict[str, Any]:
    """
    Add or remove 'chr' prefix from chromosome names in BED files.
    
    Handles both regular and gzipped BED files. Includes safety checks
    to prevent double 'chr' prefixes (e.g., 'chrchr1').
    
    Parameters
    ----------
    input_bed : str
        Path to input BED file (.bed or .bed.gz)
    output_bed : str
        Path to output BED file (.bed or .bed.gz)
    add_chr : bool, default=True
        If True, adds 'chr' prefix. If False, removes 'chr' prefix.
    
    Returns
    -------
    dict
        Dictionary with:
        - status: 'success' or 'error'
        - input: Input file path
        - output: Output file path
        - modified: Number of lines modified
        - total: Total number of lines processed
        - message: Status message
    
    Examples
    --------
    >>> # Add chr prefix
    >>> modify_chr_prefix("peaks.bed", "peaks_chr.bed", add_chr=True)
    
    >>> # Remove chr prefix
    >>> modify_chr_prefix("peaks_chr.bed", "peaks.bed", add_chr=False)
    """
    import gzip
    
    try:
        # Determine if input is gzipped
        is_input_gz = input_bed.endswith('.gz')
        is_output_gz = output_bed.endswith('.gz')
        
        # Open input file
        if is_input_gz:
            fin = gzip.open(input_bed, 'rt')
        else:
            fin = open(input_bed, 'r')
        
        # Open output file
        if is_output_gz:
            fout = gzip.open(output_bed, 'wt')
        else:
            fout = open(output_bed, 'w')
        
        modified_count = 0
        total_count = 0
        
        try:
            for line in fin:
                total_count += 1
                
                # Skip empty lines and comments
                if not line.strip() or line.startswith('#') or line.startswith('track'):
                    fout.write(line)
                    continue
                
                # Split line into fields
                fields = line.rstrip('\n').split('\t')
                chrom = fields[0]
                
                if add_chr:
                    # Add chr prefix
                    if not chrom.startswith('chr'):
                        # Safety check: ensure we're not adding to something like 'chrchr1'
                        if 'chrchr' not in chrom.lower():
                            fields[0] = 'chr' + chrom
                            modified_count += 1
                    # If it already has chr, don't modify
                    # Safety check: fix double chr if it exists
                    elif chrom.startswith('chrchr'):
                        fields[0] = chrom.replace('chrchr', 'chr', 1)
                        modified_count += 1
                else:
                    # Remove chr prefix
                    if chrom.startswith('chr'):
                        # Safety check: only remove one 'chr' at the start
                        fields[0] = chrom[3:]  # Remove first 3 characters ('chr')
                        modified_count += 1
                        
                        # Double check for any remaining chr prefix (safety)
                        if fields[0].startswith('chr'):
                            fields[0] = fields[0][3:]
                
                # Write modified line
                fout.write('\t'.join(fields) + '\n')
        
        finally:
            fin.close()
            fout.close()
        
        action = "Added" if add_chr else "Removed"
        return {
            'status': 'success',
            'input': input_bed,
            'output': output_bed,
            'modified': modified_count,
            'total': total_count,
            'message': f"{action} chr prefix for {modified_count}/{total_count} lines",
        }
    
    except Exception as e:
        return {
            'status': 'error',
            'input': input_bed,
            'output': output_bed,
            'modified': 0,
            'total': 0,
            'message': f"Error: {str(e)}",
        }


def add_chr_prefix(input_bed: str, output_bed: str) -> Dict[str, Any]:
    """
    Add 'chr' prefix to chromosome names in BED file.
    
    Convenience wrapper for modify_chr_prefix with add_chr=True.
    
    Parameters
    ----------
    input_bed : str
        Path to input BED file
    output_bed : str
        Path to output BED file
    
    Returns
    -------
    dict
        Result dictionary with status and statistics
    """
    return modify_chr_prefix(input_bed, output_bed, add_chr=True)


def remove_chr_prefix(input_bed: str, output_bed: str) -> Dict[str, Any]:
    """
    Remove 'chr' prefix from chromosome names in BED file.
    
    Convenience wrapper for modify_chr_prefix with add_chr=False.
    
    Parameters
    ----------
    input_bed : str
        Path to input BED file
    output_bed : str
        Path to output BED file
    
    Returns
    -------
    dict
        Result dictionary with status and statistics
    """
    return modify_chr_prefix(input_bed, output_bed, add_chr=False)


# =============================================================================
# BED FILE DIAGNOSTICS
# =============================================================================

def diagnose_bed(
    bed_file: str,
    name_col: int = None,
    celltype_col: int = None,
    show_plot: bool = True,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Comprehensive diagnostic analysis of a BED file.
    
    Provides statistics on peak sizes, chromosome distribution, and 
    optionally cell type distribution if the column is specified.
    
    Parameters
    ----------
    bed_file : str
        Path to BED file
    name_col : int, optional
        Column index (0-based) containing peak name/ID
    celltype_col : int, optional
        Column index (0-based) containing cell type information
    show_plot : bool
        Whether to display plots (default True)
    verbose : bool
        Print detailed statistics (default True)
    
    Returns
    -------
    dict
        Dictionary with diagnostic statistics:
        - n_peaks: total number of peaks
        - size_stats: dict with min, max, mean, median peak sizes
        - chrom_counts: dict of chromosome -> count
        - celltype_counts: dict of celltype -> count (if celltype_col provided)
        - df: pandas DataFrame with all data
    
    Example
    -------
    >>> diag = diagnose_bed("peaks.bed", celltype_col=3)
    >>> print(diag['size_stats'])
    >>> print(diag['celltype_counts'])
    """
    import numpy as np
    
    if not os.path.exists(bed_file):
        raise FileNotFoundError(f"BED file not found: {bed_file}")
    
    # Read BED file
    df = pd.read_csv(bed_file, sep='\t', header=None, comment='#')
    n_cols = len(df.columns)
    
    # Assign column names
    col_names = ['chrom', 'start', 'end']
    if n_cols > 3:
        col_names.extend([f'col{i}' for i in range(3, n_cols)])
    df.columns = col_names[:n_cols]
    
    # Calculate peak sizes
    df['size'] = df['end'] - df['start']
    
    # Basic stats
    n_peaks = len(df)
    size_stats = {
        'min': int(df['size'].min()),
        'max': int(df['size'].max()),
        'mean': float(df['size'].mean()),
        'median': float(df['size'].median()),
        'std': float(df['size'].std()),
        'q25': float(df['size'].quantile(0.25)),
        'q75': float(df['size'].quantile(0.75)),
    }
    
    # Chromosome distribution
    chrom_counts = df['chrom'].value_counts().to_dict()
    
    # Cell type distribution (if specified)
    celltype_counts = None
    if celltype_col is not None and celltype_col < n_cols:
        celltype_counts = df.iloc[:, celltype_col].value_counts().to_dict()
    
    # Name/ID distribution (if specified)
    name_counts = None
    if name_col is not None and name_col < n_cols:
        name_counts = df.iloc[:, name_col].value_counts().to_dict()
    
    if verbose:
        print("=" * 70)
        print(f"BED FILE DIAGNOSTICS: {os.path.basename(bed_file)}")
        print("=" * 70)
        
        print(f"\n📊 SUMMARY")
        print(f"   Total peaks: {n_peaks:,}")
        print(f"   Columns: {n_cols}")
        print(f"   File size: {os.path.getsize(bed_file) / 1024 / 1024:.2f} MB")
        
        print(f"\n📏 PEAK SIZE DISTRIBUTION")
        print(f"   Min:    {size_stats['min']:,} bp")
        print(f"   Max:    {size_stats['max']:,} bp")
        print(f"   Mean:   {size_stats['mean']:,.1f} bp")
        print(f"   Median: {size_stats['median']:,.1f} bp")
        print(f"   Std:    {size_stats['std']:,.1f} bp")
        print(f"   IQR:    {size_stats['q25']:,.1f} - {size_stats['q75']:,.1f} bp")
        
        print(f"\n🧬 CHROMOSOME DISTRIBUTION")
        # Sort chromosomes naturally
        sorted_chroms = sorted(chrom_counts.keys(), key=lambda x: (
            0 if x.replace('chr', '').isdigit() else 1,
            int(x.replace('chr', '')) if x.replace('chr', '').isdigit() else x
        ))
        for chrom in sorted_chroms[:25]:  # Show top 25
            count = chrom_counts[chrom]
            pct = count / n_peaks * 100
            bar = "█" * int(pct / 2)
            print(f"   {chrom:<8} {count:>8,} ({pct:>5.1f}%) {bar}")
        if len(sorted_chroms) > 25:
            print(f"   ... and {len(sorted_chroms) - 25} more chromosomes")
        
        if celltype_counts:
            print(f"\n🔬 CELL TYPE DISTRIBUTION (column {celltype_col})")
            sorted_ct = sorted(celltype_counts.items(), key=lambda x: -x[1])
            for ct, count in sorted_ct[:20]:  # Show top 20
                pct = count / n_peaks * 100
                bar = "█" * int(pct / 2)
                print(f"   {str(ct):<30} {count:>8,} ({pct:>5.1f}%) {bar}")
            if len(sorted_ct) > 20:
                print(f"   ... and {len(sorted_ct) - 20} more cell types")
        
        if name_counts and name_col != celltype_col:
            n_unique = len(name_counts)
            print(f"\n🏷️  NAME/ID DISTRIBUTION (column {name_col})")
            print(f"   Unique values: {n_unique:,}")
            if n_unique <= 20:
                for name, count in sorted(name_counts.items(), key=lambda x: -x[1]):
                    print(f"   {str(name):<30} {count:>8,}")
    
    # Plotting
    if show_plot:
        try:
            import matplotlib.pyplot as plt
            
            n_plots = 2 + (1 if celltype_counts else 0)
            fig, axes = plt.subplots(1, n_plots, figsize=(5 * n_plots, 4))
            if n_plots == 2:
                axes = [axes[0], axes[1]]
            
            # Plot 1: Peak size distribution
            ax1 = axes[0]
            ax1.hist(df['size'], bins=50, edgecolor='black', alpha=0.7)
            ax1.axvline(size_stats['median'], color='red', linestyle='--', label=f"Median: {size_stats['median']:.0f}")
            ax1.axvline(size_stats['mean'], color='orange', linestyle='--', label=f"Mean: {size_stats['mean']:.0f}")
            ax1.set_xlabel('Peak Size (bp)')
            ax1.set_ylabel('Count')
            ax1.set_title('Peak Size Distribution')
            ax1.legend()
            
            # Plot 2: Chromosome distribution (top 20)
            ax2 = axes[1]
            top_chroms = sorted_chroms[:20]
            chrom_vals = [chrom_counts[c] for c in top_chroms]
            ax2.barh(top_chroms[::-1], chrom_vals[::-1], edgecolor='black', alpha=0.7)
            ax2.set_xlabel('Count')
            ax2.set_title('Chromosome Distribution')
            
            # Plot 3: Cell type distribution (if available)
            if celltype_counts and n_plots > 2:
                ax3 = axes[2]
                sorted_ct = sorted(celltype_counts.items(), key=lambda x: -x[1])[:15]
                ct_names = [str(ct[0])[:20] for ct in sorted_ct]  # Truncate long names
                ct_vals = [ct[1] for ct in sorted_ct]
                ax3.barh(ct_names[::-1], ct_vals[::-1], edgecolor='black', alpha=0.7)
                ax3.set_xlabel('Count')
                ax3.set_title('Cell Type Distribution')
            
            plt.tight_layout()
            plt.show()
            
        except ImportError:
            if verbose:
                print("\n⚠️  matplotlib not available, skipping plots")
    
    return {
        'n_peaks': n_peaks,
        'n_cols': n_cols,
        'size_stats': size_stats,
        'chrom_counts': chrom_counts,
        'celltype_counts': celltype_counts,
        'name_counts': name_counts,
        'df': df,
    }


# Default main chromosomes for primates
DEFAULT_MAIN_CHROMS = [
    '1', '2', '2A', '2B', '3', '4', '5', '6', '7', '8', '9', '10',
    '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22',
    'X', 'Y', 'M',
    'chr1', 'chr2', 'chr2A', 'chr2B', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
    'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
    'chrX', 'chrY', 'chrM',
]


def diagnose_liftover(
    lifted_bed: str,
    unmapped_bed: str = None,
    show_plot: bool = True,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Diagnostic analysis of a liftover result.
    
    Takes a lifted-over BED file and automatically finds its unmapped mate.
    Compares chromosome distributions between mapped and unmapped peaks,
    using the source genome coordinates (last column) for the lifted file
    so both are compared in the same (source) genome coordinate space.
    
    Chromosome names are normalised (chr prefix harmonised) so that lifted
    source coordinates and unmapped coordinates are always merged correctly,
    even when the original input had a different naming convention than the
    chain file.
    
    Expected file naming convention:
        Lifted:   Consensus_Peaks_Filtered_500.hg38_Bonobo.bed
        Unmapped: Consensus_Peaks_Filtered_500.hg38_Bonobo.unmapped.bed
    
    The lifted BED file is expected to have a source coordinate column at
    the very end in the format chr:start-end (added by liftover_peaks).
    
    Parameters
    ----------
    lifted_bed : str
        Path to the lifted-over BED file
    unmapped_bed : str, optional
        Path to the unmapped BED file. If None, auto-detected by replacing
        .bed with .unmapped.bed in the lifted file path.
    show_plot : bool
        Whether to display diagnostic plots (default True)
    verbose : bool
        Print detailed statistics (default True)
    
    Returns
    -------
    dict
        Dictionary with:
        - lifted_total: number of lifted peaks
        - unmapped_total: number of unmapped peaks
        - liftover_rate: percentage of peaks successfully lifted
        - lifted_source_chroms: chromosome distribution in source genome (lifted peaks)
        - lifted_target_chroms: chromosome distribution in target genome (lifted peaks)
        - unmapped_chroms: chromosome distribution of unmapped peaks (source genome)
        - chrom_mapping: dict showing source -> target chromosome mappings
        - source_chrom_comparison: merged DataFrame comparing lifted vs unmapped per source chrom
    
    Example
    -------
    >>> diag = diagnose_liftover("peaks.hg38_Bonobo.bed")
    >>> # Automatically finds peaks.hg38_Bonobo.unmapped.bed
    >>> print(diag['liftover_rate'])
    """
    import re
    
    # Helper: normalise a chromosome name so that e.g. "1" and "chr1" become
    # the same key.  We always normalise to the "chr" form.
    def _norm_chrom(c: str) -> str:
        if not c.startswith('chr'):
            return 'chr' + c
        return c
    
    if not os.path.exists(lifted_bed):
        raise FileNotFoundError(f"Lifted BED file not found: {lifted_bed}")
    
    # Auto-detect unmapped file
    if unmapped_bed is None:
        if lifted_bed.endswith('.bed'):
            unmapped_bed = lifted_bed.replace('.bed', '.unmapped.bed')
        else:
            unmapped_bed = lifted_bed + '.unmapped'
    
    has_unmapped = os.path.exists(unmapped_bed)
    
    if verbose:
        print("=" * 70)
        print("LIFTOVER DIAGNOSTICS")
        print("=" * 70)
        print(f"📄 Lifted file:   {os.path.basename(lifted_bed)}")
        if has_unmapped:
            print(f"📄 Unmapped file: {os.path.basename(unmapped_bed)}")
        else:
            print(f"⚠️  Unmapped file not found: {os.path.basename(unmapped_bed)}")
    
    # ---- Parse lifted file ----
    lifted_target_chroms = {}   # target genome chroms (col 0) – raw names
    lifted_source_chroms = {}   # source genome chroms (normalised)
    chrom_mapping = {}          # norm_source_chrom -> set of target_chroms
    lifted_total = 0
    source_col_found = False
    
    with open(lifted_bed) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            
            lifted_total += 1
            target_chrom = parts[0]
            lifted_target_chroms[target_chrom] = lifted_target_chroms.get(target_chrom, 0) + 1
            
            # Try to parse the source coordinate from the last column
            # Format: chr:start-end or chrom:start-end
            last_col = parts[-1]
            m = re.match(r'^(.+):(\d+)-(\d+)$', last_col)
            if m:
                source_col_found = True
                src_chrom_raw = m.group(1)
                src_chrom = _norm_chrom(src_chrom_raw)
                lifted_source_chroms[src_chrom] = lifted_source_chroms.get(src_chrom, 0) + 1
                
                # Track mapping (normalised source -> raw target)
                if src_chrom not in chrom_mapping:
                    chrom_mapping[src_chrom] = set()
                chrom_mapping[src_chrom].add(target_chrom)
    
    if not source_col_found and verbose:
        print("\n⚠️  No source genome coordinates found in last column.")
        print("   Expected format: chr:start-end (added by liftover_peaks).")
        print("   Chromosome comparison will use target genome only.")
    
    # ---- Parse unmapped file ----
    # Unmapped peaks are still in (possibly chr-prefixed) source genome coords.
    # We normalise to the same namespace as the lifted source chroms.
    unmapped_chroms = {}        # normalised source chroms
    unmapped_total = 0
    
    if has_unmapped:
        with open(unmapped_bed) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    continue
                
                unmapped_total += 1
                chrom = _norm_chrom(parts[0])
                unmapped_chroms[chrom] = unmapped_chroms.get(chrom, 0) + 1
    
    # ---- Calculate statistics ----
    total_peaks = lifted_total + unmapped_total
    liftover_rate = (lifted_total / total_peaks * 100) if total_peaks > 0 else 0.0
    
    # ---- Build source chromosome comparison ----
    # Both dicts are now in the same normalised (chr-prefixed) namespace
    all_source_chroms = set()
    if source_col_found:
        all_source_chroms.update(lifted_source_chroms.keys())
    all_source_chroms.update(unmapped_chroms.keys())
    
    def _chrom_sort_key(x):
        """Sort chromosomes naturally: numeric first, then 2A/2B, then X/Y/M, then scaffolds."""
        base = x.replace('chr', '')
        if base.isdigit():
            return (0, int(base), '')
        if base in ('2A', '2B'):
            return (0, 2, base)
        if base in ('X', 'Y', 'M'):
            return (1, 0, base)
        return (2, 0, base)
    
    comparison_rows = []
    for chrom in sorted(all_source_chroms, key=_chrom_sort_key):
        n_lifted = lifted_source_chroms.get(chrom, 0) if source_col_found else 0
        n_unmapped = unmapped_chroms.get(chrom, 0)
        n_total = n_lifted + n_unmapped
        rate = (n_lifted / n_total * 100) if n_total > 0 else 0
        maps_to = ', '.join(sorted(chrom_mapping.get(chrom, set()))) if source_col_found else ''
        comparison_rows.append({
            'source_chrom': chrom,
            'lifted': n_lifted,
            'unmapped': n_unmapped,
            'total': n_total,
            'liftover_rate': round(rate, 1),
            'maps_to': maps_to,
        })
    
    comparison_df = pd.DataFrame(comparison_rows)
    
    # ---- Verbose output ----
    if verbose:
        print(f"\n📊 OVERALL STATISTICS")
        print(f"   Total input peaks:   {total_peaks:,}")
        print(f"   Successfully lifted: {lifted_total:,} ({liftover_rate:.1f}%)")
        print(f"   Unmapped:            {unmapped_total:,} ({100 - liftover_rate:.1f}%)")
        
        if source_col_found:
            print(f"\n🧬 CHROMOSOME MAPPING (source → target genome)")
            for src_chrom in sorted(chrom_mapping.keys(), key=_chrom_sort_key):
                targets = ', '.join(sorted(chrom_mapping[src_chrom]))
                n = lifted_source_chroms.get(src_chrom, 0)
                print(f"   {src_chrom:<10} → {targets:<15} ({n:,} peaks)")
        
        print(f"\n📋 PER-CHROMOSOME LIFTOVER RATES (source genome)")
        if not comparison_df.empty:
            # Filter to main chromosomes for display
            main_mask = comparison_df['source_chrom'].str.match(
                r'^(chr)?([\dXYM]+|2[AB])$', case=False
            )
            main_df = comparison_df[main_mask].copy()
            if not main_df.empty:
                print(main_df.to_string(index=False))
            
            other_df = comparison_df[~main_mask]
            if not other_df.empty:
                other_lifted = other_df['lifted'].sum()
                other_unmapped = other_df['unmapped'].sum()
                other_total = other_df['total'].sum()
                other_rate = (other_lifted / other_total * 100) if other_total > 0 else 0
                print(f"\n   Scaffold/other chromosomes: {len(other_df)} unique")
                print(f"      Total: {other_total:,}  Lifted: {other_lifted:,}  "
                      f"Unmapped: {other_unmapped:,}  Rate: {other_rate:.1f}%")
        
        # Check for noteworthy patterns
        print(f"\n🔍 NOTABLE OBSERVATIONS")
        # Check 2A/2B mapping
        has_2a_2b = any('2A' in c or '2B' in c or '2a' in c or '2b' in c 
                        for c in all_source_chroms)
        if has_2a_2b and source_col_found:
            for variant in ['chr2A', 'chr2B', '2A', '2B']:
                norm = _norm_chrom(variant)
                if norm in chrom_mapping:
                    targets = chrom_mapping[norm]
                    n_lift = lifted_source_chroms.get(norm, 0)
                    n_unmap = unmapped_chroms.get(norm, 0)
                    n_tot = n_lift + n_unmap
                    rate = (n_lift / n_tot * 100) if n_tot > 0 else 0
                    print(f"   • {norm} → {', '.join(sorted(targets))} "
                          f"({n_lift:,}/{n_tot:,} = {rate:.1f}% lifted)")
        
        # Flag chromosomes with low liftover rates
        if not comparison_df.empty:
            low_rate = comparison_df[
                (comparison_df['total'] >= 10) & 
                (comparison_df['liftover_rate'] < 50)
            ]
            if not low_rate.empty:
                print(f"   ⚠️  Low liftover rate (<50%) for:")
                for _, row in low_rate.iterrows():
                    print(f"      {row['source_chrom']}: {row['liftover_rate']:.1f}% "
                          f"({row['lifted']}/{row['total']})")
            
            # Flag chromosomes that are 100% unmapped (with >= 10 peaks)
            fully_unmapped = comparison_df[
                (comparison_df['total'] >= 10) & 
                (comparison_df['lifted'] == 0)
            ]
            if not fully_unmapped.empty:
                print(f"   ❌ Fully unmapped chromosomes:")
                for _, row in fully_unmapped.iterrows():
                    print(f"      {row['source_chrom']}: {row['total']} peaks, 0 lifted")
            
            # Highlight high-rate chromosomes
            if not main_df.empty:
                high_rate = main_df[main_df['liftover_rate'] >= 95]
                if not high_rate.empty:
                    chroms_str = ', '.join(high_rate['source_chrom'].tolist())
                    print(f"   ✅ ≥95% liftover: {chroms_str}")
    
    # ---- Plotting ----
    if show_plot:
        try:
            import matplotlib.pyplot as plt
            import numpy as np
            
            fig, axes = plt.subplots(1, 3, figsize=(18, 6))
            
            # Plot 1: Overall liftover pie chart
            ax1 = axes[0]
            ax1.pie(
                [lifted_total, unmapped_total],
                labels=[f'Lifted\n{lifted_total:,}', f'Unmapped\n{unmapped_total:,}'],
                autopct='%1.1f%%',
                colors=['#2ecc71', '#e74c3c'],
                startangle=90,
                explode=(0.03, 0.03),
            )
            ax1.set_title('Liftover Success Rate')
            
            # Plot 2: Per-chromosome stacked bar (main chroms only)
            ax2 = axes[1]
            if not comparison_df.empty:
                main_mask = comparison_df['source_chrom'].str.match(
                    r'^(chr)?([\dXYM]+|2[AB])$', case=False
                )
                plot_df = comparison_df[main_mask].copy()
                if not plot_df.empty:
                    plot_df = plot_df.sort_values('total', ascending=True)
                    y = np.arange(len(plot_df))
                    ax2.barh(y, plot_df['lifted'], color='#2ecc71', label='Lifted')
                    ax2.barh(y, plot_df['unmapped'], left=plot_df['lifted'], 
                             color='#e74c3c', label='Unmapped')
                    ax2.set_yticks(y)
                    ax2.set_yticklabels(plot_df['source_chrom'])
                    ax2.set_xlabel('Number of Peaks')
                    ax2.set_title('Per-Chromosome Breakdown\n(Source Genome)')
                    ax2.legend(loc='lower right')
            
            # Plot 3: Per-chromosome liftover rate with mapping labels
            ax3 = axes[2]
            if source_col_found and not comparison_df.empty:
                main_mask = comparison_df['source_chrom'].str.match(
                    r'^(chr)?([\dXYM]+|2[AB])$', case=False
                )
                main_df = comparison_df[main_mask].copy()
                if not main_df.empty:
                    main_df = main_df.sort_values('liftover_rate', ascending=True)
                    y = np.arange(len(main_df))
                    colors_bar = ['#e74c3c' if r < 50 else '#f39c12' if r < 90 else '#2ecc71' 
                                  for r in main_df['liftover_rate']]
                    ax3.barh(y, main_df['liftover_rate'], color=colors_bar)
                    ax3.set_yticks(y)
                    ax3.set_yticklabels(main_df['source_chrom'])
                    ax3.set_xlabel('Liftover Rate (%)')
                    ax3.set_xlim(0, 105)
                    ax3.axvline(x=90, color='green', linestyle='--', alpha=0.5, label='90%')
                    ax3.axvline(x=50, color='red', linestyle='--', alpha=0.5, label='50%')
                    ax3.set_title('Per-Chromosome Liftover Rate\n(Source Genome)')
                    ax3.legend()
                    
                    # Add mapping labels
                    for i, (_, row) in enumerate(main_df.iterrows()):
                        if row['maps_to']:
                            ax3.text(row['liftover_rate'] + 1, i, 
                                     f"→ {row['maps_to']}", va='center', fontsize=7)
            
            plt.suptitle(f'Liftover Diagnostics: {os.path.basename(lifted_bed)}', 
                         fontsize=13, y=1.02)
            plt.tight_layout()
            plt.show()
            
        except ImportError:
            if verbose:
                print("\n⚠️  matplotlib not available, skipping plots")
    
    return {
        'lifted_total': lifted_total,
        'unmapped_total': unmapped_total,
        'total_peaks': total_peaks,
        'liftover_rate': liftover_rate,
        'lifted_source_chroms': lifted_source_chroms,
        'lifted_target_chroms': lifted_target_chroms,
        'unmapped_chroms': unmapped_chroms,
        'chrom_mapping': {k: sorted(v) for k, v in chrom_mapping.items()},
        'source_chrom_comparison': comparison_df,
        'source_col_found': source_col_found,
    }


def compare_bed_files(
    bed_files: Dict[str, str],
    celltype_col: int = None,
    main_chroms: List[str] = None,
    show_plot: bool = True,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Compare multiple BED files side by side.
    
    Parameters
    ----------
    bed_files : dict
        Dictionary mapping label -> file path
    celltype_col : int, optional
        Column index for cell type information
    main_chroms : list, optional
        List of main chromosome names to highlight (default: 1-22, 2A, 2B, X, Y, M with/without chr prefix)
    show_plot : bool
        Whether to display comparison plots
    verbose : bool
        Print detailed statistics
    
    Returns
    -------
    dict
        Dictionary with:
        - summary_df: Comparison table with statistics for each file
        - chrom_df: Chromosome counts per file
        - diagnostics: Full diagnostic results per file
    """
    if main_chroms is None:
        main_chroms = DEFAULT_MAIN_CHROMS
    
    results = []
    all_chroms = {}  # {label: {chrom: count}}
    diagnostics = {}
    
    for label, filepath in bed_files.items():
        if not os.path.exists(filepath):
            print(f"⚠️  Skipping {label}: file not found")
            continue
        
        diag = diagnose_bed(filepath, celltype_col=celltype_col, show_plot=False, verbose=False)
        diagnostics[label] = diag
        all_chroms[label] = diag['chrom_counts']
        
        # Count main vs other chromosomes
        main_count = sum(count for chrom, count in diag['chrom_counts'].items() if chrom in main_chroms)
        other_count = sum(count for chrom, count in diag['chrom_counts'].items() if chrom not in main_chroms)
        
        results.append({
            'label': label,
            'n_peaks': diag['n_peaks'],
            'n_chroms': len(diag['chrom_counts']),
            'main_chroms': main_count,
            'other_chroms': other_count,
            'pct_main': main_count / diag['n_peaks'] * 100 if diag['n_peaks'] > 0 else 0,
            'size_min': diag['size_stats']['min'],
            'size_max': diag['size_stats']['max'],
            'size_mean': diag['size_stats']['mean'],
            'size_median': diag['size_stats']['median'],
            'n_celltypes': len(diag['celltype_counts']) if diag['celltype_counts'] else None,
        })
    
    summary_df = pd.DataFrame(results)
    
    # Build chromosome comparison dataframe
    all_chrom_names = set()
    for chrom_counts in all_chroms.values():
        all_chrom_names.update(chrom_counts.keys())
    
    chrom_data = []
    for chrom in sorted(all_chrom_names, key=lambda x: (
        0 if x.replace('chr', '') in ['1','2','2A','2B','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M'] else 1,
        int(x.replace('chr', '')) if x.replace('chr', '').isdigit() else 100 if x.replace('chr', '') in ['X', 'Y', 'M', '2A', '2B'] else 200,
        x
    )):
        row = {'chrom': chrom, 'is_main': chrom in main_chroms}
        for label in bed_files.keys():
            if label in all_chroms:
                row[label] = all_chroms[label].get(chrom, 0)
        chrom_data.append(row)
    
    chrom_df = pd.DataFrame(chrom_data)
    
    if verbose:
        print("=" * 80)
        print("BED FILE COMPARISON")
        print("=" * 80)
        print("\n📊 SUMMARY STATISTICS")
        print(summary_df[['label', 'n_peaks', 'main_chroms', 'other_chroms', 'pct_main', 'size_median']].to_string(index=False))
        
        print("\n🧬 CHROMOSOME COMPOSITION")
        # Show main chromosomes
        main_chrom_df = chrom_df[chrom_df['is_main']].copy()
        if not main_chrom_df.empty:
            print("\nMain chromosomes:")
            print(main_chrom_df.drop(columns=['is_main']).to_string(index=False))
        
        # Show other chromosomes summary
        other_chrom_df = chrom_df[~chrom_df['is_main']].copy()
        if not other_chrom_df.empty:
            print(f"\nOther/scaffold chromosomes: {len(other_chrom_df)} unique")
            label_cols = [col for col in other_chrom_df.columns if col not in ['chrom', 'is_main']]
            for label in label_cols:
                total = other_chrom_df[label].sum()
                print(f"   {label}: {total:,} peaks on scaffolds")
    
    if show_plot and len(results) > 1:
        try:
            import matplotlib.pyplot as plt
            import numpy as np
            
            fig, axes = plt.subplots(2, 2, figsize=(14, 10))
            
            # Plot 1: Peak counts
            ax1 = axes[0, 0]
            labels = summary_df['label']
            x = np.arange(len(labels))
            width = 0.35
            ax1.bar(x - width/2, summary_df['main_chroms'], width, label='Main chroms', alpha=0.8)
            ax1.bar(x + width/2, summary_df['other_chroms'], width, label='Other/scaffolds', alpha=0.8)
            ax1.set_xticks(x)
            ax1.set_xticklabels(labels, rotation=45, ha='right')
            ax1.set_ylabel('Number of Peaks')
            ax1.set_title('Peak Counts: Main vs Other Chromosomes')
            ax1.legend()
            
            # Plot 2: Main chromosome composition per file (stacked bar)
            ax2 = axes[0, 1]
            main_chrom_df = chrom_df[chrom_df['is_main']].copy()
            if not main_chrom_df.empty:
                label_cols = [col for col in main_chrom_df.columns if col not in ['chrom', 'is_main']]
                # Normalize to percentages
                chrom_pcts = main_chrom_df.copy()
                for col in label_cols:
                    total = chrom_pcts[col].sum()
                    if total > 0:
                        chrom_pcts[col] = chrom_pcts[col] / total * 100
                
                # Plot grouped bar for each file
                chroms_to_show = main_chrom_df['chrom'].tolist()
                x = np.arange(len(chroms_to_show))
                width = 0.8 / len(label_cols)
                
                for i, label in enumerate(label_cols):
                    offset = (i - len(label_cols)/2 + 0.5) * width
                    ax2.bar(x + offset, chrom_pcts[label], width, label=label, alpha=0.8)
                
                ax2.set_xticks(x)
                ax2.set_xticklabels(chroms_to_show, rotation=90)
                ax2.set_ylabel('Percentage of Peaks')
                ax2.set_title('Main Chromosome Composition (%)')
                ax2.legend(loc='upper right')
            
            # Plot 3: Size distributions
            ax3 = axes[1, 0]
            x = np.arange(len(summary_df))
            width = 0.35
            ax3.bar(x - width/2, summary_df['size_mean'], width, label='Mean', alpha=0.7)
            ax3.bar(x + width/2, summary_df['size_median'], width, label='Median', alpha=0.7)
            ax3.set_xticks(x)
            ax3.set_xticklabels(summary_df['label'], rotation=45, ha='right')
            ax3.set_ylabel('Peak Size (bp)')
            ax3.set_title('Peak Size Comparison')
            ax3.legend()
            
            # Plot 4: Chromosome counts per file (heatmap-style)
            ax4 = axes[1, 1]
            main_chrom_df = chrom_df[chrom_df['is_main']].copy()
            if not main_chrom_df.empty and len(label_cols) > 0:
                label_cols = [col for col in main_chrom_df.columns if col not in ['chrom', 'is_main']]
                heatmap_data = main_chrom_df[label_cols].values.T
                
                im = ax4.imshow(heatmap_data, aspect='auto', cmap='Blues')
                ax4.set_xticks(np.arange(len(main_chrom_df)))
                ax4.set_xticklabels(main_chrom_df['chrom'], rotation=90)
                ax4.set_yticks(np.arange(len(label_cols)))
                ax4.set_yticklabels(label_cols)
                ax4.set_title('Peak Counts Heatmap (Main Chromosomes)')
                plt.colorbar(im, ax=ax4, label='Count')
            
            plt.tight_layout()
            plt.show()
            
        except ImportError:
            print("⚠️  matplotlib not available, skipping plots")
    
    return {
        'summary_df': summary_df,
        'chrom_df': chrom_df,
        'diagnostics': diagnostics,
    }


# =============================================================================
# BED / PEAK FILE I/O
# =============================================================================

def load_peaks(
    peak_file: str,
    has_header: bool = False,
) -> pd.DataFrame:
    """
    Load peaks from BED file.
    
    Parameters
    ----------
    peak_file : str
        Path to BED file with peaks
    has_header : bool
        Whether the file has a header row
        
    Returns
    -------
    pd.DataFrame
        DataFrame with Chromosome, Start, End, and optionally Name columns
    """
    header = 0 if has_header else None
    df = pd.read_csv(peak_file, sep='\t', header=header, comment='#')
    
    # Standardize column names
    if df.shape[1] >= 3:
        cols = ['Chromosome', 'Start', 'End']
        if df.shape[1] >= 4:
            cols.append('Name')
        if df.shape[1] > 4:
            cols.extend([f'col_{i}' for i in range(4, df.shape[1])])
        df.columns = cols
    
    # Create region ID if no name
    if 'Name' not in df.columns:
        df['Name'] = df['Chromosome'] + ':' + df['Start'].astype(str) + '-' + df['End'].astype(str)
    
    return df


def clean_sample_name(
    filepath: str,
    pattern: Optional[str] = None,
    replacement: str = "",
    use_stem: bool = True,
) -> str:
    """
    Clean up sample name from filepath.
    
    Parameters
    ----------
    filepath : str
        Path to the input file
    pattern : str, optional
        Regex pattern to apply with re.sub
    replacement : str
        Replacement string for regex
    use_stem : bool
        If True, use only the filename without extension
        
    Returns
    -------
    str
        Cleaned sample name
    """
    import re
    
    if use_stem:
        name = Path(filepath).stem
        # Remove common double extensions
        for ext in ['.fragments', '.tn5', '.insertions', '.sorted', '.filtered']:
            if name.endswith(ext):
                name = name[:-len(ext)]
    else:
        name = Path(filepath).name
    
    if pattern:
        name = re.sub(pattern, replacement, name)
    
    return name


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


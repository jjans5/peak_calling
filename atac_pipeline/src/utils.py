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
from typing import Optional, Dict, Any, Union

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

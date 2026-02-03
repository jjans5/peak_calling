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

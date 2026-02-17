"""
Visualization Functions
=======================

Functions for plotting ATAC-seq data, genome regions, and peak distributions.
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.patheffects as PathEffects
from matplotlib.patches import ConnectionPatch
import numpy as np
import pandas as pd
import seaborn as sns
import os
from typing import Dict, List, Optional, Union, Any

# Optional imports for genome visualization
try:
    import pyBigWig
    HAS_PYBIGWIG = True
except ImportError:
    HAS_PYBIGWIG = False

try:
    import pybedtools
    HAS_PYBEDTOOLS = True
except ImportError:
    HAS_PYBEDTOOLS = False

try:
    import pyranges as pr
    HAS_PYRANGES = True
except ImportError:
    HAS_PYRANGES = False

# =============================================================================
# MATPLOTLIB CONFIGURATION
# =============================================================================

# Configuration for editable text in vector graphics
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def parse_region(region_str: str) -> tuple:
    """Parse a region string like 'chr1:1000-2000' into (chrom, start, end)."""
    region_str = region_str.replace(',', '')
    chrom, coords = region_str.split(':')
    start, end = map(int, coords.split('-'))
    return chrom, start, end


def smooth_signal(y: np.ndarray, window: int = 20) -> np.ndarray:
    """Smooth a signal using a moving average."""
    if len(y) < window:
        return y
    box = np.ones(window) / window
    return np.convolve(y, box, mode='same')


def check_chrom_format(bedtool_obj, chrom: str) -> str:
    """Check and harmonize chromosome naming format."""
    if not HAS_PYBEDTOOLS:
        return chrom
    try:
        first_interval = next(iter(bedtool_obj))
        file_has_chr = first_interval.chrom.startswith('chr')
        query_has_chr = chrom.startswith('chr')
        if file_has_chr and not query_has_chr:
            return 'chr' + chrom
        elif not file_has_chr and query_has_chr:
            return chrom.replace('chr', '')
        return chrom
    except StopIteration:
        return chrom


def get_gtf_features(gtf_file: str, chrom: str, start: int, end: int) -> pd.DataFrame:
    """Extract gene features from a GTF file for a given region."""
    if not gtf_file or not HAS_PYBEDTOOLS:
        return pd.DataFrame()
    
    import re
    gtf = pybedtools.BedTool(gtf_file)
    chrom_query = check_chrom_format(gtf, chrom)
    query_str = f'{chrom_query} {start} {end}'
    query = pybedtools.BedTool(query_str, from_string=True)
    hits = gtf.intersect(query, wa=True)
    
    if len(hits) == 0:
        return pd.DataFrame()
    
    df = hits.to_dataframe(names=[
        'chrom', 'source', 'feature', 'start', 'end', 
        'score', 'strand', 'frame', 'attributes'
    ])
    
    def parse_attr(attr_str, key):
        m = re.search(f'{key}\\s+"?([^";]+)"?', attr_str)
        return m.group(1) if m else None

    df['gene_name'] = df['attributes'].apply(lambda x: parse_attr(x, 'gene_name'))
    mask_nan = df['gene_name'].isna()
    if mask_nan.any():
        df.loc[mask_nan, 'gene_name'] = df.loc[mask_nan, 'attributes'].apply(
            lambda x: parse_attr(x, 'gene_id')
        )
    return df


# =============================================================================
# MAIN PLOTTING FUNCTIONS
# =============================================================================

def plot_genome_regions(
    bigwigs: List[str],
    selected_regions: Dict[str, str],
    gtf_file: Optional[str] = None,
    bed_tracks: Optional[List[str]] = None,
    bed_track_names: Optional[List[str]] = None,
    snp_file: Optional[str] = None,
    snp_color: str = 'red',
    plot_coordinates: bool = True,
    relative_coords: bool = False,
    colors: Optional[List[str]] = None,
    track_names: Optional[List[str]] = None,
    ymin: float = 0,
    ymax: Union[float, List[float]] = 75,
    figsize: tuple = (12, 10),
    smooth_window: int = 20,
    saveas: Optional[str] = None,
):
    """
    Plot multiple genome regions with bigWig tracks, gene annotations, and BED tracks.
    
    Parameters
    ----------
    bigwigs : list
        List of bigWig file paths
    selected_regions : dict
        Dictionary mapping region names to coordinates (e.g., {'Gene1': 'chr1:1000-2000'})
    gtf_file : str, optional
        Path to GTF file for gene annotations
    bed_tracks : list, optional
        List of BED files to display as tracks
    bed_track_names : list, optional
        Names for BED tracks
    snp_file : str, optional
        Path to BED file with SNP positions
    snp_color : str
        Color for SNP markers
    plot_coordinates : bool
        Whether to show coordinate axis
    relative_coords : bool
        Whether to show coordinates relative to region start
    colors : list, optional
        Colors for each bigWig track
    track_names : list, optional
        Names for each bigWig track
    ymin : float
        Minimum y-axis value
    ymax : float or list
        Maximum y-axis value (single value or per-region list)
    figsize : tuple
        Figure size
    smooth_window : int
        Window size for signal smoothing
    saveas : str, optional
        Path to save figure
    """
    if not HAS_PYBIGWIG:
        raise ImportError("pyBigWig is required for plot_genome_regions")
    
    # Setup data
    region_labels = list(selected_regions.keys())
    region_coords = [parse_region(r) for r in selected_regions.values()]
    n_regions = len(region_labels)
    n_bw = len(bigwigs)
    n_bed = len(bed_tracks) if bed_tracks else 0
    has_gtf = 1 if gtf_file else 0

    # Handle Y-limits
    if isinstance(ymax, (int, float)):
        ymax_list = [ymax] * n_regions
        vary_ylim = False
    else:
        if len(ymax) != n_regions:
            raise ValueError("Length of ymax list must match number of regions")
        ymax_list = ymax
        vary_ylim = True
    
    if track_names is None:
        track_names = [os.path.basename(x).split('.')[0] for x in bigwigs]
    
    if bed_tracks and bed_track_names is None:
        bed_track_names = [os.path.basename(x).split('.')[0] for x in bed_tracks]
    
    if colors is None:
        colors = ['#333333'] * n_bw
    elif len(colors) < n_bw:
        from itertools import cycle
        colors = [c for c, _ in zip(cycle(colors), range(n_bw))]

    # Pre-calculate SNPs
    region_snps = {i: [] for i in range(n_regions)}
    if snp_file and HAS_PYBEDTOOLS:
        snp_tool = pybedtools.BedTool(snp_file)
        for i, (chrom, start, end) in enumerate(region_coords):
            chrom_q = check_chrom_format(snp_tool, chrom)
            hits = snp_tool.intersect(
                pybedtools.BedTool(f'{chrom_q} {start} {end}', from_string=True), wa=True
            )
            region_snps[i] = [x.start for x in hits]

    # Setup figure
    total_rows = n_bw + n_bed + has_gtf
    height_ratios = [1.0] * n_bw + [0.2] * n_bed + ([0.5] if has_gtf else [])
    
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(
        total_rows, n_regions,
        height_ratios=height_ratios,
        wspace=0.15, hspace=0.1
    )

    col_axes_limits = {i: {'top': None, 'bottom': None} for i in range(n_regions)}

    # Plot BigWigs
    for row_idx, (bw_file, track_name, color) in enumerate(zip(bigwigs, track_names, colors)):
        bw = pyBigWig.open(bw_file)
        for col_idx, (label, (chrom, start, end)) in enumerate(zip(region_labels, region_coords)):
            ax = fig.add_subplot(gs[row_idx, col_idx])
            
            if row_idx == 0:
                col_axes_limits[col_idx]['top'] = ax
            col_axes_limits[col_idx]['bottom'] = ax
            
            current_ymax = ymax_list[col_idx]
            chrom_bw = chrom if chrom in bw.chroms() else (
                'chr' + chrom if 'chr' + chrom in bw.chroms() else chrom
            )
            
            try:
                values = bw.values(chrom_bw, start, end)
                values = np.nan_to_num(np.array(values))
                if smooth_window:
                    values = smooth_signal(values, smooth_window)
                x = np.arange(start, end)
                ax.fill_between(x, values, color=color, alpha=0.9, linewidth=0)
            except Exception as e:
                print(f"BW Error {track_name}: {e}")

            ax.set_xlim(start, end)
            ax.set_ylim(ymin, current_ymax)
            ax.set_yticks([ymin, current_ymax])
            
            if col_idx == 0:
                ax.set_ylabel(track_name, rotation=0, ha='right', va='center', fontsize=10)
                ax.tick_params(axis='y', labelsize=8)
            elif vary_ylim:
                ax.tick_params(axis='y', labelsize=8)
            else:
                ax.set_yticklabels([])

            ax.set_xticks([])
            sns.despine(ax=ax, left=False, bottom=True, right=True, top=True)
            if row_idx == 0:
                ax.set_title(label, fontsize=12)
        bw.close()

    current_row = n_bw

    # Plot BED tracks
    if bed_tracks and HAS_PYBEDTOOLS:
        for i, bed_file in enumerate(bed_tracks):
            bed = pybedtools.BedTool(bed_file)
            track_name = bed_track_names[i]
            
            for col_idx, (label, (chrom, start, end)) in enumerate(zip(region_labels, region_coords)):
                ax = fig.add_subplot(gs[current_row, col_idx])
                col_axes_limits[col_idx]['bottom'] = ax
                
                chrom_q = check_chrom_format(bed, chrom)
                hits = bed.intersect(
                    pybedtools.BedTool(f'{chrom_q} {start} {end}', from_string=True), wa=True
                )
                for hit in hits:
                    h_start, h_end = max(start, hit.start), min(end, hit.end)
                    ax.plot([h_start, h_end], [0, 0], linewidth=6, color='#555555', solid_capstyle='butt')
                ax.set_xlim(start, end)
                ax.set_ylim(-1, 1)
                ax.axis('off')
                
                if col_idx == 0:
                    ax.axis('on')
                    sns.despine(ax=ax, left=True, bottom=True, top=True, right=True)
                    ax.set_yticks([])
                    ax.set_xticks([])
                    ax.set_ylabel(track_name, rotation=0, ha='right', va='center', fontsize=10)

            current_row += 1

    # Plot GTF
    if has_gtf:
        for col_idx, (label, (chrom, start, end)) in enumerate(zip(region_labels, region_coords)):
            ax = fig.add_subplot(gs[current_row, col_idx])
            col_axes_limits[col_idx]['bottom'] = ax
            
            df = get_gtf_features(gtf_file, chrom, start, end)
            if not df.empty:
                df = df[(df['end'] > start) & (df['start'] < end)]
                genes = df['gene_name'].dropna().unique()
                if len(genes) == 0:
                    genes = df['attributes'].unique()
                for i, gene in enumerate(genes):
                    y_pos = i % 3
                    sub = df[df['gene_name'] == gene]
                    if sub.empty:
                        continue
                    g_start, g_end = sub['start'].min(), sub['end'].max()
                    ax.plot([max(start, g_start), min(end, g_end)], [y_pos, y_pos], c='k', lw=1)
                    exons = sub[sub['feature'] == 'exon']
                    for _, ex in exons.iterrows():
                        ax.plot([max(start, ex['start']), min(end, ex['end'])], [y_pos, y_pos], c='k', lw=6)
                    txt = ax.text(max(start, g_start), y_pos + 0.25, gene, fontsize=9, ha='left')
                    txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='w')])
                ax.set_ylim(-0.5, 3)
            ax.set_xlim(start, end)
            
            if plot_coordinates:
                ax.axis('on')
                sns.despine(ax=ax, left=True, bottom=False, top=True, right=True)
                ax.set_yticks([])
            else:
                ax.axis('off')

    # Draw SNP lines
    if snp_file:
        for col_idx in range(n_regions):
            snps = region_snps[col_idx]
            top_ax = col_axes_limits[col_idx]['top']
            bottom_ax = col_axes_limits[col_idx]['bottom']
            
            if top_ax and bottom_ax:
                top_ylim = top_ax.get_ylim()[1]
                bottom_ylim = bottom_ax.get_ylim()[0]
                
                for snp in snps:
                    con = ConnectionPatch(
                        xyA=(snp, top_ylim), xyB=(snp, bottom_ylim),
                        coordsA="data", coordsB="data",
                        axesA=top_ax, axesB=bottom_ax,
                        color=snp_color, linewidth=1, alpha=0.7, zorder=20
                    )
                    fig.add_artist(con)

    # Apply coordinate formatting
    if plot_coordinates:
        axes_to_label = fig.get_axes()[-n_regions:]
        for i, ax in enumerate(axes_to_label):
            chrom, start, end = region_coords[i]
            locator = ticker.MaxNLocator(nbins=3)
            ax.xaxis.set_major_locator(locator)
            
            if relative_coords:
                def rel_formatter(x, pos):
                    return f'{int(x - start):,}'
                ax.xaxis.set_major_formatter(ticker.FuncFormatter(rel_formatter))
                xlabel = f"{chrom} (relative)"
            else:
                ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
                xlabel = f"{chrom}"

            ax.set_xlabel(xlabel, fontsize=10)
            ax.tick_params(axis='x', labelsize=8)
            if not has_gtf:
                sns.despine(ax=ax, left=False, bottom=False, right=True, top=True)

    if saveas:
        plt.savefig(saveas, dpi=300, bbox_inches='tight')
    else:
        plt.show()


# =============================================================================
# SUMMARY PLOTS
# =============================================================================

def plot_peak_distribution(
    consensus_peaks,
    title: str = "Consensus Peak Width Distribution",
    target_width: int = 500,
    figsize: tuple = (10, 4),
    saveas: Optional[str] = None,
):
    """
    Plot the distribution of peak widths.
    
    Parameters
    ----------
    consensus_peaks : PyRanges
        Consensus peaks object
    title : str
        Plot title
    target_width : int
        Expected peak width (shown as vertical line)
    figsize : tuple
        Figure size
    saveas : str, optional
        Path to save figure
    """
    df = consensus_peaks.df.copy()
    df["width"] = df["End"] - df["Start"]
    
    plt.figure(figsize=figsize)
    sns.histplot(df["width"], bins=30)
    plt.title(title)
    plt.xlabel("Width (bp)")
    plt.axvline(x=target_width, color='red', linestyle='--', label=f"Target ({target_width}bp)")
    plt.legend()
    
    if saveas:
        plt.savefig(saveas, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    
    print(f"Total Peaks: {len(df)}")
    print(f"Mean Width: {df['width'].mean():.1f} bp")
    print(f"Median Width: {df['width'].median():.1f} bp")


def plot_consensus_summary(
    narrow_peaks_dict: dict,
    consensus_peaks,
    figsize: tuple = (12, 5),
    saveas: Optional[str] = None,
):
    """
    Plot contribution of each sample to consensus peaks.
    
    Parameters
    ----------
    narrow_peaks_dict : dict
        Dictionary of sample name -> PyRanges of original peaks
    consensus_peaks : PyRanges
        Final consensus peaks
    figsize : tuple
        Figure size
    saveas : str, optional
        Path to save figure
    """
    if not HAS_PYRANGES:
        raise ImportError("pyranges is required for plot_consensus_summary")
    
    counts = {}
    
    print("Calculating contribution per sample...")
    for sample, peaks_pr in narrow_peaks_dict.items():
        n_overlaps = len(peaks_pr.overlap(consensus_peaks))
        counts[sample] = n_overlaps

    counts_series = pd.Series(counts).sort_values(ascending=False)

    plt.figure(figsize=figsize)
    counts_series.plot(kind='bar', color='teal')
    plt.title("Number of Original Peaks Included in Consensus Set")
    plt.ylabel("Count of Overlapping Peaks")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    if saveas:
        plt.savefig(saveas, dpi=300, bbox_inches='tight')
    else:
        plt.show()


def plot_peak_counts_report(
    report_path: str,
    figsize: tuple = (12, 6),
    saveas: Optional[str] = None,
):
    """
    Plot peak counts from a MACS3 report file.
    
    Parameters
    ----------
    report_path : str
        Path to peak_counts_report.tsv
    figsize : tuple
        Figure size
    saveas : str, optional
        Path to save figure
    """
    df = pd.read_csv(report_path, sep='\t')
    df_sorted = df.sort_values('peak_count', ascending=False)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    colors = ['#2ecc71' if s == 'success' else '#e74c3c' for s in df_sorted['status']]
    ax.barh(df_sorted['cell_type'], df_sorted['peak_count'], color=colors)
    
    ax.set_xlabel('Peak Count')
    ax.set_ylabel('Cell Type')
    ax.set_title('Peaks Called per Cell Type')
    
    # Add summary stats
    total_peaks = df['peak_count'].sum()
    mean_peaks = df['peak_count'].mean()
    ax.axvline(mean_peaks, color='blue', linestyle='--', alpha=0.7, label=f'Mean: {mean_peaks:,.0f}')
    ax.legend()
    
    plt.tight_layout()
    
    if saveas:
        plt.savefig(saveas, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    
    print(f"Total peaks: {total_peaks:,}")
    print(f"Mean per cell type: {mean_peaks:,.0f}")


# =============================================================================
# UPSET PLOT  (pure matplotlib -- no upsetplot library)
# =============================================================================

def plot_upset(
    df: pd.DataFrame,
    set_columns: List[str],
    *,
    set_labels: Optional[List[str]] = None,
    top_n: int = 30,
    color: str = "steelblue",
    inactive_color: str = "#e0e0e0",
    title: Optional[str] = None,
    figsize: Optional[tuple] = None,
    show_counts: bool = True,
    count_fontsize: float = 7,
    count_rotation: float = 90,
    label_fontsize: float = 10,
    dot_size: float = 60,
    line_width: float = 2,
    bar_width: float = 0.7,
    saveas: Optional[str] = None,
    dpi: int = 150,
    show: bool = True,
) -> plt.Figure:
    """Create an UpSet plot showing intersection sizes for set-membership data.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame where each ``set_columns`` column is boolean / 0-1,
        indicating membership of that row in the corresponding set.
    set_columns : list of str
        Column names in *df* that represent set membership (bool or 0/1).
    set_labels : list of str, optional
        Display labels for the sets.  Defaults to *set_columns*.
    top_n : int
        Number of largest intersections to show (default 30).
    color : str
        Colour for bars and active dots (any matplotlib colour spec).
    inactive_color : str
        Colour for inactive (background) dots.
    title : str, optional
        Figure title.
    figsize : tuple, optional
        ``(width, height)`` in inches.  Auto-calculated when *None*.
    show_counts : bool
        Label each bar with its count.
    count_fontsize : float
        Font size for count labels on bars.
    count_rotation : float
        Rotation angle (degrees) for count labels.
    label_fontsize : float
        Font size for set labels on the y-axis of the dot matrix.
    dot_size : float
        Marker size for dots in the intersection matrix.
    line_width : float
        Width of vertical connector lines between active dots.
    bar_width : float
        Width of intersection-size bars (0–1).
    saveas : str, optional
        Path to save the figure.  If *None*, figure is not saved to disk.
    dpi : int
        Resolution for saved figure.
    show : bool
        Whether to call ``plt.show()``.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if set_labels is None:
        set_labels = list(set_columns)
    n_sets = len(set_columns)

    # --- build intersection counts ---
    mem = df[set_columns].astype(int)
    pattern = mem.apply(tuple, axis=1)
    combo_counts = pattern.value_counts().sort_values(ascending=False)

    combos = combo_counts.head(top_n)
    n_combos = len(combos)

    membership = np.array([list(c) for c in combos.index])   # (n_combos, n_sets)
    sizes = combos.values                                      # (n_combos,)
    set_sizes = mem.sum().values                               # (n_sets,)

    # --- figure layout ---
    if figsize is None:
        figsize = (max(14, n_combos * 0.55), 8)
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(
        2, 2,
        width_ratios=[n_combos, max(4, n_combos * 0.2)],
        height_ratios=[2.5, n_sets * 0.45],
        hspace=0.05, wspace=0.12,
    )
    ax_bar = fig.add_subplot(gs[0, 0])
    ax_mat = fig.add_subplot(gs[1, 0])
    ax_set = fig.add_subplot(gs[1, 1])
    ax_empty = fig.add_subplot(gs[0, 1])
    ax_empty.axis("off")

    x = np.arange(n_combos)

    # --- intersection-size bars ---
    ax_bar.bar(x, sizes, color=color, edgecolor="white", width=bar_width)
    if show_counts:
        for i, v in enumerate(sizes):
            ax_bar.text(
                i, v + max(sizes) * 0.01, f"{v:,}",
                ha="center", va="bottom",
                fontsize=count_fontsize, rotation=count_rotation,
            )
    ax_bar.set_ylabel("Intersection size")
    ax_bar.set_xlim(-0.6, n_combos - 0.4)
    ax_bar.set_xticks([])
    ax_bar.spines[["top", "right", "bottom"]].set_visible(False)
    if title:
        ax_bar.set_title(title, fontsize=13, pad=12)

    # --- dot matrix ---
    for j in range(n_sets):
        ax_mat.scatter(
            x, np.full(n_combos, j), s=dot_size,
            color=inactive_color, zorder=1, edgecolors="none",
        )
    for i in range(n_combos):
        active = np.where(membership[i] == 1)[0]
        if len(active) > 0:
            ax_mat.scatter(
                np.full(len(active), i), active, s=dot_size,
                color=color, zorder=2, edgecolors="none",
            )
            if len(active) > 1:
                ax_mat.plot(
                    [i, i], [active.min(), active.max()],
                    color=color, lw=line_width, zorder=1,
                )

    ax_mat.set_yticks(range(n_sets))
    ax_mat.set_yticklabels(set_labels, fontsize=label_fontsize)
    ax_mat.set_xlim(-0.6, n_combos - 0.4)
    ax_mat.set_ylim(-0.5, n_sets - 0.5)
    ax_mat.invert_yaxis()
    ax_mat.set_xticks([])
    ax_mat.spines[["top", "right", "bottom"]].set_visible(False)
    ax_mat.tick_params(left=False)

    # --- set-size bars ---
    y = np.arange(n_sets)
    ax_set.barh(y, set_sizes, color=color, edgecolor="white", height=0.6)
    for j, v in enumerate(set_sizes):
        ax_set.text(
            v + max(set_sizes) * 0.02, j, f"{v:,.0f}",
            va="center", fontsize=label_fontsize - 1,
        )
    ax_set.set_yticks([])
    ax_set.set_ylim(-0.5, n_sets - 0.5)
    ax_set.invert_yaxis()
    ax_set.set_xlabel("Set size", fontsize=label_fontsize)
    ax_set.spines[["top", "right", "left"]].set_visible(False)

    plt.tight_layout()

    if saveas:
        fig.savefig(saveas, dpi=dpi, bbox_inches="tight")
        print(f"Saved UpSet plot: {saveas}")
    if show:
        plt.show()

    return fig

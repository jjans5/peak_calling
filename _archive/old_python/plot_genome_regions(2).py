import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.patheffects as PathEffects
from matplotlib.patches import ConnectionPatch
import numpy as np
import pandas as pd
import seaborn as sns
import pyBigWig
import pybedtools
import re
import os
import matplotlib

# --- Configuration for Editable Text in Vector Graphics ---
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# --- Helper Functions ---

def parse_region(region_str):
    region_str = region_str.replace(',', '')
    chrom, coords = region_str.split(':')
    start, end = map(int, coords.split('-'))
    return chrom, start, end

def smooth_signal(y, window=20):
    if len(y) < window: return y
    box = np.ones(window) / window
    return np.convolve(y, box, mode='same')

def check_chrom_format(bedtool_obj, chrom):
    try:
        first_interval = next(iter(bedtool_obj))
        file_has_chr = first_interval.chrom.startswith('chr')
        query_has_chr = chrom.startswith('chr')
        if file_has_chr and not query_has_chr: return 'chr' + chrom
        elif not file_has_chr and query_has_chr: return chrom.replace('chr', '')
        return chrom
    except StopIteration:
        return chrom

def get_gtf_features(gtf_file, chrom, start, end):
    if not gtf_file: return pd.DataFrame()
    gtf = pybedtools.BedTool(gtf_file)
    chrom_query = check_chrom_format(gtf, chrom)
    query_str = f'{chrom_query} {start} {end}'
    query = pybedtools.BedTool(query_str, from_string=True)
    hits = gtf.intersect(query, wa=True)
    if len(hits) == 0: return pd.DataFrame()
    
    df = hits.to_dataframe(names=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'])
    
    def parse_attr(attr_str, key):
        m = re.search(f'{key}\s+"?([^";]+)"?', attr_str)
        return m.group(1) if m else None

    df['gene_name'] = df['attributes'].apply(lambda x: parse_attr(x, 'gene_name'))
    mask_nan = df['gene_name'].isna()
    if mask_nan.any():
        df.loc[mask_nan, 'gene_name'] = df.loc[mask_nan, 'attributes'].apply(lambda x: parse_attr(x, 'gene_id'))
    return df

# --- Main Plotting Function ---

def plot_genome_regions(bigwigs,
                        selected_regions,
                        gtf_file=None,
                        bed_tracks=None,
                        bed_track_names=None,  # NEW
                        snp_file=None,
                        snp_color='red',
                        plot_coordinates=True, 
                        relative_coords=False,
                        colors=None,
                        track_names=None,
                        ymin=0, 
                        ymax=75,
                        figsize=(12, 10),
                        smooth_window=20,
                        saveas=None):
    
    # 1. Setup Data
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
    
    # Handle Bed Track Names
    if bed_tracks and bed_track_names is None:
        bed_track_names = [os.path.basename(x).split('.')[0] for x in bed_tracks]
    
    if colors is None:
        colors = ['#333333'] * n_bw
    elif len(colors) < n_bw:
        from itertools import cycle
        colors = [c for c, _ in zip(cycle(colors), range(n_bw))]

    # 1.5 Pre-calculate SNPs
    region_snps = {i: [] for i in range(n_regions)}
    if snp_file:
        snp_tool = pybedtools.BedTool(snp_file)
        for i, (chrom, start, end) in enumerate(region_coords):
            chrom_q = check_chrom_format(snp_tool, chrom)
            hits = snp_tool.intersect(pybedtools.BedTool(f'{chrom_q} {start} {end}', from_string=True), wa=True)
            region_snps[i] = [x.start for x in hits]

    # 2. Setup Figure
    total_rows = n_bw + n_bed + has_gtf
    height_ratios = [1.0] * n_bw + [0.2] * n_bed + ([0.5] if has_gtf else [])
    
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(total_rows, n_regions, 
                           height_ratios=height_ratios, 
                           wspace=0.15, hspace=0.1)

    # To draw continuous lines, we need to track the top and bottom axes for each column
    col_axes_limits = {i: {'top': None, 'bottom': None} for i in range(n_regions)}

    # 3. Plot BigWigs
    for row_idx, (bw_file, track_name, color) in enumerate(zip(bigwigs, track_names, colors)):
        bw = pyBigWig.open(bw_file)
        for col_idx, (label, (chrom, start, end)) in enumerate(zip(region_labels, region_coords)):
            ax = fig.add_subplot(gs[row_idx, col_idx])
            
            # Store axes for line drawing
            if row_idx == 0: col_axes_limits[col_idx]['top'] = ax
            col_axes_limits[col_idx]['bottom'] = ax # Temporarily update bottom
            
            current_ymax = ymax_list[col_idx]
            chrom_bw = chrom if chrom in bw.chroms() else ('chr'+chrom if 'chr'+chrom in bw.chroms() else chrom)
            
            try:
                values = bw.values(chrom_bw, start, end)
                values = np.nan_to_num(np.array(values))
                if smooth_window: values = smooth_signal(values, smooth_window)
                x = np.arange(start, end)
                ax.fill_between(x, values, color=color, alpha=0.9, linewidth=0)
            except Exception as e:
                print(f"BW Error {track_name}: {e}")

            ax.set_xlim(start, end)
            ax.set_ylim(ymin, current_ymax)
            
            # Y-Axis Logic
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
            if row_idx == 0: ax.set_title(label, fontsize=12)
        bw.close()

    current_row = n_bw

    # 4. Plot Bed Tracks
    if bed_tracks:
        for i, bed_file in enumerate(bed_tracks):
            bed = pybedtools.BedTool(bed_file)
            track_name = bed_track_names[i]
            
            for col_idx, (label, (chrom, start, end)) in enumerate(zip(region_labels, region_coords)):
                ax = fig.add_subplot(gs[current_row, col_idx])
                col_axes_limits[col_idx]['bottom'] = ax # Update bottom
                
                chrom_q = check_chrom_format(bed, chrom)
                hits = bed.intersect(pybedtools.BedTool(f'{chrom_q} {start} {end}', from_string=True), wa=True)
                for hit in hits:
                    h_start, h_end = max(start, hit.start), min(end, hit.end)
                    ax.plot([h_start, h_end], [0, 0], linewidth=6, color='#555555', solid_capstyle='butt')
                ax.set_xlim(start, end)
                ax.set_ylim(-1, 1)
                ax.axis('off')
                
                # Add Name to the left-most plot
                if col_idx == 0:
                    ax.axis('on')
                    sns.despine(ax=ax, left=True, bottom=True, top=True, right=True)
                    ax.set_yticks([])
                    ax.set_xticks([])
                    ax.set_ylabel(track_name, rotation=0, ha='right', va='center', fontsize=10)

            current_row += 1

    # 5. Plot GTF
    if has_gtf:
        for col_idx, (label, (chrom, start, end)) in enumerate(zip(region_labels, region_coords)):
            ax = fig.add_subplot(gs[current_row, col_idx])
            col_axes_limits[col_idx]['bottom'] = ax # Update bottom
            
            df = get_gtf_features(gtf_file, chrom, start, end)
            if not df.empty:
                df = df[(df['end'] > start) & (df['start'] < end)]
                genes = df['gene_name'].dropna().unique()
                if len(genes) == 0: genes = df['attributes'].unique()
                for i, gene in enumerate(genes):
                    y_pos = i % 3
                    sub = df[df['gene_name'] == gene]
                    if sub.empty: continue
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

    # 6. Draw Continuous SNPs (Across subplots)
    if snp_file:
        for col_idx in range(n_regions):
            snps = region_snps[col_idx]
            top_ax = col_axes_limits[col_idx]['top']
            bottom_ax = col_axes_limits[col_idx]['bottom']
            
            if top_ax and bottom_ax:
                top_ylim = top_ax.get_ylim()[1]
                bottom_ylim = bottom_ax.get_ylim()[0]
                
                for snp in snps:
                    # ConnectionPatch allows drawing lines between points in different axes
                    # We connect (snp, max_y) in top axis to (snp, min_y) in bottom axis
                    con = ConnectionPatch(xyA=(snp, top_ylim), xyB=(snp, bottom_ylim), 
                                          coordsA="data", coordsB="data",
                                          axesA=top_ax, axesB=bottom_ax,
                                          color=snp_color, linewidth=1, alpha=0.7, zorder=20)
                    fig.add_artist(con)

    # 7. Apply Coordinate Formatting
    if plot_coordinates:
        axes_to_label = fig.get_axes()[-n_regions:]
        for i, ax in enumerate(axes_to_label):
            chrom, start, end = region_coords[i]
            locator = ticker.MaxNLocator(nbins=3)
            ax.xaxis.set_major_locator(locator)
            
            if relative_coords:
                def rel_formatter(x, pos): return f'{int(x - start):,}'
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
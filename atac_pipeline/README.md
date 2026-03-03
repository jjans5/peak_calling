# ATAC-seq Cross-Species Peak Analysis Pipeline

A reproducible pipeline for ATAC-seq peak calling, cross-species consensus peak analysis, and differential accessibility across 6 primate species.

## Overview

This pipeline provides a complete workflow for:

1. **Fragment Conversion** — Convert pseudobulk fragment files to Tn5 cut-site BED files
2. **Peak Calling** — Call peaks using MACS3 with ATAC-seq optimized parameters
3. **Liftover** — Lift peaks to human genome (hg38) for cross-species comparison
4. **Consensus Peaks** — Generate consensus peak sets per species
5. **Cross-Species Consensus** — Summit-based merge of peaks across 6 primates into a unified peak set with master annotation (`_det` / `_orth` tracking)
6. **Quantification** — BigWig-based quantification of fragment counts over the unified peak set
7. **Differential Accessibility** — DESeq2 analysis (Python & R) for cross-species comparisons
8. **Pseudobulk Matrices** — Single-cell fragment → pseudobulk count matrices

### Species

| Species | Assembly | Chain to hg38 |
|---------|----------|---------------|
| Human | hg38 | — |
| Gorilla | gorGor4 | gorGor4ToHg38 |
| Chimpanzee | panTro5 | panTro5ToHg38 |
| Bonobo | panPan2 | panPan2ToHg38 |
| Macaque | rheMac10 | rheMac10ToHg38 |
| Marmoset | calJac1 | calJac1→calJac4→hg38 (two-step) |

## Project Structure

```
atac_pipeline/
├── config/
│   ├── config.yaml              # Main configuration
│   └── assembly_mapping.json    # Assembly ↔ species mapping
├── notebooks/
│   ├── 01_peak_calling_workflow.ipynb    # Config-driven peak calling
│   ├── 02_flexible_workflow.ipynb        # Manual-path peak calling
│   ├── 03_liftover_consensus.ipynb       # Liftover consensus BEDs to hg38
│   ├── 04_cross_species_consensus.ipynb  # Cross-species unified peak set
│   ├── 05_quantification.ipynb           # BigWig quantification
│   ├── 06_deseq2_cross_species.ipynb     # DESeq2 (Python/PyDESeq2)
│   ├── 07_fragment_matrices_pseudobulk.ipynb  # sc fragments → pseudobulk
│   └── 08_deseq2_R.ipynb                # DESeq2 (R, pairwise comparisons)
├── scripts/
│   ├── liftover_fragments_par.sh         # Parallel fragment liftover
│   └── subsample_fragments.sh            # Fragment subsampling
├── src/
│   ├── __init__.py              # Package exports
│   ├── peak_calling.py          # MACS3 peak calling
│   ├── consensus.py             # Per-species consensus peaks
│   ├── liftover.py              # Liftover functions + chain files
│   ├── cross_species.py         # Cross-species consensus pipeline
│   ├── quantification.py        # Fragment/BigWig quantification
│   ├── bigwig.py                # BigWig generation
│   ├── visualization.py         # Plotting functions
│   └── utils.py                 # Utility functions
└── README.md
```

## Installation

### Dependencies

```bash
# Create conda environment
conda create -n atac_pipeline python=3.10
conda activate atac_pipeline

# Core dependencies
pip install numpy pandas pyranges pybedtools pyBigWig polars numba
pip install matplotlib seaborn

# Peak calling
conda install -c bioconda macs3

# Liftover (UCSC tools)
conda install -c bioconda ucsc-liftover
```

### Setup

```bash
# Clone or copy the pipeline
cd /path/to/your/analysis
cp -r atac_pipeline .

# Make scripts executable
chmod +x atac_pipeline/scripts/*.sh
```

## Usage

### Quick Start

1. Edit `config/config.yaml` with your paths and settings
2. Open `notebooks/01_peak_calling_workflow.ipynb`
3. Run the notebook cells sequentially

### Using as a Python Package

```python
import sys
sys.path.insert(0, '/path/to/atac_pipeline')

from src import (
    run_peak_calling,
    get_consensus_peaks,
    liftover_peaks,
    plot_genome_regions,
)

# Run peak calling
results = run_peak_calling(
    species="Gorilla",
    frag_dir="/path/to/fragments",
    out_dir="/path/to/output",
    max_workers=15,
)

# Generate consensus peaks
consensus = get_consensus_peaks(
    narrow_peaks_dict=peaks_dict,
    peak_half_width=250,
    chromsizes=chromsizes,
)
```

### Command-Line Scripts

```bash
# Liftover fragments to human genome
bash scripts/liftover_fragments_par.sh \
    --i input.fragments.tsv.gz \
    --c chain.over.chain.gz \
    --o output.hg38.fragments.tsv.gz \
    --ncpu 16

# Subsample fragments for testing
bash scripts/subsample_fragments.sh \
    --input full.fragments.tsv.gz \
    --output subsampled.fragments.tsv.gz \
    --ncells 1000 \
    --min-frags 500
```

## Workflow Steps

### Step 1: Fragment to Cut-site Conversion

Converts paired-end fragment files to single-nucleotide Tn5 cut-site BED files.

```python
from src import process_all_fragments

results = process_all_fragments(
    input_dir="/path/to/fragments",
    output_dir="/path/to/cutsites",
    max_workers=8,
)
```

### Step 2: MACS3 Peak Calling

Calls peaks using MACS3 with parameters optimized for ATAC-seq:

- `shift=-73`: Centers reads on Tn5 cut site
- `extsize=146`: Extension for nucleosome-free regions
- `nomodel`: Uses fixed shift/extsize
- `nolambda`: Fixed background estimation

```python
from src import run_peak_calling

results = run_peak_calling(
    species="Gorilla",
    frag_dir="/path/to/cutsites",
    out_dir="/path/to/peaks",
    qvalue=0.01,
    min_length=200,
)
```

### Step 3: Liftover to Human Genome

Lifts peaks to hg38 for cross-species comparison:

```python
from src import liftover_peaks

result = liftover_peaks(
    input_bed="peaks.bed",
    output_bed="peaks.hg38.bed",
    chain_file="gorGor4ToHg38.over.chain.gz",
)
```

### Step 4: Consensus Peak Calling

Generates consensus peaks by:
1. Extending peaks from summit
2. Normalizing scores (CPM)
3. Iteratively resolving overlaps by score

```python
from src import get_consensus_peaks, load_narrowpeaks

# Load peaks
peaks_dict = load_narrowpeaks(
    peak_dir="/path/to/peaks",
    q_value_threshold=0.05,
)

# Generate consensus
consensus = get_consensus_peaks(
    narrow_peaks_dict=peaks_dict,
    peak_half_width=250,  # 500bp total width
    chromsizes=chromsizes,
)

# Save
consensus.to_bed("consensus_peaks.bed")
```

### Step 5: Generate BigWig Files

Create genome coverage bigWig files from fragment files (inspired by [scatac_fragment_tools](https://github.com/aertslab/scatac_fragment_tools)):

```python
from src import create_bigwig, process_all_fragments_to_bigwig

# Single file
result = create_bigwig(
    fragments="sample.fragments.tsv.gz",
    chromsizes="hg38.chrom.sizes",
    output="sample.bw",
    cut_sites=True,  # Use Tn5 cut sites (recommended for ATAC)
    normalize=True,  # Normalize by RPM
)

# Batch processing
results = process_all_fragments_to_bigwig(
    input_dir="/path/to/fragments",
    output_dir="/path/to/bigwigs",
    chrom_sizes_file="hg38.chrom.sizes",
    pattern="*.tsv.gz",
    max_workers=4,
)
```

**Note**: BigWig generation requires `polars` and either `pyBigWig` or `pybigtools`. For best performance, also install `numba`.

## Configuration

Edit `config/config.yaml` to customize:

- **Species**: Target species for analysis
- **Paths**: Input/output directories
- **MACS3 parameters**: Peak calling settings
- **Consensus parameters**: Peak width, filters
- **Parallel workers**: CPU/thread allocation

## Cross-Species Approach

Our pipeline uses **pairwise UCSC liftOver chains** to map peaks between assemblies. The cross-species consensus is built via `summit_based_merge()`:

1. All species' consensus peaks are lifted to hg38
2. Peaks are clustered using bedtools (`-d 250`) on hg38 coordinates
3. Within each cluster, a score-weighted center is computed
4. Fixed 500bp windows are placed around each center
5. The unified set is lifted back to each native assembly for species-specific analysis

The master annotation tracks two distinct concepts per species:
- **`{Species}_det`** — A peak was *called* (detected) in that species near this location
- **`{Species}_orth`** — The unified peak could be *lifted back* to that species (ortholog exists)

This mirrors the orthology-vs-activity distinction used by the Zoonomia consortium (Kaplow et al. Science 2023), though they use Cactus multi-genome alignment + HALPER across 222 mammals, while we use pairwise liftOver for 6 closely-related primates with direct ATAC-seq measurement in each.

## Key Output Files

### Per Species
- `{sample}_peaks.narrowPeak` — Peak calls (BED6+4)
- `{sample}_summits.bed` — Peak summit positions
- `Consensus_Peaks_Filtered_500.bed` — Per-species consensus (500bp fixed-width)

### Cross-Species
- `unified_peaks.bed` — Merged peak set on hg38 (all species)
- `master_annotation.tsv` — Full annotation with `_det`, `_orth`, `n_species_det`, `n_species_orth`
- `peak_matrix.tsv` — Binary presence/absence matrix
- `{Species}_specific_peaks.bed` — Peaks unique to one species

## Visualization

```python
from src import plot_genome_regions, plot_peak_distribution

# Plot genome browser view
plot_genome_regions(
    bigwigs=["sample1.bw", "sample2.bw"],
    selected_regions={"Gene1": "chr1:1000000-1100000"},
    gtf_file="genes.gtf",
)

# Plot peak width distribution
plot_peak_distribution(consensus_peaks, target_width=500)
```

## License

MIT License

## References

- MACS3: Zhang et al., 2008
- pyranges: Stovner & Sætrom, 2020
- Zoonomia TACIT: Kaplow et al., Science 2023
- HALPER: Zhang et al., Bioinformatics 2020
- UCSC liftOver: Hinrichs et al., 2006

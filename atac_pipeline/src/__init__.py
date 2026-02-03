"""
ATAC-seq Peak Calling Pipeline
==============================

A reproducible pipeline for:
1. Fragment file processing (conversion to cut-sites)
2. Peak calling with MACS3
3. Liftover to human genome (hg38)
4. Consensus peak calling
5. Visualization

Author: J. Janssens
"""

from .peak_calling import (
    convert_fragments_to_cutsites,
    process_all_fragments,
    build_macs3_command,
    run_macs3_worker,
    run_peak_calling,
    EFFECTIVE_GENOME_SIZES,
    DEFAULT_MACS3_PARAMS,
)

from .consensus import (
    get_consensus_peaks,
    calculate_peaks_and_extend,
    iterative_peak_filtering,
    harmonize_chromosomes,
    load_narrowpeaks,
)

from .liftover import (
    liftover_peaks,
    liftover_two_step,
    liftover_fragments_parallel,
    print_chain_info,
    get_chain_file,
    CHAIN_FILES,
    DEFAULT_CHAIN_DIR,
)

from .visualization import (
    plot_genome_regions,
    plot_peak_distribution,
    plot_consensus_summary,
)

from .bigwig import (
    fragments_to_bigwig,
    create_bigwig,
    process_all_fragments_to_bigwig,
)

from .utils import (
    get_chromsizes,
    load_config,
    save_parameters,
    modify_chr_prefix,
    add_chr_prefix,
    remove_chr_prefix,
)

from .cross_species import (
    cross_species_consensus_pipeline,
    merge_bed_files,
    add_peak_ids,
    liftback_peaks,
    create_peak_matrix,
    get_reverse_chain_file,
    REVERSE_CHAIN_FILES,
)

__version__ = "1.0.0"
__all__ = [
    # Peak calling
    "convert_fragments_to_cutsites",
    "process_all_fragments", 
    "build_macs3_command",
    "run_macs3_worker",
    "run_peak_calling",
    "EFFECTIVE_GENOME_SIZES",
    "DEFAULT_MACS3_PARAMS",
    # Consensus
    "get_consensus_peaks",
    "calculate_peaks_and_extend",
    "iterative_peak_filtering",
    "harmonize_chromosomes",
    "load_narrowpeaks",
    # Liftover
    "liftover_peaks",
    "liftover_two_step",
    "liftover_fragments_parallel",
    "print_chain_info",
    "get_chain_file",
    "CHAIN_FILES",
    "DEFAULT_CHAIN_DIR",
    # Cross-species
    "cross_species_consensus_pipeline",
    "merge_bed_files",
    "add_peak_ids",
    "liftback_peaks",
    "create_peak_matrix",
    "get_reverse_chain_file",
    "REVERSE_CHAIN_FILES",
    # BigWig
    "fragments_to_bigwig",
    "create_bigwig",
    "process_all_fragments_to_bigwig",
    # Visualization
    "plot_genome_regions",
    "plot_peak_distribution",
    "plot_consensus_summary",
    # Utils
    "get_chromsizes",
    "load_config",
    "save_parameters",
    "modify_chr_prefix",
    "add_chr_prefix",
    "remove_chr_prefix",
]

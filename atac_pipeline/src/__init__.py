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
    plot_upset,
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
    diagnose_bed,
    diagnose_liftover,
    compare_bed_files,
    load_peaks,
    clean_sample_name,
)

from .cross_species import (
    cross_species_consensus_pipeline,
    merge_with_species_tracking,
    merge_bed_files,
    add_peak_ids,
    liftback_peaks,
    create_peak_matrix,
    build_master_annotation,
    create_peak_annotation,
    extract_gene_bed_from_gtf,
    annotate_with_closest_gene,
    classify_peak_distance,
    find_species_specific_peaks,
    cross_map_species_specific_peaks,
    get_reverse_chain_file,
    summit_based_merge,
    liftover_summits_with_conservation,
    filter_liftback_by_size,
    update_master_after_filtering,
    filter_and_update_master,
    resize_peaks,
    reciprocal_liftover_check,
    finalize_peak_sets,
    REVERSE_CHAIN_FILES,
    CROSS_SPECIES_ROUTES,
    DEFAULT_GTF_FILES,
)

from .pipeline_steps import (
    lift_to_human,
    merge_consensus,
    lift_back_to_species,
    filter_liftback,
    extract_species_specific,
    generate_master_annotation,
    rescue_unmapped_peaks,
    classify_distances,
    export_final_beds,
)

from .quantification import (
    quantify,
    quantify_matrix,
    quantify_bigwig,
    quantify_bigwig_matrix,
    fragments_to_bigwigs,
    save_matrix,
    load_matrix,
)

from .fragment_matrices import (
    reindex_nhp,
    reindex_human,
    load_species_data,
    load_regions_as_polars,
    harmonize_chroms,
    build_fragment_matrix,
    create_pseudobulk,
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
    "merge_with_species_tracking",
    "merge_bed_files",  # deprecated
    "add_peak_ids",
    "liftback_peaks",
    "create_peak_matrix",  # deprecated
    "build_master_annotation",
    "create_peak_annotation",
    "extract_gene_bed_from_gtf",
    "annotate_with_closest_gene",
    "classify_peak_distance",
    "find_species_specific_peaks",
    "cross_map_species_specific_peaks",
    "get_reverse_chain_file",
    "summit_based_merge",
    "liftover_summits_with_conservation",
    "filter_liftback_by_size",
    "update_master_after_filtering",
    "filter_and_update_master",
    "resize_peaks",
    "reciprocal_liftover_check",
    "finalize_peak_sets",
    "REVERSE_CHAIN_FILES",
    "CROSS_SPECIES_ROUTES",
    "DEFAULT_GTF_FILES",
    # Quantification
    "quantify",
    "quantify_matrix",
    "quantify_bigwig",
    "quantify_bigwig_matrix",
    "fragments_to_bigwigs",
    "save_matrix",
    "load_matrix",
    # Fragment matrices & pseudobulk
    "reindex_nhp",
    "reindex_human",
    "load_species_data",
    "load_regions_as_polars",
    "harmonize_chroms",
    "build_fragment_matrix",
    "create_pseudobulk",
    # Pipeline steps (high-level wrappers)
    "lift_to_human",
    "merge_consensus",
    "lift_back_to_species",
    "filter_liftback",
    "extract_species_specific",
    "generate_master_annotation",
    "rescue_unmapped_peaks",
    "classify_distances",
    "export_final_beds",
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
    "diagnose_bed",
    "diagnose_liftover",
    "compare_bed_files",
    "load_peaks",
    "clean_sample_name",
]

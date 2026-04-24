#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.publication_figures import (
    DEFAULT_CELL_TYPES,
    DEFAULT_CONTRASTS,
    build_volcano_suite,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build publication-ready cross-species differential accessibility volcano panels."
    )
    parser.add_argument(
        "--results-root",
        required=True,
        help="Directory containing per-cell-type evolutionary branch DESeq2 CSVs.",
    )
    parser.add_argument(
        "--annotation-file",
        required=True,
        help="Peak annotation TSV used for nearest-gene labeling.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for the manifest and volcano figure files.",
    )
    parser.add_argument(
        "--cell-types",
        nargs="+",
        default=list(DEFAULT_CELL_TYPES),
        help="Cell types to include in each multi-panel figure.",
    )
    parser.add_argument(
        "--contrasts",
        nargs="+",
        default=list(DEFAULT_CONTRASTS),
        help="Contrast names to render, matching branch result CSV stems.",
    )
    parser.add_argument("--padj-thresh", type=float, default=0.05)
    parser.add_argument("--lfc-thresh", type=float, default=1.0)
    parser.add_argument("--labels-per-direction", type=int, default=8)
    args = parser.parse_args()

    outputs = build_volcano_suite(
        results_root=args.results_root,
        annotation_path=args.annotation_file,
        output_dir=args.output_dir,
        cell_types=args.cell_types,
        contrasts=args.contrasts,
        padj_thresh=args.padj_thresh,
        lfc_thresh=args.lfc_thresh,
        n_per_direction=args.labels_per_direction,
    )

    print(f"Manifest: {outputs['manifest']}")
    for item in outputs["outputs"]:
        print(f"[{item['contrast']}] pdf={item['pdf']}")
        print(f"[{item['contrast']}] png={item['png']}")
        print(f"[{item['contrast']}] labels={item['labels']}")


if __name__ == "__main__":
    main()
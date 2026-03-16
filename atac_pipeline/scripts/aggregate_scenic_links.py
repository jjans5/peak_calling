#!/usr/bin/env python3
import argparse
import os

from src.scenicplus_prep import aggregate_links_for_species


def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate SCENIC+ link tables across seeds")
    parser.add_argument("--scenicplus-dir", required=True, help="Base scenicplus directory")
    parser.add_argument("--species", required=True, nargs="+", help="Species names")
    parser.add_argument("--output-dir", required=True, help="Output folder for consensus links")
    parser.add_argument("--min-seed-support", type=int, default=2)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    for species in args.species:
        species_out = os.path.join(args.output_dir, species)
        outputs = aggregate_links_for_species(
            scenicplus_dir=args.scenicplus_dir,
            species=species,
            output_dir=species_out,
            min_seed_support=args.min_seed_support,
        )
        print(f"[{species}] tf_to_gene: {outputs['tf_to_gene']}")
        print(f"[{species}] region_to_gene: {outputs['region_to_gene']}")


if __name__ == "__main__":
    main()

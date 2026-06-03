#!/usr/bin/env python3

import os
import sys
import argparse
from ete3 import NCBITaxa
import contextlib

bulk_reference_dir = os.environ.get("BULK_REFERENCES_DIR")

if bulk_reference_dir is None:
    raise ValueError("BULK_REFERENCES_DIR not set")


dbfile = os.path.join(bulk_reference_dir, "taxonomy", "taxa.sqlite")
os.makedirs(os.path.dirname(dbfile), exist_ok=True)

try:
    if not os.path.exists(dbfile):
        # Suppress any output produced internally by ete3 during DB download/update
        with open(os.devnull, "w") as devnull, \
             contextlib.redirect_stdout(devnull), \
             contextlib.redirect_stderr(devnull):
            ncbi = NCBITaxa(dbfile=dbfile)
            ncbi.update_taxonomy_database()
    else:
        ncbi = NCBITaxa(dbfile=dbfile)
except Exception as e:
    print(f"Failed to initialize NCBI taxonomy database: {e}", file=sys.stderr)
    sys.exit(1)


# Reference taxids
GROUPS = {
    "bryophytes_plants": {3208, 3195, 3209},
    "monocot_plants": {4447},
    "dicot_plants": {71240, 91827, 91835},
    "yeast": {4895, 4894, 4890, 4930},
    "protist": {554915, 33630, 33634, 543769}
}

def classify_species(taxid):
    try:
        lineage = set(ncbi.get_lineage(int(taxid)))

        for group, taxids in GROUPS.items():
            if lineage & taxids:
                return group

        return "default"

    except Exception:
        return "default"


def main():
    parser = argparse.ArgumentParser(description="Classify a single TaxID")
    parser.add_argument("taxid", type=int, help="NCBI TaxID (e.g. 3702)")
    args = parser.parse_args()

    category = classify_species(args.taxid)
    print(category)


if __name__ == "__main__":
    main()

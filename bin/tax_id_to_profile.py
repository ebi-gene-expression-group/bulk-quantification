#!/usr/bin/env python3

import os
import argparse
from ete3 import NCBITaxa

dbfile = "./taxonomy/taxa.sqlite"
os.makedirs(os.path.dirname(dbfile), exist_ok=True)

try:
    if not os.path.exists(dbfile):
        print(f"Downloading taxonomy database to {dbfile}...", file=sys.stderr)
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
    "yeast": {4892, 147537, 4893},
    "protist": {554915, 33630, 33634, 543769}
}

def classify_species(taxid):
    try:
        lineage = set(ncbi.get_lineage(int(taxid)))

        for group, taxids in GROUPS.items():
            if lineage & taxids:
                return group

        return "other"

    except Exception:
        return "other"


def main():
    parser = argparse.ArgumentParser(description="Classify a single TaxID")
    parser.add_argument("taxid", help="NCBI TaxID (e.g. 3702)")
    args = parser.parse_args()

    category = classify_species(args.taxid)
    print(category)


if __name__ == "__main__":
    main()

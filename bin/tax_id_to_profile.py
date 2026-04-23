#!/usr/bin/env python3

import argparse
from ete3 import NCBITaxa

ncbi = NCBITaxa()

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

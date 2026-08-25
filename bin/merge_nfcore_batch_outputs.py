#!/usr/bin/env python3

import argparse
import json
import os
import shutil
import sys
from pathlib import Path

import pandas as pd


MATRIX_PATTERNS = [
    "salmon.merged.gene_counts.tsv",
    "salmon.merged.gene_counts_scaled.tsv",
    "salmon.merged.gene_counts_length_scaled.tsv",
    "salmon.merged.gene_tpm.tsv",
    "salmon.merged.gene_lengths.tsv",
    "salmon.merged.transcript_counts.tsv",
    "salmon.merged.transcript_tpm.tsv",
    "salmon.merged.transcript_lengths.tsv",
    "kallisto.merged.gene_counts.tsv",
    "kallisto.merged.gene_tpm.tsv",
    "kallisto.merged.transcript_counts.tsv",
    "kallisto.merged.transcript_tpm.tsv",
    "rsem.merged.gene_counts.tsv",
    "rsem.merged.gene_tpm.tsv",
    "rsem.merged.transcript_counts.tsv",
    "rsem.merged.transcript_tpm.tsv",
    "featureCounts.merged.counts.tsv",
    "featureCounts.merged.counts.txt",
    "featureCounts.merged.counts.tsv.summary",
    "featureCounts.merged.counts.txt.summary",
]


QUANT_DIRS = [
    "star_salmon",
    "salmon",
    "star_rsem",
    "rsem",
    "star_featurecounts",
    "hisat2",
    "star",
]


METADATA_PATHS = [
    "pipeline_info",
    "multiqc",
    "multiqc_report.html",
    "multiqc_data",
    "samplesheet.csv",
    "samplesheet.valid.csv",
    "params.json",
]


def info(msg):
    print(f"[INFO] {msg}", file=sys.stderr)


def warn(msg):
    print(f"[WARN] {msg}", file=sys.stderr)


def fail(msg):
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(1)


def read_table(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.shape[1] < 2:
        raise ValueError(f"{path} has fewer than 2 columns")

    first_col = df.columns[0]
    df = df.set_index(first_col)

    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="ignore")

    return df, first_col


def merge_matrices(paths, out_path, strict_features=True):
    merged = None
    feature_col = None
    seen_samples = set()

    for path in paths:
        df, this_feature_col = read_table(path)

        if feature_col is None:
            feature_col = this_feature_col

        dup_samples = seen_samples.intersection(df.columns)
        if dup_samples:
            raise ValueError(
                f"Duplicate sample columns found while merging {path}: "
                f"{', '.join(sorted(dup_samples))}"
            )

        seen_samples.update(df.columns)

        if merged is None:
            merged = df
        else:
            if strict_features:
                if not merged.index.equals(df.index):
                    missing_left = df.index.difference(merged.index)
                    missing_right = merged.index.difference(df.index)
                    raise ValueError(
                        f"Feature IDs/order do not match in {path}\n"
                        f"Features only in this file: {len(missing_left)}\n"
                        f"Features missing from this file: {len(missing_right)}\n"
                        f"Use --relaxed-features to outer-join instead."
                    )
                merged = pd.concat([merged, df], axis=1)
            else:
                merged = merged.join(df, how="outer")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.reset_index(names=feature_col).to_csv(out_path, sep="\t", index=False)
    return merged.shape


def safe_symlink_or_copy(src, dst, copy=False):
    dst.parent.mkdir(parents=True, exist_ok=True)

    if dst.exists() or dst.is_symlink():
        raise FileExistsError(f"Destination already exists: {dst}")

    if copy:
        if src.is_dir():
            shutil.copytree(src, dst)
        else:
            shutil.copy2(src, dst)
    else:
        os.symlink(src.resolve(), dst)


def copy_metadata(batch_dirs, out_dir):
    meta_root = out_dir / "_batch_inputs"
    meta_root.mkdir(parents=True, exist_ok=True)

    for batch in batch_dirs:
        batch_name = batch.name.rstrip("/")
        dest = meta_root / batch_name
        dest.mkdir(parents=True, exist_ok=True)

        for rel in METADATA_PATHS:
            src = batch / rel
            if src.exists():
                target = dest / rel
                if target.exists():
                    continue
                if src.is_dir():
                    shutil.copytree(src, target)
                else:
                    shutil.copy2(src, target)


def find_matrix_files(batch_dirs):
    found = {}

    for batch in batch_dirs:
        for pattern in MATRIX_PATTERNS:
            for path in batch.rglob(pattern):
                rel = path.relative_to(batch)
                found.setdefault(str(rel), []).append(path)

    return found


def link_sample_dirs(batch_dirs, out_dir, copy=False):
    linked = []

    for batch in batch_dirs:
        for qdir in QUANT_DIRS:
            src_root = batch / qdir
            if not src_root.is_dir():
                continue

            out_root = out_dir / qdir
            out_root.mkdir(parents=True, exist_ok=True)

            for child in src_root.iterdir():
                if not child.is_dir():
                    continue

                skip_names = {
                    "logs",
                    "log",
                    "bigwig",
                    "samtools_stats",
                    "rseqc",
                    "deseq2_qc",
                    "rsem_merge_counts",
                }
                if child.name in skip_names:
                    continue

                dst = out_root / child.name
                if dst.exists() or dst.is_symlink():
                    raise FileExistsError(
                        f"Duplicate sample or output directory name: {dst}\n"
                        f"Conflicting source: {child}\n"
                        f"Rename samples or use unique sample IDs before merging."
                    )

                safe_symlink_or_copy(child, dst, copy=copy)
                linked.append((str(child), str(dst)))

    return linked


def copy_reference_like_files(batch_dirs, out_dir):
    candidates = [
        "star_salmon/tx2gene.tsv",
        "salmon/tx2gene.tsv",
        "kallisto/tx2gene.tsv",
        "genome",
        "igenomes.config",
    ]

    copied = []

    for rel in candidates:
        for batch in batch_dirs:
            src = batch / rel
            if src.exists():
                dst = out_dir / rel
                dst.parent.mkdir(parents=True, exist_ok=True)
                if not dst.exists():
                    if src.is_dir():
                        shutil.copytree(src, dst)
                    else:
                        shutil.copy2(src, dst)
                    copied.append((str(src), str(dst)))
                break

    return copied


def write_manifest(out_dir, batch_dirs, merged_records, linked_records, copied_records):
    manifest = {
        "batch_dirs": [str(p.resolve()) for p in batch_dirs],
        "merged_matrices": merged_records,
        "linked_sample_dirs": linked_records,
        "copied_reference_or_shared_files": copied_records,
        "notes": [
            "This is a reconstructed nf-core/rnaseq-like combined output.",
            "RDS/SummarizedExperiment objects are not regenerated by this script.",
            "Check pipeline_info/params.json and software_versions.yml across batches before trusting merged results.",
        ],
    }

    with open(out_dir / "batch_merge_manifest.json", "w") as handle:
        json.dump(manifest, handle, indent=2)


def main():
    parser = argparse.ArgumentParser(
        description="Merge batch nf-core/rnaseq outputs into one nf-core-like combined directory."
    )

    parser.add_argument(
        "-o",
        "--outdir",
        default="batch_all",
        help="Output directory. Default: batch_all",
    )

    parser.add_argument(
        "--copy",
        action="store_true",
        help="Copy sample directories instead of symlinking them.",
    )

    parser.add_argument(
        "--relaxed-features",
        action="store_true",
        help="Outer-join feature IDs instead of requiring identical feature rows/order.",
    )

    parser.add_argument(
        "batch_dirs",
        nargs="+",
        help="nf-core/rnaseq batch result directories.",
    )

    args = parser.parse_args()

    batch_dirs = [Path(p).resolve() for p in args.batch_dirs]
    out_dir = Path(args.outdir).resolve()

    for batch in batch_dirs:
        if not batch.is_dir():
            fail(f"Input batch directory does not exist or is not a directory: {batch}")

    if out_dir.exists():
        fail(f"Output directory already exists: {out_dir}")

    out_dir.mkdir(parents=True)

    info(f"Creating combined output: {out_dir}")
    info(f"Number of input batches: {len(batch_dirs)}")

    copy_metadata(batch_dirs, out_dir)
    copied_records = copy_reference_like_files(batch_dirs, out_dir)

    linked_records = link_sample_dirs(batch_dirs, out_dir, copy=args.copy)

    matrix_map = find_matrix_files(batch_dirs)
    merged_records = []

    if not matrix_map:
        warn("No known merged matrix files found. Sample directories/metadata were still linked/copied.")

    for rel, paths in sorted(matrix_map.items()):
        if len(paths) < 1:
            continue

        out_path = out_dir / rel

        info(f"Merging {rel} from {len(paths)} batch file(s)")
        shape = merge_matrices(
            paths,
            out_path,
            strict_features=not args.relaxed_features,
        )

        merged_records.append(
            {
                "relative_output": rel,
                "input_files": [str(p) for p in paths],
                "rows": int(shape[0]),
                "columns": int(shape[1]),
            }
        )

    write_manifest(out_dir, batch_dirs, merged_records, linked_records, copied_records)

    info("Done.")
    info(f"Manifest: {out_dir / 'batch_merge_manifest.json'}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import pandas as pd


TARGETS = [
    "star_salmon/salmon.merged.gene_counts.tsv",
    "star_salmon/salmon.merged.transcript_counts.tsv",
    "star_salmon/salmon.merged.gene_tpm.tsv",
    "star_salmon/salmon.merged.transcript_tpm.tsv",
]


def info(msg):
    print(f"[INFO] {msg}", file=sys.stderr)


def warn(msg):
    print(f"[WARN] {msg}", file=sys.stderr)


def fail(msg):
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(1)


def detect_annotation_columns(df):
    known_annotation_cols = {
        "tx",
        "gene_id",
        "gene_name",
        "gene_biotype",
        "gene_type",
        "transcript_id",
        "transcript_name",
        "transcript_biotype",
        "transcript_type",
        "Name",
        "target_id",
    }

    annotation_cols = []
    for col in df.columns:
        if col in known_annotation_cols:
            annotation_cols.append(col)
        else:
            break

    if not annotation_cols:
        annotation_cols = [df.columns[0]]

    sample_cols = [c for c in df.columns if c not in annotation_cols]
    if not sample_cols:
        fail(f"Could not detect sample columns. Columns found: {list(df.columns)}")

    return annotation_cols, sample_cols


def make_feature_key(df, annotation_cols):
    key = df[annotation_cols].astype(str).agg("||".join, axis=1)
    if key.duplicated().any():
        examples = key[key.duplicated()].head(10).tolist()
        fail(
            f"Duplicate feature keys detected using annotation columns {annotation_cols}. "
            f"Examples: {examples}"
        )
    return key


def read_batch_matrix(path):
    df = pd.read_csv(path, sep="\t", dtype=str)

    if df.shape[1] < 2:
        fail(f"Matrix has fewer than 2 columns: {path}")

    annotation_cols, sample_cols = detect_annotation_columns(df)
    feature_key = make_feature_key(df, annotation_cols)

    annotation = df[annotation_cols].copy()
    annotation.index = feature_key

    samples = df[sample_cols].copy()
    samples.index = feature_key

    for col in samples.columns:
        samples[col] = pd.to_numeric(samples[col], errors="ignore")

    return annotation, samples, annotation_cols


def merge_one_target(batch_dirs, rel_path, outdir, relaxed=False):
    input_files = []
    for batch_dir in batch_dirs:
        path = batch_dir / rel_path
        if path.is_file():
            input_files.append(path)
        else:
            warn(f"Missing {rel_path} in {batch_dir}")

    if not input_files:
        fail(f"No input files found for {rel_path}")

    info(f"Merging {rel_path} from {len(input_files)} batch files")

    merged_samples = None
    reference_annotation = None
    reference_index = None
    reference_annotation_cols = None
    seen_samples = set()

    for path in input_files:
        annotation, samples, annotation_cols = read_batch_matrix(path)

        duplicate_samples = seen_samples.intersection(samples.columns)
        if duplicate_samples:
            fail(
                f"Duplicate sample columns while merging {rel_path}:\n"
                f"{', '.join(sorted(duplicate_samples))}\n"
                f"Problem file: {path}"
            )

        seen_samples.update(samples.columns)

        if reference_annotation is None:
            reference_annotation = annotation
            reference_index = samples.index
            reference_annotation_cols = annotation_cols
            merged_samples = samples
        else:
            if annotation_cols != reference_annotation_cols:
                fail(
                    f"Annotation columns differ in {path}\n"
                    f"Expected: {reference_annotation_cols}\n"
                    f"Found:    {annotation_cols}"
                )

            if not relaxed:
                if not reference_index.equals(samples.index):
                    fail(
                        f"Feature rows/order differ in {path}\n"
                        f"Use --relaxed to outer-join by feature IDs/annotation."
                    )
                merged_samples = pd.concat([merged_samples, samples], axis=1)
            else:
                merged_samples = merged_samples.join(samples, how="outer")
                missing_keys = annotation.index.difference(reference_annotation.index)
                if len(missing_keys) > 0:
                    reference_annotation = pd.concat(
                        [reference_annotation, annotation.loc[missing_keys]],
                        axis=0,
                    )

    if relaxed:
        reference_annotation = reference_annotation.loc[merged_samples.index]

    output = pd.concat(
        [
            reference_annotation.reset_index(drop=True),
            merged_samples.reset_index(drop=True),
        ],
        axis=1,
    )

    out_path = outdir / rel_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(out_path, sep="\t", index=False)

    info(
        f"Wrote {out_path} "
        f"features={output.shape[0]} samples={merged_samples.shape[1]}"
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Merge STAR-Salmon count matrices from multiple batch result directories."
        )
    )

    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory for merged matrices.",
    )

    parser.add_argument(
        "--relaxed",
        action="store_true",
        help="Outer-join features instead of requiring identical rows/order.",
    )

    parser.add_argument(
        "batch_dirs",
        nargs="+",
        help="Batch result directories, e.g. outdir/batches/batch_000001",
    )

    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    batch_dirs = [Path(x).resolve() for x in args.batch_dirs]

    for batch_dir in batch_dirs:
        if not batch_dir.is_dir():
            fail(f"Batch directory does not exist: {batch_dir}")

    outdir.mkdir(parents=True, exist_ok=True)

    for rel_path in TARGETS:
        merge_one_target(
            batch_dirs=batch_dirs,
            rel_path=rel_path,
            outdir=outdir,
            relaxed=args.relaxed,
        )

    info("Done")


if __name__ == "__main__":
    main()

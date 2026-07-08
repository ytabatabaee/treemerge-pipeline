#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd
from ete3 import Tree


def parse_time_to_seconds(time_str):
    time_str = time_str.strip()

    days = 0
    if "-" in time_str:
        day_part, time_str = time_str.split("-", 1)
        days = int(day_part)

    parts = time_str.split(":")

    if len(parts) == 2:
        hours = 0
        minutes, seconds = parts
    elif len(parts) == 3:
        hours, minutes, seconds = parts
    else:
        raise ValueError(f"Cannot parse elapsed time: {time_str}")

    return (
        days * 86400
        + int(hours) * 3600
        + int(minutes) * 60
        + float(seconds)
    )


def parse_time_memory_file(path):
    text = path.read_text(errors="replace")

    elapsed_match = re.search(
        r"Elapsed \(wall clock\) time .*?:\s*([0-9:\.\-]+)",
        text,
    )
    mem_match = re.search(
        r"Maximum resident set size \(kbytes\):\s*(\d+)",
        text,
    )

    elapsed_seconds = (
        parse_time_to_seconds(elapsed_match.group(1))
        if elapsed_match
        else None
    )

    max_rss_kb = int(mem_match.group(1)) if mem_match else None

    return elapsed_seconds, max_rss_kb


def read_first_tree(path):
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                return line

    raise RuntimeError(f"No tree found in {path}")


def normalized_rf(true_tree_path, inferred_tree_path):
    true_tree_str = read_first_tree(true_tree_path)
    inferred_tree_str = read_first_tree(inferred_tree_path)

    true_tree = Tree(true_tree_str, format=1)
    inferred_tree = Tree(inferred_tree_str, format=1)

    common_taxa = sorted(
        set(true_tree.get_leaf_names())
        & set(inferred_tree.get_leaf_names())
    )

    if len(common_taxa) < 4:
        return None

    true_tree.prune(common_taxa, preserve_branch_length=True)
    inferred_tree.prune(common_taxa, preserve_branch_length=True)

    rf, max_rf, *_ = true_tree.robinson_foulds(
        inferred_tree,
        unrooted_trees=True,
    )

    if max_rf == 0:
        return 0.0

    return rf / max_rf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base_dir",
        default="/u/syt3/scratch/treemerge-data",
    )
    parser.add_argument(
        "--output",
        default="runtime_summary_by_N.tsv",
    )
    args = parser.parse_args()

    base = Path(args.base_dir)
    rows = []

    for time_file in base.glob("*/*/treemerge_estimatedgenetre/time_memory.txt"):
        n = int(time_file.parts[-4])
        rep_id = time_file.parts[-3]

        rep_dir = base / str(n) / rep_id
        true_tree = rep_dir / "s_tree.trees"

        outdir = time_file.parent
        treemerge_tree = outdir / "treemerge_species_tree.nwk"

        elapsed_seconds, max_rss_kb = parse_time_memory_file(time_file)

        treemerge_success = (
            treemerge_tree.exists()
            and treemerge_tree.stat().st_size > 0
        )

        rf_distance = None
        if treemerge_success and true_tree.exists() and true_tree.stat().st_size > 0:
            try:
                rf_distance = normalized_rf(true_tree, treemerge_tree)
            except Exception as e:
                print(
                    f"Warning: could not compute RF for "
                    f"N={n}, replicate={rep_id}: {e}"
                )

        rows.append(
            {
                "N": n,
                "replicate": rep_id,
                "elapsed_minutes": elapsed_seconds / 60
                if elapsed_seconds is not None
                else None,
                "max_rss_mb": max_rss_kb / 1024
                if max_rss_kb is not None
                else None,
                "normalized_rf": rf_distance,
                "treemerge_success": treemerge_success,
            }
        )

    df = pd.DataFrame(rows)

    if df.empty:
        raise RuntimeError(f"No time_memory.txt files found under {base}")

    successful = df[df["treemerge_success"]].copy()
    failed = df[~df["treemerge_success"]].copy()

    summary = (
        df.groupby("N")
        .agg(
            total_replicates=("replicate", "count"),
        )
        .reset_index()
    )

    successful_summary = (
        successful.groupby("N")
        .agg(
            average_elapsed_minutes=("elapsed_minutes", "mean"),
            average_max_rss_mb=("max_rss_mb", "mean"),
            average_normalized_rf=("normalized_rf", "mean"),
        )
        .reset_index()
    )

    failed_summary = (
        failed.groupby("N")
        .agg(
            failed_replicates=("replicate", "count"),
            failed_replicate_ids=("replicate", lambda reps: ",".join(sorted(reps))),
        )
        .reset_index()
    )

    summary = summary.merge(successful_summary, on="N", how="left")
    summary = summary.merge(failed_summary, on="N", how="left")

    summary["failed_replicates"] = summary["failed_replicates"].fillna(0).astype(int)
    summary["failed_replicate_ids"] = summary["failed_replicate_ids"].fillna("")

    summary = summary[
        [
            "N",
            "total_replicates",
            "average_elapsed_minutes",
            "average_max_rss_mb",
            "average_normalized_rf",
            "failed_replicates",
            "failed_replicate_ids",
        ]
    ]

    summary = summary.sort_values("N")
    summary.to_csv(args.output, sep="\t", index=False)

    print(summary.to_string(index=False))
    print(f"\nWrote summary to: {args.output}")


if __name__ == "__main__":
    main()

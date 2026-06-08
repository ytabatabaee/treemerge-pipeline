#!/usr/bin/env python3

import argparse
import multiprocessing as mp
import os
import subprocess
import sys
from pathlib import Path

from ete3 import Tree


def run(cmd, log_file=None):
    cmd = [str(x) for x in cmd]
    print("\n[RUN]")
    print(" ".join(cmd))
    sys.stdout.flush()

    if log_file:
        with open(log_file, "w") as log:
            result = subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
    else:
        result = subprocess.run(cmd)

    if result.returncode != 0:
        raise RuntimeError(f"Command failed:\n{' '.join(cmd)}")


def mkdir(path):
    os.makedirs(path, exist_ok=True)


def run_astrid(astrid_bin, gene_trees, agid_matrix, log_file):
    run([astrid_bin, "-i", gene_trees, "-c", agid_matrix], log_file)


def infer_starting_tree(fastme_bin, agid_matrix, starting_tree, log_file):
    run([fastme_bin, "-i", agid_matrix, "-o", starting_tree], log_file)


def centroid_bisect_taxa(component_taxa, starting_tree_file, min_subset_size=5):
    """
    PASTA-style centroid bisection:
    choose an edge whose deletion gives the most balanced valid split.
    """

    t = Tree(starting_tree_file, format=1)
    t.prune(list(component_taxa), preserve_branch_length=True)

    all_taxa = set(t.get_leaf_names())

    best_split = None
    best_score = None

    for node in t.traverse("postorder"):
        if node.is_root():
            continue

        side1 = set(node.get_leaf_names())
        side2 = all_taxa - side1

        if len(side1) < min_subset_size:
            continue
        if len(side2) < min_subset_size:
            continue

        score = abs(len(side1) - len(side2))

        if best_score is None or score < best_score:
            best_score = score
            best_split = (side1, side2)

    return best_split


def decompose_like_original_helper(
    starting_tree_file,
    max_subset_size,
    min_subset_size=5,
):
    """
    Reimplementation of tools/build_subsets_from_tree.py behavior:

    - repeatedly bisect using centroid-style edge cuts
    - keep splitting subsets larger than max_subset_size
    - require each final subset to have at least 5 taxa
    """

    t = Tree(starting_tree_file, format=1)
    all_taxa = set(t.get_leaf_names())

    next_components = [all_taxa]
    done_components = []

    while next_components:
        current_components = next_components
        next_components = []

        for component in current_components:
            split = centroid_bisect_taxa(
                component,
                starting_tree_file,
                min_subset_size,
            )

            if split is None:
                if len(component) <= max_subset_size:
                    if len(component) >= min_subset_size:
                        done_components.append(component)
                    else:
                        raise RuntimeError(
                            f"Component has fewer than {min_subset_size} taxa."
                        )
                else:
                    raise RuntimeError(
                        f"Could not split component of size {len(component)}. "
                        f"Try increasing --max_subset_size."
                    )

                continue

            side1, side2 = split

            for side_name, side in [("T1", side1), ("T2", side2)]:
                if len(side) > max_subset_size:
                    next_components.append(side)
                else:
                    if len(side) >= min_subset_size:
                        done_components.append(side)
                    else:
                        raise RuntimeError(
                            f"{side_name} has fewer than {min_subset_size} leaves."
                        )

    return done_components


def write_subset_files(subsets, subset_dir):
    mkdir(subset_dir)

    n = len(subsets)
    pad = len(str(n))
    subset_records = []

    for i, taxa in enumerate(subsets, start=1):
        j = str(i).zfill(pad)
        subset_file = Path(subset_dir) / f"subset_{j}-outof-{n}.txt"

        with open(subset_file, "w") as out:
            out.write("\n".join(sorted(taxa)) + "\n")

        subset_records.append((str(subset_file), taxa))

    return subset_records


def prune_tree_to_subset(tree_str, subset_taxa, min_taxa):
    t = Tree(tree_str, format=1)

    keep = [
        taxon for taxon in t.get_leaf_names()
        if taxon in subset_taxa
    ]

    if len(keep) < min_taxa:
        return None

    t.prune(keep, preserve_branch_length=True)

    return t.write(format=1)


def build_subset_gene_trees(
    gene_trees,
    subset_taxa,
    outfile,
    min_taxa,
):
    count = 0

    with open(gene_trees) as infile, open(outfile, "w") as out:
        for line in infile:
            tree_str = line.strip()

            if not tree_str:
                continue

            pruned = prune_tree_to_subset(
                tree_str,
                subset_taxa,
                min_taxa,
            )

            if pruned is not None:
                out.write(pruned + "\n")
                count += 1

    return count


def run_astral4_subset(job):
    astral4_bin, subset_gene_trees, output_tree, log_file = job

    run(
        [
            astral4_bin,
            "-i",
            subset_gene_trees,
            "-o",
            output_tree,
        ],
        log_file,
    )

    return output_tree


def run_treemerge(
    treemerge_script,
    paup,
    starting_tree,
    subset_species_trees,
    agid_matrix,
    taxlist,
    output_tree,
    workdir,
    log_file,
):
    cmd = [
        sys.executable,
        treemerge_script,
        "-s",
        starting_tree,
        "-t",
    ]

    cmd.extend(subset_species_trees)

    cmd.extend(
        [
            "-m",
            agid_matrix,
            "-x",
            taxlist,
            "-o",
            output_tree,
            "-w",
            workdir,
            "-p",
            paup,
        ]
    )

    run(cmd, log_file)


def score_species_tree(
    astral4_bin,
    species_tree,
    gene_trees,
    output_tree,
    log_file,
):
    run(
        [
            astral4_bin,
            "-q",
            species_tree,
            "-i",
            gene_trees,
            "-o",
            output_tree,
        ],
        log_file,
    )


def main():
    parser = argparse.ArgumentParser(
        description="TreeMerge divide-and-conquer species tree pipeline."
    )

    parser.add_argument("--gene_trees", required=True)
    parser.add_argument("--outdir", required=True)

    parser.add_argument("--astrid_bin", required=True)
    parser.add_argument("--fastme_bin", required=True)
    parser.add_argument("--astral4_bin", required=True)
    parser.add_argument("--treemerge_script", required=True)
    parser.add_argument("--paup", required=True)

    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--max_subset_size", type=int, default=100)

    parser.add_argument(
        "--min_subset_size",
        type=int,
        default=5,
        help="Minimum subset size. Original helper requires at least 5 leaves.",
    )

    parser.add_argument(
        "--min_taxa_per_gene_tree",
        type=int,
        default=4,
        help="Minimum taxa required after pruning a gene tree.",
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    intermediate = outdir / "intermediate"
    logs = outdir / "logs"

    mkdir(outdir)
    mkdir(intermediate)
    mkdir(logs)

    agid_matrix = intermediate / "agid_matrix.txt"
    taxlist = str(agid_matrix) + "_taxlist"
    starting_tree = intermediate / "starting_tree.nwk"

    print("\n========== Step 1: ASTRID AGID matrix ==========")

    run_astrid(
        args.astrid_bin,
        args.gene_trees,
        str(agid_matrix),
        logs / "astrid.log",
    )

    print("\n========== Step 2: FastME starting tree ==========")

    infer_starting_tree(
        args.fastme_bin,
        str(agid_matrix),
        str(starting_tree),
        logs / "fastme.log",
    )

    print("\n========== Step 3: Centroid subset decomposition ==========")

    subsets = decompose_like_original_helper(
        str(starting_tree),
        args.max_subset_size,
        args.min_subset_size,
    )

    subset_dir = intermediate / "subsets"
    subset_records = write_subset_files(subsets, subset_dir)

    print(f"Generated {len(subset_records)} subsets")

    subset_sizes_file = intermediate / "subset_sizes.tsv"
    with open(subset_sizes_file, "w") as out:
        out.write("subset_file\tsize\n")
        for subset_file, taxa in subset_records:
            print(f"{Path(subset_file).name}: {len(taxa)} taxa")
            out.write(f"{subset_file}\t{len(taxa)}\n")

    print("\n========== Step 4: Prune gene trees to subsets ==========")

    subset_gene_tree_dir = intermediate / "subset_gene_trees"
    mkdir(subset_gene_tree_dir)

    subset_gene_tree_files = []

    for i, (subset_file, subset_taxa) in enumerate(subset_records, start=1):
        outfile = subset_gene_tree_dir / f"subset_{i}.trees"

        count = build_subset_gene_trees(
            args.gene_trees,
            subset_taxa,
            str(outfile),
            args.min_taxa_per_gene_tree,
        )

        print(f"[Subset {i}] {count} pruned gene trees")

        if count == 0:
            raise RuntimeError(
                f"Subset {i} produced zero usable gene trees. "
                f"Subset taxa file: {subset_file}"
            )

        subset_gene_tree_files.append(str(outfile))

    print("\n========== Step 5: ASTRAL4 subset trees ==========")

    subset_species_tree_dir = intermediate / "subset_species_trees"
    mkdir(subset_species_tree_dir)

    jobs = []
    subset_species_trees = []

    for i, subset_gene_trees in enumerate(subset_gene_tree_files, start=1):
        output_tree = subset_species_tree_dir / f"subset_{i}.nwk"
        log_file = logs / f"astral4_subset_{i}.log"

        jobs.append(
            (
                args.astral4_bin,
                subset_gene_trees,
                str(output_tree),
                str(log_file),
            )
        )

        subset_species_trees.append(str(output_tree))

    with mp.Pool(args.threads) as pool:
        pool.map(run_astral4_subset, jobs)

    print("\n========== Step 6: TreeMerge ==========")

    merged_tree = outdir / "treemerge_species_tree.nwk"

    run_treemerge(
        args.treemerge_script,
        args.paup,
        str(starting_tree),
        subset_species_trees,
        str(agid_matrix),
        taxlist,
        str(merged_tree),
        str(intermediate),
        logs / "treemerge.log",
    )

    print("\n========== Step 7: Final ASTRAL4 scoring ==========")

    final_tree = outdir / "final_species_tree_scored.nwk"

    score_species_tree(
        args.astral4_bin,
        str(merged_tree),
        args.gene_trees,
        str(final_tree),
        logs / "astral4_final_scoring.log",
    )

    print("\n===================================")
    print("PIPELINE COMPLETED")
    print("===================================")
    print(f"\nTreeMerge tree:\n{merged_tree}")
    print(f"\nFinal scored species tree:\n{final_tree}")


if __name__ == "__main__":
    main()

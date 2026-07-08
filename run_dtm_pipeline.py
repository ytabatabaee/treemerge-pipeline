#!/usr/bin/env python3

import argparse
import multiprocessing as mp
import os
import platform
import subprocess
import sys
import tempfile
from pathlib import Path

from ete3 import Tree


REPO_ROOT = Path(__file__).resolve().parent
SOFTWARE_DIR = REPO_ROOT / "software"

LINUX_X86_64_DEFAULTS = {
    "astrid_bin": SOFTWARE_DIR / "ASTRID-linux",
    "fastme_bin": SOFTWARE_DIR / "fastme-2.1.5-linux64",
    "astral4_bin": SOFTWARE_DIR / "astral4-linux",
    "treemerge_script": SOFTWARE_DIR / "treemerge.py",
    "paup": [
        SOFTWARE_DIR / "paup4a169_ubuntu64",
        SOFTWARE_DIR / "paup4a168_centos64",
    ],
}

DARWIN_COMMON_DEFAULTS = {
    "astrid_bin": SOFTWARE_DIR / "ASTRID-osx",
    "fastme_bin": SOFTWARE_DIR / "fastme-2.1.5-osx",
    "treemerge_script": SOFTWARE_DIR / "treemerge.py",
    "paup": SOFTWARE_DIR / "paup4a168_osx",
}


def bundled_defaults():
    system = platform.system()
    machine = platform.machine().lower()

    if system == "Linux" and machine in {"x86_64", "amd64"}:
        return LINUX_X86_64_DEFAULTS

    if system == "Darwin":
        defaults = dict(DARWIN_COMMON_DEFAULTS)
        if machine == "arm64":
            defaults["astral4_bin"] = SOFTWARE_DIR / "astral4-osx"
        return defaults

    return {}


def executable_label(name):
    return name.replace("_", "-")


def paup_smoke_test(path):
    with tempfile.NamedTemporaryFile(
        "w",
        suffix=".nex",
        prefix="treemerge-paup-check-",
        delete=False,
    ) as handle:
        handle.write("#NEXUS\nbegin paup;\nq;\nend;\n")
        check_file = handle.name

    try:
        result = subprocess.run(
            [str(path), "-n", check_file],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=15,
        )
        return result.returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False
    finally:
        try:
            os.unlink(check_file)
        except OSError:
            pass


def resolve_paup_candidate(candidates):
    if not isinstance(candidates, list):
        candidates = [candidates]

    existing = []
    for candidate in candidates:
        path = Path(candidate).expanduser()
        if path.exists() and os.access(path, os.X_OK):
            existing.append(path)
            if paup_smoke_test(path):
                return path

    if existing:
        return existing[0]

    return Path(candidates[0]).expanduser()


def resolve_tool_paths(args):
    defaults = bundled_defaults()

    for name in [
        "astrid_bin",
        "fastme_bin",
        "astral4_bin",
        "treemerge_script",
        "paup",
    ]:
        value = getattr(args, name)
        if value is None and name in defaults:
            if name == "paup":
                value = resolve_paup_candidate(defaults[name])
            else:
                value = defaults[name]

        if value is not None:
            setattr(args, name, str(Path(value).expanduser()))

    missing = [
        name
        for name in [
            "astrid_bin",
            "fastme_bin",
            "astral4_bin",
            "treemerge_script",
            "paup",
        ]
        if getattr(args, name) is None
    ]

    if missing:
        labels = ", ".join(f"--{executable_label(name)}" for name in missing)
        system = platform.system() or "unknown"
        machine = platform.machine() or "unknown"
        raise RuntimeError(
            f"No bundled defaults are available for {system} {machine}: {labels}. "
            "Pass these paths explicitly."
        )

    validate_tool_paths(args)


def validate_tool_paths(args):
    for name in [
        "astrid_bin",
        "fastme_bin",
        "astral4_bin",
        "paup",
    ]:
        path = Path(getattr(args, name))
        if not path.exists():
            raise RuntimeError(f"--{executable_label(name)} does not exist: {path}")
        if not os.access(path, os.X_OK):
            raise RuntimeError(
                f"--{executable_label(name)} is not executable: {path}\n"
                f"Run: chmod +x {path}"
            )

    treemerge_script = Path(args.treemerge_script)
    if not treemerge_script.exists():
        raise RuntimeError(f"--treemerge-script does not exist: {treemerge_script}")
    if not treemerge_script.is_file():
        raise RuntimeError(f"--treemerge-script is not a file: {treemerge_script}")


def run(cmd, log_file=None):
    cmd = [str(x) for x in cmd]

    print("\n[RUN]")
    print(" ".join(cmd))
    sys.stdout.flush()

    if log_file:
        log_file = str(log_file)

        with open(log_file, "w") as log:
            result = subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
            )

        if result.returncode != 0:
            try:
                with open(log_file) as f:
                    log_contents = f.read()
            except OSError:
                log_contents = "(Could not read log file.)"

            raise RuntimeError(
                f"Command failed with exit code {result.returncode}:\n"
                f"{' '.join(cmd)}\n\n"
                f"Log file: {log_file}\n\n"
                f"Log contents:\n{log_contents}"
            )

    else:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"Command failed with exit code {result.returncode}:\n"
                f"{' '.join(cmd)}\n\n"
                f"STDOUT:\n{result.stdout}\n\n"
                f"STDERR:\n{result.stderr}"
            )


def mkdir(path):
    os.makedirs(path, exist_ok=True)


def require_nonempty_file(path, label):
    path = Path(path)
    if not path.exists():
        raise RuntimeError(f"{label} was not created: {path}")
    if path.stat().st_size == 0:
        raise RuntimeError(f"{label} is empty: {path}")


def run_astrid(astrid_bin, gene_trees, agid_matrix, log_file):
    run([astrid_bin, "-i", gene_trees, "-c", agid_matrix], log_file)
    require_nonempty_file(agid_matrix, "ASTRID AGID matrix")


def infer_starting_tree(fastme_bin, agid_matrix, starting_tree, log_file):
    run([fastme_bin, "-i", agid_matrix, "-o", starting_tree, "-mN"], log_file)
    require_nonempty_file(starting_tree, "FastME starting tree")


def centroid_bisect_taxa(component_taxa, starting_tree_file, min_subset_size=5):
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

    require_nonempty_file(subset_gene_trees, "ASTRAL4 subset input gene tree file")

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

    require_nonempty_file(output_tree, "ASTRAL4 subset output tree")

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
    for tree in subset_species_trees:
        require_nonempty_file(tree, "TreeMerge input subset species tree")

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
    require_nonempty_file(output_tree, "TreeMerge output tree")


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
            "-C -c",
            species_tree,
            "-i",
            gene_trees,
            "-o",
            output_tree,
        ],
        log_file,
    )

    require_nonempty_file(output_tree, "Final ASTRAL4 scored tree")


def main():
    parser = argparse.ArgumentParser(
        description="TreeMerge divide-and-conquer species tree pipeline."
    )

    parser.add_argument("--gene_trees", required=True)
    parser.add_argument("--outdir", required=True)

    parser.add_argument(
        "--astrid-bin",
        "--astrid_bin",
        dest="astrid_bin",
        help="Path to ASTRID.",
    )
    parser.add_argument(
        "--fastme-bin",
        "--fastme_bin",
        dest="fastme_bin",
        help="Path to FastME.",
    )
    parser.add_argument(
        "--astral4-bin",
        "--astral4_bin",
        dest="astral4_bin",
        help="Path to ASTRAL4.",
    )
    parser.add_argument(
        "--treemerge-script",
        "--treemerge_script",
        dest="treemerge_script",
        help="Path to treemerge.py.",
    )
    parser.add_argument(
        "--paup",
        help="Path to PAUP*.",
    )

    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--max_subset_size", type=int, default=10)

    parser.add_argument(
        "--parallel_astral_subsets",
        action="store_true",
        help="Run ASTRAL4 subset jobs in parallel. Default is sequential for safer debugging.",
    )

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
    resolve_tool_paths(args)

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
                f"Subset {i} produced zero usable gene trees.\n"
                f"Subset taxa file: {subset_file}\n"
                f"Output file: {outfile}"
            )

        require_nonempty_file(outfile, f"Subset {i} pruned gene tree file")
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

    if args.parallel_astral_subsets:
        print(f"Running ASTRAL4 subset jobs in parallel with {args.threads} workers")
        with mp.Pool(args.threads) as pool:
            pool.map(run_astral4_subset, jobs)
    else:
        print("Running ASTRAL4 subset jobs sequentially")
        for job in jobs:
            run_astral4_subset(job)

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

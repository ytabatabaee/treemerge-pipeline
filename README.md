Divide-and-conquer species tree estimation pipeline using TreeMerge
=======
This repository provides a fully automated implementation of a divide-and-conquer species tree estimation pipeline based on [TreeMerge](https://github.com/ekmolloy/treemerge) introduced in [Molloy and Warnow, Bioinformatics 2019](https://academic.oup.com/bioinformatics/article/35/14/i417/5529167). This pipeline is based on the [Trees in the desert tutorial](https://github.com/ekmolloy/trees-in-the-desert-tutorial). 

The pipeline automates the full workflow using a single command-line interface.

1. AGID matrix estimation from gene trees using ASTRID
2. Starting tree estimation using Neighbor Joining (as implemented in FastME)
3. Recursive taxon decomposition 
4. Subset species tree estimation using ASTRAL4
5. TreeMerge for merging subset trees
6. Final branch support and branch length estimation using ASTRAL4

Centroid decomposition provenance
---------------------------------
The subset decomposition follows the TreeMerge experimental centroid
decomposer used by Molloy and Warnow (2019). The reference routine is the
DendroPy implementation in the cited experimental scripts: it first treats
the guide tree as unrooted, resolves polytomies, collapses the basal
bifurcation, and updates bipartitions. While a component has more than the
requested maximum size, it removes the first minimum-balance edge encountered
in postorder and recursively processes the two resulting components. The
components are physically disjoint, and recursion order determines output
order. There is no minimum subset-size filter and no support filter.

The previous implementation used ete3, filtered both sides of every candidate
edge by `min_subset_size=5`, performed breadth-first component processing, and
raised `Could not split component of size 11` when a valid original split had a
side smaller than five taxa. The current `--min_subset_size` option remains
only for command-line compatibility and does not affect decomposition.

The centroid core depends on DendroPy and its bipartition bitmasks; the
original command-line wrapper additionally used Click and Biopython for
alignment partitioning. The current pipeline adapts the centroid core to
DendroPy 4.6.1 (the original scripts used DendroPy 4.3-era APIs). The
checked-in regression reference preserves the
original edge scan, mutation-based split, stopping condition, and recursion
order. It also tests the 11-taxon case that previously failed. The original
TreeMerge merger and ASTRAL4 settings are unchanged.

The local validation interpreter exposed DendroPy 4.5.2; the project
requirements remain pinned to 4.6.1, and the tested API calls are shared by
both versions.

The decomposition source was compared with the [Illinois Data Bank TreeMerge
dataset](https://databank.illinois.edu/datasets/IDB-9570561) and the [original
TreeMerge repository](https://github.com/ekmolloy/treemerge). The TreeMerge
repository contains the merger; the experimental centroid routine is in the
cited experimental scripts and is preserved in the original centroid mode of
QuickScripts. A later QuickScripts revision added a support threshold mode;
that later option is deliberately excluded here because it is not part of the
original centroid mode used for these experiments. The Data Bank page and file
listing were available for inspection, but its direct file endpoint returned
HTTP 403 in this environment, so no Data Bank archive is vendored in this
repository.

For the 100-taxon ASTRAL-III-SU experiments, keep `--max_subset_size 10`.
The supplied SLURM script already uses 10 for `N=2000` and does not change the
task mapping or runtime measurement.


Installation
------------
Recommended installation with Python `venv` and `pip`:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Alternatively, use can microconda locally in your project directory.
```bash
cd ~/treemerge-pipeline

# download micromamba
mkdir -p ~/bin
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C ~/bin --strip-components=1 bin/micromamba

# initialize shell
~/bin/micromamba shell init -s bash -r ~/micromamba
source ~/.bashrc
```

Then create the environment 

```bash
micromamba create -n treemerge python=3.9 -y
micromamba activate treemerge

python --version
python -m pip install --upgrade pip setuptools wheel
python -m pip install -r requirements.txt
```

Check that the Python packages and bundled external tools are correctly installed:

```bash
python check_install.py
```

This code reports which bundled tools are available and which paths need to be supplied manually.

Usage
-----------
On Linux x86_64 or Apple Silicon macOS with the bundled tools:

```bash
python run_dtm_pipeline.py \
    --gene_trees example/estimated_gene_trees.txt \
    --outdir results \
    --threads 4 \
    --max_subset_size 10
```

If a bundled binary is not available for your platform, or if you want to use a
different installation of a tool, override the path explicitly:

```bash
python run_dtm_pipeline.py \
    --gene_trees example/estimated_gene_trees.txt \
    --outdir results \
    --astrid-bin /path/to/ASTRID \
    --fastme-bin /path/to/fastme \
    --astral4-bin /path/to/astral4 \
    --treemerge-script software/treemerge.py \
    --paup /path/to/paup \
    --threads 4 \
    --max_subset_size 10
```

Run the decomposition regression tests with:

```bash
python -m unittest discover -s tests -v
```

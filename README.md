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
    --max_subset_size 100
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
    --max_subset_size 100
```

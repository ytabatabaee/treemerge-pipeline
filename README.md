Divide-and-conquer species tree estimation pipeline using TreeMerge
=======
This repository provides a fully automated implementation of a divide-and-conquer species tree estimation pipeline based on [TreeMerge](https://github.com/ekmolloy/treemerge) introduced in [Molloy and Warnow, Bioinformatics 2019](https://academic.oup.com/bioinformatics/article/35/14/i417/5529167). This pipeline is based on the [Trees in the desert tutorial](https://github.com/ekmolloy/trees-in-the-desert-tutorial). 

The pipeline automates the full workflow using a single command-line interface.

1. AGID matrix estimation from gene trees using ASTRID
2. Starting tree estimation using FastME
3. Recursive taxon decomposition 
4. Subset species tree estimation using ASTRAL4
5. TreeMerge for merging subset trees
6. Final branch support and branch length estimation using ASTRAL4


INSTALLATION
------------
Use a fresh environment. Do not run the pipeline from the conda `base`
environment.

Recommended installation with Python `venv` and `pip`:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

If you use conda, create a fresh environment instead of updating an existing one:

```bash
conda config --set channel_priority strict
conda env remove -n treemerge -y
conda env create -f environment.yml
conda activate treemerge
```

Avoid these commands for initial installation:

```bash
conda env update -f environment.yml
conda install -n base -c conda-forge mamba -y
```

Both commands can make conda spend a long time reconciling old packages in an
existing environment. If conda solving is slow on your machine, use the `venv`
instructions above.

If conda is not available, use microconda locally in your project directory.
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

Check that the Python packages and bundled external tools are available:

```bash
python check_install.py
```

On Linux x86_64 and Apple Silicon macOS, the pipeline detects the bundled ASTRID,
FastME, ASTRAL4, TreeMerge, and PAUP* tools automatically. For PAUP* on Linux, it tries
`software/paup4a169_ubuntu64` first and falls back to
`software/paup4a168_centos64` if the Ubuntu binary does not pass a smoke test.
The bundled `software/astral4-osx` binary is arm64, so Intel macOS users need to
provide a compatible ASTRAL4 binary with `--astral4-bin` or run on Linux.

On other platforms, `check_install.py` reports which bundled tools are available
and which paths need to be supplied manually.

USAGE
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

TreeMerge Divide-and-Conquer Species Tree Estimation Pipeline
=======
This repository provides a fully automated implementation of a divide-and-conquer species tree estimation pipeline based on [TreeMerge](https://github.com/ekmolloy/treemerge) introduced in [Molloy and Warnow, Bioinformatics 2019](https://academic.oup.com/bioinformatics/article/35/14/i417/5529167). This pipeline is based on the [Trees in the desert tutorial](https://github.com/ekmolloy/trees-in-the-desert-tutorial). 

The pipeline automates the full workflow using a single command-line interface.

1. AGID matrix estimation from gene trees using ASTRID
2. Starting tree estimation using FastME
3. Recursive taxon decomposition 
4. Subset species tree estimation using ASTRAL4
5. TreeMerge for merging subset trees
6. Final branch support and branch length estimation using ASTRAL4


INSTALATION
------------
```
conda create -n treemerge python=3.9 -y
conda activate treemerge
conda env update -f environment.yml
```

USAGE
-----------
```bash
python run_dtm_pipeline.py \
    --gene_trees example/estimated_gene_trees.txt \
    --outdir results \
    --astrid_bin software/ASTRID-linux \
    --astral4_bin software/astral4-linux \
    --fastme_bin software/fastme-2.1.5-linux64 \
    --treemerge_script software/treemerge.py \
    --paup "$(pwd)/software/paup4a168_centos64" \
    --threads 4 \
    --max_subset_size 100
```


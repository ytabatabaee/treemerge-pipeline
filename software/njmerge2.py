"""
This file is a python prototype of NJMerge-2 from the paper:

Molloy, E.K., Warnow, T. (2019). TreeMerge: A new method for improving the
scalability of species tree estimation methods.

Copyright (c) 2018-19 Erin K. Molloy
All rights reserved.

License: 3-Clause BSD,
see https://opensource.org/licenses/BSD-3-Clause
"""
import njmergepair

import argparse
from copy import deepcopy
import dendropy
import networkx
import numpy
import os
import os.path
import sys
import time


def tree_to_mst(strefile, treefiles):
    """
    Derives the spanning tree on subsets from a starting tree using
    topological distances between one representative taxon per subset.

    One representative leaf is selected from each subset tree. The starting
    tree is pruned to contain only those representatives, all branch lengths
    are set to 1 (so distances are purely topological), and a phylogenetic
    distance matrix is computed. The minimum spanning tree of that matrix
    gives the correct adjacency structure between subsets as implied by the
    guide tree topology.

    Using a single representative per subset (rather than labelling all
    internal nodes via Fitch's algorithm) ensures that edge weights reflect
    true inter-subset distances, yielding a unique and topologically correct
    MST.  When all edge weights are equal (as in the Fitch approach) the MST
    is not unique and NetworkX may return a star or other structure that does
    not match the guide tree, causing incorrect merge ordering and degraded
    RF distances.

    Parameters
    ----------
    strefile : file containing starting tree on full taxon set (newick format)
    treefiles : a list of files containing subset trees (newick format)

    Returns
    -------
    adjacency matrix for minimum spanning tree (MST) connecting subset tree files
    """
    taxa = dendropy.TaxonNamespace()
    stre = dendropy.Tree.get(path=strefile,
                             schema="newick",
                             taxon_namespace=taxa)

    # Pick one representative leaf label per subset and record its subset index
    ltos = {}
    for i, f in enumerate(treefiles):
        ti = dendropy.Tree.get(path=f,
                               schema="newick",
                               taxon_namespace=taxa)
        li = ti.leaf_nodes()[0].taxon.label
        ltos[li] = i

    labs = list(ltos.keys())

    # Prune starting tree to representatives only
    stre.retain_taxa_with_labels(labs)
    stre.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    # Set all branch lengths to 1 so distances are purely topological
    edges = [e for e in stre.preorder_edge_iter()]
    for e in edges[1:]:
        e.length = 1.0

    # Build pairwise distance matrix between representatives
    pdm = stre.phylogenetic_distance_matrix()
    n = len(labs)
    dmat = numpy.zeros((n, n))
    for li in labs:
        i = ltos[li]
        xi = pdm.taxon_namespace.get_taxon(li)
        for lj in labs:
            j = ltos[lj]
            xj = pdm.taxon_namespace.get_taxon(lj)
            dmat[i, j] = pdm.patristic_distance(xi, xj)

    # Turn distance matrix into graph and find MST
    graph = networkx.from_numpy_array(dmat)
    mst = networkx.minimum_spanning_tree(graph)

    # Return MST as an adjacency matrix
    amat = networkx.to_numpy_array(mst)
    amat[amat > 0] = 1

    return amat


def get_leaf_list(subtree):
    """Return list of leaf labels

    Parameters
    ----------
    subtree : dendropy tree or node object

    Returns
    -------
    list of str

    """
    return [l.taxon.label for l in subtree.leaf_nodes()]


def njmerge2(dmatfile, taxafile, treefiles, mstmat):
    """
    Parameters
    ----------
    dmatfile : file containing distance matrix on full leaf set (phylip format)
    taxafile : file containing taxon labels for rows in distance matrix
    treefiles : a list of files containing subset trees (newick format)
    mstmat : adjacency matrix for minimum spanning tree on subsetes

    Results
    -------
    tree on full leaf set
    """
    # Turn matrix into graph
    graph = networkx.Graph(mstmat)

    # Merge trees
    final_tree = None
    added = []
    nodes = list(deepcopy(graph.nodes()))
    root = None
    next_roots = [nodes[0]]
    while len(nodes) > 0:
        if root is None:
            root = next_roots[0]
        neighbors = graph.neighbors(root)
        next_roots = list(set(next_roots).union(set(neighbors)))
        if len(neighbors) == 0:
            nodes.remove(root)
            next_roots.remove(root)
            root = None
        else:
            sys.stdout.write("Combining %d and %d...\n" % (root, neighbors[0]))

            if root < neighbors[0]:
                i = root
                j = neighbors[0]
            else:
                i = neighbors[0]
                j = root

            if final_tree is None:
                ti = dendropy.Tree.get(path=treefiles[i], schema="newick")
                tj = dendropy.Tree.get(path=treefiles[j], schema="newick")
                added.append(i)
                added.append(j)
            else:
                if i in added:
                    ti = final_tree
                    tj = dendropy.Tree.get(path=treefiles[j], schema="newick")
                    added.append(j)
                else:
                    ti = dendropy.Tree.get(path=treefiles[i], schema="newick")
                    tj = final_tree
                    added.append(i)

            [dij, tij] = njmergepair.run(dmatfile, taxafile, ti, tj)
            final_tree = tij

            sys.stdout.write("...combined tree has %d leaves!\n"
                             % len(get_leaf_list(final_tree)))

            graph.remove_edge(root, neighbors[0])
    return final_tree


def main(args):
    output = args.output
    if os.path.isfile(output):
        sys.exit(output + " already exists!\n")

    sys.stdout.write("Running NJMerge2...\n")
    total_ts = time.time()

    # Input
    strefile = args.start
    treefiles = args.trees
    dmatfile = args.matrix
    taxafile = args.taxa

    # Build MST from starting tree
    ts = time.time()
    mst = tree_to_mst(strefile, treefiles)
    print(mst)
    rt = time.time() - ts
    sys.stdout.write("...computed MST in %d seconds.\n" % rt)

    # Combine pairs of trees using NJMerge in order
    # indicated by the MST
    ts = time.time()
    final_tree = njmerge2(dmatfile, taxafile, treefiles, mst)
    rt = time.time() - ts
    sys.stdout.write("...combined all trees in %d seconds.\n" % rt)

    with open (output, 'w') as f:
        f.write(final_tree.as_string(schema="newick")[5:])

    total = time.time() - total_ts
    sys.stdout.write("Total execution time: %d seconds\n" % total)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--start", type=str,
                        help="Input starting tree file", required=True)

    parser.add_argument("-t", "--trees", type=str, nargs="+",
                        help="Input subset tree files", required=True)

    parser.add_argument("-m", "--matrix", type=str,
                        help="Input distance matrix", required=True)

    parser.add_argument("-x", "--taxa", type=str,
                        help="Taxon names corresponding to rows in distance matrix, necessary when using output from FastME",
                        required=False)

    parser.add_argument("-o", "--output", type=str,
                        help="Output tree",
                        required=True)

    main(parser.parse_args())

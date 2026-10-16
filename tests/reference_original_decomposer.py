"""Minimal copy of the original TreeMerge-era centroid decomposer."""

import dendropy
from dendropy.utility import bitprocessing


def _centroid_edge(tree):
    n = bitprocessing.num_set_bits(tree.seed_node.tree_leafset_bitmask)
    best_balance = float("inf")
    for edge in tree.postorder_edge_iter():
        if edge.tail_node is None:
            continue
        balance = abs(
            n / 2.0
            - bitprocessing.num_set_bits(edge.bipartition.leafset_bitmask)
        )
        if balance < best_balance:
            best_balance = balance
            best_edge = edge
    return best_edge


def _split(tree, edge):
    new_root = edge.head_node
    edge.tail_node.remove_child(new_root)
    new_tree = dendropy.Tree(
        seed_node=new_root,
        taxon_namespace=tree.taxon_namespace,
    )
    tree.update_bipartitions()
    new_tree.update_bipartitions()
    return tree, new_tree


def decompose_tree(tree, max_subset_size):
    n = len(tree.leaf_nodes())
    if n <= max_subset_size:
        return [[leaf.taxon.label for leaf in tree.leaf_node_iter()]]
    tree1, tree2 = _split(tree, _centroid_edge(tree))
    return decompose_tree(tree1, max_subset_size) + decompose_tree(
        tree2, max_subset_size
    )


def read_and_decompose(path, max_subset_size):
    tree = dendropy.Tree.get(path=path, schema="newick")
    tree.is_rooted = False
    tree.resolve_polytomies(limit=2)
    tree.collapse_basal_bifurcation()
    tree.update_bipartitions()
    return decompose_tree(tree, max_subset_size)

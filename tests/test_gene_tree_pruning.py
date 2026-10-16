import unittest

import dendropy

from run_dtm_pipeline import prune_tree_to_subset


class GeneTreePruningTests(unittest.TestCase):
    def test_numeric_taxa_and_internal_support_are_distinct(self):
        tree = "((1:0.1,2:0.2)0.95:0.3,3:0.4);"

        pruned = prune_tree_to_subset(tree, {"1", "2"}, min_taxa=2)
        parsed = dendropy.Tree.get(data=pruned, schema="newick")

        self.assertEqual(
            {leaf.taxon.label for leaf in parsed.leaf_node_iter()}, {"1", "2"}
        )
        self.assertIn(
            "0.95",
            [node.label for node in parsed.internal_nodes() if node.label],
        )
        self.assertEqual(
            sorted(leaf.edge.length for leaf in parsed.leaf_node_iter()),
            [0.1, 0.2],
        )
        self.assertIn(":0.3", pruned)

    def test_duplicate_leaf_labels_are_rejected(self):
        tree = "((1:0.1,1:0.2)0.90:0.3,2:0.4);"

        with self.assertRaisesRegex(ValueError, "duplicate leaf labels: 1"):
            prune_tree_to_subset(tree, {"1", "2"}, min_taxa=1)

    def test_min_taxa_still_filters_incomplete_gene_trees(self):
        tree = "((1:0.1,2:0.2)0.95:0.3,3:0.4);"

        self.assertIsNone(prune_tree_to_subset(tree, {"1"}, min_taxa=2))


if __name__ == "__main__":
    unittest.main()

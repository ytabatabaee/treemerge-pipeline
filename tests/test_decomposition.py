import tempfile
import unittest
from pathlib import Path

import dendropy

import run_dtm_pipeline
from tests.reference_original_decomposer import read_and_decompose


TREE_11 = "((((a,b),(c,d)),((e,f),(g,h))),((i,j),k));\n"


class CentroidDecompositionTests(unittest.TestCase):
    def setUp(self):
        self.tree_file = tempfile.NamedTemporaryFile(
            mode="w", suffix=".nwk", delete=False
        )
        self.tree_file.write(TREE_11)
        self.tree_file.close()

    def tearDown(self):
        Path(self.tree_file.name).unlink(missing_ok=True)

    def test_matches_original_order_and_membership(self):
        expected = read_and_decompose(self.tree_file.name, 10)
        actual = run_dtm_pipeline.decompose_like_original_helper(
            self.tree_file.name,
            max_subset_size=10,
            min_subset_size=5,
        )
        self.assertEqual(actual, expected)

    def test_eleven_taxa_no_longer_hits_unsplittable_component(self):
        subsets = run_dtm_pipeline.decompose_like_original_helper(
            self.tree_file.name, max_subset_size=10
        )
        flattened = [taxon for subset in subsets for taxon in subset]

        self.assertEqual(sorted(map(len, subsets)), [4, 7])
        self.assertEqual(len(flattened), 11)
        self.assertEqual(len(set(flattened)), 11)
        self.assertLessEqual(max(map(len, subsets)), 10)

    def test_minimum_size_is_not_a_decomposition_constraint(self):
        tree = dendropy.Tree.get(data="(a,b,c,d);", schema="newick")
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".nwk", delete=False
        ) as handle:
            handle.write(tree.as_string(schema="newick"))
            path = handle.name
        try:
            subsets = run_dtm_pipeline.decompose_like_original_helper(
                path, max_subset_size=10, min_subset_size=5
            )
        finally:
            Path(path).unlink(missing_ok=True)
        self.assertEqual(len(subsets), 1)
        self.assertEqual(set(subsets[0]), {"a", "b", "c", "d"})


if __name__ == "__main__":
    unittest.main()

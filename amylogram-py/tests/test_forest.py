from __future__ import annotations

import unittest

from amylogram_py.forest import RangerForest, RangerTree, forest_from_export, tree_from_export


class RangerForestTests(unittest.TestCase):
    def test_single_tree_traversal_uses_ranger_numeric_split_semantics(self) -> None:
        tree = RangerTree(
            left_children=[1, 0, 0],
            right_children=[2, 0, 0],
            split_var_ids=[0, 0, 0],
            split_values=[0.0, 0.0, 0.0],
            terminal_class_counts=[
                [0, 0],
                [9, 1],
                [2, 8],
            ],
        )
        self.assertAlmostEqual(tree.predict_probability([0]), 0.1)
        self.assertAlmostEqual(tree.predict_probability([1]), 0.8)

    def test_forest_averages_tree_probabilities(self) -> None:
        low = RangerTree([0], [0], [0], [0.0], [[9, 1]])
        high = RangerTree([0], [0], [0], [0.0], [[1, 9]])
        forest = RangerForest([low, high])
        self.assertAlmostEqual(forest.predict_probability([0]), 0.5)

    def test_tree_from_export_accepts_ranger_child_pair_layout(self) -> None:
        tree = tree_from_export(
            child_node_ids=[[1, 0, 0], [2, 0, 0]],
            split_var_ids=[0, 0, 0],
            split_values=[0.0, 0.0, 0.0],
            terminal_class_counts=[[0, 0], [9, 1], [2, 8]],
        )
        self.assertAlmostEqual(tree.predict_probability([1]), 0.8)

    def test_forest_from_export_normalizes_one_based_r_split_ids(self) -> None:
        forest = forest_from_export(
            {
                "num_trees": 1,
                "child_node_ids": [[[1, 0, 0], [2, 0, 0]]],
                "split_var_ids": [[1, 0, 0]],
                "split_values": [[0.0, 0.0, 0.0]],
                "terminal_class_counts": [[[0, 0], [9, 1], [2, 8]]],
            },
            feature_count=1,
        )
        self.assertAlmostEqual(forest.predict_probability([1]), 0.8)


if __name__ == "__main__":
    unittest.main()

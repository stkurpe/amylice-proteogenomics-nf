"""Ranger-style forest traversal for exported AmyloGram models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence


@dataclass(frozen=True)
class RangerTree:
    """One binary tree exported from `ranger`.

    Node IDs are zero-based, matching the IDs stored by ranger's R object.
    Terminal nodes are represented by `left_child == right_child == 0`.
    """

    left_children: Sequence[int]
    right_children: Sequence[int]
    split_var_ids: Sequence[int]
    split_values: Sequence[float]
    terminal_class_counts: Sequence[Sequence[float]]

    def predict_probability(self, features: Sequence[float]) -> float:
        node = 0
        while True:
            left = int(self.left_children[node])
            right = int(self.right_children[node])
            if left == 0 and right == 0:
                counts = self.terminal_class_counts[node]
                if len(counts) < 2:
                    raise ValueError("terminal_class_counts must contain two classes")
                total = float(counts[0]) + float(counts[1])
                if total <= 0:
                    return 0.0
                return float(counts[1]) / total

            feature_index = int(self.split_var_ids[node])
            split_value = float(self.split_values[node])
            node = left if float(features[feature_index]) <= split_value else right


@dataclass(frozen=True)
class RangerForest:
    trees: Sequence[RangerTree]

    def predict_probability(self, features: Sequence[float]) -> float:
        if not self.trees:
            raise ValueError("forest has no trees")
        return sum(tree.predict_probability(features) for tree in self.trees) / len(self.trees)


def tree_from_export(
    child_node_ids: Sequence[Sequence[int]],
    split_var_ids: Sequence[int],
    split_values: Sequence[float],
    terminal_class_counts: Sequence[Sequence[float]],
    split_var_offset: int = 0,
) -> RangerTree:
    """Create a RangerTree from one R-exported tree object."""
    if len(child_node_ids) != 2:
        raise ValueError("child_node_ids must contain [left_children, right_children]")
    return RangerTree(
        left_children=child_node_ids[0],
        right_children=child_node_ids[1],
        split_var_ids=[int(value) + split_var_offset for value in split_var_ids],
        split_values=split_values,
        terminal_class_counts=terminal_class_counts,
    )


def forest_from_export(payload: dict[str, Any], feature_count: int) -> RangerForest:
    """Create a RangerForest from the JSON payload exported from R AmyloGram."""
    child_node_ids = payload["child_node_ids"]
    split_var_ids = payload["split_var_ids"]
    split_values = payload["split_values"]
    terminal_class_counts = payload["terminal_class_counts"]

    if not (
        len(child_node_ids)
        == len(split_var_ids)
        == len(split_values)
        == len(terminal_class_counts)
        == int(payload["num_trees"])
    ):
        raise ValueError("Inconsistent ranger forest export lengths")

    nonterminal_split_ids: list[int] = []
    for children, tree_split_ids in zip(child_node_ids, split_var_ids):
        for left, right, split_id in zip(children[0], children[1], tree_split_ids):
            if int(left) != 0 or int(right) != 0:
                nonterminal_split_ids.append(int(split_id))

    if not nonterminal_split_ids:
        raise ValueError("Ranger forest export contains no non-terminal splits")

    min_split_id = min(nonterminal_split_ids)
    max_split_id = max(nonterminal_split_ids)
    if min_split_id == 0 and max_split_id < feature_count:
        split_var_offset = 0
    elif min_split_id >= 1 and max_split_id <= feature_count:
        split_var_offset = -1
    else:
        raise ValueError(
            f"Unexpected split variable ID range {min_split_id}..{max_split_id} "
            f"for {feature_count} features"
        )

    return RangerForest(
        [
            tree_from_export(
                child_node_ids=children,
                split_var_ids=tree_split_ids,
                split_values=tree_split_values,
                terminal_class_counts=tree_terminal_counts,
                split_var_offset=split_var_offset,
            )
            for children, tree_split_ids, tree_split_values, tree_terminal_counts in zip(
                child_node_ids,
                split_var_ids,
                split_values,
                terminal_class_counts,
            )
        ]
    )

import pytest
import tskit
import numpy.testing as nt
import numpy as np

import model
import utils


def single_tree_example_ts():
    # 2.00┊    6    ┊
    #     ┊  ┏━┻━┓  ┊
    # 1.00┊  4   5  ┊
    #     ┊ ┏┻┓ ┏┻┓ ┊
    # 0.00┊ 0 1 2 3 ┊
    #     0         10
    ts = tskit.Tree.generate_balanced(4, span=10).tree_sequence
    tables = ts.dump_tables()
    for j in range(6):
        tables.sites.add_row(position=j + 1, ancestral_state="A")
        tables.mutations.add_row(site=j, derived_state="T", node=j)
    return tables.tree_sequence()


def single_tree_recurrent_mutation_example_ts():
    # 2.00 ┊                    6                    ┊
    #      ┊            ┏━━━━━━━┻━━━━━━━┓            ┊
    #      ┊      4:A→T x               x 5:A→T      ┊
    #      ┊            |               x 6:A→G      ┊
    # 1.00 ┊            4               5            ┊
    #      ┊       ┏━━━━┻━━━━┓     ┏━━━━┻━━━━┓       ┊
    #      ┊ 0:A→T x   1:A→T x     x 2:A→T   x 3:A→T ┊
    #      ┊       |         |     |         |       ┊
    # 0.00 ┊       0         1     2         3       ┊
    #      0                                        10
    ts = tskit.Tree.generate_balanced(4, span=10).tree_sequence
    tables = ts.dump_tables()
    for j in range(6):
        tables.sites.add_row(position=j + 1, ancestral_state="A")
        tables.mutations.add_row(site=j, derived_state="T", node=j)
    tables.mutations.add_row(site=j, derived_state="G", node=j, parent=j)
    ts = tables.tree_sequence()
    return tables.tree_sequence()

def multiple_trees_example_ts():
    # 2.00┊   4   ┊   4   ┊
    #     ┊ ┏━┻┓  ┊  ┏┻━┓ ┊
    # 1.00┊ ┃  3  ┊  3  ┃ ┊
    #     ┊ ┃ ┏┻┓ ┊ ┏┻┓ ┃ ┊
    # 0.00┊ 0 1 2 ┊ 0 1 2 ┊
    #     0       5      10
    ts = tskit.Tree.generate_balanced(3, span=10).tree_sequence
    tables = ts.dump_tables()
    tables.edges[1] = tables.edges[1].replace(right=5)
    tables.edges[2] = tables.edges[2].replace(right=5)
    tables.edges.add_row(5, 10, 3, 0)
    tables.edges.add_row(5, 10, 4, 2)
    tables.sort()
    return tables.tree_sequence()


def single_tree_with_polytomies_example_ts():
    # 3.00┊         8         ┊
    #     ┊  ┏━━━━━━╋━━━━━━━┓ ┊
    # 2.00┊  ┃      7       ┃ ┊
    #     ┊  ┃  ┏━━━╋━━━━┓  ┃ ┊
    # 1.00┊  5  ┃   6    ┃  ┃ ┊
    #     ┊ ┏┻┓ ┃ ┏━╋━━┓ ┃  ┃ ┊
    # 0.00┊ 0 1 2 3 4 11 9 10 ┊
    #     0                  10
    ts = tskit.Tree.generate_balanced(5, span=10).tree_sequence
    tables = ts.dump_tables()
    tables.nodes.add_row(flags=1, time=0)
    tables.edges.add_row(0, 10, 7, 9)
    tables.nodes.add_row(flags=1, time=0)
    tables.edges.add_row(0, 10, 8, 10)
    tables.nodes.add_row(flags=1, time=0)
    tables.edges.add_row(0, 10, 6, 11)
    tables.sort()
    return tables.tree_sequence()


def multi_tree_with_polytomies_example_ts():
    # 3.00┊     8       ┊     8       ┊
    #     ┊  ┏━━┻━┓     ┊  ┏━━┻━━┓    ┊
    # 2.00┊  ┃    7     ┊  ┃     7    ┊
    #     ┊  ┃  ┏━┻━┓   ┊  ┃  ┏━━╋━━┓ ┊
    # 1.00┊  5  ┃   6   ┊  5  ┃  6  ┃ ┊
    #     ┊ ┏┻┓ ┃ ┏━╋━┓ ┊ ┏┻┓ ┃ ┏┻┓ ┃ ┊
    # 0.00┊ 0 1 2 3 4 9 ┊ 0 1 2 3 4 9 ┊
    #     0             5            10
    ts = tskit.Tree.generate_balanced(5, span=10).tree_sequence
    tables = ts.dump_tables()
    tables.nodes.add_row(flags=1, time=0)
    tables.edges.add_row(0, 5, 6, 9)
    tables.edges.add_row(5, 10, 7, 9)
    tables.sort()
    return tables.tree_sequence()


class TestMutationDataTable:
    def test_single_tree_example(self):
        ts = single_tree_example_ts()
        tsm = model.TSModel(ts)
        df = tsm.mutations_df
        assert len(df) == 6
        nt.assert_array_equal(df.node, list(range(6)))
        nt.assert_array_equal(df.position, list(range(1, 7)))
        nt.assert_array_equal(df.time, [0, 0, 0, 0, 1, 1])
        nt.assert_array_equal(df.derived_state, ["T"] * 6)
        nt.assert_array_equal(df.inherited_state, ["A"] * 6)
        nt.assert_array_equal(df.num_parents, [0] * 6)
        nt.assert_array_equal(df.num_descendants, [1] * 4 + [2] * 2)
        nt.assert_array_equal(df.num_inheritors, [1] * 4 + [2] * 2)

    def test_single_tree_recurrent_mutation_example(self):
        ts = single_tree_recurrent_mutation_example_ts()
        tsm = model.TSModel(ts)
        df = tsm.mutations_df
        assert len(df) == 7
        nt.assert_array_equal(df.node, [0, 1, 2, 3, 4, 5, 5])
        nt.assert_array_equal(df.position, [1, 2, 3, 4, 5, 6, 6])
        nt.assert_array_equal(df.time, [0, 0, 0, 0, 1, 1, 1])
        nt.assert_array_equal(df.derived_state, ["T"] * 6 + ["G"])
        nt.assert_array_equal(df.inherited_state, ["A"] * 6 + ["T"])
        nt.assert_array_equal(df.num_parents, [0] * 6 + [1])
        nt.assert_array_equal(df.num_descendants, [1] * 4 + [2] * 3)
        nt.assert_array_equal(df.num_inheritors, [1] * 4 + [2, 0, 2])


class TestEdgeDataTable:
    def test_single_tree_example(self):
        ts = single_tree_example_ts()
        tsm = model.TSModel(ts)
        df = tsm.edges_df
        assert len(df) == 6
        nt.assert_array_equal(df.left, [0, 0, 0, 0, 0, 0])
        nt.assert_array_equal(df.right, [10, 10, 10, 10, 10, 10])
        nt.assert_array_equal(df.parent, [4, 4, 5, 5, 6, 6])
        nt.assert_array_equal(df.child, [0, 1, 2, 3, 4, 5])
        nt.assert_array_equal(df.child_time, [0, 0, 0, 0, 1, 1])
        nt.assert_array_equal(df.parent_time, [1, 1, 1, 1, 2, 2])

    def test_multiple_trees_example(self):
        ts = multiple_trees_example_ts()
        tsm = model.TSModel(ts)
        df = tsm.edges_df
        assert len(df) == 6
        nt.assert_array_equal(df.left, [5, 0, 0, 0, 5, 0])
        nt.assert_array_equal(df.right, [10, 10, 5, 5, 10, 10])
        nt.assert_array_equal(df.parent, [3, 3, 3, 4, 4, 4])
        nt.assert_array_equal(df.child, [0, 1, 2, 0, 2, 3])
        nt.assert_array_equal(df.child_time, [0, 0, 0, 0, 0, 1])
        nt.assert_array_equal(df.parent_time, [1, 1, 1, 2, 2, 2])


class TestNodeDataTable:
    def test_single_tree_example(self):
        ts = single_tree_example_ts()
        tsm = model.TSModel(ts)
        df = tsm.nodes_df
        assert len(df) == 7
        nt.assert_array_equal(df.time, [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0])
        nt.assert_array_equal(df.num_mutations, [1, 1, 1, 1, 1, 1, 0])
        nt.assert_array_equal(df.ancestors_span, [
                              10, 10, 10, 10, 10, 10, -np.inf])

    def test_multiple_tree_example(self):
        ts = multiple_trees_example_ts()
        tsm = model.TSModel(ts)
        df = tsm.nodes_df
        assert len(df) == 5
        nt.assert_array_equal(df.time, [0.0, 0.0, 0.0, 1.0, 2.0])
        nt.assert_array_equal(df.num_mutations, [0, 0, 0, 0, 0])
        nt.assert_array_equal(df.ancestors_span, [10, 10, 10, 10, -np.inf])


class TestTreesDataTable:
    def test_single_tree_example(self):
        ts = single_tree_example_ts()
        tsm = model.TSModel(ts)
        df = tsm.trees_df
        assert len(df) == 1
        nt.assert_array_equal(df.left, 0)
        nt.assert_array_equal(df.right, 10)
        nt.assert_array_equal(df.total_branch_length, 6.0)
        # nt.assert_array_equal(df.mean_internal_arity, 2.0)
        nt.assert_array_equal(df.max_internal_arity, 2.0)

    def test_single_tree_with_polytomies_example(self):
        ts = single_tree_with_polytomies_example_ts()
        tsm = model.TSModel(ts)
        df = tsm.trees_df
        assert len(df) == 1
        nt.assert_array_equal(df.left, 0)
        nt.assert_array_equal(df.right, 10)
        nt.assert_array_equal(df.total_branch_length, 16.0)
        # nt.assert_array_equal(df.mean_internal_arity, 2.75)
        nt.assert_array_equal(df.max_internal_arity, 3.0)

    def test_multi_tree_with_polytomies_example(self):
        ts = multi_tree_with_polytomies_example_ts()
        tsm = model.TSModel(ts)
        df = tsm.trees_df
        assert len(df) == 2
        nt.assert_array_equal(df.left, [0, 5])
        nt.assert_array_equal(df.right, [5, 10])
        nt.assert_array_equal(df.total_branch_length, [11.0, 12.0])
        # nt.assert_array_equal(df.mean_internal_arity, [2.25, 2.25])
        nt.assert_array_equal(df.max_internal_arity, [3.0, 3.0])

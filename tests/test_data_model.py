import pytest
import tskit
import numpy.testing as nt

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


class TestMutationDataTable:
    def test_single_tree_example(self):
        ts = single_tree_example_ts()
        ti = utils.TreeInfo(ts, 0)
        df = ti.mutations_data()
        assert len(df) == 6
        nt.assert_array_equal(df.node, list(range(6)))
        nt.assert_array_equal(df.position, list(range(1, 7)))
        nt.assert_array_equal(df.time, [0, 0, 0, 0, 1, 1])


class TestEdgeDataTable:
    def test_single_tree_example(self):
        ts = single_tree_example_ts()
        ti = utils.TreeInfo(ts, 0)
        df = ti.edges_data()
        assert len(df) == 6
        nt.assert_array_equal(df.left, [0, 0, 0, 0, 0, 0])
        nt.assert_array_equal(df.right, [10, 10, 10, 10, 10, 10])
        nt.assert_array_equal(df.parent, [4, 4, 5, 5, 6, 6])
        nt.assert_array_equal(df.child, [0, 1, 2, 3, 4, 5])
        nt.assert_array_equal(df.child_time, [0, 0, 0, 0, 1, 1])
        nt.assert_array_equal(df.parent_time, [1, 1, 1, 1, 2, 2])

    def test_multiple_trees_example(self):
        ts = multiple_trees_example_ts()
        ti = utils.TreeInfo(ts, 0)
        df = ti.edges_data()
        assert len(df) == 6
        nt.assert_array_equal(df.left, [5, 0, 0, 0, 5, 0])
        nt.assert_array_equal(df.right, [10, 10, 5, 5, 10, 10])
        nt.assert_array_equal(df.parent, [3, 3, 3, 4, 4, 4])
        nt.assert_array_equal(df.child, [0, 1, 2, 0, 2, 3])
        nt.assert_array_equal(df.child_time, [0, 0, 0, 0, 0, 1])
        nt.assert_array_equal(df.parent_time, [1, 1, 1, 2, 2, 2])


class TestNodesDataTable:
    def test_single_tree_example(self):
        ts = single_tree_example_ts()
        ti = utils.TreeInfo(ts, 0)
        df = ti.nodes_data()
        assert len(df) == 7
        nt.assert_array_equal(df.time, [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0])
        nt.assert_array_equal(df.num_mutations, [1, 1, 1, 1, 1, 1, 0])
        nt.assert_array_equal(
            df.ancestors_span, [10, 10, 10, 10, 10, 10, -1_000_000_000]
        )

    def test_single_tree_example(self):
        ts = multiple_trees_example_ts()
        ti = utils.TreeInfo(ts, 0)
        df = ti.nodes_data()
        assert len(df) == 5
        nt.assert_array_equal(df.time, [0.0, 0.0, 0.0, 1.0, 2.0])
        nt.assert_array_equal(df.num_mutations, [0, 0, 0, 0, 0])
        nt.assert_array_equal(df.ancestors_span, [
                              10, 10, 10, 10, -1_000_000_000])

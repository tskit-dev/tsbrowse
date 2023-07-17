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


class TestMutationDataTable:
    def test_single_tree_example(self):
        ts = single_tree_example_ts()
        ti = utils.TreeInfo(ts, 0)
        df = ti.mutations_data()
        assert len(df) == 6
        nt.assert_array_equal(df.node, list(range(6)))
        nt.assert_array_equal(df.position, list(range(1, 7)))
        nt.assert_array_equal(df.time, [0, 0, 0, 0, 1, 1])

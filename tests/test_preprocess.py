import os

import msprime
import numpy as np
import numpy.testing as nt
import pytest
import tskit
import tszip
import zarr

from tsbrowse import preprocess
from tsbrowse import TSBROWSE_DATA_VERSION


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
    tables.sites.add_row(position=7, ancestral_state="FOOBAR")
    tables.mutations.add_row(site=6, derived_state="FOOBARD", node=6)
    tables.compute_mutation_times()
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
        m = preprocess.mutations(ts)
        for array in m.values():
            assert len(array) == 7
        nt.assert_array_equal(m["inherited_state"], ["A"] * 6 + ["FOOBAR"])
        nt.assert_array_equal(m["num_parents"], [0] * 7)
        nt.assert_array_equal(m["num_descendants"], [1] * 4 + [2] * 2 + [4])
        nt.assert_array_equal(m["num_inheritors"], [1] * 4 + [2] * 2 + [4])

    def test_single_tree_recurrent_mutation_example(self):
        ts = single_tree_recurrent_mutation_example_ts()
        m = preprocess.mutations(ts)
        for array in m.values():
            assert len(array) == 7
        nt.assert_array_equal(m["inherited_state"], ["A"] * 6 + ["T"])
        nt.assert_array_equal(m["num_parents"], [0] * 6 + [1])
        nt.assert_array_equal(m["num_descendants"], [1] * 4 + [2] * 3)
        nt.assert_array_equal(m["num_inheritors"], [1] * 4 + [2, 0, 2])


class TestNodeDataTable:
    def test_single_tree_example(self):
        ts = single_tree_example_ts()
        n = preprocess.nodes(ts)
        for array in n.values():
            assert len(array) == 7
        nt.assert_array_equal(n["num_mutations"], [1, 1, 1, 1, 1, 1, 1])
        nt.assert_array_equal(n["ancestors_span"], [10, 10, 10, 10, 10, 10, -np.inf])

    def test_multiple_tree_example(self):
        ts = multiple_trees_example_ts()
        n = preprocess.nodes(ts)
        for array in n.values():
            assert len(array) == 5
        nt.assert_array_equal(n["num_mutations"], [0, 0, 0, 0, 0])
        nt.assert_array_equal(n["ancestors_span"], [10, 10, 10, 10, -np.inf])


class TestEdgeDataTable:
    def test_single_tree_example(self):
        ts = single_tree_example_ts()
        e = preprocess.edges(ts)
        for array in e.values():
            assert len(array) == 6
        nt.assert_array_equal(e["child_time"], [0, 0, 0, 0, 1, 1])
        nt.assert_array_equal(e["parent_time"], [1, 1, 1, 1, 2, 2])
        nt.assert_array_equal(e["branch_length"], [1, 1, 1, 1, 1, 1])
        nt.assert_array_equal(e["span"], [10, 10, 10, 10, 10, 10])

    def test_multiple_trees_example(self):
        ts = multiple_trees_example_ts()
        e = preprocess.edges(ts)
        for array in e.values():
            assert len(array) == 6
        nt.assert_array_equal(e["child_time"], [0, 0, 0, 0, 0, 1])
        nt.assert_array_equal(e["parent_time"], [1, 1, 1, 2, 2, 2])
        nt.assert_array_equal(e["branch_length"], [1, 1, 1, 2, 2, 1])
        nt.assert_array_equal(e["span"], [5, 10, 5, 5, 5, 10])


class TestSiteDataTable:
    def test_single_tree_example(self):
        ts = single_tree_example_ts()
        s = preprocess.sites(ts)
        for array in s.values():
            assert len(array) == 7
        nt.assert_array_equal(s["num_mutations"], [1, 1, 1, 1, 1, 1, 1])

    def test_single_tree_recurrent_mutation_example(self):
        ts = single_tree_recurrent_mutation_example_ts()
        s = preprocess.sites(ts)
        for array in s.values():
            assert len(array) == 6
        nt.assert_array_equal(s["num_mutations"], [1, 1, 1, 1, 1, 2])


class TestTreesDataTable:
    def test_single_tree_example(self):
        ts = single_tree_example_ts()
        t = preprocess.trees(ts)
        for array in t.values():
            assert len(array) == 1
        nt.assert_array_equal(t["left"], 0)
        nt.assert_array_equal(t["right"], 10)
        nt.assert_array_equal(t["total_branch_length"], 6.0)
        # nt.assert_array_equal(t['mean_internal_arity'], 2.0)
        nt.assert_array_equal(t["max_internal_arity"], 2.0)

    def test_single_tree_with_polytomies_example(self):
        ts = single_tree_with_polytomies_example_ts()
        t = preprocess.trees(ts)
        for array in t.values():
            assert len(array) == 1
        nt.assert_array_equal(t["left"], 0)
        nt.assert_array_equal(t["right"], 10)
        nt.assert_array_equal(t["total_branch_length"], 16.0)
        # nt.assert_array_equal(t['mean_internal_arity'], 2.75)
        nt.assert_array_equal(t["max_internal_arity"], 3.0)

    def test_multi_tree_with_polytomies_example(self):
        ts = multi_tree_with_polytomies_example_ts()
        t = preprocess.trees(ts)
        for array in t.values():
            assert len(array) == 2
        nt.assert_array_equal(t["left"], [0, 5])
        nt.assert_array_equal(t["right"], [5, 10])
        nt.assert_array_equal(t["total_branch_length"], [11.0, 12.0])
        # nt.assert_array_equal(t['mean_internal_arity'], [2.25, 2.25])
        nt.assert_array_equal(t["max_internal_arity"], [3.0, 3.0])


class TestMutationFrequencies:
    def example_ts(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=10_000)
        demography.add_population(name="B", initial_size=5_000)
        demography.add_population(name="C", initial_size=1_000)
        demography.add_population_split(time=1000, derived=["A", "B"], ancestral="C")
        return msprime.sim_ancestry(
            samples={"A": 1, "B": 1},
            demography=demography,
            random_seed=12,
            sequence_length=10_000,
        )

    def compute_mutation_counts(self, ts):
        pop_mutation_count = np.zeros((ts.num_populations, ts.num_mutations), dtype=int)
        for pop in ts.populations():
            for tree in ts.trees(tracked_samples=ts.samples(population=pop.id)):
                for mut in tree.mutations():
                    count = tree.num_tracked_samples(mut.node)
                    pop_mutation_count[pop.id, mut.id] = count
        return pop_mutation_count

    def check_ts(self, ts):
        C1 = self.compute_mutation_counts(ts)
        C2 = preprocess.compute_population_mutation_counts(ts)
        nt.assert_array_equal(C1, C2)
        m = preprocess.mutations(ts)
        nt.assert_array_equal(m["pop_A_freq"], C1[0] / ts.num_samples)
        nt.assert_array_equal(m["pop_B_freq"], C1[1] / ts.num_samples)
        nt.assert_array_equal(m["pop_C_freq"], C1[2] / ts.num_samples)

    def test_all_nodes(self):
        ts = self.example_ts()
        tables = ts.dump_tables()
        for u in range(ts.num_nodes - 1):
            site_id = tables.sites.add_row(u, "A")
            tables.mutations.add_row(site=site_id, node=u, derived_state="T")
        ts = tables.tree_sequence()
        self.check_ts(ts)

    @pytest.mark.parametrize("seed", range(1, 7))
    def test_simulated_mutations(self, seed):
        ts = msprime.sim_mutations(self.example_ts(), rate=1e-6, random_seed=seed)
        assert ts.num_mutations > 0
        self.check_ts(ts)

    def test_no_metadata_schema(self):
        ts = msprime.sim_mutations(self.example_ts(), rate=1e-6, random_seed=43)
        assert ts.num_mutations > 0
        tables = ts.dump_tables()
        tables.populations.metadata_schema = tskit.MetadataSchema(None)
        self.check_ts(tables.tree_sequence())

    def test_no_populations(self):
        tables = single_tree_example_ts().dump_tables()
        tables.populations.add_row(b"{}")
        with pytest.raises(ValueError, match="must be assigned to populations"):
            preprocess.mutations(tables.tree_sequence())


class TestNodeIsSample:
    def test_simple_example(self):
        ts = single_tree_example_ts()
        is_sample = preprocess.node_is_sample(ts)
        for node in ts.nodes():
            assert node.is_sample() == is_sample[node.id]

    @pytest.mark.parametrize("bit", [1, 2, 17, 31])
    def test_sample_and_other_flags(self, bit):
        tables = single_tree_example_ts().dump_tables()
        flags = tables.nodes.flags
        tables.nodes.flags = flags | (1 << bit)
        ts = tables.tree_sequence()
        is_sample = preprocess.node_is_sample(ts)
        for node in ts.nodes():
            assert node.is_sample() == is_sample[node.id]
            assert (node.flags & (1 << bit)) != 0


def test_preprocess_calculate_mutation_times(tmpdir):
    ts = msprime.sim_ancestry(5, sequence_length=1e4, random_seed=42)
    ts = msprime.sim_mutations(ts, rate=0.1, random_seed=43)
    assert ts.num_mutations > 0
    # Wipe out mutation times
    tables = ts.dump_tables()
    tables.mutations.time = np.full_like(tables.mutations.time, tskit.UNKNOWN_TIME)
    ts = tables.tree_sequence()
    input_path = os.path.join(tmpdir, "test.trees")
    ts.dump(input_path)
    output_path = os.path.join(tmpdir, "test.trees")
    with pytest.warns(UserWarning, match="All mutation times are unknown"):
        preprocess.preprocess(input_path, output_path)
    ts = tszip.load(output_path)
    assert not np.any(tskit.is_unknown_time(ts.tables.mutations.time))


@pytest.mark.parametrize("use_tszip", [True, False])
def test_preprocess(tmpdir, use_tszip):
    input_path = os.path.join(tmpdir, "test_input.tszip")
    output_path = os.path.join(tmpdir, "test_output.tsbrowse")

    ts = single_tree_example_ts()
    if use_tszip:
        tszip.compress(ts, input_path)
    else:
        ts.dump(input_path)

    preprocess.preprocess(input_path, output_path)

    assert os.path.exists(output_path)
    # Check that the output is still a valid tszip
    tszip.load(output_path).tables.assert_equals(ts.tables)

    # Check that the file contains the expected arrays
    with zarr.ZipStore(output_path, mode="r") as zarr_store:
        root = zarr.group(store=zarr_store)
        assert root.attrs["tsbrowse"]["data_version"] == TSBROWSE_DATA_VERSION
        for array_name in [
            "mutations/position",
            "mutations/inherited_state",
            "mutations/num_descendants",
            "mutations/num_inheritors",
            "mutations/num_parents",
            "nodes/num_mutations",
            "nodes/ancestors_span",
            "trees/left",
            "trees/right",
            "trees/total_branch_length",
            "trees/mean_internal_arity",
            "trees/max_internal_arity",
            "trees/num_sites",
            "trees/num_mutations",
            "edges/parent_time",
            "edges/child_time",
            "edges/branch_length",
            "edges/span",
            "sites/num_mutations",
        ]:
            assert array_name in root

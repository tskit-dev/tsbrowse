import msprime
import numpy as np
import pytest
import tskit
import tszip
import zarr

from tsbrowse import model
from tsbrowse import preprocess


def test_model(tmpdir):
    # Generate a tree sequence with populations and migrations
    N = 1000
    demography = msprime.Demography()
    demography.add_population(name="pop1", initial_size=N)
    demography.add_population(name="pop2", initial_size=N)
    demography.add_population(name="ancestral", initial_size=N)
    demography.set_symmetric_migration_rate(["pop1", "pop2"], 0.01)
    demography.add_population_split(
        time=1000, derived=["pop1", "pop2"], ancestral="ancestral"
    )
    ts = msprime.sim_ancestry(
        samples={"pop1": 5, "pop2": 5},
        demography=demography,
        sequence_length=1e4,
        record_migrations=True,
        random_seed=42,
    )
    ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=42)
    assert ts.num_populations > 0
    assert ts.num_sites > 0
    assert ts.num_migrations > 0
    assert ts.num_mutations > 0

    tables = ts.tables
    tables.nodes.metadata_schema = tskit.MetadataSchema({"codec": "json"})

    # Give each individual a location
    indiv_copy = tables.individuals.copy()
    tables.individuals.clear()
    for i, ind in enumerate(indiv_copy):
        tables.individuals.append(ind.replace(location=[i / 2, i + 1]))

    ts = tables.tree_sequence()

    tszip.compress(ts, tmpdir / "test.tszip")
    preprocess.preprocess(tmpdir / "test.tszip", tmpdir / "test.tsbrowse")
    tsm = model.TSModel(tmpdir / "test.tsbrowse")

    assert tsm.ts == ts
    assert tsm.name == "test"
    assert tsm.file_uuid == ts.file_uuid
    assert len(tsm.summary_df) == 9
    assert len(tsm.trees_df) == ts.num_trees

    assert len(tsm.edges_df) == ts.num_edges
    for col in ["left", "right", "parent", "child"]:
        assert np.array_equal(tsm.edges_df[col].values, getattr(ts.tables.edges, col))

    assert len(tsm.mutations_df) == ts.num_mutations
    for m1, m2 in zip(ts.mutations(), tsm.mutations_df.to_dict("records")):
        assert m1.derived_state == m2["derived_state"]
        assert m1.site == m2["site"]
        assert m1.node == m2["node"]
        assert m1.parent == m2["parent"]
        assert m1.time == m2["time"]

    assert len(tsm.nodes_df) == ts.num_nodes
    for col in ["time", "flags", "population", "individual"]:
        assert np.array_equal(tsm.nodes_df[col].values, getattr(ts.tables.nodes, col))

    assert len(tsm.sites_df) == ts.num_sites
    for m1, m2 in zip(ts.sites(), tsm.sites_df.to_dict("records")):
        assert m1.ancestral_state == m2["ancestral_state"]
        assert m1.position == m2["position"]

    assert len(tsm.individuals_df) == ts.num_individuals
    for m1, m2 in zip(ts.individuals(), tsm.individuals_df.to_dict("records")):
        assert m1.flags == m2["flags"]
        assert np.array_equal(m1.location, m2["location"])
        assert np.array_equal(m1.parents, m2["parents"])

    assert len(tsm.populations_df) == ts.num_populations
    for m1, m2 in zip(ts.populations(), tsm.populations_df.to_dict("records")):
        assert str(m1.metadata) == m2["metadata"]

    assert len(tsm.migrations_df) == ts.num_migrations
    for col in ["left", "right", "node", "source", "dest", "time"]:
        assert np.array_equal(
            tsm.migrations_df[col].values, getattr(ts.tables.migrations, col)
        )


def test_model_errors(tmpdir):
    # Write an empty zarr ZipStore
    with zarr.ZipStore(tmpdir / "test.tsbrowse", mode="w") as z:
        g = zarr.group(store=z)
        g.attrs["foo"] = "bar"
    with pytest.raises(ValueError, match="File is not a tsbrowse file"):
        model.TSModel(tmpdir / "test.tsbrowse")

    ts = msprime.sim_ancestry(
        recombination_rate=1e-3, samples=2, sequence_length=1000, random_seed=42
    )
    tszip.compress(ts, tmpdir / "test.tszip")
    preprocess.preprocess(tmpdir / "test.tszip", tmpdir / "test.tsbrowse")
    with zarr.ZipStore(tmpdir / "test.tsbrowse", mode="w") as z:
        g = zarr.group(store=z)
        g.attrs["tsbrowse"] = {"data_version": 0}
    with pytest.raises(ValueError, match="File .* has version .*"):
        model.TSModel(tmpdir / "test.tsbrowse")

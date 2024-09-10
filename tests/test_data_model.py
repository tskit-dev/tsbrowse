import msprime
import pytest
import tskit
import tszip
import zarr

from tsbrowse import model
from tsbrowse import preprocess


def test_model(tmpdir):
    ts = msprime.sim_ancestry(
        recombination_rate=1e-3, samples=10, sequence_length=1_000, random_seed=42
    )
    ts = msprime.sim_mutations(ts, rate=1e-2, random_seed=43)
    tables = ts.tables
    tables.nodes.metadata_schema = tskit.MetadataSchema({"codec": "json"})
    ts = tables.tree_sequence()

    tszip.compress(ts, tmpdir / "test.tszip")
    preprocess.preprocess(tmpdir / "test.tszip", tmpdir / "test.tsbrowse")
    tsm = model.TSModel(tmpdir / "test.tsbrowse")

    assert tsm.ts == ts
    assert tsm.name == "test"
    assert tsm.file_uuid == ts.file_uuid
    assert len(tsm.summary_df) == 9
    assert len(tsm.edges_df) == ts.num_edges
    assert len(tsm.trees_df) == ts.num_trees
    assert len(tsm.mutations_df) == ts.num_mutations
    assert len(tsm.nodes_df) == ts.num_nodes
    assert len(tsm.sites_df) == ts.num_sites


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

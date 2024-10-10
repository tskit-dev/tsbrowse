import pathlib
from functools import cached_property

import daiquiri
import numpy as np
import pandas as pd
import tszip
import zarr

from . import TSBROWSE_DATA_VERSION

logger = daiquiri.getLogger("tsbrowse")


class TSModel:
    """
    A wrapper around a tskit.TreeSequence object that provides some
    convenience methods for analysing the tree sequence.
    """

    def __init__(self, tsbrowse_path):
        tsbrowse_path = pathlib.Path(tsbrowse_path)
        root = zarr.open(zarr.ZipStore(tsbrowse_path, mode="r"))
        if "tsbrowse" not in root.attrs or "data_version" not in root.attrs["tsbrowse"]:
            raise ValueError("File is not a tsbrowse file, run tsbrowse preprocess")
        if root.attrs["tsbrowse"]["data_version"] != TSBROWSE_DATA_VERSION:
            raise ValueError(
                f"File {tsbrowse_path} has version "
                f"{root.attrs['tsbrowse']['data_version']}, "
                f"but this version of tsbrowse expects version "
                f"{TSBROWSE_DATA_VERSION} rerun tsbrowse preprocess"
            )
        self.ts = tszip.load(tsbrowse_path)
        self.name = tsbrowse_path.stem
        for table_name in [
            "edges",
            "trees",
            "mutations",
            "nodes",
            "sites",
            "individuals",
            "populations",
            "migrations",
            "provenances",
        ]:
            # filter out ragged arrays with offset
            array_names = set(root[table_name].keys())
            ragged_array_names = {
                "_".join(name.split("_")[:-1])
                for name in array_names
                if "offset" in name
            }
            array_names -= set(ragged_array_names)
            array_names -= {"metadata_schema"}
            array_names -= {f"{name}_offset" for name in ragged_array_names}
            arrays = {name: root[table_name][name][:] for name in array_names}
            ragged_array_names -= {"metadata"}
            for name in ragged_array_names:
                array = root[table_name][name][:]
                offsets = root[table_name][f"{name}_offset"][:]
                arrays[name] = np.array(
                    [
                        array[s].tobytes().decode("utf-8")
                        for s in (
                            slice(start, end)
                            for start, end in zip(offsets[:-1], offsets[1:])
                        )
                    ]
                )
            df = pd.DataFrame(arrays)
            df["id"] = df.index
            setattr(self, f"{table_name}_df", df)

    @property
    def file_uuid(self):
        return self.ts.file_uuid

    @cached_property
    def summary_df(self):
        data = [
            ("samples", self.ts.num_samples),
            ("nodes", self.ts.num_nodes),
            ("mutations", self.ts.num_mutations),
            ("nodes_with_zero_muts", np.sum(self.nodes_df.num_mutations == 0)),
            ("sites_with_zero_muts", np.sum(self.sites_df.num_mutations == 0)),
            ("max_mutations_per_site", np.max(self.sites_df.num_mutations)),
            ("mean_mutations_per_site", np.mean(self.sites_df.num_mutations)),
            ("median_mutations_per_site", np.median(self.sites_df.num_mutations)),
            ("max_mutations_per_node", np.max(self.nodes_df.num_mutations)),
        ]
        df = pd.DataFrame(
            {"property": [d[0] for d in data], "value": [d[1] for d in data]}
        )
        logger.info("Computed summary dataframe")
        return df.set_index("property")

    def _repr_html_(self):
        return self.summary_df._repr_html_()

    def genes_df(self, genes_file):
        genes_df = pd.read_csv(genes_file, sep=";")
        # TODO file checks!
        genes_df.columns = ["chr", "position", "end", "strand", "id", "name"]
        genes_df = genes_df[
            (genes_df["position"] >= self.ts.first().interval.left)
            & (genes_df["end"] <= self.ts.last().interval.right)
        ]
        logger.info("Computed genes dataframe")
        return genes_df

    def calc_polytomy_fractions(self):
        """
        Calculates the fraction of polytomies for each tree in the
        tree sequence
        """
        assert self.ts.num_samples > 2
        polytomy_fractions = []
        for tree in self.ts.trees():
            if tree.num_edges == 0:
                polytomy_fractions.append(None)
            else:
                polytomy_fractions.append(
                    float(
                        (tree.num_edges - self.ts.num_samples)
                        / (self.ts.num_samples - 2)
                    )
                )
        return polytomy_fractions

    def map_stats_to_genome(self, to_map):
        """
        Converts a list of tree-based stats to genomic coordinates
        """
        mapped = np.zeros(int(self.ts.sequence_length))
        for i, tree in enumerate(self.ts.trees()):
            left, right = map(int, tree.interval)
            mapped[left:right] = to_map[i]
        return mapped

    def make_sliding_windows(self, iterable, size, overlap=0):
        start = 0
        assert overlap < size, "overlap must be smaller then window size"
        end = size
        step = size - overlap

        length = len(iterable)
        while end < length:
            yield iterable[start:end]
            start += step
            end += step
        yield iterable[start:]

    def calc_mean_node_arity(self):
        span_sums = np.bincount(
            self.ts.edges_parent,
            weights=self.ts.edges_right - self.ts.edges_left,
            minlength=self.ts.num_nodes,
        )
        node_spans = self.ts.sample_count_stat(
            [self.ts.samples()],
            lambda x: (x > 0),
            1,
            polarised=True,
            span_normalise=False,
            strict=False,
            mode="node",
        )[:, 0]
        return span_sums / node_spans

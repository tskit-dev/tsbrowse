import dataclasses
import json
import pathlib
import warnings

import daiquiri
import numba
import numpy as np
import tskit
import tszip
import zarr
from tqdm import tqdm

from . import jit
from tsbrowse import TSBROWSE_DATA_VERSION

logger = daiquiri.getLogger("tsbrowse")


def node_is_sample(ts):
    sample_flag = np.full_like(ts.nodes_flags, tskit.NODE_IS_SAMPLE)
    return np.bitwise_and(ts.nodes_flags, sample_flag) != 0


spec = [
    ("num_edges", numba.int64),
    ("sequence_length", numba.float64),
    ("edges_left", numba.float64[:]),
    ("edges_right", numba.float64[:]),
    ("edge_insertion_order", numba.int32[:]),
    ("edge_removal_order", numba.int32[:]),
    ("edge_insertion_index", numba.int64),
    ("edge_removal_index", numba.int64),
    ("interval", numba.float64[:]),
    ("in_range", numba.int64[:]),
    ("out_range", numba.int64[:]),
]


@jit.numba_jitclass(spec)
class TreePosition:
    def __init__(
        self,
        num_edges,
        sequence_length,
        edges_left,
        edges_right,
        edge_insertion_order,
        edge_removal_order,
    ):
        self.num_edges = num_edges
        self.sequence_length = sequence_length
        self.edges_left = edges_left
        self.edges_right = edges_right
        self.edge_insertion_order = edge_insertion_order
        self.edge_removal_order = edge_removal_order
        self.edge_insertion_index = 0
        self.edge_removal_index = 0
        self.interval = np.zeros(2)
        self.in_range = np.zeros(2, dtype=np.int64)
        self.out_range = np.zeros(2, dtype=np.int64)

    def next(self):  # noqa
        left = self.interval[1]
        j = self.in_range[1]
        k = self.out_range[1]
        self.in_range[0] = j
        self.out_range[0] = k
        M = self.num_edges
        edges_left = self.edges_left
        edges_right = self.edges_right
        out_order = self.edge_removal_order
        in_order = self.edge_insertion_order

        while k < M and edges_right[out_order[k]] == left:
            k += 1
        while j < M and edges_left[in_order[j]] == left:
            j += 1
        self.out_range[1] = k
        self.in_range[1] = j

        right = self.sequence_length
        if j < M:
            right = min(right, edges_left[in_order[j]])
        if k < M:
            right = min(right, edges_right[out_order[k]])
        self.interval[:] = [left, right]
        return j < M or left < self.sequence_length


# Helper function to make it easier to communicate with the numba class
def alloc_tree_position(ts):
    return TreePosition(
        num_edges=ts.num_edges,
        sequence_length=ts.sequence_length,
        edges_left=ts.edges_left,
        edges_right=ts.edges_right,
        edge_insertion_order=ts.indexes_edge_insertion_order,
        edge_removal_order=ts.indexes_edge_removal_order,
    )


@jit.numba_jit()
def _compute_population_mutation_counts(
    tree_pos,
    num_nodes,
    num_mutations,
    num_populations,
    edges_parent,
    edges_child,
    nodes_is_sample,
    nodes_population,
    mutations_position,
    mutations_node,
    mutations_parent,
):
    num_pop_samples = np.zeros((num_nodes, num_populations), dtype=np.int32)

    pop_mutation_count = np.zeros((num_populations, num_mutations), dtype=np.int32)
    parent = np.zeros(num_nodes, dtype=np.int32) - 1

    for u in range(num_nodes):
        if nodes_is_sample[u]:
            num_pop_samples[u, nodes_population[u]] = 1

    mut_id = 0
    while tree_pos.next():
        for j in range(tree_pos.out_range[0], tree_pos.out_range[1]):
            e = tree_pos.edge_removal_order[j]
            c = edges_child[e]
            p = edges_parent[e]
            parent[c] = -1
            u = p
            while u != -1:
                for k in range(num_populations):
                    num_pop_samples[u, k] -= num_pop_samples[c, k]
                u = parent[u]

        for j in range(tree_pos.in_range[0], tree_pos.in_range[1]):
            e = tree_pos.edge_insertion_order[j]
            p = edges_parent[e]
            c = edges_child[e]
            parent[c] = p
            u = p
            while u != -1:
                for k in range(num_populations):
                    num_pop_samples[u, k] += num_pop_samples[c, k]
                u = parent[u]

        left, right = tree_pos.interval
        while mut_id < num_mutations and mutations_position[mut_id] < right:
            assert mutations_position[mut_id] >= left
            mutation_node = mutations_node[mut_id]
            for pop in range(num_populations):
                pop_mutation_count[pop, mut_id] = num_pop_samples[mutation_node, pop]
            mut_id += 1

    return pop_mutation_count


def compute_population_mutation_counts(ts):
    """
    Return a (num_populations, num_mutations) array that gives the frequency
    of each mutation in each of the populations in the specified tree sequence.
    """
    logger.info(
        f"Computing mutation frequencies within {ts.num_populations} populations"
    )
    mutations_position = ts.sites_position[ts.mutations_site].astype(int)

    if np.any(ts.nodes_population[ts.samples()] == -1):
        raise ValueError("Sample nodes must be assigned to populations")

    return _compute_population_mutation_counts(
        alloc_tree_position(ts),
        ts.num_nodes,
        ts.num_mutations,
        ts.num_populations,
        ts.edges_parent,
        ts.edges_child,
        node_is_sample(ts),
        ts.nodes_population,
        mutations_position,
        ts.mutations_node,
        ts.mutations_parent,
    )


@dataclasses.dataclass
class MutationCounts:
    num_parents: np.ndarray
    num_inheritors: np.ndarray
    num_descendants: np.ndarray


def compute_mutation_counts(ts):
    logger.info("Computing mutation inheritance counts")
    tree_pos = alloc_tree_position(ts)
    mutations_position = ts.sites_position[ts.mutations_site].astype(int)
    num_descendants, num_inheritors = _compute_mutation_inheritance_counts(
        tree_pos,
        ts.num_nodes,
        ts.num_mutations,
        ts.edges_parent,
        ts.edges_child,
        ts.samples(),
        mutations_position,
        ts.mutations_node,
        ts.mutations_parent,
    )
    num_parents = _compute_mutation_parent_counts(ts.mutations_parent)
    return MutationCounts(num_parents, num_inheritors, num_descendants)


@jit.numba_jit()
def _compute_mutation_parent_counts(mutations_parent):
    N = mutations_parent.shape[0]
    num_parents = np.zeros(N, dtype=np.int32)

    for j in range(N):
        u = j
        while mutations_parent[u] != -1:
            num_parents[j] += 1
            u = mutations_parent[u]
    return num_parents


@jit.numba_jit()
def _compute_mutation_inheritance_counts(
    tree_pos,
    num_nodes,
    num_mutations,
    edges_parent,
    edges_child,
    samples,
    mutations_position,
    mutations_node,
    mutations_parent,
):
    parent = np.zeros(num_nodes, dtype=np.int32) - 1
    num_samples = np.zeros(num_nodes, dtype=np.int32)
    num_samples[samples] = 1
    mutations_num_descendants = np.zeros(num_mutations, dtype=np.int32)
    mutations_num_inheritors = np.zeros(num_mutations, dtype=np.int32)

    mut_id = 0

    while tree_pos.next():
        for j in range(tree_pos.out_range[0], tree_pos.out_range[1]):
            e = tree_pos.edge_removal_order[j]
            c = edges_child[e]
            p = edges_parent[e]
            parent[c] = -1
            u = p
            while u != -1:
                num_samples[u] -= num_samples[c]
                u = parent[u]

        for j in range(tree_pos.in_range[0], tree_pos.in_range[1]):
            e = tree_pos.edge_insertion_order[j]
            p = edges_parent[e]
            c = edges_child[e]
            parent[c] = p
            u = p
            while u != -1:
                num_samples[u] += num_samples[c]
                u = parent[u]
        left, right = tree_pos.interval
        while mut_id < num_mutations and mutations_position[mut_id] < right:
            assert mutations_position[mut_id] >= left
            mutation_node = mutations_node[mut_id]
            descendants = num_samples[mutation_node]
            mutations_num_descendants[mut_id] = descendants
            mutations_num_inheritors[mut_id] = descendants
            # Subtract this number of descendants from the parent mutation. We are
            # guaranteed to list parents mutations before their children
            mut_parent = mutations_parent[mut_id]
            if mut_parent != -1:
                mutations_num_inheritors[mut_parent] -= descendants
            mut_id += 1

    return mutations_num_descendants, mutations_num_inheritors


def mutations(ts):
    # FIXME use tskit's impute mutations time
    mutations_time = ts.mutations_time.copy()
    mutations_node = ts.mutations_node.copy()
    unknown = tskit.is_unknown_time(mutations_time)
    mutations_time[unknown] = ts.nodes_time[mutations_node[unknown]]

    position = ts.sites_position[ts.mutations_site]

    tables = ts.tables
    derived_state = tables.mutations.derived_state
    offsets = tables.mutations.derived_state_offset
    derived_state = np.array(
        [
            derived_state[s].tobytes().decode("utf-8")
            for s in (
                slice(start, end) for start, end in zip(offsets[:-1], offsets[1:])
            )
        ]
    )
    ancestral_state = tables.sites.ancestral_state
    offsets = tables.sites.ancestral_state_offset
    ancestral_state = np.array(
        [
            ancestral_state[s].tobytes().decode("utf-8")
            for s in (
                slice(start, end) for start, end in zip(offsets[:-1], offsets[1:])
            )
        ]
    )
    del tables
    inherited_state = ancestral_state[ts.mutations_site]

    mutations_with_parent = ts.mutations_parent != -1

    parent = ts.mutations_parent[mutations_with_parent]
    assert np.all(parent >= 0)
    inherited_state[mutations_with_parent] = derived_state[parent]
    mutations_inherited_state = inherited_state

    population_data = {}
    if ts.num_populations > 0:
        pop_mutation_count = compute_population_mutation_counts(ts)
        for pop in ts.populations():
            name = f"pop{pop.id}"
            if isinstance(pop.metadata, bytes):
                try:
                    metadata_dict = json.loads(pop.metadata.decode("utf-8"))
                except (json.JSONDecodeError, UnicodeDecodeError):
                    metadata_dict = {}
            else:
                metadata_dict = pop.metadata
            if "name" in metadata_dict:
                name = metadata_dict["name"]
            col_name = f"pop_{name}_freq"
            population_data[col_name] = pop_mutation_count[pop.id] / ts.num_samples

    counts = compute_mutation_counts(ts)
    logger.info("Preprocessed mutations")
    return {
        "position": position,
        "inherited_state": mutations_inherited_state,
        "num_descendants": counts.num_descendants,
        "num_inheritors": counts.num_inheritors,
        "num_parents": counts.num_parents,
        **population_data,
    }


def edges(ts):
    parent_time = ts.nodes_time[ts.edges_parent]
    child_time = ts.nodes_time[ts.edges_child]
    branch_length = parent_time - child_time
    span = ts.edges_right - ts.edges_left

    logger.info("Preprocessed edges")
    return {
        "parent_time": parent_time,
        "child_time": child_time,
        "branch_length": branch_length,
        "span": span,
    }


@jit.numba_jit()
def child_bounds(num_nodes, edges_left, edges_right, edges_child):
    num_edges = edges_left.shape[0]
    child_left = np.zeros(num_nodes, dtype=np.float64) + np.inf
    child_right = np.zeros(num_nodes, dtype=np.float64)

    for e in range(num_edges):
        u = edges_child[e]
        if edges_left[e] < child_left[u]:
            child_left[u] = edges_left[e]
        if edges_right[e] > child_right[u]:
            child_right[u] = edges_right[e]
    return child_left, child_right


def nodes(ts):
    child_left, child_right = child_bounds(
        ts.num_nodes, ts.edges_left, ts.edges_right, ts.edges_child
    )
    nodes_num_mutations = np.bincount(ts.mutations_node, minlength=ts.num_nodes)
    logger.info("Preprocessed nodes")
    return {
        "num_mutations": nodes_num_mutations,
        "ancestors_span": child_right - child_left,
    }


def sites(ts):
    sites_num_mutations = np.bincount(ts.mutations_site, minlength=ts.num_sites)
    logger.info("Preprocessed sites")
    return {
        "num_mutations": sites_num_mutations,
    }


@jit.numba_jit()
def _compute_per_tree_stats(
    tree_pos, num_trees, num_nodes, nodes_time, edges_parent, edges_child
):
    tbl = np.zeros(num_trees)
    num_internal_nodes = np.zeros(num_trees)
    max_arity = np.zeros(num_trees, dtype=np.int32)
    num_children = np.zeros(num_nodes, dtype=np.int32)
    nodes_with_arity = np.zeros(num_nodes, dtype=np.int32)

    current_tbl = 0
    tree_index = 0
    current_num_internal_nodes = 0
    current_max_arity = 0
    while tree_pos.next():
        for j in range(tree_pos.out_range[0], tree_pos.out_range[1]):
            e = tree_pos.edge_removal_order[j]
            p = edges_parent[e]
            nodes_with_arity[num_children[p]] -= 1
            if (
                num_children[p] == current_max_arity
                and nodes_with_arity[num_children[p]] == 1
            ):
                current_max_arity -= 1

            num_children[p] -= 1
            if num_children[p] == 0:
                current_num_internal_nodes -= 1
            else:
                nodes_with_arity[num_children[p]] += 1
            c = edges_child[e]
            branch_length = nodes_time[p] - nodes_time[c]
            current_tbl -= branch_length

        for j in range(tree_pos.in_range[0], tree_pos.in_range[1]):
            e = tree_pos.edge_insertion_order[j]
            p = edges_parent[e]
            if num_children[p] == 0:
                current_num_internal_nodes += 1
            else:
                nodes_with_arity[num_children[p]] -= 1
            num_children[p] += 1
            nodes_with_arity[num_children[p]] += 1
            if num_children[p] > current_max_arity:
                current_max_arity = num_children[p]
            c = edges_child[e]
            branch_length = nodes_time[p] - nodes_time[c]
            current_tbl += branch_length
        tbl[tree_index] = current_tbl
        num_internal_nodes[tree_index] = current_num_internal_nodes
        max_arity[tree_index] = current_max_arity
        tree_index += 1
        # print("tree", tree_index, nodes_with_arity)

    return tbl, num_internal_nodes, max_arity


def compute_per_tree_stats(ts):
    """
    Returns the per-tree statistics
    """
    tree_pos = alloc_tree_position(ts)
    return _compute_per_tree_stats(
        tree_pos,
        ts.num_trees,
        ts.num_nodes,
        ts.nodes_time,
        ts.edges_parent,
        ts.edges_child,
    )


def trees(ts):
    num_trees = ts.num_trees
    total_branch_length, num_internal_nodes, max_arity = compute_per_tree_stats(ts)

    # FIXME - need to add this to the computation above
    mean_internal_arity = np.zeros(num_trees)

    site_tree_index = (
        np.searchsorted(ts.breakpoints(as_array=True), ts.sites_position, side="right")
        - 1
    )
    unique_values, counts = np.unique(site_tree_index, return_counts=True)
    sites_per_tree = np.zeros(ts.num_trees, dtype=np.int64)
    sites_per_tree[unique_values] = counts

    mutation_tree_index = site_tree_index[ts.mutations_site]
    unique_mutation_values, mutation_counts = np.unique(
        mutation_tree_index, return_counts=True
    )
    mutations_per_tree = np.zeros(ts.num_trees, dtype=np.int64)
    mutations_per_tree[unique_mutation_values] = mutation_counts

    breakpoints = ts.breakpoints(as_array=True)
    logger.info("Pre processed trees")
    return {
        "left": breakpoints[:-1],
        "right": breakpoints[1:],
        "total_branch_length": total_branch_length,
        "mean_internal_arity": mean_internal_arity,
        "max_internal_arity": max_arity,
        "num_sites": sites_per_tree,
        "num_mutations": mutations_per_tree,
    }


def preprocess(tszip_path, output_path, show_progress=False):
    tszip_path = pathlib.Path(tszip_path)
    preprocessors = [mutations, nodes, trees, edges, sites]
    with tqdm(
        total=2 + len(preprocessors), desc="Processing", disable=not show_progress
    ) as pbar:
        logger.info(f"Loading tree sequence from {tszip_path}")
        ts = tszip.load(tszip_path)

        # Check if all mutation times are unknown, calculate them if so
        if np.all(tskit.is_unknown_time(ts.mutations_time)):
            warnings.warn(
                "All mutation times are unknown. Calculating mutation times"
                " from tree sequence",
                stacklevel=1,
            )
            tables = ts.dump_tables()
            tables.compute_mutation_times()
            ts = tables.tree_sequence()

        pbar.update(1)

        logger.info(f"Compressing tree sequence to {output_path}")
        tszip.compress(ts, output_path)
        pbar.update(1)

        # Preprocess the data first so we can error out before writing to the file
        data = {}
        for preprocessor in preprocessors:
            group_name = preprocessor.__name__.split(".")[-1]
            logger.info(f"Processing {group_name}")
            pbar.set_description(f"Processing {group_name}")
            data[group_name] = preprocessor(ts)
            pbar.update(1)

    logger.info(f"Writing preprocessed data to {output_path}")
    with zarr.ZipStore(output_path, mode="a") as zarr_store:
        root = zarr.group(store=zarr_store)
        total_arrays = sum(len(arrays) for arrays in data.values())
        with tqdm(
            total=total_arrays, desc="Writing", disable=not show_progress
        ) as pbar:
            for table_name, arrays in data.items():
                for array_name, array in arrays.items():
                    logger.info(f"Writing {table_name}/{array_name}")
                    root[f"{table_name}/{array_name}"] = array
                    pbar.update(1)
        root.attrs["tsbrowse"] = {"data_version": TSBROWSE_DATA_VERSION}
    logger.info("Finished writing preprocessed data")

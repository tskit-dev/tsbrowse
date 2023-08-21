"""
QC functions for tsinfer trees
"""
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def plot_polytomy_fractions(
    tsm, region_start=None, region_end=None, window_size=100_000, overlap=0
):
    """
    Plots the fraction of polytomies in windows actoss the genomic sequence
    """
    if region_start is None:
        region_start = max(0, tsm.ts.tables.sites.position[0] - 50_000)
    if region_end is None:
        region_end = tsm.ts.tables.sites.position[-1] + 50_000
    fig, ax = plt.subplots(figsize=(20, 5))
    polytomy_fractions = tsm.calc_polytomy_fractions()
    poly_fracs_by_pos = tsm.map_stats_to_genome(polytomy_fractions)
    poly_fracs_means = []
    poly_fracs_sd = []
    genomic_positions = []
    for poly_win in tsm.make_sliding_windows(
        poly_fracs_by_pos, window_size, overlap
    ):
        poly_fracs_means.append(np.mean(poly_win))
        poly_fracs_sd.append(np.std(poly_win))
    for gen_win in tsm.make_sliding_windows(
        np.arange(1, tsm.ts.sequence_length), window_size, overlap
    ):
        genomic_positions.append(gen_win[0] / 1_000_000)
    ax.plot(
        genomic_positions,
        poly_fracs_means,
        label="mean",
        linewidth=0.5,
    )
    ax.fill_between(
        genomic_positions,
        np.array(poly_fracs_means) - np.array(poly_fracs_sd),
        np.array(poly_fracs_means) + np.array(poly_fracs_sd),
        alpha=0.3,
        label="mean +/- std",
    )
    missing_vals = np.take(genomic_positions, np.where(np.isnan(poly_fracs_means)))
    ax.plot(
        missing_vals,
        np.zeros(len(missing_vals)),
        color="red",
        marker="o",
        label="missing data",
    )
    ax.set_xlabel(f"Position (Mb)", fontsize=10)
    ax.set_ylabel("Window mean", fontsize=10)
    ax.set_title("Polytomy score", fontsize=10)
    ax.set_ylim(0, 1)
    ax.set_xlim(region_start / 1_000_000, region_end / 1_000_000)
    handles, labels = ax.get_legend_handles_labels()
    unique = [
        (h, l)
        for i, (h, l) in enumerate(zip(handles, labels))
        if l not in labels[:i]
    ]
    ax.legend(*zip(*unique))
    plt.show()

def plot_mutations_per_site(tsm, max_num_muts=None, show_counts=False):
    fig, ax = plt.subplots()
    bins = None
    plt.xlabel("Number of mutations")
    if max_num_muts is not None:
        bins = range(max_num_muts + 1)
        sites_with_many_muts = np.sum(tsm.sites_num_mutations > max_num_muts)
        plt.xlabel(
            f"Number of mutations\n\n\nThere are {sites_with_many_muts:,} sites with more than {max_num_muts:,} mutations"
        )
    counts, edges, bars = plt.hist(
        tsm.sites_num_mutations, bins=bins, edgecolor="black"
    )
    ax.set_xticks(edges)
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, pos: "{:,}".format(int(x)))
    )
    plt.ylabel("Number of sites")
    plt.title("Mutations-per-site distribution")
    if show_counts:
        plt.bar_label(bars, fmt="{:,.0f}")

def plot_mutations_per_site_along_seq(
    tsm, region_start=None, region_end=None, hist_bins=1000
):
    count = tsm.sites_num_mutations
    pos = tsm.ts.sites_position
    if region_start is None:
        region_start = pos[0]
    if region_end is None:
        region_end = pos[-1]
    grid = sns.jointplot(
        x=pos / 1_000_000,
        y=count,
        kind="scatter",
        marginal_ticks=True,
        alpha=0.5,
        marginal_kws=dict(bins=hist_bins),
        xlim=(region_start / 1_000_000, region_end / 1_000_000),
    )
    grid.ax_marg_y.remove()
    grid.fig.set_figwidth(20)
    grid.fig.set_figheight(8)
    grid.ax_joint.set_xlabel("Position on genome (Mb)")
    grid.ax_joint.set_ylabel("Number of mutations")

def plot_mutations_per_node(tsm, max_num_muts=None, show_counts=False):
    fig, ax = plt.subplots()
    bins = None
    plt.xlabel(f"Number of mutations")
    if max_num_muts is not None:
        bins = range(max_num_muts + 1)
        nodes_with_many_muts = np.sum(tsm.nodes_num_mutations > max_num_muts)
        plt.xlabel(
            f"Number of mutations \n\n\nThere are {nodes_with_many_muts:,} nodes with more than {max_num_muts:,} mutations"
        )

    counts, edges, bars = plt.hist(
        tsm.nodes_num_mutations, bins=bins, edgecolor="black"
    )
    ax.set_xticks(edges)
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, pos: "{:,}".format(int(x)))
    )
    plt.ylabel("Number of nodes")
    plt.title("Mutations-per-node distribution")
    if show_counts:
        plt.bar_label(bars, fmt="{:,.0f}")

def plot_tree_spans(
    tsm, log_transform=True, region_start=None, region_end=None, show_counts=False
):
    fig, ax = plt.subplots()
    bins = None
    breakpoints = tsm.ts.breakpoints(as_array=True)
    start_idx = 2
    end_idx = len(breakpoints) - 1

    if region_start is not None:
        start_idx = max(start_idx, np.argmax(breakpoints > region_start))
    if region_end is not None:
        end_idx = min(np.argmax(breakpoints >= region_end), end_idx)

    spans = (
        breakpoints[start_idx:end_idx] - breakpoints[start_idx - 1 : end_idx - 1]
    )
    xlabel = "span"
    if log_transform:
        spans = np.log10(spans)
        xlabel = "span (log10)"
        bins = range(int(np.min(spans)), int(np.max(spans)) + 2)

    counts, edges, bars = plt.hist(spans, edgecolor="black", bins=bins)
    ax.set_xticks(edges)
    if show_counts:
        plt.bar_label(bars, fmt="{:,.0f}")
    ax.set_xlabel(xlabel)
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, pos: "{:,}".format(int(x)))
    )
    plt.title(f"Distribution of {len(spans):,} tree spans")


def plot_mean_node_arity(tsm, show_counts=False):
    fig, ax = plt.subplots()
    mean_arity = tsm.calc_mean_node_arity()
    counts, edges, bars = plt.hist(mean_arity, bins=None, edgecolor="black")
    ax.set_xlabel("Mean node arity")
    ax.set_ylabel("Number of nodes")
    ax.set_title("Mean-node-arity distribution")
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, pos: "{:,}".format(int(x)))
    )
    if show_counts:
        plt.bar_label(bars, fmt="{:,.0f}")

def plot_mutations_per_tree(tsm, max_num_muts=None, show_counts=False):
    fig, ax = plt.subplots()
    tree_mutations = tsm.calc_mutations_per_tree()
    bins = max(100, int(np.sqrt(tsm.ts.num_trees)))
    plt.xlabel(f"Number of mutations")
    if max_num_muts is not None:
        bins = range(max_num_muts + 1)
        trees_with_many_muts = np.sum(tree_mutations > max_num_muts)
        plt.xlabel(
            f"Number of mutations\n\n\nThere are {trees_with_many_muts:,} trees with more than {max_num_muts:,} mutations"
        )

    counts, edges, bars = plt.hist(
        tsm.calc_mutations_per_tree(), bins=bins, edgecolor="black"
    )
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, pos: "{:,}".format(int(x)))
    )
    plt.ylabel("Number of trees")
    plt.title("Mutations-per-tree distribution")
    if show_counts:
        plt.bar_label(bars, fmt="{:,.0f}")

def plot_mutations_per_tree_along_seq(
    tsm, region_start=None, region_end=None, hist_bins=1000
):
    tree_mutations = tsm.calc_mutations_per_tree()
    tree_mutations = tree_mutations[1:-1]
    breakpoints = tsm.ts.breakpoints(as_array=True)
    tree_mids = breakpoints[1:] - ((breakpoints[1:] - breakpoints[:-1]) / 2)
    tree_mids = tree_mids[1:-1]
    if region_start is None or region_start < tree_mids[0]:
        region_start = tree_mids[0]
    if region_end is None or region_end > tree_mids[-1]:
        region_end = tree_mids[-1]

    grid = sns.jointplot(
        x=tree_mids / 1_000_000,
        y=tree_mutations,
        kind="scatter",
        marginal_ticks=True,
        alpha=0.5,
        marginal_kws=dict(bins=hist_bins),
        xlim=(region_start / 1_000_000, region_end / 1_000_000),
        # set ylim to the max number of sites in a tree in the region
        ylim=(
            0,
            np.max(
                tree_mutations[
                    (tree_mids >= region_start) & (tree_mids <= region_end)
                ]
            ),
        ),
    )
    grid.ax_marg_y.remove()
    grid.fig.set_figwidth(20)
    grid.fig.set_figheight(8)
    grid.ax_joint.set_xlabel("Position on genome (Mb)")
    grid.ax_joint.set_ylabel("Number of mutations per tree")

def plot_sites_per_tree(tsm, max_num_sites=None, show_counts=False):
    fig, ax = plt.subplots()
    bins = max(100, int(np.sqrt(tsm.ts.num_trees)))
    plt.xlabel(f"Number of sites")
    if max_num_sites is not None:
        bins = range(max_num_sites + 1)
        trees_with_many_sites = np.sum(tsm.calc_sites_per_tree() > max_num_sites)
        plt.xlabel(
            f"Number of sites\n\n\nThere are {trees_with_many_sites:,} trees with more than {max_num_sites:,} sites"
        )

    counts, edges, bars = plt.hist(
        tsm.calc_sites_per_tree(), bins=bins, edgecolor="black"
    )
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, pos: "{:,}".format(int(x)))
    )

    plt.ylabel("Number of trees")
    plt.title("Sites-per-tree distribution")
    if show_counts:
        plt.bar_label(bars, fmt="{:,.0f}")

def plot_sites_per_tree_along_seq(
    tsm, region_start=None, region_end=None, hist_bins=500
):
    tree_sites = tsm.calc_sites_per_tree()
    tree_sites = tree_sites[1:-1]
    breakpoints = tsm.ts.breakpoints(as_array=True)
    tree_mids = breakpoints[1:] - ((breakpoints[1:] - breakpoints[:-1]) / 2)
    tree_mids = tree_mids[1:-1]
    if region_start is None or region_start < tree_mids[0]:
        region_start = tree_mids[0]
    if region_end is None or region_end > tree_mids[-1]:
        region_end = tree_mids[-1]

    grid = sns.jointplot(
        x=tree_mids / 1_000_000,
        y=tree_sites,
        kind="scatter",
        marginal_ticks=True,
        alpha=0.5,
        marginal_kws=dict(bins=hist_bins),
        xlim=(region_start / 1_000_000, region_end / 1_000_000),
        ylim=(
            0,
            np.max(
                tree_sites[(tree_mids >= region_start) & (tree_mids <= region_end)]
            ),
        ),
    )
    grid.ax_marg_y.remove()
    grid.fig.set_figwidth(20)
    grid.fig.set_figheight(8)
    grid.ax_joint.set_xlabel("Position on genome (Mb)")
    grid.ax_joint.set_ylabel("Number of sites per tree")

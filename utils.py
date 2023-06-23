import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from bokeh.models import RangeSlider
from bokeh.layouts import layout
from bokeh.plotting import figure, output_notebook, show

"""
QC functions for tsinfer trees
"""


class TreeInfo:
    """
    Class for storing tree information
    """

    def __init__(self, ts, chr):
        self.ts = ts
        self.chr = chr

        self.sites_num_mutations = np.bincount(
            self.ts.mutations_site, minlength=self.ts.num_sites
        )
        self.nodes_num_mutations = np.bincount(
            self.ts.mutations_node, minlength=self.ts.num_nodes
        )

    def summary(self):
        nodes_with_zero_muts = np.sum(self.nodes_num_mutations == 0)
        sites_with_zero_muts = np.sum(self.sites_num_mutations == 0)

        data = [
            ("samples", self.ts.num_samples),
            ("nodes", self.ts.num_nodes),
            ("mutations", self.ts.num_mutations),
            ("nodes_with_zero_muts", nodes_with_zero_muts),
            ("sites_with_zero_muts", sites_with_zero_muts),
            ("max_mutations_per_site", np.max(self.sites_num_mutations)),
            ("mean_mutations_per_site", np.mean(self.sites_num_mutations)),
            ("median_mutations_per_site", np.median(self.sites_num_mutations)),
            ("max_mutations_per_node", np.max(self.nodes_num_mutations)),
        ]
        df = pd.DataFrame(
            {"property": [d[0] for d in data], "value": [d[1] for d in data]}
        )
        return df.set_index("property")

    def _repr_html_(self):
        return self.summary()._repr_html_()

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

    def plot_polytomy_fractions(
        self, region_start=None, region_end=None, window_size=100_000, overlap=0
    ):
        """
        Plots the fraction of polytomies in windows actoss the genomic sequence
        """
        if region_start is None:
            region_start = max(0, self.ts.tables.sites.position[0] - 50_000)
        if region_end is None:
            region_end = self.ts.tables.sites.position[-1] + 50_000
        fig, ax = plt.subplots(figsize=(20, 5))
        polytomy_fractions = self.calc_polytomy_fractions()
        poly_fracs_by_pos = self.map_stats_to_genome(polytomy_fractions)
        poly_fracs_means = []
        poly_fracs_sd = []
        genomic_positions = []
        for poly_win in self.make_sliding_windows(
            poly_fracs_by_pos, window_size, overlap
        ):
            poly_fracs_means.append(np.mean(poly_win))
            poly_fracs_sd.append(np.std(poly_win))
        for gen_win in self.make_sliding_windows(
            np.arange(1, self.ts.sequence_length), window_size, overlap
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
        missing_vals = np.take(
            genomic_positions, np.where(np.isnan(poly_fracs_means)))
        ax.plot(
            missing_vals,
            np.zeros(len(missing_vals)),
            color="red",
            marker="o",
            label="missing data",
        )
        ax.set_xlabel(f"Position on chr {self.chr}(Mb)", fontsize=10)
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

    def plot_mutations_per_site(self, max_num_muts=None, show_counts=False):
        fig, ax = plt.subplots()
        bins = None
        if max_num_muts is not None:
            bins = range(max_num_muts + 1)
            sites_with_many_muts = np.sum(
                self.sites_num_mutations > max_num_muts)
            ax.text(
                0.5,
                0.9,
                f"there are {sites_with_many_muts:,} sites\nwith more than {max_num_muts:,} mutations",
                transform=ax.transAxes,
            )
        counts, edges, bars = plt.hist(
            self.sites_num_mutations, bins=bins, edgecolor="black"
        )
        ax.set_xticks(edges)
        ax.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, pos: "{:,}".format(int(x)))
        )
        plt.xlabel("Number of mutations")
        plt.ylabel("Number of sites")
        plt.title("Mutations-per-site distribution")
        if show_counts:
            plt.bar_label(bars, fmt="{:,.0f}")

    def plot_mutations_per_site_along_seq(
        self, region_start=None, region_end=None, hist_bins=1000
    ):
        count = self.sites_num_mutations
        pos = self.ts.sites_position
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

    def plot_mutations_per_node(self, max_num_muts=None, show_counts=False):
        fig, ax = plt.subplots()
        bins = None
        if max_num_muts is not None:
            bins = range(max_num_muts + 1)
            nodes_with_many_muts = np.sum(
                self.nodes_num_mutations > max_num_muts)
            ax.text(
                0.5,
                0.9,
                f"there are {nodes_with_many_muts:,} nodes\nwith more than {max_num_muts:,} mutations",
                transform=ax.transAxes,
            )
        counts, edges, bars = plt.hist(
            self.nodes_num_mutations, bins=bins, edgecolor="black"
        )
        ax.set_xticks(edges)
        ax.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, pos: "{:,}".format(int(x)))
        )
        plt.xlabel("Number of mutations")
        plt.ylabel("Number of nodes")
        plt.title("Mutations-per-node distribution")
        if show_counts:
            plt.bar_label(bars, fmt="{:,.0f}")

    def plot_tree_spans(
        self, log_transform=True, region_start=None, region_end=None, show_counts=False
    ):
        fig, ax = plt.subplots()
        bins = None
        breakpoints = self.ts.breakpoints(as_array=True)
        start_idx = 2
        end_idx = len(breakpoints) - 1

        if region_start is not None:
            start_idx = max(start_idx, np.argmax(breakpoints > region_start))
        if region_end is not None:
            end_idx = min(np.argmax(breakpoints >= region_end), end_idx)

        spans = (
            breakpoints[start_idx:end_idx] -
            breakpoints[start_idx - 1: end_idx - 1]
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

    def plot_mean_node_arity(self, show_counts=False):
        fig, ax = plt.subplots()
        mean_arity = self.calc_mean_node_arity()
        counts, edges, bars = plt.hist(
            mean_arity, bins=None, edgecolor="black")
        ax.set_xlabel("Mean node arity")
        ax.set_ylabel("Number of nodes")
        ax.set_title("Mean-node-arity distribution")
        ax.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, pos: "{:,}".format(int(x)))
        )
        if show_counts:
            plt.bar_label(bars, fmt="{:,.0f}")

    def plot_mutations_per_site_along_seq_bokeh(self, num_sites=100_000):
        counts = self.sites_num_mutations
        pos = self.ts.sites_position
        pos = pos[:num_sites]
        counts = counts[:num_sites]
        p = figure(
            title="mutations per site along the genome",
            x_axis_label="position (bp)",
            y_axis_label="mutations per site",
        )
        p.height = 400
        p.width = 1000
        p.circle(pos, counts, size=5, alpha=0.5)
        range_slider = RangeSlider(
            start=pos[0],
            end=pos[-1],
            value=(pos[0], pos[-1]),
            step=100_000,
            title="start position (bp)",
        )
        range_slider.js_link("value", p.x_range, "start", attr_selector=0)
        range_slider.js_link("value", p.x_range, "end", attr_selector=1)
        lo = layout([range_slider], [p])
        output_notebook()
        show(lo)
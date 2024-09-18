import numpy as np
import panel as pn

from .. import config
from ..plot_helpers import make_hist


class TreesPage:
    key = "trees"
    title = "Trees"

    def __init__(self, tsm):
        df_trees = tsm.trees_df
        bins = min(50, int(np.sqrt(len(df_trees))))

        spans = df_trees.right - df_trees.left

        log_y_checkbox = pn.widgets.Checkbox(name="log y-axis", value=True)

        plot_options = pn.Column(
            pn.pane.Markdown("# Plot Options"),
            log_y_checkbox,
        )

        def make_tree_hist_panel(tsm, log_y):
            sites_hist = make_hist(
                df_trees.num_sites,
                "Sites per tree",
                bins,
                log_y=log_y,
                plot_width=config.PLOT_WIDTH,
                xlabel="number of sites",
                ylabel="number of trees",
            )

            spans_hist = make_hist(
                spans,
                "Genomic span per tree",
                bins,
                log_y=log_y,
                plot_width=config.PLOT_WIDTH,
                xlabel="genomic span",
                ylabel="number of trees",
            )

            muts_hist = make_hist(
                df_trees.num_mutations,
                "Mutations per tree",
                bins,
                log_y=log_y,
                plot_width=config.PLOT_WIDTH,
                xlabel="number of mutations",
                ylabel="number of trees",
            )

            tbl_hist = make_hist(
                df_trees.total_branch_length,
                "Total branch length per tree",
                bins,
                log_y=log_y,
                plot_width=config.PLOT_WIDTH,
                xlabel="total branch length",
                ylabel="number of trees",
            )

            mean_arity_hist = make_hist(
                df_trees.mean_internal_arity,
                "Mean arity per tree \n(not yet implemented)",
                bins,
                log_y=log_y,
                plot_width=config.PLOT_WIDTH,
                xlabel="mean arity",
                ylabel="number of trees",
            )

            max_arity_hist = make_hist(
                df_trees.max_internal_arity,
                "Max arity per tree",
                bins,
                log_y=log_y,
                plot_width=config.PLOT_WIDTH,
                xlabel="max arity",
                ylabel="number of trees",
            )
            return pn.Column(
                pn.Row(
                    sites_hist,
                    spans_hist,
                ),
                pn.Row(
                    muts_hist,
                    tbl_hist,
                ),
                pn.Row(mean_arity_hist, max_arity_hist),
            )

        hist_panel = pn.bind(make_tree_hist_panel, log_y=log_y_checkbox, tsm=tsm)

        self.content = pn.Column(hist_panel)
        self.sidebar = plot_options

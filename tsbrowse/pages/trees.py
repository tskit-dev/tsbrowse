import numpy as np
import panel as pn

from ..plot_helpers import make_hist


class TreesPage:
    key = "trees"
    title = "Trees"

    def __init__(self, tsm):
        df_trees = tsm.trees_df
        bins = min(50, int(np.sqrt(len(df_trees))))

        spans = df_trees.right - df_trees.left

        log_y_checkbox = pn.widgets.Checkbox(name="Log Y-axis", value=True)

        plot_options = pn.Column(
            pn.pane.Markdown("# Trees"),
            log_y_checkbox,
        )

        def make_tree_hist_panel(tsm, log_y):
            gspec = pn.GridSpec(nrows=3, ncols=2, sizing_mode="stretch_both")

            sites_hist = make_hist(
                df_trees.num_sites,
                "Sites per tree",
                bins,
                log_y=log_y,
                xlabel="Number of Sites",
                ylabel="Number of Trees",
            )

            spans_hist = make_hist(
                spans,
                "Genomic span per tree",
                bins,
                log_y=log_y,
                xlabel="Genomic Span",
                ylabel="Number of Trees",
            )

            muts_hist = make_hist(
                df_trees.num_mutations,
                "Mutations per tree",
                bins,
                log_y=log_y,
                xlabel="Number of Mutations",
                ylabel="Number of Trees",
            )

            tbl_hist = make_hist(
                df_trees.total_branch_length,
                "Total branch length per tree",
                bins,
                log_y=log_y,
                xlabel="Total Banch Length",
                ylabel="Number of Trees",
            )

            mean_arity_hist = make_hist(
                df_trees.mean_internal_arity,
                "Mean arity per tree \n(not yet implemented)",
                bins,
                log_y=log_y,
                xlabel="Mean Arity",
                ylabel="Number of Trees",
            )

            max_arity_hist = make_hist(
                df_trees.max_internal_arity,
                "Max arity per tree",
                bins,
                log_y=log_y,
                xlabel="Max Arity",
                ylabel="Number of Trees",
            )

            # Assign histograms to grid positions
            gspec[0, 0] = sites_hist
            gspec[0, 1] = spans_hist
            gspec[1, 0] = muts_hist
            gspec[1, 1] = tbl_hist
            gspec[2, 0] = mean_arity_hist
            gspec[2, 1] = max_arity_hist

            return gspec

        hist_panel = pn.bind(make_tree_hist_panel, log_y=log_y_checkbox, tsm=tsm)

        self.content = pn.Column(hist_panel, sizing_mode="stretch_both")
        self.sidebar = plot_options

import panel as pn
import holoviews as hv
import numpy as np
import functools

from plot_helpers import make_hist_matplotlib


def page(tsm):
    hv.extension("matplotlib")
    df_trees = tsm.trees_df
    bins = min(50, int(np.sqrt(len(df_trees))))

    sites_hist_func = functools.partial(
        make_hist_matplotlib,
        df_trees.num_sites,
        "Sites per tree",
        num_bins=bins,
        log_y=True,
    )

    log_y_checkbox = pn.widgets.Checkbox(name="log y-axis", value=True)

    sites_hist_panel = pn.bind(
        sites_hist_func,
        log_y=log_y_checkbox,
    )

    spans = df_trees.right - df_trees.left

    spans_hist_func = functools.partial(
        make_hist_matplotlib,
        spans,
        "Genomic span per tree",
        num_bins=bins,
        log_y=True,
    )

    spans_hist_panel = pn.bind(
        spans_hist_func,
        log_y=log_y_checkbox,
    )

    muts_hist_func = functools.partial(
        make_hist_matplotlib,
        df_trees.num_mutations,
        "Mutations per tree",
        num_bins=bins,
        log_y=True,
    )

    muts_hist_panel = pn.bind(
        muts_hist_func,
        log_y=log_y_checkbox,
    )

    tbl_hist_func = functools.partial(
        make_hist_matplotlib,
        df_trees.total_branch_length,
        "Total branch length per tree",
        num_bins=bins,
        log_y=True,
    )

    tbl_hist_panel = pn.bind(
        tbl_hist_func,
        log_y=log_y_checkbox,
    )

    mean_arity_hist_func = functools.partial(
        make_hist_matplotlib,
        df_trees.mean_internal_arity,
        f"Mean arity per tree \n(not yet implemented)",
        num_bins=bins,
        log_y=True,
    )

    mean_arity_hist_panel = pn.bind(
        mean_arity_hist_func,
        log_y=log_y_checkbox,
    )

    max_arity_hist_func = functools.partial(
        make_hist_matplotlib,
        df_trees.max_internal_arity,
        "Max arity per tree",
        num_bins=bins,
        log_y=True,
    )

    max_arity_hist_panel = pn.bind(
        max_arity_hist_func,
        log_y=log_y_checkbox,
    )

    plot_options = pn.Column(
        pn.pane.Markdown("# Plot Options"),
        log_y_checkbox,
    )

    hist_panel = pn.Column(
        pn.Row(
            sites_hist_panel,
            spans_hist_panel,
            muts_hist_panel,
        ),
        pn.Row(tbl_hist_panel, mean_arity_hist_panel, max_arity_hist_panel),
    )

    return pn.Column(hist_panel, plot_options)

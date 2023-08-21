import panel as pn
import holoviews as hv
import numpy as np
import functools

import plot_helpers


def page(tsm):
    hv.extension("matplotlib")
    trees_df = tsm.trees_df
    max_num_sites = int(trees_df.num_sites.max())
    min_num_sites = int(trees_df.num_sites.min())
    sites_default_num_bins = min(int(np.sqrt(max_num_sites)), 50)

    spans = trees_df.right - trees_df.left
    min_span = int(np.min(spans))
    max_span = int(np.max(spans))
    spans_default_num_bins = min(50, int(np.sqrt(max_span)))

    sites_hist_func = functools.partial(
        plot_helpers.make_hist_matplotlib,
        trees_df.num_sites,
        "Sites per tree",
        num_bins=sites_default_num_bins,
        log_y=True,
        xlim=(min_num_sites, max_num_sites),
    )

    spans_hist_func = functools.partial(
        plot_helpers.make_hist_matplotlib,
        spans,
        "Genomic span per tree",
        num_bins=spans_default_num_bins,
        log_y=True,
        xlim=(min_span, max_span),
    )

    if len(trees_df) > 1:
        sites_log_checkbox = pn.widgets.Checkbox(
            name="log y-axis", value=True, align=("center")
        )

        sites_xlim_slider = pn.widgets.IntRangeSlider(
            name="x limits",
            start=min_num_sites,
            value=(min_num_sites, max_num_sites),
            end=max_num_sites,
            step=10,  # TODO set a sensible step size
            align=("center"),
        )

        sites_bins_slider = pn.widgets.IntSlider(
            name="Number of bins",
            start=1,
            end=max(500, int(max_num_sites)),  # TODO set a sensible end value
            value=sites_default_num_bins,
            step=1,  # TODO set a sensible step size
            align=("center"),
        )

        sites_reset_checkbox = pn.widgets.Checkbox(
            name="reset", value=True, align=("center")
        )

        def reset_sites(reset=False):
            if reset:
                sites_log_checkbox.value = True
                sites_xlim_slider.value = (min_num_sites, max_num_sites)
                sites_bins_slider.value = sites_default_num_bins
                sites_reset_checkbox.value = False

        sites_reset_panel = pn.bind(reset_sites, reset=sites_reset_checkbox)

        sites_hist_panel = pn.bind(
            sites_hist_func,
            log_y=sites_log_checkbox,
            xlim=sites_xlim_slider,
            num_bins=sites_bins_slider,
        )

        spans_log_checkbox = pn.widgets.Checkbox(
            name="log y-axis", value=True, align=("center")
        )

        spans_xlim_slider = pn.widgets.IntRangeSlider(
            name="Maximum tree span",
            start=min_span,
            value=(min_span, max_span),
            end=max_span,
            step=1_000,  # TODO set a sensible step size
            align=("center"),
        )

        spans_bins_slider = pn.widgets.IntSlider(
            name="Number of bins",
            start=1,
            end=min(int(max_span), 500),  # TODO set a sensible end value
            value=spans_default_num_bins,
            step=1,  # TODO set a sensible step size
            align=("center"),
        )

        spans_reset_checkbox = pn.widgets.Checkbox(
            name="reset", value=True, align=("center")
        )

        def reset_spans(reset=False):
            if reset:
                spans_log_checkbox.value = True
                spans_xlim_slider.value = (min_span, max_span)
                spans_bins_slider.value = spans_default_num_bins
                spans_reset_checkbox.value = False

        spans_reset_panel = pn.bind(reset_spans, reset=spans_reset_checkbox)

        spans_hist_panel = pn.bind(
            spans_hist_func,
            log_y=spans_log_checkbox,
            xlim=spans_xlim_slider,
            num_bins=spans_bins_slider,
        )

        hist_panel = pn.Row(
            pn.Column(
                sites_reset_panel,
                sites_hist_panel,
                pn.Column(
                    pn.pane.Markdown("Plot Options:", align=("center")),
                    sites_reset_checkbox,
                    sites_log_checkbox,
                    sites_xlim_slider,
                    sites_bins_slider,
                ),
            ),
            pn.Column(
                spans_reset_panel,
                spans_hist_panel,
                pn.Column(
                    pn.pane.Markdown("Plot Options:", align=("center")),
                    spans_reset_checkbox,
                    spans_log_checkbox,
                    spans_xlim_slider,
                    spans_bins_slider,
                ),
            ),
        )
    else:
        sites_hist_single_tree = plot_helpers.make_hist_matplotlib(
            trees_df.num_sites, "Sites per tree", num_bins="auto", log_y=False
        )
        sites_hist_single_tree.opts(xticks=[trees_df.num_sites[0]], yticks=[1])
        spans_hist_single_tree = plot_helpers.make_hist_matplotlib(
            spans, "Genomic span per tree", num_bins="auto", log_y=False
        )
        spans_hist_single_tree.opts(xticks=[spans[0]], yticks=[1])
        hist_panel = pn.Row(
            pn.Column(sites_hist_single_tree), pn.Column(spans_hist_single_tree)
        )

    return hist_panel

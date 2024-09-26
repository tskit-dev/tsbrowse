import holoviews as hv
import numpy as np
import panel as pn
from holoviews import opts


def make_one_way_stats(ts, windows_trees, windows_count, span_normalise, statistic):
    windows_count = windows_count
    if windows_trees:
        windows = "trees"
        x_axis = ts.breakpoints(as_array=True)[:-1]
    else:
        windows = np.concatenate(
            (
                np.arange(
                    0, ts.sequence_length, max(1, ts.sequence_length // windows_count)
                ),
                [ts.sequence_length],
            )
        )
        x_axis = windows[:-1]

    if statistic == "Diversity":
        site_statistic = ts.diversity(
            span_normalise=span_normalise, windows=windows, mode="site"
        )
        branch_statistic = ts.diversity(
            span_normalise=span_normalise, windows=windows, mode="branch"
        )
    else:
        site_statistic = ts.segregating_sites(
            span_normalise=span_normalise, windows=windows, mode="site"
        )
        branch_statistic = ts.segregating_sites(
            span_normalise=span_normalise, windows=windows, mode="branch"
        )

    ratio = np.divide(branch_statistic, site_statistic, where=(site_statistic != 0))
    modes = {"Site": site_statistic, "Branch": branch_statistic, "Branch/Site": ratio}
    stat_curves = [
        hv.Curve((x_axis, v), "Genomic position", k).opts(interpolation="steps-post")
        for k, v in modes.items()
    ]

    layout = hv.Layout(stat_curves).cols(1)
    layout.opts(opts.Curve(height=200, responsive=True))
    return layout


class PopgenPage:
    key = "popgen"
    title = "Pop Gen"

    def __init__(self, tsm):
        ts = tsm.ts
        windows_trees = pn.widgets.Checkbox(name="Tree-based windows", value=False)
        windows_count = pn.widgets.IntSlider(
            name="Window count", start=1, end=100_000, value=1000
        )
        span_normalise = pn.widgets.Checkbox(name="Span normalise", value=True)
        stat_radiobox = pn.widgets.RadioButtonGroup(
            name="Statistic",
            options=["Diversity", "Segregating Sites"],
            value="Segregating Sites",
        )
        plot_options = pn.Column(
            pn.Row(pn.Column(windows_trees, span_normalise), stat_radiobox),
            windows_count,
        )
        one_way_panel = pn.bind(
            make_one_way_stats,
            ts=ts,
            windows_trees=windows_trees,
            windows_count=windows_count,
            span_normalise=span_normalise,
            statistic=stat_radiobox,
        )
        windows_trees.jslink(windows_count, value="disabled")
        self.content = pn.Column(one_way_panel)
        self.sidebar = pn.Column(pn.pane.Markdown("# Pop Gen"), plot_options)

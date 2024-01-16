import holoviews as hv
import numpy as np
import panel as pn
from holoviews import opts


def make_one_way_stats(ts, windows_trees, windows_count, span_normalise, mode):
    mode = mode.lower().split()[0]
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
    diversity = ts.diversity(span_normalise=span_normalise, windows=windows, mode=mode)
    tajimas = ts.Tajimas_D(windows=windows, mode=mode)
    diversity_curve = hv.Curve(
        (x_axis, diversity), "Genomic position", "Diversity"
    ).opts(interpolation="steps-post")
    tajimas_curve = hv.Curve((x_axis, tajimas), "Genomic position", "Tajimas D").opts(
        interpolation="steps-post"
    )

    layout = hv.Layout([diversity_curve, tajimas_curve]).cols(1)
    layout.opts(opts.Curve(height=200, responsive=True))
    return layout


def page(tsm):
    hv.extension("bokeh")
    ts = tsm.ts
    windows_trees = pn.widgets.Checkbox(name="Tree-based windows", value=False)
    windows_count = pn.widgets.IntSlider(
        name="Window count", start=1, end=100_000, value=1000
    )
    span_normalise = pn.widgets.Checkbox(name="Span normalise", value=True)
    mode = pn.widgets.RadioButtonGroup(
        name="Mode", options=["Site mode", "Branch mode"], value="Branch mode"
    )
    plot_options = pn.Column(
        pn.Row(pn.Column(windows_trees, span_normalise), mode),
        windows_count,
    )
    one_way_panel = pn.bind(
        make_one_way_stats,
        ts=ts,
        windows_trees=windows_trees,
        windows_count=windows_count,
        span_normalise=span_normalise,
        mode=mode,
    )
    windows_trees.jslink(windows_count, value="disabled")
    return pn.Column(plot_options, one_way_panel)
